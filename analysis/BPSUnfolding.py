import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.stats import poisson
from scipy.signal import convolve
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import ROOT
import uproot
import matplotlib
matplotlib.use('Agg')  # or 'Agg' for non-interactive
import os
import json
import time
ROOT.gSystem.Load("libRooUnfold")
from concurrent.futures import ProcessPoolExecutor

def numpy_to_th1(name, title, values, edges):
    values = np.ascontiguousarray(values, dtype=np.float64)
    edges = np.ascontiguousarray(edges, dtype=np.float64)

    n_bins = len(values)

    h = ROOT.TH1D(name, title, n_bins, edges)
    for i, v in enumerate(values):
        h.SetBinContent(i+1, v)
        h.SetBinError(i+1, np.sqrt(v))  # Poisson errors
    return h

def numpy_to_th2(name, title, R, xedges, yedges):
    
    R = np.ascontiguousarray(R, dtype=np.float64)
    xedges = np.ascontiguousarray(xedges, dtype=np.float64)
    yedges = np.ascontiguousarray(yedges, dtype=np.float64)    
    
    h = ROOT.TH2D(
        name, title,
        len(xedges)-1, xedges,
        len(yedges)-1, yedges
    )
    for ix in range(R.shape[0]):      # measured
        for iy in range(R.shape[1]):  # true
            h.SetBinContent(ix+1, iy+1, R[ix, iy])
    return h

def calcBinEdges(centers):
    centers = np.array(centers)
    edges = np.zeros(len(centers) + 1)
    edges[1:-1] = (centers[1:] + centers[:-1]) / 2
    edges[0] = centers[0] - (centers[1] - centers[0]) / 2
    edges[-1] = centers[-1] + (centers[-1] - centers[-2]) / 2
    return edges

def th1_to_numpy(h):
    """Return bin centers and bin contents"""
    bins = np.array([h.GetBinCenter(i+1) for i in range(h.GetNbinsX())])
    contents = np.array([h.GetBinContent(i+1) for i in range(h.GetNbinsX())])
    return bins, contents

class SiPM:
    """
    Simple container for SiPM single-photoelectron parameters.
    All values are assumed to be in consistent units (e.g. ADC counts or charge).
    """

    def __init__(self, gain, gain_err, gain_width, gain_width_err, elec_width):
        """
        Parameters
        ----------
        gain : float
            Mean gain (1 p.e. peak position)
        gain_width : float
            Width (sigma) of the 1 p.e. gain distribution
        elec_width : float
            Width (sigma) of the electronic noise (pedestal)
        """
        self.gain = float(gain)
        self.gain_err = float(gain_width)
        self.gain_width = float(gain_width)
        self.gain_width_err = float(gain_width)
        self.elec_width = float(elec_width)

    def __repr__(self):
        return (
            f"SiPM(gain={self.gain}, "
            f"gain_width={self.gain_width}, "
            f"elec_width={self.elec_width})"
        )

class ChannelCalibration:
    """
    Represents the energy calibration of one detector channel using two points.
    """

    def __init__(self, a, b, a_err, b_err, channelWidth, simVar, simStatVar, ADCVar, ADCStatVar, covAB, quenched, ADCchannel):
        self.kB = 12.68 * 0.001
        self.kbVariance = abs(self.kB*0.3)**2
        self.a = a
        self.b = b
        self.quenched = quenched
        self.ADCchannel = ADCchannel
        self.a_err = a_err
        self.b_err = b_err
        self.channelWidth = channelWidth
        self.simVar = simVar
        self.simStatVar = simStatVar
        self.ADCVar = ADCVar
        self.ADCStatVar = ADCStatVar
        self.covAB = covAB
        
    def linear_adc_to_energy(self, adc_value):
        return self.a * adc_value + self.b

    def adc_to_energy(self, adc_value):
        linearenergy = self.a * adc_value + self.b
        return 1/(1/linearenergy - self.kB/self.channelWidth)

    def energy_to_adc(self, energy):
        energy = energy/(1 + self.kB * energy/self.channelWidth)
        return (energy - self.b) / self.a

    def __repr__(self):
        return (
            f"ChannelCalibration(\n"
            f"  slope = {self.a:.6f}, offset = {self.b:.6f}\n"
        )

class Detector:
    """
    Represents a detector with multiple calibrated channels.
    """

    def __init__(self, n_channels=32):
        self.n_channels = n_channels
        self.channels = {}
        self.sipms = {}

    def add_channel(self, channel_id, calibration: ChannelCalibration, sipm : SiPM = None):
        if channel_id < 0 or channel_id >= self.n_channels:
            raise ValueError(f"Channel must be between 0 and {self.n_channels-1}")
        self.channels[channel_id] = calibration
        self.sipms[channel_id] = sipm

    def linear_adc_to_energy(self, channel_id, adc_value):
        return self.channels[channel_id].linear_adc_to_energy(adc_value)

    def adc_to_energy(self, channel_id, adc_value):
        return self.channels[channel_id].adc_to_energy(adc_value)

    def energy_to_adc(self, channel_id, energy):
        return self.channels[channel_id].energy_to_adc(energy)

    def sipm_gain(self, channel_id):
        sipm = self.sipms.get(channel_id)
        if sipm is None:
            raise ValueError(f"No SiPM parameters for channel {channel_id}")
        return sipm.gain

    def sipm_gain_width(self, channel_id):
        sipm = self.sipms.get(channel_id)
        if sipm is None:
            raise ValueError(f"No SiPM parameters for channel {channel_id}")
        return sipm.gain_width
    
    def sipm_elec_width(self, channel_id):
        sipm = self.sipms.get(channel_id)
        if sipm is None:
            raise ValueError(f"No SiPM parameters for channel {channel_id}")
        return sipm.elec_width
    
    def sipm_gain_err(self, channel_id):
        sipm = self.sipms.get(channel_id)
        if sipm is None:
            raise ValueError(f"No SiPM parameters for channel {channel_id}")
        return sipm.gain_err
    
    def sipm_gain_width_err(self, channel_id):
        sipm = self.sipms.get(channel_id)
        if sipm is None:
            raise ValueError(f"No SiPM parameters for channel {channel_id}")
        return sipm.gain_width_err

    def __repr__(self):
        return f"DetectorCalibration({len(self.channels)}/{self.n_channels} channels)"

class BPSUnfolding:
    """
    Bayesian Unfolding with a functional Detector Response Function (DRF)
    """
    
    def __init__(self, R, true_energies, measured_energies, n_iterations=10, channel = 0):
        self.R = R
        self.measured_energies = measured_energies
        self.true_energies = true_energies
        self.n_iterations = n_iterations
        self.channel = channel
            
    def unfold(self, measured_spectrum, prior, efficiency=None):
        """
        Unfold measured spectrum using Bayesian method
        
        Parameters:
        measured_spectrum: 1D array of measured counts (binned according to measured_edges)
        prior: initial guess for true spectrum
        efficiency: detection efficiency for each true bin
        
        Returns:
        unfolded_spectrum: estimated true spectrum
        """
        if not hasattr(self, 'R'):
            raise ValueError("Must build response matrix first using build_response_matrix()")
        
        g = measured_spectrum.astype(float)
        f = prior.copy()
        
        print(f"Sum g: {np.sum(g)}")
        if np.sum(g) != 0:                
            f = f / np.sum(f) * np.sum(g)
        else:
            print("No data -> No Unfolding")
            self.history = [f.copy()*0]
            self.chi2 = [0]
            self.closure_errors = [0]
            self.unfolded = f*0
            return self.unfolded
        
        if efficiency is None:
            efficiency = np.sum(self.R, axis=0)
        
        max_index = np.argmax(measured_spectrum)

        self.history = [f.copy()]
        self.chi2 = [0]
        self.closure_errors = [0]
        print(f"\nChannel: {self.channel} Starting Bayesian unfolding...")
        for iteration in range(self.n_iterations):
            g_expected = np.dot(self.R, f)
            
            g_expected[g_expected == 0] = 1e-10
            
            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency
            f_unsmoothed = f_new.copy()
            f_new = self._smooth_spectrum(f_new, strength=0.1)
            # f_new = self._tikhonov_regularize(f_new, f, 0.01)
            f_new = np.maximum(f_new, 0)
            
            # mask = (self.history[0] > 0) & (f_new > 0)
            # chi2 = np.sum((f_new[mask] - self.history[0][mask])**2 / self.history[0][mask])

            predicted_measured = np.dot(self.R, f_new)
            mask = (measured_spectrum > 0) & (predicted_measured > 0)
            # mask[:max_index+130] = False
            n_dof = np.sum(mask) - 1
            chi2 = np.sum((predicted_measured[mask] - measured_spectrum[mask])**2 / measured_spectrum[mask])/n_dof
            closure_error = np.sqrt(np.mean(((measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask])**2))
            self.chi2.append(chi2)
            self.closure_errors.append(closure_error)
            
            if(self.chi2[-2]-self.chi2[-1]<0.1 and iteration>3):
                print(f"Converged after {iteration+1} iterations")
                f = f_unsmoothed
                self.history.append(f.copy())
                
                break

            f = f_new
            self.history.append(f.copy())
                                
        self.unfolded = f
        return f

    def unfold_regularized(self, measured_spectrum, prior=None, regularization_strength=0.01, early_stopping=True):
        """
        Unfold with Tikhonov regularization for better convergence
        """
        if not hasattr(self, 'R'):
            raise ValueError("Must build response matrix first")

        g = measured_spectrum.astype(float)
        f = prior.copy() if prior is not None else np.ones(self.n_true_bins)

        # Normalize prior
        f = f / np.sum(f) * np.sum(g)

        efficiency = np.sum(self.R, axis=0)

        self.history = [f.copy()]
        self.closure_errors = []

        print("\nStarting regularized Bayesian unfolding...")
        for iteration in range(self.n_iterations):
            g_expected = np.dot(self.R, f)
            g_expected[g_expected == 0] = 1e-10

            correction_factor = np.dot(self.R.T, g / g_expected)
            f_new = f * correction_factor / efficiency

            f_new = self._tikhonov_regularize(f_new, f, regularization_strength)
            f_new = np.maximum(f_new, 0)

            relative_change = np.sum(np.abs(f_new - f)) / np.sum(f)

            f = f_new
            self.history.append(f.copy())

            pred_meas = np.dot(self.R, f)
            closure_err = np.sqrt(np.mean(((g - pred_meas) / (g + 1e-10))**2))
            self.closure_errors.append(closure_err)

            print(f"Iteration {iteration+1}: Total counts = {np.sum(f):.1f}, "
                  f"Change = {relative_change:.2e}, Closure = {closure_err:.4f}")

            # Early stopping
            if early_stopping and relative_change < 1e-4 and iteration > 3:
                print(f"Converged after {iteration+1} iterations")
                break
            
        self.unfolded = f
        return f

    def _tikhonov_regularize(self, f_new, f_old, strength):
        """Apply Tikhonov regularization to reduce oscillations"""
        # Simple gradient penalty: minimize curvature
        from scipy import sparse

        # Create second derivative matrix
        n = len(f_new)
        D2 = sparse.diags([1, -2, 1], [0, 1, 2], shape=(n-2, n)).toarray()

        # Calculate regularization term
        reg_term = strength * np.dot(D2.T, np.dot(D2, f_new))

        # Apply regularization (subtract gradient)
        f_regularized = f_new - reg_term

        return np.maximum(f_regularized, 0)

    def _smooth_spectrum(self, spectrum, strength=0.1):
        """Apply mild smoothing to reduce oscillations"""
        from scipy.ndimage import gaussian_filter1d
        return (1 - strength) * spectrum + strength * gaussian_filter1d(spectrum, sigma=10.0)
    
    def calculate_closure(self, measured_spectrum):
        """Calculate how well the unfolded spectrum reproduces the measurement"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)
            
            # Avoid division by zero
            mask = measured_spectrum > 0
            
            # Calculate normalized chi-squared
            chi2 = np.sum((measured_spectrum[mask] - predicted_measured[mask])**2 / measured_spectrum[mask])
            n_dof = np.sum(mask) - 1  # degrees of freedom
            
            # Calculate relative closure error
            closure_error = np.sqrt(np.mean(((measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask])**2))
            
            return closure_error, chi2/n_dof, predicted_measured
        else:
            raise ValueError("Must run unfold() first")
    
    def calculate_closure_fixed(self, measured_spectrum):
        """Proper closure calculation that handles zero bins correctly"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)

            # Method 1: Use only bins with sufficient counts
            threshold = np.max(measured_spectrum) * 0.01  # 1% of peak
            mask = measured_spectrum > threshold

            if np.sum(mask) < 5:  # If too few bins, use non-zero bins
                mask = measured_spectrum > 0

            # Calculate metrics only on meaningful bins
            if np.sum(mask) > 0:
                # Relative RMS error on meaningful bins
                rel_errors = np.abs(measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask]
                closure_rms = np.sqrt(np.mean(rel_errors**2))

                # Chi-squared per DOF
                chi2 = np.sum((measured_spectrum[mask] - predicted_measured[mask])**2 / measured_spectrum[mask])
                chi2_per_dof = chi2 / len(mask)

                # Correlation
                correlation = np.corrcoef(measured_spectrum[mask], predicted_measured[mask])[0, 1]
            else:
                closure_rms = float('inf')
                chi2_per_dof = float('inf')
                correlation = 0

            # Method 2: Overall shape metrics (less sensitive to zeros)
            total_measured = np.sum(measured_spectrum)
            total_predicted = np.sum(predicted_measured)
            total_error = np.abs(total_measured - total_predicted) / total_measured

            # Method 3: Wasserstein distance (shape similarity)
            from scipy.stats import wasserstein_distance
            wasserstein = wasserstein_distance(measured_spectrum, predicted_measured)

            print(f"\n=== ROBUST CLOSURE METRICS ===")
            print(f"Bins used for metrics: {np.sum(mask)}/{len(measured_spectrum)}")
            print(f"Relative RMS Error: {closure_rms:.4f} ({closure_rms*100:.1f}%)")
            print(f"Chi-squared per DOF: {chi2_per_dof:.2f}")
            print(f"Correlation: {correlation:.4f}")
            print(f"Total Counts Error: {total_error:.4f} ({total_error*100:.1f}%)")
            print(f"Wasserstein Distance: {wasserstein:.4f}")

            return {
                'rms_error': closure_rms,
                'chi2_per_dof': chi2_per_dof,
                'correlation': correlation,
                'total_error': total_error,
                'wasserstein': wasserstein,
                'predicted': predicted_measured,
                'mask': mask
            }

    def check_response_matrix(self):
        """Diagnose potential issues with the response matrix"""
        print("\nResponse Matrix Diagnostics:")
        print(f"Shape: {self.R.shape}")
        print(f"Min value: {np.min(self.R):.2e}")
        print(f"Max value: {np.max(self.R):.2e}")
        print(f"Matrix condition number: {np.linalg.cond(self.R):.2e}")

        # Check for columns with very small sums (ill-conditioned)
        column_sums = np.sum(self.R, axis=0)
        zero_columns = np.sum(column_sums < 1e-10)
        print(f"Columns with near-zero sum: {zero_columns}")

        if zero_columns > 0:
            print("WARNING: Response matrix has columns with near-zero sum!")
            print("This can cause convergence issues.")

    def calculate_closure_robust(self, measured_spectrum):
        """Robust closure calculation that handles edge cases"""
        if hasattr(self, 'unfolded'):
            predicted_measured = np.dot(self.R, self.unfolded)

            # Method 1: Normalized RMS (handles zeros)
            diff = measured_spectrum - predicted_measured
            rms_error = np.sqrt(np.mean(diff**2))
            normalization = np.sqrt(np.mean(measured_spectrum**2))
            closure_rms = rms_error / normalization if normalization > 0 else float('inf')

            # Method 2: Relative error in total counts
            total_measured = np.sum(measured_spectrum)
            total_predicted = np.sum(predicted_measured)
            closure_total = np.abs(total_measured - total_predicted) / total_measured

            # Method 3: Correlation coefficient
            correlation = np.corrcoef(measured_spectrum, predicted_measured)[0, 1]

            # Method 4: Exclude very low count bins
            mask = measured_spectrum > np.max(measured_spectrum) * 0.01  # Exclude lowest 1%
            if np.sum(mask) > 10:  # Only use if we have enough bins
                rel_error_masked = np.sqrt(np.mean(((measured_spectrum[mask] - predicted_measured[mask]) / measured_spectrum[mask])**2))
            else:
                rel_error_masked = float('inf')

            print(f"\nClosure Diagnostics:")
            print(f"RMS Error: {closure_rms:.4f} ({closure_rms*100:.1f}%)")
            print(f"Total Counts Error: {closure_total:.4f} ({closure_total*100:.1f}%)")
            print(f"Correlation: {correlation:.4f}")
            print(f"Masked Relative Error: {rel_error_masked:.4f} ({rel_error_masked*100:.1f}%)")

            return {
                'rms_error': closure_rms,
                'total_error': closure_total, 
                'correlation': correlation,
                'masked_error': rel_error_masked,
                'predicted': predicted_measured
            }
        else:
            raise ValueError("Must run unfold() first")

    def debug_closure_issue(self, measured_spectrum, predicted_measured):
        """Debug why closure error is huge but plots look good"""
        print("\n=== CLOSURE DEBUGGING ===")
        print(f"Measured spectrum sum: {np.sum(measured_spectrum):.1f}")
        print(f"Predicted spectrum sum: {np.sum(predicted_measured):.1f}")
        print(f"Ratio (predicted/measured): {np.sum(predicted_measured)/np.sum(measured_spectrum):.3f}")

        # Check individual bins
        diff = measured_spectrum - predicted_measured
        rel_error = np.abs(diff) / (measured_spectrum + 1e-10)  # Avoid division by zero

        print(f"\nTop 5 largest relative errors:")
        worst_indices = np.argsort(rel_error)[-5:][::-1]
        for idx in worst_indices:
            print(f"  Bin {idx}: Measured={measured_spectrum[idx]:.1f}, "
                  f"Predicted={predicted_measured[idx]:.1f}, "
                  f"Error={rel_error[idx]*100:.1f}%")

        # Check for zeros in measured spectrum
        zero_bins = np.sum(measured_spectrum == 0)
        print(f"\nBins with zero measured counts: {zero_bins}")

        if zero_bins > 0:
            print("WARNING: Zero bins in measured spectrum can inflate closure error!")
    
    def saveClosurePlot(self):
        plt.plot(range(len(unfolder.chi2))[1:], unfolder.chi2[1:], 'bo-')
        plt.xlabel('Iteration')
        plt.ylabel('Chi-squared/N_DOF')
        plt.title('Convergence of Unfolding')
        plt.grid(True)

        plt.tight_layout()
        plt.show()
        plt.close()
    
    def saveUnfoldedComparison(self, true_energies, measured_spectrum, unfolded_spectrum, show=False):
        plt.step(true_energies, measured_spectrum, color=targetColorMap[0], linewidth=2, label='Energy converted Spectrum')
        plt.plot(true_energies, unfolded_spectrum, color=targetColorMap[1], linewidth=4, label=f'Unfolded Energy Spectrum')

        plt.xlabel('Energy / MeV')
        plt.ylabel('Counts')
        plt.title('Measured vs Unfolded Spectrum')
        plt.legend()
        plt.grid(True)
        # plt.savefig(output_filename, dpi=300)
        if(show):
            plt.show()
        plt.close()
    
    def savePredictedComparison(self, measured_energies, measured_spectrum, predicted_measured, show=False):
        plt.step(measured_energies, measured_spectrum, color=targetColorMap[0], linewidth=2, label='Measured Charge Spectrum')
        plt.plot(measured_energies, predicted_measured, color=targetColorMap[1], linewidth=4, label=f'Predicted Charge Spectrum from Unfolded')
        plt.xlabel('Charge / ADC Count')
        plt.ylabel('Counts')
        plt.title('Measured vs Predicted Spectrum')
        plt.legend()
        plt.grid(True)
        # plt.savefig(output_filename, dpi=300)
        if(show):
            plt.show()
        plt.close()

    def plotResponseMatrix(self, R, trueEnergies, measuredEnergies):
        im = plt.imshow(R, aspect='auto', origin='lower',
                    extent=[trueEnergies[0], trueEnergies[-1],
                            measuredEnergies[0], measuredEnergies[-1]],
                    norm=LogNorm(vmin=1e-4, vmax=1))

        plt.xlabel('True Photon Count')
        plt.ylabel('Measured Photon Count')
        plt.title('Detector Response Matrix')
        plt.colorbar(im)
        plt.show()

def rebin_measured_axis(R, factor):
    """
    Rebin response matrix along measured axis only.
    R shape: (n_meas, n_true)
    """
    n_meas, n_true = R.shape
    assert n_meas % factor == 0

    R_rebinned = R.reshape(
        n_meas // factor, factor, n_true
    ).sum(axis=1)

    return R_rebinned

def rebin_hist(hist, factor):
    assert len(hist) % factor == 0
    return hist.reshape(-1, factor).sum(axis=1)


targetColorMap = ["#000000","#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]

if __name__ == "__main__":
    def UnfoldChannel(ch, hist_ch, ehist_ch, hist_edges_ch, eedges_ch, detector, datapath, detectorpath, depth, deptherr):
        # print(f"Channel: {ch}")
        # print(f"Depth: {depth} ± {deptherr} mm")
        # print(f"Histogram sum: {hist_ch}")
        # print(f"Energy Histogram sum: {ehist_ch}")
        # print(f"Histogram edges: {hist_edges_ch}")
        # print(f"detectorpath {detectorpath}")
        # print(f"datapath {datapath}")
        n_iter = 15
        if(np.sum(hist_ch)<50):
            n_iter = 1
        start = time.perf_counter()
        data = np.load(f"{detectorpath}/responseMatrix_CH{ch}.npz")

        R  = data["response_matrix"]
        Ecov = data["Energy_covariance_matrix"]
        E_true = data["true_energy"]
        initial_prior = np.ones_like(E_true)
        measured = hist_ch
        Emeasured = ehist_ch
        Emeas_edges = eedges_ch
        true_edges = np.array(calcBinEdges(E_true))

        h_meas  = numpy_to_th1(f"h_meas{ch}", "Measured", measured, hist_edges_ch)
        h_Emeas = numpy_to_th1(f"h_Emeas{ch}", "Energy Measured", Emeasured, Emeas_edges)

        h_true  = numpy_to_th1(f"h_true{ch}", "Truth", initial_prior, true_edges)
        
        h_resp  = numpy_to_th2(f"h_resp{ch}", "Response", R, hist_edges_ch, true_edges)
        response = ROOT.RooUnfoldResponse(h_meas, h_true, h_resp)
        unfold = ROOT.RooUnfoldBayes(response , h_meas, n_iter)
        h_unfold = unfold.Hunfold()
        print(f"Unfolding complete of channel {ch}.")
        cov_matrix = unfold.Eunfold()
        print("Covariance matrix retrieved.")
        
        x_meas, y_meas = th1_to_numpy(h_meas)
        x_Emeas, y_Emeas = th1_to_numpy(h_Emeas)
        x_true, y_true = th1_to_numpy(h_true)
        x_unf, y_unf = th1_to_numpy(h_unfold)

        n_bins = h_unfold.GetNbinsX()
        cov = np.zeros((n_bins, n_bins))

        for i in range(n_bins):
            for j in range(n_bins):
                cov[i, j] = cov_matrix[i][j]
        bin_errs = np.sqrt(np.diag(cov))

        end = time.perf_counter()
        print(f"Time taken: {end - start:.6f} s")

        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        plt.step(x_Emeas, y_Emeas, where='mid', label='Energy Converted Spectrum', color=targetColorMap[0], linewidth=2)
        plt.plot(x_unf, y_unf, label='Unfolded Energy Spectrum', color=targetColorMap[1], linewidth=4)
        # plt.fill_between(x_unf, y_unf - bin_errs, y_unf + bin_errs, color='tab:red', alpha=0.5)
        plt.xlabel('Deposited Energy / MeV')
        plt.ylabel('Counts')
        # plt.title(f'Channel {ch}: Unfolded Spectrum vs Measured Energy Spectrum')
        plt.legend()
        plt.grid(True)
        # plt.savefig(f"{datapath}/Redunfold/img/UnfoldedMeasured_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/UnfoldedMeasured_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/UnfoldedMeasured_CH{ch}.svg", bbox_inches='tight')
        plt.close()
        
        h_predicted = response.ApplyToTruth(h_unfold)# or .Hmeas()
        x_pred, y_pred = th1_to_numpy(h_predicted)
        
        y_meas_err = np.array([h_meas.GetBinError(i+1) for i in range(len(y_meas))])

        # Calculate χ² between measurement and refolded distribution
        # Using measurement errors in denominator
        chi2_numerator = 0.0
        ndof = 0

        for i in range(len(y_meas)):
            # Skip bins with zero or very small measurement error
            if y_meas_err[i] > 0 and y_meas[i] > 0:
                residual = y_meas[i] - y_pred[i]
                chi2_contribution = (residual ** 2) / (y_meas_err[i] ** 2)
                chi2_numerator += chi2_contribution
                ndof += 1

        # For Bayesian unfolding, adjust degrees of freedom for regularization
        # A simple approximation: ndf = number of bins * (1 - 1/n_iter)
        regularization_factor = 1.0 / n_iter
        effective_ndof = ndof * (1.0 - regularization_factor)

        # Ensure we don't divide by zero
        effective_ndof = max(effective_ndof, 1)

        reduced_chi2 = chi2_numerator / effective_ndof
        print(f"Reduced Chi2: {reduced_chi2}")
        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        plt.step(x_meas, y_meas, where='mid', label='Measured Charge', color=targetColorMap[0], linewidth=2)
        plt.plot(x_pred, y_pred, label='Predicted Charge from Unfolded', color=targetColorMap[1], linewidth=4)
        plt.xlabel('Charge / ADC Counts')
        plt.ylabel('Counts')
        # plt.title(f'Channel {ch}: Measured vs Predicted')
        plt.legend()
        plt.grid(True)
        # plt.savefig(f"{datapath}/Redunfold/img/PredictedMeasured_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/PredictedMeasured_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/PredictedMeasured_CH{ch}.svg", bbox_inches='tight')
        plt.close()

        conv = 1.60218e-13

        volume = detector.channels[ch].channelWidth*0.1*3*3 # cm^3
        mass = volume*8.28*0.001 #kg

        dose = 0
        unfoldedDose=0
        for idx, value in enumerate(x_unf):
            unfoldedDose += value * y_unf[idx]
        
        for idx, value in enumerate(x_Emeas):
            dose += value * y_Emeas[idx]

        unfoldedDose = unfoldedDose / mass * 1e6 * conv  # μGy
        dose = dose / mass * 1e6 * conv  # μGy

        Ncov = np.zeros((n_bins, n_bins))
        for i in range(n_bins):
            Ncov[i, i] = y_Emeas[i]

        DoseVariance = np.dot(x_Emeas, np.dot(cov, x_Emeas))+np.dot(y_Emeas, np.dot(Ecov, y_Emeas))
        DoseVariance *= conv**2 / mass**2 * 1e6**2 


        E = np.array(x_unf)
        N = np.array(y_unf) 
        
        unfoldedDoseVariance = np.dot(E, np.dot(cov, E))+np.dot(N, np.dot(Ecov, N))
        # doseVariance = np.sum(E[:, None] * E[None, :] * cov)+np.sum(N[:, None] * N[None, :] * Ecov)
        unfoldedDosePartial = np.dot(N, np.dot(Ecov, N))
        unfoldedDosePartial *= conv**2 / mass**2 * 1e6**2 

        unfoldedDoseVariance *= conv**2 / mass**2 * 1e6**2 

        print(f"Unfolded dose for channel {ch}: {unfoldedDose:.2f} ± {np.sqrt(unfoldedDoseVariance):.2f} μGy")
        print(f"Normally calculated dose for channel {ch}: {dose:.2f} ± {np.sqrt(DoseVariance):.2f} μGy")

        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        im = plt.imshow(Ecov, aspect='auto', origin='lower',
                    extent=[x_true[0], x_true[-1],
                            x_true[0], x_true[-1]],
                    norm=SymLogNorm(linthresh=1e-6,vmin=Ecov.min(), vmax=Ecov.max()))

        plt.xlabel('Energy Deposition / MeV')
        plt.ylabel('Energy Deposition / MeV')
        plt.title('Energy Covariance Matrix')
        plt.colorbar(im)
        plt.savefig(f"{datapath}/unfold/img/EnergyCovariance_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/EnergyCovariance_CH{ch}.svg", bbox_inches='tight')
        # plt.savefig(f"{datapath}/Redunfold/img/EnergyCovariance_CH{ch}.pdf", bbox_inches='tight')
        plt.close()
        
        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        im = plt.imshow(cov, aspect='auto', origin='lower',
                    extent=[x_meas[0], x_meas[-1],
                            x_meas[0], x_meas[-1]],
                    norm=SymLogNorm(linthresh=1e-6,vmin=cov.min(), vmax=cov.max()))

        plt.xlabel('Charge / ADC Count')
        plt.ylabel('Charge / ADC Count')
        plt.title('Unfold Covariance Matrix')
        plt.colorbar(im)
        plt.savefig(f"{datapath}/unfold/img/UnfoldCovariance_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/UnfoldCovariance_CH{ch}.svg", bbox_inches='tight')
        # plt.savefig(f"{datapath}/Redunfold/img/UnfoldCovariance_CH{ch}.pdf", bbox_inches='tight')
        plt.close()
        
        plt.rcParams.update({'font.size': 32})
        plt.figure(figsize=(16,12))
        im = plt.imshow(R, aspect='auto', origin='lower',
            extent=[x_true[0], x_true[-1],
                    x_meas[0], x_meas[-1]],
            norm=SymLogNorm(linthresh=1e-6,vmin=R.min(), vmax=R.max()))

        plt.xlabel('Energy Deposition / MeV')
        plt.ylabel('Charge / ADC Count')
        plt.colorbar(im)
        plt.savefig(f"{datapath}/unfold/img/ResponseMatrix_CH{ch}.pdf", bbox_inches='tight')
        plt.savefig(f"{datapath}/unfold/img/ResponseMatrix_CH{ch}.svg", bbox_inches='tight')
        # plt.savefig(f"{datapath}/Redunfold/img/ResponseMatrix_CH{ch}.pdf", bbox_inches='tight')
        plt.close()
    
        np.savez(
            f"{datapath}/unfold/data/Unfolded{ch}.npz",
            # f"{datapath}/Redunfold/data/Unfolded{ch}.npz",
            # covarianceMatrix=cov,
            # EcovarianceMatrix=Ecov,
            # counts_variance = Ncov,
            Phi0 = np.sum(hist),
            true_energy=[x_true, y_true],
            measured_energy=[x_meas,y_meas],
            unfolded=[x_unf, y_unf],
            predicted_measured=[x_pred, y_pred],
            unfoldedDose =[unfoldedDose, unfoldedDoseVariance],
            dose =[dose, DoseVariance],
            depth = [depth, deptherr],
            n_iter_used = n_iter,
            reduced_chi2 = reduced_chi2
        )

        h_meas.Delete()
        h_true.Delete()
        h_resp.Delete()
        h_unfold.Delete()
        cov_matrix.Delete()
        
        del Ncov
        del response
        del unfold
        del Ecov
        del cov
        return ch
    
    np.random.seed(42)
    
    with open("config.json", "r") as fileT:
        fullConfig = json.load(fileT)

    detectorSelect = fullConfig["detectorSelect"]
    targetSelect = fullConfig["targetSelect"]
    plotEnable = fullConfig["plotEnable"]

    config = fullConfig["detectors"][detectorSelect]

    datasetSelect = config["datasetSelect"]
    
    nLayers = config["nLayers"]
    
    targetThickness = config["targetThickness"]

    datasets = ["MIT_05_2024", "simulation", "paperBeamtime"]
    in_data = ["notarget", "homotarget", "heterotarget"]
    in_title = ["without a target", "with the homogeneous target", "with the heterogeneous target"]

    dataset = datasets[datasetSelect]
    print(f"Dataset: {dataset}")
    filename = in_data[targetSelect]
    if(targetSelect == 1):
        filename = f"{filename}{targetThickness}"

    datapath = f"../data/paperBeamtime/{filename}/output"
    detectorpath = "../data/paperBeamtime/detector"
    
    detector = Detector()
    speParams = np.load(f"{detectorpath}/SPEparams.npy")
    calibParamsFront = np.load(f"{detectorpath}/calibParamsFront.npy")
    calibParamsEnd = np.load(f"{detectorpath}/calibParamsEnd.npy")
    channelWidth = 3
    
    print(f"BPS Unfolding of {filename}")
    print("=" * 50)
    
    print("Loading measured spectrum...")   

    for crystal in range(32):
        if crystal > 10:
            channelWidth = 2
            
        if crystal < 29:
            calib = ChannelCalibration(
            a = calibParamsFront[crystal][0],
            b = 0,
            a_err = 0,
            b_err = 0,
            ADCchannel=calibParamsFront[crystal][1],
            quenched=calibParamsFront[crystal][2],
            simVar=calibParamsFront[crystal][3],
            simStatVar=calibParamsFront[crystal][4],
            ADCVar=calibParamsFront[crystal][5],
            ADCStatVar=calibParamsFront[crystal][6],
            covAB=0, 
            channelWidth=channelWidth
        )

        else:
            calib = ChannelCalibration(
            a = calibParamsEnd[crystal-29][0],
            b = calibParamsEnd[crystal-29][2],
            a_err = calibParamsEnd[crystal-29][1],
            b_err = calibParamsEnd[crystal-29][3],
            covAB=calibParamsEnd[crystal-29][4],
            ADCchannel=0,
            quenched=0,
            simVar=0,
            simStatVar=0,
            ADCVar=0,
            ADCStatVar=0,
            channelWidth=channelWidth
        )

        sipm = SiPM(
            gain = speParams[crystal][0],
            gain_err = speParams[crystal][1],
            gain_width = speParams[crystal][2],
            gain_width_err = speParams[crystal][3],
            elec_width = speParams[crystal][4],
        )
    
        detector.add_channel(crystal, calibration=calib, sipm=sipm)
    
    file = uproot.open(f"{datapath}/coincidence_output.root")
    tree = file["events"]

    charge = tree["Charge"].array(library="np")   # shape (N_events, 32)
    
    hist = []
    ehist = []
    hist_edges = []
    eedges = []

    energy = []
    selectedCharge = []
    idx = 0

    for event in charge:
        idx += 1
        energyEvent = []
        for layer, layerCharge in enumerate(event):
            if(layerCharge == 0):# or layer == 8):
                energyEvent.append(0)
            else:
                energyEvent.append(detector.adc_to_energy(layer, layerCharge))
        energy.append(np.sum(energyEvent))
    
    TEnergyCounts, TEnergyEdges = np.histogram(energy, bins=2000)
    TEnergyCenters = 0.5 * (TEnergyEdges[:-1] + TEnergyEdges[1:])
    peak_idx = np.argmax(TEnergyCounts)
    peak_x = TEnergyCenters[peak_idx]
    fit_mask = TEnergyCenters > peak_x - 30
    x_fit = TEnergyCenters[fit_mask]
    y_fit = TEnergyCounts[fit_mask]

    def gauss(x, A, mu, sigma):
        return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

    p0 = [TEnergyCounts[peak_idx], peak_x, np.std(energy) / 5]

    popt, pcov = curve_fit(gauss, x_fit, y_fit, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    energy = []
    A, mu, sigma = popt
    FWHM = 2.355 * sigma
    
    energyCutOff = mu - (FWHM) / 2
    if(targetThickness == 0):
        energyCutOff = 195
    elif (targetThickness == 52):
        energyCutOff = mu - 2.5*sigma
    elif (targetThickness == 53):
        energyCutOff = mu - 2.5*sigma   
    elif (targetThickness == 54):
        energyCutOff = mu - 2.5*sigma   
    elif (targetThickness == 55):
        energyCutOff = mu - 2.5*sigma   
    elif (targetThickness == 56):
        energyCutOff = mu - 2.5*sigma   
    elif (targetThickness == 57):
        energyCutOff = mu - 2.5*sigma   
    elif (targetThickness == 200):
        energyCutOff = mu - 2.5*sigma   
        
    idx = 0
    for event in charge:
        idx += 1
        energyEvent = []
        for layer, layerCharge in enumerate(event):
            if(layerCharge == 0):
                energyEvent.append(0)
            else:
                energyEvent.append(detector.adc_to_energy(layer, layerCharge))
        if(np.sum(energyEvent)>energyCutOff):
            selectedCharge.append(event)
            energy.append(np.sum(energyEvent))
    
    mask = TEnergyCenters >= energyCutOff
    TCutEnergyCounts = TEnergyCounts[mask]
    TCutEnergyCenters = TEnergyCenters[mask]

    x_smooth = np.linspace(TEnergyCenters[0], TEnergyCenters[-1], 1000)
    y_smooth = gauss(x_smooth, *popt)

    print(f"Peak at {mu+8.1} MeV +. {perr[1]}: Energycut at {energyCutOff} MeV")
    print(f"Difference: {207.7-(mu+8.1)}")
    print(f"{221.6*(mu+8.1)/207.7}")
    plt.rcParams.update({'font.size': 32})
    plt.figure(figsize=(16,12))
    plt.step(TEnergyCenters, TEnergyCounts, where="mid", color=targetColorMap[0])
    # plt.step(TCutEnergyCenters, TCutEnergyCounts, where="mid", color='green')
    # plt.plot(x_smooth, y_smooth, label="Gaussian fit")
    plt.axvline(
    x=221.6,
    linestyle="--",
    linewidth=5,
    color=targetColorMap[1],
    label="Nominal energy (221.6 MeV)")
    plt.grid()
    plt.xlabel("Total Deposited Energy / MeV")
    plt.ylabel("Counts")
    plt.tight_layout()
    plt.savefig(f"{datapath}/unfold/img/TotalEnergyDeposition.pdf", format="pdf", bbox_inches="tight")
    plt.savefig(f"{datapath}/unfold/img/TotalEnergyDeposition.svg", format="svg", bbox_inches="tight")
    # plt.savefig(f"{datapath}/Redunfold/img/TotalEnergyDeposition.pdf", format="pdf", bbox_inches="tight")
    # plt.show()
    plt.close()

    # exit()
    selectedCharge = np.array(selectedCharge)
    
    for i in range(0,32):
        q = selectedCharge[:, i]
        q_nz = q[q != 0]
        c, e = np.histogram(q_nz, bins=2048, range=(0, 65536/2))
        hist.append(c)
        hist_edges.append(e)
        # c, e = np.histogram(detector.linear_adc_to_energy(i, q_nz), bins=2048, range=(0, detector.linear_adc_to_energy(i, 65536/2)))
        
        c, e = np.histogram(detector.adc_to_energy(i, q_nz), bins=2048, range=(0, detector.adc_to_energy(i, 65536/2)))
        ehist.append(c)
        eedges.append(e)
    
    hist = np.array(hist) 
    ehist = np.array(ehist) 
    eCenters = []
    Centers = []
    for value in hist_edges:
        Centers.append((value[1:] + value[:-1]) / 2)
    for value in eedges:
        eCenters.append((value[1:] + value[:-1]) / 2)  

    print(f"Getting depth information from ../data/{dataset}/{filename}/output/{filename}Means.csv")
    data = np.loadtxt(f"../data/{dataset}/{filename}/output/{filename}Means.csv", delimiter=",")
    depth, deptherr, mean, error = data.T  # transpose to get separate vectors
    
    depth = depth/10
    deptherr = deptherr/10
    


    print("\nPerforming Bayesian unfolding...")
    
    channels = list(range(0, 1))
    # UnfoldChannel(0, hist[0], ehist[0], hist_edges[0], eedges[0], detector, datapath, detectorpath, depth[0], deptherr[0])
    # exit()
    args = [(hist[ch], ehist[ch], hist_edges[ch], eedges[ch], detector, datapath, detectorpath, depth, deptherr) for ch in channels]


    with ProcessPoolExecutor(max_workers=1) as executor:
        for ch in executor.map(UnfoldChannel,
        channels,
        [hist[ch] for ch in channels],
        [ehist[ch] for ch in channels],
        [hist_edges[ch] for ch in channels],
        [eedges[ch] for ch in channels],
        [detector]* len(channels),
        [datapath]* len(channels),
        [detectorpath]* len(channels),
        [depth[ch] for ch in channels],
        [deptherr[ch] for ch in channels],
    ):
            print(f"Finished channel {ch}")

    
    unfoldedDoses = []
    Doses = []
    old = False
    for ch in range(0,32):
        if(old):
            data = np.load(f"{detectorpath}/responseMatrix_CH{ch}.npz")

            R  = data["response_matrix"]
            Ecov = data["Energy_covariance_matrix"]
            E_true = data["true_energy"]
            E_meas = data["measured_energy"]
            measured = hist[ch]
            initial_prior = np.ones_like(E_true)
            unfolder = BPSUnfolding(R, E_true, E_meas, n_iterations=100, channel = ch)


            unfolded_spectrum = unfolder.unfold(measured , prior=initial_prior)
            #unfolded_spectrum = unfolder.unfold_regularized(measured_spectrum, prior=measured_spectrum)
            predicted_measured = np.dot(unfolder.R, unfolded_spectrum)

            # unfolder.saveUnfoldedComparison(E_true, ehist[ch], unfolded_spectrum, True)

            fig, axs = plt.subplots(1, 2, figsize=(20, 10))

            axs[0].plot(E_true, unfolded_spectrum, color=targetColorMap[0], linewidth=4, label='Unfolded Spectrum')
            axs[0].step(eCenters[ch], ehist[ch], color=targetColorMap[1], linewidth=2, label=f'Measured Spectrum')

            axs[0].set_xlabel('Energy / MeV')
            axs[0].set_ylabel('Counts')
            axs[0].set_title('Measured vs Unfolded Spectrum')
            axs[0].legend()
            axs[0].grid(True)

            axs[1].step(Centers[ch], measured, color=targetColorMap[0], linewidth=2, label=f'Measured Spectrum')
            axs[1].plot(Centers[ch], predicted_measured, color=targetColorMap[1], linewidth=4, label='Predicted Spectrum')

            axs[1].set_xlabel('Energy / MeV')
            axs[1].set_ylabel('Counts')
            axs[1].set_title('Measured vs Unfolded Spectrum')
            axs[1].legend()
            axs[1].grid(True)

            plt.savefig(f"{detectorpath}/preoutput/unfoldedSpectrumCH{ch}.pdf", format="pdf", bbox_inches='tight')
            plt.close()

            conv = 1.60218e-13

            unfoldedDose=0
            Dose=0
            volume = detector.channels[ch].channelWidth*0.1*3*3 # cm^3
            mass = volume*8.28*0.001 #kg

            for idx, value in enumerate(E_true):
                unfoldedDose += value*unfolded_spectrum[idx]
            unfoldedDoses.append(unfoldedDose/mass*1e6*conv)

            for idx, value in enumerate(eCenters[ch]):
                dose += value*ehist[ch][idx]
            Doses.append(dose/mass*1e6*conv)

