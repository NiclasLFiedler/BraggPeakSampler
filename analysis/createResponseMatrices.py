# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import curve_fit

# %%
from scipy.stats import poisson, binom, norm
from scipy.signal import convolve
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import ROOT
import uproot
import matplotlib
matplotlib.use('TkAgg')  # or 'Agg' for non-interactive
import os


class ChannelCalibration:
    """
    Represents the energy calibration of one detector channel using two points.
    """

    def __init__(self,
                 E1, E2,
                 adc1, adc2,
                 std1=None, std2=None,
                 err1=None, err2=None):
        """
        Parameters:
        E1, E2     : known energies [MeV or keV]
        adc1, adc2 : corresponding ADC means
        std1, std2 : standard deviations of ADC peaks
        err1, err2 : error of ADC means
        """

        self.E1 = E1
        self.E2 = E2
        self.adc1 = adc1
        self.adc2 = adc2
        self.std1 = std1
        self.std2 = std2
        self.err1 = err1
        self.err2 = err2

        # Calculate linear calibration: E = m * ADC + b
        self.m, self.b = self._calculate_calibration()

    def _calculate_calibration(self):
        m = (self.E2 - self.E1) / (self.adc2 - self.adc1)
        b = self.E1 - m * self.adc1
        return m, b

    def adc_to_energy(self, adc_value):
        return self.m * adc_value + self.b

    def energy_to_adc(self, energy):
        return (energy - self.b) / self.m

    def __repr__(self):
        return (
            f"ChannelCalibration(\n"
            f"  slope = {self.m:.6f}, offset = {self.b:.6f}\n"
            f"  Points: ({self.adc1} → {self.E1}), ({self.adc2} → {self.E2})\n)"
        )

class DetectorCalibration:
    """
    Represents a detector with multiple calibrated channels.
    """

    def __init__(self, n_channels=1):
        self.n_channels = n_channels
        self.channels = {}

    def add_channel(self, channel_id, calibration: ChannelCalibration):
        if channel_id < 0 or channel_id >= self.n_channels:
            raise ValueError(f"Channel must be between 0 and {self.n_channels-1}")
        self.channels[channel_id] = calibration

    def adc_to_energy(self, channel_id, adc_value):
        return self.channels[channel_id].adc_to_energy(adc_value)

    def energy_to_adc(self, channel_id, energy):
        return self.channels[channel_id].energy_to_adc(energy)

    def __repr__(self):
        return f"DetectorCalibration({len(self.channels)}/{self.n_channels} channels)"
    
class ResponseMatrixBuilder:
    """
    ResponseMatrixBuilder
    """
    
    def __init__(self, detector_calibration: DetectorCalibration, n_measured_bins=4001, n_iterations=10):
        self.n_measured_bins = n_measured_bins
        self.n_iterations = n_iterations
        self.detector_calibration = detector_calibration

        trueBins = []
        for i in range(self.detector_calibration.n_channels):
            trueBins.append(self.detector_calibration.adc_to_energy(i, n_measured_bins))
        
        self.n_true_bins = trueBins
    
    def drf_function(self, measured_energy, true_energy, resolution=0.15, tail_strength=0.2):
        """
        EXAMPLE DRF function - REPLACE THIS WITH YOUR ACTUAL DRF!
        
        This is a realistic DRF for a scintillator-SiPM system:
        - Gaussian main peak with energy-dependent resolution
        - Low-energy tail from incomplete light collection
        
        Parameters:
        measured_energy: the observed energy (MeV)
        true_energy: the true deposited energy (MeV)
        resolution: energy resolution (σ/E) at 1 MeV
        tail_strength: relative strength of low-energy tail
        
        Returns:
        probability density for measuring 'measured_energy' given 'true_energy'
        """
        # Energy resolution typically scales as 1/sqrt(E)
        sigma = resolution * np.sqrt(true_energy) * true_energy
        
        # Main Gaussian peak
        gaussian = np.exp(-0.5 * ((measured_energy - true_energy) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))
        
        # Low-energy tail (common in scintillators due to light collection effects)
        if measured_energy < true_energy:
            tail = tail_strength * np.exp(-(true_energy - measured_energy) / (0.1 * true_energy)) / true_energy
        else:
            tail = 0
            
        return gaussian + tail
    
    def build_response_matrix(self, custom_drf=None):
        """
        Build response matrix from DRF function
        
        Parameters:
        true_energy_range: tuple (min, max) for true energies
        measured_energy_range: tuple (min, max) for measured energies  
        custom_drf: your custom DRF function. If None, uses the example above
        
        Returns:
        response_matrix: R[i,j] = P(measured_bin i | true_bin j)
        true_energies: center of true energy bins
        measured_energies: center of measured energy bins
        """
        # Use custom DRF if provided, otherwise use example
        drf = custom_drf if custom_drf is not None else self.drf_function
        
        # Create energy bins
        self.true_edges = np.linspace(0, true_energy_range[1], self.n_true_bins + 1)

    
        self.measured_edges = np.linspace(0, self.n_measured_bins, self.n_measured_bins + 1)
        
        self.true_centers = (self.true_edges[1:] + self.true_edges[:-1]) / 2
        self.measured_centers = (self.measured_edges[1:] + self.measured_edges[:-1]) / 2
        bin_widths = self.measured_edges[1:] - self.measured_edges[:-1]
        
        response_matrix = np.zeros((self.n_measured_bins, self.n_true_bins))
        
        print("Building response matrix...")
        for j, true_energy in enumerate(self.true_centers):
            if j % 100 == 0:
                print(f"  Processing true energy bin {j+1}/{self.n_true_bins}")
            
            # For each true energy, calculate probability in each measured bin
            # for i in range(self.n_measured_bins):
            #     # Integrate DRF over the measured energy bin
            #     low = self.measured_edges[i]
            #     high = self.measured_edges[i + 1]
                
            #     # Numerical integration of DRF over the bin
            #     integral, error = quad(drf, low, high, args=(true_energy))
            #     response_matrix[i, j] = integral
            
            energyResponse = drf(true_energy, self.measured_centers)
            
            response_matrix[:, j] = energyResponse * bin_widths
        
        # Normalize each column (true energy) to sum to 1
        response_matrix = response_matrix / np.sum(response_matrix, axis=0, keepdims=True)

        self.R = response_matrix
        return response_matrix, self.true_centers, self.measured_centers
    
    def load_response_matrix(self, true_energy_range, measured_energy_range, filename):
        self.true_edges = np.linspace(true_energy_range[0], true_energy_range[1], self.n_true_bins + 1)
        self.measured_edges = np.linspace(measured_energy_range[0], measured_energy_range[1], self.n_measured_bins + 1)
        
        self.true_centers = (self.true_edges[1:] + self.true_edges[:-1]) / 2
        self.measured_centers = (self.measured_edges[1:] + self.measured_edges[:-1]) / 2
        file = uproot.open(filename)
        hist = file["h_responseMatrix"]
        response_matrix, _, _ = hist.to_numpy()
        
        self.R = response_matrix
        return response_matrix, self.true_centers, self.measured_centers
        
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

def save_response_matrix_th2d(response_matrix, true_edges, measured_edges, 
                             filename="response_matrix.root", histname="response_matrix"):
    """
    Save response matrix as TH2D in ROOT file
    
    Parameters:
    response_matrix: 2D numpy array (n_measured_bins, n_true_bins)
    true_edges: 1D array of true energy bin edges
    measured_edges: 1D array of measured energy bin edges  
    filename: output ROOT filename
    histname: name of the TH2D histogram
    """
    
    # Convert edges to ROOT arrays
    true_edges_root = np.array('d', true_edges)
    measured_edges_root = np.array('d', measured_edges)
    
    n_true_bins = len(true_edges) - 1
    n_measured_bins = len(measured_edges) - 1
    
    # Create TH2D
    h_response = ROOT.TH2D(histname, "Detector Response Matrix",
                          n_true_bins, true_edges_root,
                          n_measured_bins, measured_edges_root)
    
    # Set axis titles
    h_response.GetXaxis().SetTitle("True Energy")
    h_response.GetYaxis().SetTitle("Measured Energy")
    h_response.GetZaxis().SetTitle("Probability")
    
    # Fill the histogram
    for i in range(n_measured_bins):      # y-axis (measured)
        for j in range(n_true_bins):      # x-axis (true)
            # TH2D bin indexing: (x_bin, y_bin) where x=true, y=measured
            h_response.SetBinContent(j+1, i+1, response_matrix[i, j])
    
    # Save to file
    output_file = ROOT.TFile(filename, "RECREATE")
    h_response.Write()
    output_file.Close()
    
    print(f"Response matrix saved to {filename}")
    print(f"Histogram name: {histname}")
    print(f"Dimensions: {n_true_bins} true bins × {n_measured_bins} measured bins")
    
    return h_response

def save_comprehensive_response_data(response_matrix, true_edges, measured_edges, 
                                   true_centers, measured_centers, efficiency,
                                   filename="response_data.root"):
    """
    Save comprehensive response matrix data including 1D projections
    """
    # Convert edges to ROOT arrays
    true_edges_root = np.array('d', true_edges)
    measured_edges_root = np.array('d', measured_edges)
    true_centers_root = np.array('d', true_centers)
    measured_centers_root = np.array('d', measured_centers)
    
    n_true_bins = len(true_edges) - 1
    n_measured_bins = len(measured_edges) - 1
    
    # Create output file
    output_file = ROOT.TFile(filename, "RECREATE")
    
    # 1. Main response matrix (TH2D)
    h_response = ROOT.TH2D("response_matrix", "Detector Response Matrix R(i,j)=P(measured i|true j)",
                          n_true_bins, true_edges_root,
                          n_measured_bins, measured_edges_root)
    h_response.GetXaxis().SetTitle("True Energy [MeV]")
    h_response.GetYaxis().SetTitle("Measured Energy [MeV]")
    h_response.GetZaxis().SetTitle("Probability Density")
    
    # 2. Efficiency curve (TH1D)
    h_efficiency = ROOT.TH1D("efficiency", "Detection Efficiency",
                            n_true_bins, true_edges_root)
    h_efficiency.GetXaxis().SetTitle("True Energy [MeV]")
    h_efficiency.GetYaxis().SetTitle("Efficiency")
    
    # 3. Example response slices (TH1D for different true energies)
    h_response_slices = []
    
    # Fill histograms
    for i in range(n_measured_bins):
        for j in range(n_true_bins):
            h_response.SetBinContent(j+1, i+1, response_matrix[i, j])
    
    for j in range(n_true_bins):
        h_efficiency.SetBinContent(j+1, efficiency[j])
        h_efficiency.SetBinError(j+1, 0.01)  # Small error for plotting
    
    # Create example slices at different true energies
    slice_positions = [0.25, 0.5, 0.75]  # Fractions of energy range
    for frac in slice_positions:
        j = int(n_true_bins * frac)
        if j >= n_true_bins:
            j = n_true_bins - 1
        
        slice_name = f"response_slice_true_{true_centers[j]:.2f}MeV"
        h_slice = ROOT.TH1D(slice_name, f"Response at {true_centers[j]:.2f} MeV",
                           n_measured_bins, measured_edges_root)
        h_slice.GetXaxis().SetTitle("Measured Energy [MeV]")
        h_slice.GetYaxis().SetTitle("Probability Density")
        
        for i in range(n_measured_bins):
            h_slice.SetBinContent(i+1, response_matrix[i, j])
        
        h_slice.Write()
        h_response_slices.append(h_slice)
    
    # Write all histograms
    h_response.Write()
    h_efficiency.Write()
    
    # Create a directory for projections
    output_file.mkdir("projections")
    output_file.cd("projections")
    
    # Save projections
    h_projX = h_response.ProjectionX("projection_true")
    h_projY = h_response.ProjectionY("projection_measured")
    h_projX.Write()
    h_projY.Write()
    
    output_file.cd()  # Back to main directory
    
    # Save metadata as TTree
    tree = ROOT.TTree("metadata", "Response Matrix Metadata")
    
    # Create variables for the tree
    n_true = np.array([n_true_bins], dtype='i')
    n_measured = np.array([n_measured_bins], dtype='i')
    true_min = np.array([true_edges[0]], dtype='f')
    true_max = np.array([true_edges[-1]], dtype='f')
    measured_min = np.array([measured_edges[0]], dtype='f')
    measured_max = np.array([measured_edges[-1]], dtype='f')
    
    # Create branches
    tree.Branch("n_true_bins", n_true, "n_true_bins/I")
    tree.Branch("n_measured_bins", n_measured, "n_measured_bins/I")
    tree.Branch("true_energy_min", true_min, "true_energy_min/F")
    tree.Branch("true_energy_max", true_max, "true_energy_max/F")
    tree.Branch("measured_energy_min", measured_min, "measured_energy_min/F")
    tree.Branch("measured_energy_max", measured_max, "measured_energy_max/F")
    
    tree.Fill()
    tree.Write()
    
    output_file.Close()
    
    print(f"Comprehensive response data saved to {filename}")
    print("Contents:")
    print("  - response_matrix: TH2D main response matrix")
    print("  - efficiency: TH1D detection efficiency")
    print("  - response_slice_*: TH1D example response slices")
    print("  - projections/: directory with 1D projections")
    print("  - metadata: TTree with binning information")

# %%
def save_response_matrix(self, filename="response_matrix.root"):
    """Save the current response matrix to ROOT file"""
    if not hasattr(self, 'R'):
        print("Error: Build response matrix first using build_response_matrix()")
        return
    
    save_comprehensive_response_data(
        self.R, self.true_edges, self.measured_edges,
        self.true_centers, self.measured_centers,
        np.sum(self.R, axis=0),  # efficiency
        filename
    )

# %%
# Example usage with your data
if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    print("Bayesian Unfolding with Functional DRF")
    print("=" * 50)

    
    def gaussian(x, amplitude, mean, sigma):
        return amplitude * np.exp(-0.5 * ((x - mean) / sigma)**2)

    def EnergytoADC(energy):
        a = 0.1
        b = 1
        return a*energy/(1+12.68*0.001+energy/3)+b

    def EnergytoCell(energy):
        adc = EnergytoADC(energy)
        spe = 100
        return adc/spe
    def sipm_response(energy, pe_axis, n_max=1000):
        """
        Calculate SiPM response function for a given energy deposition
        Parameters:
        - energy: True energy value
        - mean_SPE: mean of single photoelectron peak
        - sigma_SPE: Standard deviation of single photoelectron peak
        - mean_ped: mean of pedastal
        - sigma_ped: Standard deviation of pedastal
        - pe_range: Tuple (min, max, n_points) for pe axis
        - n_max: Maximum number of photoelectrons to consider
        Returns:
        - pe_axis: Array of pe values
        - response: Probability density function
        - components: Individual photoelectron peaks for plotting
        """

        Ncells = 8334
        mean_SPE = 110.342
        sigma_SPE= 10.0004
        mean_ped= 85.5295
        sigma_ped= 0.48595
        PDE = 0.2316
        sigma_SPE = (sigma_SPE) / (mean_SPE-mean_ped)#/0.2316
        sigma_ped = (sigma_ped) / (mean_SPE-mean_ped)#/0.2316
        response = np.zeros_like(pe_axis)
        bin_width = builder.true_centers[1]-builder.true_centers[0]

        underflow_axis = np.arange(-(bin_width)*400+builder.true_centers[0], -builder.true_centers[0], builder.true_centers[1]-builder.true_centers[0])
        overflow_axis = np.arange(builder.true_centers[-1]+builder.true_centers[0], (bin_width)*400+builder.true_centers[-1], builder.true_centers[1]-builder.true_centers[0])

        underflow_response = np.zeros_like(underflow_axis)
        overflow_response = np.zeros_like(overflow_axis)
        components = []

        for n in range(0, n_max + 1):   
            cell_prob = 1 - np.exp(-EnergytoCell(energy) / Ncells)
            binomial_prob = binom.pmf(Ncells, n, cell_prob) 
            if binomial_prob > 1e-10:  
                sigma_n = np.sqrt((np.sqrt(n) * sigma_SPE)**2) if n > 0 else np.sqrt(sigma_SPE**2)
                if sigma_n < 1e-10:
                    sigma_n = sigma_SPE
                gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                          np.exp(-0.5 * ((pe_axis*PDE - n) / sigma_n) ** 2)
                component = binomial_prob * gaussian

                if(component[0]>1e-6):
                    underflow_gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                              np.exp(-0.5 * ((underflow_axis*PDE - n) / sigma_n) ** 2)
                    underflow_component = binomial_prob * underflow_gaussian

                    underflow_response += underflow_component
                    component[0] += np.sum(underflow_component)
                
                if(component[-1]>1e-6):
                    overflow_gaussian = (1 / (np.sqrt(2 * np.pi) * sigma_n)) * \
                              np.exp(-0.5 * ((overflow_axis*PDE - n) / sigma_n) ** 2)
                    overflow_component = binomial_prob * overflow_gaussian

                    overflow_response += overflow_component
                    component[-1] += np.sum(overflow_component)
                components.append((n, n, sigma_n, component))
                response += component
                
        return response

    
    def saveGif():
        import imageio.v2 as imageio
        os.makedirs("frames", exist_ok=True)
        filenames = []
        for i, E in enumerate(true_energies):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(R[:,i], color='C0')
            ax.set_title(f"Response for Energy = {E:.1f} MeV")
            ax.set_xlabel("Detector Channel")
            ax.set_ylabel("Signal (a.u.)")
            ax.set_ylim(0, R[:,i].max()*1.1)

            fname = f"gif/frame_{i:03d}.png"
            plt.savefig(fname)
            plt.close(fig)
            filenames.append(fname)

        # Create GIF
        with imageio.get_writer('gif/response_matrix.gif', mode='I', duration=0.1) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
        for filename in filenames:
            os.remove(filename)

    def plotSingleResponse(R, true_energies):
        _, ax = plt.subplots(figsize=(8, 5))

        indices = [30, len(true_energies)//2, len(true_energies)-30]  # first, middle, last
        labels = ['Low edge', 'Mid range', 'High edge']

        for idx, label in zip(indices, labels):
            ax.plot(true_energies, R[:, idx], label=f"{label} ({true_energies[idx]:.2f} ph)")

        ax.set_xlabel("Photon Count")
        ax.set_ylabel("Probability Density")
        ax.set_title("Response Matrix Edge Cases")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.show()
    
    def saveResponseMatrix(R, output_filename="response_matrix.pdf", show=False):
        im = plt.imshow(R, aspect='auto', origin='lower',
                    extent=[true_energies[0], true_energies[-1],
                            measured_energies[0], measured_energies[-1]],
                    norm=LogNorm(vmin=1e-4, vmax=1))

        plt.xlabel('True Photon Count')
        plt.ylabel('Measured Photon Count')
        plt.title('Detector Response Matrix')
        plt.colorbar(im)
        plt.savefig(output_filename)
        if(show):
            plt.show()
        plt.close()
        
    print("Loading measured spectrum...")   

    maxMeasured = 16000

    mean_SPE = 110.342
    sigma_SPE= 10.0004
    mean_ped= 85.5295
    sigma_ped= 0.48595

    saveMatrix = False
    debug = False
    lateral = False

    detector = DetectorCalibration()
    
    for crystal in range(1):
        cal = ChannelCalibration(
            E1=0, E2=5.7038,
            adc1=0, adc2=5419.79,  
            std1=0, std2=1278.59,
            err1=0, err2=10
        )

        detector.add_channel(crystal, cal)

    builder = ResponseMatrixBuilder(detector, n_measured_bins=maxMeasured, n_iterations=40)

    R, true_energies, measured_energies = builder.build_response_matrix(custom_drf=sipm_response)

    maxE = (4001-pedastal_idx-mean_ped)/(mean_SPE-mean_ped)/PDE
    
    useExistingResponseMatrix = True
    if(not useExistingResponseMatrix):
        print("Building response matrix from DRF function...")
        
    else:
        print("Loading response matrix from DRF function...")
        R, true_energies, measured_energies = builder.load_response_matrix(
            true_energy_range=(0, maxE), 
            measured_energy_range=(0, maxE),
            filename="responseMatrix.root")
    
    if (saveMatrix):
        responseShape  = R.shape
        hist1d = ROOT.TH1D("h_measurement", "Measurement", responseShape[0], 0, maxE)
        hist2d = ROOT.TH2D("h_responseMatrix", "Response Matrix", responseShape[0], 0, maxE, responseShape[1], 0, maxE)

        for i in range(responseShape[0]):
            for j in range(responseShape[1]):
                hist2d.SetBinContent(i+1, j+1, R[i, j])  # ROOT bins start at 1
        for i in range(responseShape[0]):
            hist1d.SetBinContent(i+1, measured_spectrum[i])  # ROOT bins start at 1

        root_file = ROOT.TFile("responseMatrix.root", "RECREATE")
        hist2d.Write()
        hist1d.Write()
        root_file.Close()

    saveResponseMatrix(R, output_filename=f"{datapath}/output/responseMatrix.pdf")
    

