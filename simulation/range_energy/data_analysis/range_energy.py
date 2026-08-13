import ROOT
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uproot  # uproot is a great library for reading ROOT files in Python
matplotlib.use('TkAgg')  # or 'Qt5Agg'
from copy import deepcopy
from analysisFunctions import EnergyRangeData, range_energy_relationship, range_energy_sum, load_range_data, gaussian

def load_range_data(file_folder, name=None, colors=None, UseSumFit=False):
    data = np.load(f"{file_folder}/ranges.npz")

    return EnergyRangeData(
        name=name,
        colors=colors,
        energies=data["energies"],
        ranges=data["ranges"]/10,
        range_errors=data["range_errors"]/10,
        sigmas=data["sigmas"],
        sigma_errors=data["sigma_errors"],
        useSumFit=UseSumFit
    )

def fit_range_energy(data: EnergyRangeData, output=False):

    if not data.useSumFit:
        popt, pcov = curve_fit(
            range_energy_relationship,
            data.energies,
            data.ranges,
            p0=[0.002, 1.75],
            maxfev=10000
        )

        data.alpha, data.p = [popt[0]], popt[1]
        data.alpha_error, data.p_error = [np.sqrt(np.diag(pcov))[0]], np.sqrt(np.diag(pcov))[1]

        if output:
            print(f"Fitted parameters for {data.name}:")
            print(f"alpha = {data.alpha[0]:.6g} ± {data.alpha_error[0]:.6g}")
            print(f"p     = {data.p:.6f} ± {data.p_error:.6f}")
    else:
        popt, pcov = curve_fit(
            range_energy_sum,
            data.energies,
            data.ranges,
            p0=[6.94656e-3, 8.13116e-4, -1.21068e-6, 1.053e-9],
            maxfev=10000
        )
    
        data.alpha = popt
        data.alpha_error = np.sqrt(np.diag(pcov))
    
        if output:
            print(f"Fitted parameters for {data.name}:")
            print(f"alpha1 = {data.alpha[0]:.6g} ± {data.alpha_error[0]:.6g}")
            print(f"alpha2 = {data.alpha[1]:.6g} ± {data.alpha_error[1]:.6g}")
            print(f"alpha3 = {data.alpha[2]:.6g} ± {data.alpha_error[2]:.6g}")
            print(f"alpha4 = {data.alpha[3]:.6g} ± {data.alpha_error[3]:.6g}")

    calculate_residuals(data)

    if output:
        print(f"chi2          = {data.chi2:.4f}")
        print(f"reduced chi2  = {data.reduced_chi2:.4f}")
        print(f"RMS residual  = {data.rms_residual:.6g} cm")
        
    return data

def calculate_residuals(data: EnergyRangeData):
    if data.useSumFit:
        fitted_ranges = range_energy_sum(
            data.energies,
            data.alpha[0],
            data.alpha[1],
            data.alpha[2],
            data.alpha[3]
        )
        n_parameters = 4
    else:
        fitted_ranges = range_energy_relationship(
            data.energies,
            data.alpha[0],
            data.p
        )
        n_parameters = 2
    
    data.residuals = data.ranges - fitted_ranges
    
    # Normalized residuals: divide by measurement uncertainty
    # Assumes range_errors are 1-sigma uncertainties on each measurement
    data.normalized_residuals = data.residuals / data.range_errors
    
    # Chi-squared: sum of squared normalized residuals
    data.chi2 = np.sum(data.normalized_residuals**2)
    
    dof = len(data.energies) - n_parameters
    data.reduced_chi2 = data.chi2 / dof
    
    # RMS residual in physical units
    data.rms_residual = np.sqrt(np.mean(data.residuals**2))
    
    return data

def plotInit():
    plt.rcParams.update({'font.size': 26})
    plt.figure(figsize=(20, 15))
    plt.xlabel('Energy / MeV')
    plt.ylabel('Range / cm')
    plt.grid(True)
    return

def plotEnd(path = None):
    plt.legend()
    if path is not None:
        plt.savefig(f"{path}/rangeenergy.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    return

def plot_range_energy(data: EnergyRangeData):

    fitEnergies = np.linspace(
        data.energies.min(),
        data.energies.max() + 20,
        200
    )

    if data.useSumFit:
        fitRanges = range_energy_sum(
            fitEnergies,
            data.alpha[0],
            data.alpha[1],
            data.alpha[2],
            data.alpha[3]
        )
    else:
        fitRanges = range_energy_relationship(
            fitEnergies,
            data.alpha[0],
            data.p
        )

    plt.errorbar(
        data.energies,
        data.ranges,
        yerr=data.range_errors,
        fmt='o',
        color=data.colors[0],
        capsize=5
    )

    if data.useSumFit:
        plt.plot(
            fitEnergies,
            fitRanges,
            color=data.colors[1],
            linewidth=2,
            label=(
                f'{data.name} Fit: '
                f'$\\alpha_1 = {data.alpha[0]:.3e}$ '
                f'$\\alpha_2 = {data.alpha[1]:.3e}$ '
                f'$\\alpha_3 = {data.alpha[2]:.3e}$ '
                f'$\\alpha_4 = {data.alpha[3]:.3e}$ '
            )
        )
    else:
        plt.plot(
            fitEnergies,
            fitRanges,
            color=data.colors[1],
            linewidth=2,
            label=(
                f'{data.name} Fit: '
                f'$\\alpha_{{{data.name}}} = {data.alpha[0]:.3e}$ '
                f'$\\frac{{cm}}{{MeV^p}}$; '
                f'$p_{{{data.name}}}$ = {data.p:.3e}'
            )
        )
    return

def plot_residuals(data: EnergyRangeData):
    fit_type = "Polynomial" if data.useSumFit else "Power law"
    plt.errorbar(
        data.energies,
        data.residuals,
        yerr=data.range_errors,
        fmt='o--',
        color=data.colors[0],
        capsize=5
    )

    plt.axhline(
        0,
        color='black',
        linestyle='--',
        linewidth=2,
        label=(
            f'{data.name} ({fit_type}), '
            fr'$\chi^2_\nu={data.reduced_chi2:.2f}$, '
            fr'RMS={data.rms_residual:.2e}\,\mathrm{{cm}}'
        )
    )
    plt.ylabel(r'Residual $R_\mathrm{data}-R_\mathrm{fit}$ / cm')

def plot_normalized_residuals(data: EnergyRangeData):

    fit_type = "Polynomial" if data.useSumFit else "Power law"

    plt.plot(
        data.energies,
        data.normalized_residuals,
        'o--',
        color=data.colors[0],
        label=(
            f'{data.name} ({fit_type}), '
            fr'$\chi^2_\nu={data.reduced_chi2:.2f}$'
        )
    )

    plt.axhline(
        0,
        color='black',
        linestyle='--',
        linewidth=2
    )

    plt.axhline(
        1,
        color='gray',
        linestyle=':',
        linewidth=1.5
    )

    plt.axhline(
        -1,
        color='gray',
        linestyle=':',
        linewidth=1.5
    )

    plt.ylabel(
        r'Normalized residual '
        r'$(R_\mathrm{data}-R_\mathrm{fit})/\sigma_R$'
    )

def save_range_data(data: EnergyRangeData, filename):

    np.savez(
        filename,
        name=data.name,
        energies=data.energies,
        ranges=data.ranges,
        range_errors=data.range_errors,
        sigmas=data.sigmas,
        sigma_errors=data.sigma_errors,

        alpha=data.alpha,
        alpha_error=data.alpha_error,
        p=data.p,
        p_error=data.p_error,

        useSumFit=data.useSumFit,

        residuals=data.residuals,
        normalized_residuals=data.normalized_residuals,

        chi2=data.chi2,
        reduced_chi2=data.reduced_chi2,
        rms_residual=data.rms_residual
    )

    print(f"Saved {data.name} to {filename}")

def main():
    UseSumFit = True

    targetColorMap = ["#000000","#1f77b4", "#4e79a7", "#76b7b2", "#bab0ac", "#f28e2b", "#e15759", "#9c755f"]

    colors = {
        "PbWO4":    (targetColorMap[0], targetColorMap[0]),
        "PTFE":  (targetColorMap[1], targetColorMap[1]),
        "H2O":   (targetColorMap[2], targetColorMap[2]),
        "ICRU_ALU": (targetColorMap[3], targetColorMap[3]),
        "ICRU_AIR": (targetColorMap[4], targetColorMap[4]),
        "ICRU_H2O": (targetColorMap[5], targetColorMap[5]),
    }

    rho_alu = 2.69890
    rho_h2o_icru = 1.000
    rho_h2o_air=1.2048e-3
    
    ICRUenergies = [3, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250]
    CSDA_ranges_h20 = [1.499e-2, 3.623e-2, 1.23e-1, 2.539e-1, 4.26e-1, 8.853e-1, 1.489, 2.227, 3.093, 4.080, 5.184, 6.398, 7.718, 1.146e1, 1.577e1, 2.062e1, 2.596e1, 3.174e1, 3.794e1] #g/cm^2
    CSDA_ranges_air= [1.737e-2, 4.173e-2, 1.408e-1, 2.899e-1, 4.855e-1, 1.007, 1.691, 2.528, 3.509, 4.628, 5.876, 7.250, 8.744, 1.297e1, 1.786e1, 2.334e1, 2.937e1, 3.590e1, 4.290e1]
    CSDA_ranges_alu = [2.193e-02, 5.157e-02, 1.697e-01, 3.448e-01, 5.726e-01, 1.175e+00, 1.961e+00, 2.918e+00, 4.037e+00, 5.309e+00, 6.727e+00, 8.284e+00, 9.976e+00, 1.476e+01, 2.026e+01, 2.644e+01, 3.322e+01, 4.057e+01, 4.843e+01]
    
    ICRU_H2O_data = EnergyRangeData(
        name="ICRU_H2O",
        colors = colors["ICRU_H2O"],
        energies=np.array(ICRUenergies),
        ranges=np.array(CSDA_ranges_h20) / rho_h2o_icru,
        range_errors=np.array(CSDA_ranges_h20) / rho_h2o_icru * 0.015,
        sigmas=np.full(len(ICRUenergies), np.nan),
        sigma_errors=np.full(len(ICRUenergies), np.nan),
        useSumFit=UseSumFit
    )

    ICRU_AIR_data = EnergyRangeData(
        name="ICRU_AIR",
        colors = colors["ICRU_AIR"],
        energies=np.array(ICRUenergies),
        ranges=np.array(CSDA_ranges_air) / rho_h2o_air,
        range_errors=np.array(CSDA_ranges_air) / rho_h2o_air * 0.015,
        sigmas=np.full(len(ICRUenergies), np.nan),
        sigma_errors=np.full(len(ICRUenergies), np.nan),
        useSumFit=UseSumFit
    )

    ICRU_ALU_data = EnergyRangeData(
        name="ICRU_ALU",
        colors = colors["ICRU_ALU"],
        energies=np.array(ICRUenergies),
        ranges=np.array(CSDA_ranges_alu) / rho_alu,
        range_errors=np.array(CSDA_ranges_alu) / rho_alu * 0.015,
        sigmas=np.full(len(ICRUenergies), np.nan),
        sigma_errors=np.full(len(ICRUenergies), np.nan),
        useSumFit=UseSumFit
    )

    h2o_data = load_range_data("h2o", name="H2O", colors=colors["H2O"], UseSumFit=UseSumFit)
    pbwo4_data = load_range_data("pbwo4", name="PbWO4", colors=colors["PbWO4"], UseSumFit=UseSumFit)

    h2o_data = fit_range_energy(h2o_data)
    pbwo4_data = fit_range_energy(pbwo4_data)

    ICRU_H2O_data = fit_range_energy(ICRU_H2O_data, output=True)
    # ICRU_AIR_data = fit_range_energy(ICRU_AIR_data)
    ICRU_ALU_data = fit_range_energy(ICRU_ALU_data)

    # ICRU_H2O_data.alpha = [6.94656e-3, 8.13116e-4, -1.21068e-6, 1.053e-9] #paper fit paramas
    plotInit()
    plot_range_energy(pbwo4_data)
    plot_range_energy(h2o_data)
    plot_range_energy(ICRU_H2O_data)

    pbwo4_data_alt = deepcopy(pbwo4_data)
    pbwo4_data_alt.useSumFit = not pbwo4_data_alt.useSumFit
    pbwo4_data_alt.colors = colors["ICRU_AIR"]
    pbwo4_data_alt = fit_range_energy(pbwo4_data_alt, output=True)

    h2o_data_alt = deepcopy(h2o_data)
    h2o_data_alt.useSumFit = not h2o_data_alt.useSumFit
    h2o_data_alt.colors = colors["ICRU_AIR"]
    h2o_data_alt = fit_range_energy(h2o_data_alt, output=True)
    plot_range_energy(pbwo4_data_alt)
    plot_range_energy(h2o_data_alt)
    plotEnd("h2o")

    plotInit()
    plot_residuals(pbwo4_data)
    plot_residuals(pbwo4_data_alt)
    plotEnd()

    save_range_data(h2o_data, "h2o_range_energy.npz")
    save_range_data(h2o_data_alt, "h2o_alt_range_energy.npz")
    save_range_data(pbwo4_data, "pbwo4_range_energy.npz")
    save_range_data(ICRU_H2O_data, "ICRU_H2O_range_energy.npz")
    save_range_data(ICRU_ALU_data, "ICRU_ALU_range_energy.npz")

if __name__ == "__main__":
    main()