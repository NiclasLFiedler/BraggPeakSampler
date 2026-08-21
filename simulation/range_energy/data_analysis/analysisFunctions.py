from dataclasses import dataclass
import numpy as np

@dataclass
class EnergyRangeData:
    name: str
    colors : dict

    energies: np.ndarray
    ranges: np.ndarray
    range_errors: np.ndarray
    sigmas: np.ndarray
    sigma_errors: np.ndarray

    alpha: np.ndarray = None
    alpha_error: np.ndarray = None
    p: float = None
    p_error: float = None

    useSumFit : bool = None

    residuals: np.ndarray = None
    normalized_residuals: np.ndarray = None

    chi2: float = None
    reduced_chi2: float = None
    rms_residual: float = None

def range_energy_relationship(energy, alpha, p):
    return alpha * energy**p

def range_energy_sum(energy, alpha1, alpha2, alpha3, alpha4):
    term1 = alpha1 * energy
    term2 = alpha2 * energy**2
    term3 = alpha3 * energy**3 
    term4 = alpha4 * energy**4

    return term1 + term2 + term3 + term4

def range_energy(data : EnergyRangeData, energy):
    if data.useSumFit:
        return range_energy_sum(energy, *data.alpha)
    else:
        return range_energy_relationship(energy, data.alpha[0], data.p)

def stoppingPower(data : EnergyRangeData, energy):
    return 1/(data.p*data.alpha[0])*energy**(1-data.p)

# Define a Gaussian function
def gaussian(x, mean, sigma, amplitude):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

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


def load_EnergyRange(filename):
    data = np.load(filename, allow_pickle=True)

    return EnergyRangeData(
        name=str(data["name"]),
        colors={},  # colors are currently not saved

        energies=data["energies"],
        ranges=data["ranges"],
        range_errors=data["range_errors"],
        sigmas=data["sigmas"],
        sigma_errors=data["sigma_errors"],

        alpha=data["alpha"],
        alpha_error=data["alpha_error"],
        p=data["p"],
        p_error=data["p_error"],

        useSumFit=bool(data["useSumFit"]),

        residuals=data["residuals"],
        normalized_residuals=data["normalized_residuals"],

        chi2=float(data["chi2"]),
        reduced_chi2=float(data["reduced_chi2"]),
        rms_residual=float(data["rms_residual"])
    )