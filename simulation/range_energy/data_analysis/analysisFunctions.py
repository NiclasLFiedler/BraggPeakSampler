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