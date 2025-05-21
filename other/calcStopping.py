import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Constants
K = 0.307075 # MeV·cm²/mol
me_c2 = 0.510998950  # Electron mass [MeV]

def beta_gamma(E, M):
    """Compute beta and gamma from kinetic energy and mass."""
    E_total = E + M
    gamma = E_total / M
    beta = np.sqrt(1 - 1 / gamma**2)
    return beta, gamma

def Wmax(beta, gamma, M):
    """Maximum energy transfer to an electron."""
    return (2 * me_c2 * beta**2 * gamma**2) / (1 + 2 * gamma * me_c2 / M + (me_c2 / M)**2)

def bethe_interpolate(energy_data, stopping_power_data):
    logE = np.log10(energy_data)  # energy_data from PSTAR
    S = stopping_power_data       # corresponding stopping powers

    return CubicSpline(logE, S)

def bethe_bloch(E, z, M, Z, A, I):
    """
    Calculate stopping power -dE/dx using the Bethe-Bloch equation.
    E: kinetic energy [MeV]
    z: charge of incident particle
    M: mass of incident particle [MeV/c²]
    Z: atomic number of absorber
    A: atomic mass of absorber [g/mol]
    I: mean excitation energy [MeV]
    """
    beta, gamma = beta_gamma(E, M)
    W_max = Wmax(beta, gamma, M)
    coeff = K * z**2 * Z / A
    log_term = np.log(2 * me_c2 * beta**2 * gamma**2 * W_max / I**2)
    dEdx = coeff * (1 / beta**2) * (1/2*log_term - beta**2)
    return dEdx  # MeV/cm

elements = {
    'H': {'Z': 1, 'A': 1.00794, 'I': 19.2e-6},
    'C': {'Z': 6, 'A': 12.011, 'I': 81.0e-6},
    'N': {'Z': 7, 'A': 14.007, 'I': 82.0e-6},
    'O': {'Z': 8, 'A': 15.999, 'I': 95.0e-6},
    'Na': {'Z': 11, 'A': 22.990, 'I': 149.0e-6},
    'Mg': {'Z': 12, 'A': 24.305, 'I': 156.0e-6},
    'P': {'Z': 15, 'A': 30.974, 'I': 173.0e-6},
    'S': {'Z': 16, 'A': 32.06, 'I': 180.0e-6},
    'Cl': {'Z': 17, 'A': 35.45, 'I': 174.0e-6},
    'K': {'Z': 19, 'A': 39.098, 'I': 190.0e-6},
    'Ca': {'Z': 20, 'A': 40.078, 'I': 191.0e-6},
    'Fe': {'Z': 26, 'A': 55.845, 'I': 286.0e-6},
    'Zn': {'Z': 30, 'A': 65.38, 'I': 322.0e-6},
}

stern_params = {
    'H': {'C': 2.2, 'x0': 1.75, 'x1': 4, 'a': 0.12, 'm': 3},
    'O': {'C': 5.4, 'x0': 1.754, 'x1': 4, 'a': 0.2666, 'm': 2.825},
    'H2O': {'C': 3.502, 'x0': 0.2400, 'x1': 2.5, 'a': 0.2065, 'm': 3.007},
}

# Weight fractions
total_mass = 2 * elements['H']['A'] + elements['O']['A']
weights = {
    'H': (2 * elements['H']['A']) / total_mass,
    'O': elements['O']['A'] / total_mass
}

weights_lung = {
    'H': 0.101278,
    'C':  0.102310,
    'N':  0.02865,
    'O':  0.757072,
    'Na': 0.001840,
    'Mg': 0.000730,
    'P':  0.0008,
    'S':  0.002250,
    'Cl': 0.002660,
    'K':  0.001940,
    'Ca': 0.000090,
    'Fe': 0.000370,
    'Zn': 0.000010,
}

A_h2o= 6.005
Z_h2o= 7.42

def stopping_power_water(E):
    dEdx_total = 0
    for el in ['H', 'O']:
        Z = elements[el]['Z']
        A = elements[el]['A']
        I = elements[el]['I']
        w = weights[el]
        dEdx = bethe_bloch(E, z, M, Z, A, I)  # MeV/cm
        dEdx_total += w * dEdx
    coeffh2o = K * z**2 * Z_h2o / A_h2o
    return dEdx_total/density
    return (dEdx_total-coeffh2o*density_correction(*beta_gamma(E,M), 'H2O')) / density  # Convert to MeV·cm²/g

def density_correction(beta, gamma, element):
    x = np.log10(beta * gamma)
    C = stern_params[element]['C']
    x0 = stern_params[element]['x0']
    x1 = stern_params[element]['x1']
    a = stern_params[element]['a']
    m = stern_params[element]['m']
    if x < x0:
        return 0
    elif x < x1:
        return 2 * np.log(10) * x - C + a * (x1 - x) ** m
    else:
        return 2 * np.log(10) * x - C



# Example: Proton in water
z = 1            # Proton charge
M = 938.272      # Proton mass [MeV/c²]
Z = 1         # Effective Z for water (H2O: 2*1 + 8)
A = 1.0078      # g/mol
I = 19.2e-6        # MeV (mean excitation energy for water)
density = 1
#rho = 8.37480E-05
# Energy range
energies = np.linspace(1, 5000, 500)  # MeV
dEdx_mass = [stopping_power_water(E) for E in energies]

ICRUstoppingTotal_H = [9.730E+02, 1.087E+03,
 1.197E+03, 1.300E+03, 1.398E+03, 1.578E+03, 1.741E+03, 1.890E+03, 2.030E+03, 2.161E+03, 2.284E+03, 2.402E+03, 2.621E+03, 2.807E+03, 2.968E+03, 3.107E+03, 3.229E+03, 3.335E+03, 3.427E+03, 3.506E+03, 3.633E+03, 3.723E+03, 3.783E+03, 3.818E+03, 3.833E+03, 3.831E+03, 3.816E+03, 3.789E+03, 3.753E+03, 3.710E+03, 3.661E+03, 3.608E+03, 3.552E+03, 3.493E+03, 3.188E+03, 2.895E+03, 2.633E+03, 2.406E+03, 2.207E+03, 2.034E+03, 1.884E+03, 1.755E+03, 1.546E+03, 1.385E+03, 1.260E+03, 1.160E+03, 1.078E+03, 1.009E+03, 9.481E+02, 8.952E+02, 8.485E+02, 8.068E+02, 7.695E+02, 7.357E+02, 7.051E+02, 6.771E+02, 5.673E+02, 4.902E+02, 4.329E+02, 3.885E+02, 3.530E+02, 3.238E+02, 2.995E+02, 2.788E+02, 2.455E+02, 2.197E+02, 1.992E+02, 1.825E+02, 1.685E+02, 1.567E+02, 1.465E+02, 1.377E+02, 1.299E+02, 1.230E+02, 1.169E+02, 1.114E+02, 1.065E+02, 1.019E+02, 8.444E+01, 7.239E+01, 6.356E+01, 5.679E+01, 4.707E+01, 4.346E+01, 4.041E+01, 3.554E+01, 3.182E+01, 2.887E+01, 2.649E+01, 2.451E+01, 2.285E+01, 2.143E+01, 2.020E+01, 1.913E+01, 1.818E+01, 1.734E+01, 1.659E+01, 1.591E+01, 1.530E+01, 1.295E+01, 1.135E+01, 1.020E+01, 9.328E+00, 8.645E+00, 8.096E+00, 7.645E+00, 7.269E+00, 6.679E+00, 6.238E+00, 5.897E+00, 5.627E+00, 5.409E+00, 5.230E+00, 5.081E+00, 4.955E+00, 4.848E+00, 4.757E+00, 4.678E+00, 4.609E+00, 4.549E+00, 4.497E+00, 4.213E+00, 4.125E+00, 4.107E+00, 4.118E+00, 4.172E+00, 4.240E+00, 4.307E+00, 4.372E+00, 4.431E+00, 4.487E+00, 4.539E+00]

ICRUstoppingTotal_H2o = [1.769E+02 , 1.984E+02 , 2.184E+02 , 2.370E+02 , 2.544E+02 , 2.864E+02 , 3.153E+02 , 3.420E+02 , 3.667E+02 , 3.900E+02 , 4.120E+02 , 4.329E+02 , 4.745E+02 , 5.110E+02 , 5.437E+02 , 5.733E+02 , 6.001E+02 , 6.245E+02 , 6.467E+02 , 6.671E+02 , 7.028E+02 , 7.324E+02 , 7.569E+02 , 7.768E+02 , 7.927E+02 , 8.050E+02 , 8.142E+02 , 8.205E+02 , 8.243E+02 , 8.260E+02 , 8.258E+02 , 8.239E+02 , 8.206E+02 , 8.161E+02 , 7.814E+02 , 7.371E+02 , 6.969E+02 , 6.613E+02 , 6.294E+02 , 6.006E+02 , 5.744E+02 , 5.504E+02 , 5.080E+02 , 4.719E+02 , 4.406E+02 , 4.132E+02 , 3.891E+02 , 3.680E+02 , 3.492E+02 , 3.325E+02 , 3.175E+02 , 3.039E+02 , 2.917E+02 , 2.805E+02 , 2.702E+02 , 2.608E+02 , 2.229E+02 , 1.957E+02 , 1.749E+02 , 1.586E+02 , 1.454E+02 , 1.344E+02 , 1.251E+02 , 1.172E+02 , 1.042E+02 , 9.404E+01 , 8.586E+01 , 7.911E+01 , 7.343E+01 , 6.858E+01 , 6.438E+01 , 6.071E+01 , 5.747E+01 , 5.460E+01 , 5.202E+01 , 4.969E+01 , 4.759E+01 , 4.567E+01 , 3.815E+01 , 3.292E+01 , 2.905E+01 , 2.607E+01 , 2.175E+01 , 2.013E+01 , 1.876E+01 , 1.656E+01 , 1.488E+01 , 1.354E+01 , 1.245E+01 , 1.154E+01 , 1.078E+01 , 1.013E+01 , 9.559E+00 , 9.063E+00 , 8.625E+00 , 8.236E+00 , 7.888E+00 , 7.573E+00 , 7.289E+00 , 6.192E+00 , 5.445E+00 , 4.903E+00 , 4.492E+00 , 4.170E+00 , 3.911E+00 , 3.698E+00 , 3.520E+00 , 3.241E+00 , 3.032E+00 , 2.871E+00 , 2.743E+00 , 2.640E+00 , 2.556E+00 , 2.485E+00 , 2.426E+00 , 2.376E+00 , 2.333E+00 , 2.296E+00 , 2.264E+00 , 2.236E+00 , 2.211E+00 , 2.070E+00 , 2.021E+00 , 2.004E+00 , 2.001E+00 , 2.012E+00 , 2.031E+00 , 2.052E+00 , 2.072E+00 , 2.091E+00 , 2.109E+00 , 2.126E+00]

ICRU_E = [1.000E-03, 1.500E-03, 2.000E-03, 2.500E-03, 3.000E-03, 4.000E-03, 5.000E-03, 6.000E-03, 7.000E-03, 8.000E-03, 9.000E-03, 1.000E-02, 1.250E-02, 1.500E-02, 1.750E-02, 2.000E-02, 2.250E-02, 2.500E-02, 2.750E-02, 3.000E-02, 3.500E-02, 4.000E-02, 4.500E-02, 5.000E-02, 5.500E-02, 6.000E-02, 6.500E-02, 7.000E-02, 7.500E-02, 8.000E-02, 8.500E-02, 9.000E-02, 9.500E-02, 1.000E-01, 1.250E-01, 1.500E-01, 1.750E-01, 2.000E-01, 2.250E-01, 2.500E-01, 2.750E-01, 3.000E-01, 3.500E-01, 4.000E-01, 4.500E-01, 5.000E-01, 5.500E-01, 6.000E-01, 6.500E-01, 7.000E-01, 7.500E-01, 8.000E-01, 8.500E-01, 9.000E-01, 9.500E-01, 1.000E+00, 1.250E+00, 1.500E+00, 1.750E+00, 2.000E+00, 2.250E+00, 2.500E+00, 2.750E+00, 3.000E+00, 3.500E+00, 4.000E+00, 4.500E+00, 5.000E+00, 5.500E+00, 6.000E+00, 6.500E+00, 7.000E+00, 7.500E+00, 8.000E+00, 8.500E+00, 9.000E+00, 9.500E+00, 1.000E+01, 1.250E+01, 1.500E+01, 1.750E+01, 2.000E+01, 2.500E+01, 2.750E+01, 3.000E+01, 3.500E+01, 4.000E+01, 4.500E+01, 5.000E+01, 5.500E+01, 6.000E+01, 6.500E+01, 7.000E+01, 7.500E+01, 8.000E+01, 8.500E+01, 9.000E+01, 9.500E+01, 1.000E+02, 1.250E+02, 1.500E+02, 1.750E+02, 2.000E+02, 2.250E+02, 2.500E+02, 2.750E+02, 3.000E+02, 3.500E+02, 4.000E+02, 4.500E+02, 5.000E+02, 5.500E+02, 6.000E+02, 6.500E+02, 7.000E+02, 7.500E+02, 8.000E+02, 8.500E+02, 9.000E+02, 9.500E+02, 1.000E+03, 1.500E+03, 2.000E+03, 2.500E+03, 3.000E+03, 4.000E+03, 5.000E+03, 6.000E+03, 7.000E+03, 8.000E+03, 9.000E+03, 1.000E+04]
# Plot
plt.figure(figsize=(8, 5))
plt.plot(energies, dEdx_mass, label='Bethe-Bloch Stopping Power')
spline = bethe_interpolate(ICRU_E, ICRUstoppingTotal_H2o)
plt.plot(energies, spline(np.log10(energies)), label='ICRU Stopping Power (H2O)')
plt.plot(ICRU_E, ICRUstoppingTotal_H2o, 'o', label='ICRU Stopping Power')
#plt.plot(ICRU_E, [rho*Stopping for Stopping in ICRUstoppingTotal], 'o', label='ICRU Stopping Power')
plt.xlabel('Kinetic Energy [MeV]')
plt.ylabel('Stopping Power [-dE/dx] [MeV/cm]')
plt.title('Stopping Power of Protons in Water (Bethe-Bloch)')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
