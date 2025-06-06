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
    'N': {'Z': 7, 'A': 14.007, 'I': 82.0e-6},
    'H': {'Z': 1, 'A': 1.00794, 'I': 19.2e-6},
    'C': {'Z': 6, 'A': 12.011, 'I': 81.0e-6},
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
    return (dEdx_total-coeffh2o*density_correction(*beta_gamma(E,M), 'H2O'))

def stopping_power_lung(E):
    dEdx_total = 0
    for el in elements:
        Z = elements[el]['Z']
        A = elements[el]['A']
        I = elements[el]['I']
        w = weights_lung[el]
        dEdx = bethe_bloch(E, z, M, Z, A, I)  # MeV/cm
        dEdx_total += w * dEdx
    coeffh2o = K * z**2 * Z_h2o / A_h2o
    #if density_correction(*beta_gamma(E,M), 'H2O') != 0:
    #    print(f'Density correction is not zero at E={E} MeV')
    return (dEdx_total-coeffh2o*density_correction(*beta_gamma(E,M), 'H2O'))

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
density_lung = 0.26
density_air = 0.00120479
density_pmma = 1.18
#rho = 8.37480E-05
# Energy range
energies = np.linspace(1, 1000, 50000)  # MeV
dEdx_mass_lung = [stopping_power_lung(E) for E in energies]

ICRU_E = [1.000E-03, 1.500E-03, 2.000E-03, 2.500E-03, 3.000E-03, 4.000E-03, 5.000E-03, 6.000E-03, 7.000E-03, 8.000E-03, 9.000E-03, 1.000E-02, 1.250E-02, 1.500E-02, 1.750E-02, 2.000E-02, 2.250E-02, 2.500E-02, 2.750E-02, 3.000E-02, 3.500E-02, 4.000E-02, 4.500E-02, 5.000E-02, 5.500E-02, 6.000E-02, 6.500E-02, 7.000E-02, 7.500E-02, 8.000E-02, 8.500E-02, 9.000E-02, 9.500E-02, 1.000E-01, 1.250E-01, 1.500E-01, 1.750E-01, 2.000E-01, 2.250E-01, 2.500E-01, 2.750E-01, 3.000E-01, 3.500E-01, 4.000E-01, 4.500E-01, 5.000E-01, 5.500E-01, 6.000E-01, 6.500E-01, 7.000E-01, 7.500E-01, 8.000E-01, 8.500E-01, 9.000E-01, 9.500E-01, 1.000E+00, 1.250E+00, 1.500E+00, 1.750E+00, 2.000E+00, 2.250E+00, 2.500E+00, 2.750E+00, 3.000E+00, 3.500E+00, 4.000E+00, 4.500E+00, 5.000E+00, 5.500E+00, 6.000E+00, 6.500E+00, 7.000E+00, 7.500E+00, 8.000E+00, 8.500E+00, 9.000E+00, 9.500E+00, 1.000E+01, 1.250E+01, 1.500E+01, 1.750E+01, 2.000E+01, 2.500E+01, 2.750E+01, 3.000E+01, 3.500E+01, 4.000E+01, 4.500E+01, 5.000E+01, 5.500E+01, 6.000E+01, 6.500E+01, 7.000E+01, 7.500E+01, 8.000E+01, 8.500E+01, 9.000E+01, 9.500E+01, 1.000E+02, 1.250E+02, 1.500E+02, 1.750E+02, 2.000E+02, 2.250E+02, 2.500E+02, 2.750E+02, 3.000E+02, 3.500E+02, 4.000E+02, 4.500E+02, 5.000E+02, 5.500E+02, 6.000E+02, 6.500E+02, 7.000E+02, 7.500E+02, 8.000E+02, 8.500E+02, 9.000E+02, 9.500E+02, 1.000E+03, 1.500E+03, 2.000E+03, 2.500E+03, 3.000E+03, 4.000E+03, 5.000E+03, 6.000E+03, 7.000E+03, 8.000E+03, 9.000E+03, 1.000E+04]


ICRUstoppingTotal_H = [9.730E+02, 1.087E+03,1.197E+03, 1.300E+03, 1.398E+03, 1.578E+03, 1.741E+03, 1.890E+03, 2.030E+03, 2.161E+03, 2.284E+03, 2.402E+03, 2.621E+03, 2.807E+03, 2.968E+03, 3.107E+03, 3.229E+03, 3.335E+03, 3.427E+03, 3.506E+03, 3.633E+03, 3.723E+03, 3.783E+03, 3.818E+03, 3.833E+03, 3.831E+03, 3.816E+03, 3.789E+03, 3.753E+03, 3.710E+03, 3.661E+03, 3.608E+03, 3.552E+03, 3.493E+03, 3.188E+03, 2.895E+03, 2.633E+03, 2.406E+03, 2.207E+03, 2.034E+03, 1.884E+03, 1.755E+03, 1.546E+03, 1.385E+03, 1.260E+03, 1.160E+03, 1.078E+03, 1.009E+03, 9.481E+02, 8.952E+02, 8.485E+02, 8.068E+02, 7.695E+02, 7.357E+02, 7.051E+02, 6.771E+02, 5.673E+02, 4.902E+02, 4.329E+02, 3.885E+02, 3.530E+02, 3.238E+02, 2.995E+02, 2.788E+02, 2.455E+02, 2.197E+02, 1.992E+02, 1.825E+02, 1.685E+02, 1.567E+02, 1.465E+02, 1.377E+02, 1.299E+02, 1.230E+02, 1.169E+02, 1.114E+02, 1.065E+02, 1.019E+02, 8.444E+01, 7.239E+01, 6.356E+01, 5.679E+01, 4.707E+01, 4.346E+01, 4.041E+01, 3.554E+01, 3.182E+01, 2.887E+01, 2.649E+01, 2.451E+01, 2.285E+01, 2.143E+01, 2.020E+01, 1.913E+01, 1.818E+01, 1.734E+01, 1.659E+01, 1.591E+01, 1.530E+01, 1.295E+01, 1.135E+01, 1.020E+01, 9.328E+00, 8.645E+00, 8.096E+00, 7.645E+00, 7.269E+00, 6.679E+00, 6.238E+00, 5.897E+00, 5.627E+00, 5.409E+00, 5.230E+00, 5.081E+00, 4.955E+00, 4.848E+00, 4.757E+00, 4.678E+00, 4.609E+00, 4.549E+00, 4.497E+00, 4.213E+00, 4.125E+00, 4.107E+00, 4.118E+00, 4.172E+00, 4.240E+00, 4.307E+00, 4.372E+00, 4.431E+00, 4.487E+00, 4.539E+00]

ICRUstoppingTotal_H2o = [1.769E+02 , 1.984E+02 , 2.184E+02 , 2.370E+02 , 2.544E+02 , 2.864E+02 , 3.153E+02 , 3.420E+02 , 3.667E+02 , 3.900E+02 , 4.120E+02 , 4.329E+02 , 4.745E+02 , 5.110E+02 , 5.437E+02 , 5.733E+02 , 6.001E+02 , 6.245E+02 , 6.467E+02 , 6.671E+02 , 7.028E+02 , 7.324E+02 , 7.569E+02 , 7.768E+02 , 7.927E+02 , 8.050E+02 , 8.142E+02 , 8.205E+02 , 8.243E+02 , 8.260E+02 , 8.258E+02 , 8.239E+02 , 8.206E+02 , 8.161E+02 , 7.814E+02 , 7.371E+02 , 6.969E+02 , 6.613E+02 , 6.294E+02 , 6.006E+02 , 5.744E+02 , 5.504E+02 , 5.080E+02 , 4.719E+02 , 4.406E+02 , 4.132E+02 , 3.891E+02 , 3.680E+02 , 3.492E+02 , 3.325E+02 , 3.175E+02 , 3.039E+02 , 2.917E+02 , 2.805E+02 , 2.702E+02 , 2.608E+02 , 2.229E+02 , 1.957E+02 , 1.749E+02 , 1.586E+02 , 1.454E+02 , 1.344E+02 , 1.251E+02 , 1.172E+02 , 1.042E+02 , 9.404E+01 , 8.586E+01 , 7.911E+01 , 7.343E+01 , 6.858E+01 , 6.438E+01 , 6.071E+01 , 5.747E+01 , 5.460E+01 , 5.202E+01 , 4.969E+01 , 4.759E+01 , 4.567E+01 , 3.815E+01 , 3.292E+01 , 2.905E+01 , 2.607E+01 , 2.175E+01 , 2.013E+01 , 1.876E+01 , 1.656E+01 , 1.488E+01 , 1.354E+01 , 1.245E+01 , 1.154E+01 , 1.078E+01 , 1.013E+01 , 9.559E+00 , 9.063E+00 , 8.625E+00 , 8.236E+00 , 7.888E+00 , 7.573E+00 , 7.289E+00 , 6.192E+00 , 5.445E+00 , 4.903E+00 , 4.492E+00 , 4.170E+00 , 3.911E+00 , 3.698E+00 , 3.520E+00 , 3.241E+00 , 3.032E+00 , 2.871E+00 , 2.743E+00 , 2.640E+00 , 2.556E+00 , 2.485E+00 , 2.426E+00 , 2.376E+00 , 2.333E+00 , 2.296E+00 , 2.264E+00 , 2.236E+00 , 2.211E+00 , 2.070E+00 , 2.021E+00 , 2.004E+00 , 2.001E+00 , 2.012E+00 , 2.031E+00 , 2.052E+00 , 2.072E+00 , 2.091E+00 , 2.109E+00 , 2.126E+00]

ICRUstoppingTotal_PMMA = [2.147E+02 ,2.463E+02 ,2.746E+02 ,3.004E+02 ,3.242E+02 ,3.635E+02 ,3.994E+02 ,4.322E+02 ,4.624E+02 ,4.899E+02 ,5.153E+02 ,5.391E+02 ,5.870E+02 ,6.277E+02 ,6.632E+02 ,6.948E+02 ,7.229E+02 ,7.479E+02 ,7.700E+02 ,7.898E+02 ,8.244E+02 ,8.537E+02 ,8.779E+02 ,8.972E+02 ,9.118E+02 ,9.224E+02 ,9.296E+02 ,9.340E+02 ,9.358E+02 ,9.355E+02 ,9.334E+02 ,9.296E+02 ,9.245E+02 ,9.183E+02 ,8.757E+02 ,8.243E+02 ,7.723E+02 ,7.232E+02 ,6.769E+02 ,6.346E+02 ,5.969E+02 ,5.634E+02 ,5.077E+02 ,4.639E+02 ,4.289E+02 ,4.006E+02 ,3.770E+02 ,3.564E+02 ,3.383E+02 ,3.222E+02 ,3.078E+02 ,2.948E+02 ,2.830E+02 ,2.722E+02 ,2.623E+02 ,2.532E+02 ,2.168E+02 ,1.905E+02 ,1.705E+02 ,1.546E+02 ,1.418E+02 ,1.311E+02 ,1.221E+02 ,1.143E+02 ,1.017E+02 ,9.179E+01 ,8.379E+01 ,7.719E+01 ,7.164E+01 ,6.690E+01 ,6.280E+01 ,5.921E+01 ,5.605E+01 ,5.324E+01 ,5.073E+01 ,4.845E+01 ,4.640E+01 ,4.452E+01 ,3.719E+01 ,3.208E+01 ,2.831E+01 ,2.539E+01 ,2.118E+01 ,1.961E+01 ,1.827E+01 ,1.613E+01 ,1.449E+01 ,1.318E+01 ,1.212E+01 ,1.124E+01 ,1.050E+01 ,9.858E+00 ,9.306E+00 ,8.823E+00 ,8.397E+00 ,8.018E+00 ,7.678E+00 ,7.372E+00 ,7.095E+00 ,6.027E+00 ,5.300E+00 ,4.772E+00 ,4.372E+00 ,4.058E+00 ,3.806E+00 ,3.599E+00 ,3.426E+00 ,3.154E+00 ,2.951E+00 ,2.794E+00 ,2.670E+00 ,2.569E+00 ,2.487E+00 ,2.418E+00 ,2.361E+00 ,2.312E+00 ,2.269E+00 ,2.232E+00 ,2.199E+00 ,2.170E+00 ,2.145E+00 ,2.004E+00 ,1.954E+00 ,1.937E+00 ,1.934E+00 ,1.945E+00 ,1.964E+00 ,1.985E+00 ,2.005E+00 ,2.024E+00 ,2.042E+00 ,2.059E+00]

ICRUstoppingTotal_Air = [1.414E+02, 1.651E+02, 1.855E+02, 2.038E+02, 2.206E+02, 2.507E+02, 2.776E+02, 3.021E+02, 3.248E+02, 3.460E+02, 3.660E+02, 3.850E+02, 4.224E+02, 4.552E+02, 4.843E+02, 5.106E+02, 5.343E+02, 5.558E+02, 5.755E+02, 5.934E+02, 6.246E+02, 6.506E+02, 6.721E+02, 6.897E+02, 7.038E+02, 7.149E+02, 7.233E+02, 7.293E+02, 7.333E+02, 7.355E+02, 7.360E+02, 7.352E+02, 7.332E+02, 7.301E+02, 7.038E+02, 6.680E+02, 6.298E+02, 5.928E+02, 5.589E+02, 5.284E+02, 5.011E+02, 4.767E+02, 4.353E+02, 4.015E+02, 3.736E+02, 3.501E+02, 3.300E+02, 3.123E+02, 2.967E+02, 2.826E+02, 2.701E+02, 2.589E+02, 2.486E+02, 2.393E+02, 2.308E+02, 2.229E+02, 1.912E+02, 1.683E+02, 1.509E+02, 1.371E+02, 1.258E+02, 1.165E+02, 1.086E+02, 1.018E+02, 9.068E+01, 8.197E+01, 7.492E+01, 6.909E+01, 6.417E+01, 5.997E+01, 5.633E+01, 5.315E+01, 5.033E+01, 4.783E+01, 4.559E+01, 4.357E+01, 4.173E+01, 4.006E+01, 3.351E+01, 2.894E+01, 2.555E+01, 2.294E+01, 1.915E+01, 1.773E+01, 1.653E+01, 1.460E+01, 1.312E+01, 1.194E+01, 1.099E+01, 1.019E+01, 9.517E+00, 8.942E+00, 8.443E+00, 8.006E+00, 7.620E+00, 7.277E+00, 6.970E+00, 6.693E+00, 6.443E+00, 5.475E+00, 4.816E+00, 4.338E+00, 3.976E+00, 3.691E+00, 3.462E+00, 3.275E+00, 3.118E+00, 2.871E+00, 2.687E+00, 2.544E+00, 2.431E+00, 2.340E+00, 2.266E+00, 2.203E+00, 2.151E+00, 2.107E+00, 2.069E+00, 2.037E+00, 2.008E+00, 1.984E+00, 1.963E+00, 1.850E+00, 1.820E+00, 1.818E+00, 1.828E+00, 1.861E+00, 1.898E+00, 1.934E+00, 1.967E+00, 1.998E+00, 2.026E+00, 2.052E+00]

Geant4Lung = [207.317, 236.938, 262.307, 286.276, 308.601, 349.292, 385.741, 419.294, 450.471, 479.523, 507.014, 533.162, 584.034, 628.145, 667.156, 701.991, 733.253, 761.366, 786.757, 809.671, 849.11, 881.203, 905.603, 924.495, 938.794, 949.126, 956.046, 960.024, 961.469, 960.754, 958.205, 954.076, 948.62, 942.045, 898.031, 845.477, 792.3, 742.123, 696.341, 655.237, 618.598, 585.991, 530.978, 486.684, 450.376, 420.068, 394.349, 372.2, 352.88, 335.841, 320.67, 307.052, 294.744, 283.547, 273.307, 263.896, 225.987, 198.734, 177.943, 161.475, 147.636, 136.203, 126.376, 117.985, 104.378, 94.168, 85.9249, 79.1048, 73.3741, 68.4849, 64.2606, 60.5722, 57.3219, 54.4276, 51.8372, 49.5043, 47.3913, 45.468, 37.9436, 32.715, 28.8533, 25.8772, 21.5757, 19.9665, 18.6046, 16.4218, 14.7465, 13.4177, 12.3366, 11.439, 10.6811, 10.0324, 9.47048, 8.97886, 8.54495, 8.15905, 7.81353, 7.5023, 7.22044, 6.13329, 5.39319, 4.85649, 4.44951, 4.13044, 3.87377, 3.66302, 3.48705, 3.21047, 3.00374, 2.84409, 2.71767, 2.61556, 2.53176, 2.46209, 2.40352, 2.35385, 2.31139, 2.27485, 2.24325, 2.21504, 2.18996, 2.05142, 2.00513, 1.99145, 1.99159, 2.00765, 2.03004, 2.05293, 2.07467, 2.09488, 2.11345, 2.13055]

Geant4PWO = [30.6453, 36.3445, 41.1871, 45.4843, 49.3927, 56.3776, 62.5641, 68.1969, 73.3795, 78.2047, 82.7814, 87.1078, 95.8031, 103.442, 110.345, 116.636, 122.415, 127.752, 132.702, 137.31, 145.631, 152.918, 159.323, 164.958, 169.911, 174.256, 178.048, 181.339, 184.175, 186.599, 188.65, 190.357, 191.75, 192.858, 194.939, 192.934, 188.544, 182.87, 176.621, 170.239, 163.987, 158.014, 147.152, 137.787, 129.776, 122.907, 116.973, 111.803, 107.252, 103.21, 99.5876, 96.3166, 93.3418, 90.6184, 88.1111, 85.7909, 76.2303, 69.1035, 63.461, 58.8365, 55.163, 51.998, 48.8877, 46.1436, 41.5807, 38.5634, 35.9551, 33.679, 31.683, 29.9395, 28.4142, 27.054, 25.8426, 24.7448, 23.7486, 22.8366, 22.0015, 21.2338, 18.1463, 15.9268, 14.242, 12.9174, 10.9585, 10.2118, 9.57361, 8.53863, 7.72878, 7.07857, 6.54213, 6.09747, 5.71979, 5.39466, 5.1136, 4.86604, 4.64642, 4.45015, 4.27359, 4.11391, 3.96883, 3.40435, 3.01569, 2.73141, 2.51462, 2.34358, 2.20541, 2.09166, 1.99637, 1.8463, 1.73443, 1.648, 1.57961, 1.52471, 1.47979, 1.44256, 1.4114, 1.38511, 1.36277, 1.34369, 1.32731, 1.3132, 1.30102, 1.24097, 1.22999, 1.2347, 1.24543, 1.27219, 1.2992, 1.32421, 1.34687, 1.36779, 1.38649, 1.40354]

# Plot
plt.figure(figsize=(8, 5))
plt.plot(energies, dEdx_mass_lung, label='Bethe-Bloch Stopping Power')
splineH2O = bethe_interpolate(ICRU_E, ICRUstoppingTotal_H2o)
splinePMMA = bethe_interpolate(ICRU_E, ICRUstoppingTotal_PMMA)
splineAir = bethe_interpolate(ICRU_E, ICRUstoppingTotal_Air)
splineLung = bethe_interpolate(ICRU_E, Geant4Lung)
splinePWO = bethe_interpolate(ICRU_E, Geant4PWO)

sLung = [Geant4Lung[i]*0.248 + (1-0.248)*ICRUstoppingTotal_Air[i] for i in range(len(Geant4Lung))]
splineLungInf = bethe_interpolate(ICRU_E, sLung)

plt.plot(energies, splinePWO(np.log10(energies)), label='ICRU Stopping Power (PWO)')
plt.plot(energies, splineH2O(np.log10(energies)), label='ICRU Stopping Power (H2O)')
plt.plot(energies, splinePMMA(np.log10(energies)), label='ICRU Stopping Power (PMMA)')
plt.plot(energies, splineAir(np.log10(energies)), label='ICRU Stopping Power (Air)')
plt.plot(energies, splineLung(np.log10(energies)), label='ICRU Stopping Power (Lung)')
plt.plot(energies, splineLungInf(np.log10(energies)), label='ICRU Stopping Power (Lung Inf)')
#plt.plot(ICRU_E, ICRUstoppingTotal_H2o, 'o', label='ICRU Stopping Power')
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
plt.close()

def FillProb(energy):
    #S_l = stopping_power_lung(energy)*density_lung
    S_l = splineLung(np.log10(energy))
    S_l = splineLungInf(np.log10(energy))
    S_Air = splineAir(np.log10(energy))
    S_pmma = splinePMMA(np.log10(energy))
    S_h2o = splineH2O(np.log10(energy))
    #S_pmma = S_h2o
    return (S_l-S_Air)/(S_h2o-S_Air)

def CalcPmod(energy):
    #S_l = stopping_power_lung(energy)*density_lung
    S_l = splineLung(np.log10(energy))*density_lung
    S_l = splineLungInf(np.log10(energy))*0.26
    S_h2o = splineH2O(np.log10(energy))*density
    S_Air = splineAir(np.log10(energy))*density_air
    S_pmma = splinePMMA(np.log10(energy))*density_pmma
    #S_pmma = S_h2o
    w = FillProb(energy)
    return w*(1-w)*(S_h2o)**2/(S_h2o*S_l)

def CalcD(energy, pmod):
    #S_l = stopping_power_lung(energy)*density_lung
    S_l = splineLung(np.log10(energy))*density_lung
    S_l = splineLungInf(np.log10(energy))*0.26
    S_h2o = splineH2O(np.log10(energy))*density
    S_Air = splineAir(np.log10(energy))*density_air
    S_pmma = splinePMMA(np.log10(energy))*density_pmma
    #S_pmma = S_h2o
    w = FillProb(energy)
    return pmod/(w*(1-w)*(S_h2o)**2/(S_h2o*S_l))

plt.plot(energies, [FillProb(energ) for index, energ in enumerate(energies)], 'o', label='ICRU Stopping Power')
plt.xlabel('Kinetic Energy [MeV]')
plt.ylabel('Fill probability')
plt.ylim(0,1)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.plot(energies, [CalcD(energ, 300) for index, energ in enumerate(energies)], 'o', label='ICRU Stopping Power')
plt.xlabel('Kinetic Energy [MeV]')
plt.ylabel('d / um')
#plt.ylim(0,1)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()