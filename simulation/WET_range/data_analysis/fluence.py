import ROOT
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uproot  # uproot is a great library for reading ROOT files in Python
matplotlib.use('TkAgg')  # or 'Qt5Agg'

def range_energy_relationship(energy, alpha, p):
    return alpha * energy**p

def linear(x, slope, intercept):
    return x*slope+intercept

def main():
    target = "fluence"
    fluenceFile = uproot.open(f"{target}.root")
    fluenceTree = fluenceFile["fluence"]
    depth = fluenceTree["depth"].array().to_numpy()
    counts = fluenceTree["counts"].array().to_numpy()
    print("depth ", depth)
    print("counts ", counts)
    
    energy = 250
    a_h2o=2.585e-3
    a_pwo=7.275e-4
    p_h2o=1.738
    p_pwo=1.690
    
    R0 = range_energy_relationship(energy, a_h2o, p_h2o)    
    print(f"R0: {R0}")
    depth = [R0-(x * 0.1) for x in depth]
    popt, pcov = curve_fit(linear, depth[:-4], counts[:-4])
    counts = [x/linear(0, popt[0], popt[1]) for x in counts]
    popt, pcov = curve_fit(linear, depth[:-4], counts[:-4])
    print(f"Beta: {popt[0]} 1/cm")
    plt.rcParams["text.usetex"] = True
    plt.plot(depth, counts, marker='o', linestyle='',label="Fluence simulation")
    plt.plot(np.linspace(0, 40, 50), linear(np.linspace(0, 40, 50), popt[0], popt[1]))
    plt.xlabel(r"$R_0-z$ / mm")
    plt.ylabel(r"$\Phi_0$ / arb. units")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()