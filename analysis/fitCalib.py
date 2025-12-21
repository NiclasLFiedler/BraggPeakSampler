import ROOT
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# -------- Known gamma energies (keV) ----------
E1 = 1.17323
E2 = 1.33249
# E1 = (E1+E2)/2
# E1 = 0.963
# E2 = 1.117
# -------- Define the composite model ----------
def two_gauss_calib(x, a, b, A1, A2, sigma1, sigma2):
    # Convert energy to channel
    c1 = (E1 - b) / a
    c2 = (E2 - b) / a
    g1 = A1 * np.exp(-0.5 * ((x - c1) / sigma1)**2)
    g2 = A2 * np.exp(-0.5 * ((x - c2) / sigma2)**2)
    return g1 + g2

def two_gauss_calib2(x, a, b, A1, A2, sigma):
    # Convert energy to channel
    c1 = (E1 - b) / a
    c2 = (E2 - b) / a
    g1 = A1 * np.exp(-0.5 * ((x - c1) / sigma)**2)
    g2 = A2 * np.exp(-0.5 * ((x - c2) / sigma)**2)
    return g1 + g2

def two_gauss_calib3(x, a, b, A, sigma):
    # Convert energy to channel
    c1 = (E1 - b) / a
    c2 = (E2 - b) / a
    g1 = A * np.exp(-0.5 * ((x - c1) / sigma)**2)
    g2 = A * np.exp(-0.5 * ((x - c2) / sigma)**2)
    return g1 + g2

def gauss_calib1(x, a, b, A1, sigma1):
    # Convert energy to channel
    c1 = (E1 - b) / a

    g1 = A1 * np.exp(-0.5 * ((x - c1) / sigma1)**2)
    return g1

def gauss_calib2(x, a, b, A1, sigma1):
    # Convert energy to channel
    c1 = (E2 - b) / a

    g1 = A1 * np.exp(-0.5 * ((x - c1) / sigma1)**2)
    return g1

# -------- Function to extract x,y from TH1D ----------
def hist_to_numpy(hist):
    """
    Converts a ROOT TH1D into numpy arrays of bin centers and bin contents.
    """
    nb = hist.GetNbinsX()
    x = np.zeros(nb)
    y = np.zeros(nb)

    for i in range(1, nb+1):
        x[i-1] = hist.GetBinCenter(i)
        y[i-1] = hist.GetBinContent(i)

    return x, y


# -------- Load your file and histogram ----------
# Example:
# f = ROOT.TFile.Open("yourfile.root")
# h = f.Get("your_hist_name")

# Replace these two lines with your actual code:

useCutoff = False
show = True
if(useCutoff):
    f = ROOT.TFile.Open("../data/paperBeamtime/notarget/output/Histograms.root")

cutoff = [4000, 3500, 4500, 4000, 3500, 4000, 3500, 3500, 5000, 1000, 3500, 4000, 500, 500, 500, 800, 4000, 4000, 4000, 4000, 5000, 5000, 5500, 5000, 5000, 5000, 5000, 5500, 5000, 5000, 5000, 5000]

highgain = [0,4,5]
highgain = [0]
params = []
cutoffch = [12,13,14,15]
for channel in range(18,32):
    f = ROOT.TFile.Open(f"../data/paperBeamtime/notarget/output/coincHistogram{channel}.root")
    if(channel in cutoffch ):
        useCutoff = True
    else:
        useCutoff = False

    h = f.Get(f"h_coincCharge_{channel}")

    xdata, ydata = hist_to_numpy(h)
    a_guess = 1/500
    xdata = xdata/4
    if channel not in highgain:
        xdata = xdata/4
        a_guess = 1/50
    b_guess = 0.0005
    A1_guess = 30
    A2_guess = A1_guess
    sigma_guess = 300

    # p0 = [a_guess, b_guess,
    #       A1_guess, A2_guess,
    #       sigma_guess, sigma_guess+50]
    
    # p0 = [a_guess, b_guess,
    #       A1_guess, A2_guess,
    #       sigma_guess]

    p0 = [a_guess, b_guess,
          A1_guess,
          sigma_guess]

    # Bounds (to keep model physical)
    # lower = [0, 0, 0, 0, 0.1, 0.1]
    # upper = [1,  0.05,  2000, 2000, 20000, 20000]

    # lower = [0, 0, 15, 15, 0.1]
    # upper = [1,  0.05,  2000, 2000, 20000]

    lower = [0, 0, 0, 0.1]
    upper = [1,  0.01,  20000, 20000]

    # -------- Perform fit ----------
    p02 = [a_guess, b_guess, A1_guess, sigma_guess]
    lower2 = [0, 0, 0, 0.1]
    upper2 = [1, 0.1, 20000, 20000]
    
    prefitchannels = [8200, 7110, 8900, 8750, 8600, 8750, 8250, 8521, 11830, 8550, 8240, 4230, 4300, 4000, 4250, 4500, 3668, 4470, 3871, 4200, 4260, 3540, 3700, 4870, 4900, 7725, 8381, 4700, 4472, 5500, 5400]
    
    popt2, pcov2 = curve_fit(gauss_calib1, xdata, ydata, p0=p02, bounds=(lower2, upper2), maxfev=1000000)

    a_fit2, b_fit2, A1_fit2, sigma1_fit2 = popt2[0], popt2[1], popt2[2], popt2[3]

    c3 = ((E1+E2)/2 - b_fit2) / a_fit2

    print(f"länge: {len(xdata)}" )
    # popt, pcov = curve_fit(two_gauss_calib2, xdata, ydata, p0=p0, bounds=(lower, upper), maxfev=1000000)

    print(xdata[int(cutoff[channel]/160)])
    diff = xdata[1]-xdata[0]
    if useCutoff:
        popt, pcov = curve_fit(two_gauss_calib3, xdata[int(cutoff[channel]/diff):], ydata[int(cutoff[channel]/diff):], p0=p0, bounds=(lower, upper), maxfev=1000000)
    else:
        popt, pcov = curve_fit(two_gauss_calib3, xdata, ydata, p0=p0, bounds=(lower, upper), maxfev=1000000)

    perr = np.sqrt(np.diag(pcov))
    # a_fit, b_fit, A1_fit, A2_fit, sigma1_fit, sigma2_fit = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]
    # a_err, b_err, A1_error, A2_error, sigma1_error, sigma2_error = perr[0], perr[1], perr[2], perr[3], perr[4], perr[5]

    # a_fit, b_fit, A1_fit, A2_fit, sigma1_fit = popt[0], popt[1], popt[2], popt[3], popt[4]
    # a_err, b_err, A1_error, A2_error, sigma1_error = perr[0], perr[1], perr[2], perr[3], perr[4]

    a_fit, b_fit, A1_fit, sigma1_fit = popt[0], popt[1], popt[2], popt[3]
    a_err, b_err, A1_error, sigma1_error = perr[0], perr[1], perr[2], perr[3]

    params.append([a_fit, a_err, b_fit, b_err])

    # -------- Print results ----------
    print(f"\n=== Energy Calibration Fit Results {channel} ===")
    print(f"a = {a_fit:.6f} ± {a_err:.6f}   (keV/channel)")
    print(f"b = {b_fit:.3f} ± {b_err:.3f}   (keV)")
    print(f"A1 = {A1_fit:.2f} ± {A1_error:.2f}   (counts)")
    # print(f"A2 = {A2_fit:.2f} ± {A2_error:.2f}   (counts)")
    print(f"sigma1 = {sigma1_fit:.2f} ± {sigma1_error:.2f}   (channels)")
    # print(f"sigma2 = {sigma2_fit:.2f} ± {sigma2_error:.2f}   (channels)")

    c1 = (E1 - b_fit) / a_fit
    c2 = (E2 - b_fit) / a_fit
    print(f"Peak centers (channels): c1 = {c1:.2f}, c2 = {c2:.2f}")
    print(f"Separation: {c2 - c1:.2f} channels")
    print(f"BeamPo: {(a_fit*prefitchannels[channel]+b_fit):.2f} MeV")

    # ---------- Plot ----------
    plt.figure(figsize=(10,6))
    # plt.step(xdata, two_gauss_calib(xdata, *p0), where="mid", label="Fit", color="red", linewidth=2)
    plt.step(xdata, ydata, where="mid", label="Data")
    # plt.plot(xdata, gauss_calib1(xdata, a_fit, b_fit, A1_fit, sigma1_fit), label="1173.23 keV", color="green", linewidth=2)
    plt.plot(xdata, gauss_calib1(xdata, a_fit, b_fit, A1_fit, sigma1_fit), label="1173.23 keV", color="green", linewidth=2)
    # plt.plot(xdata, gauss_calib2(xdata, a_fit, b_fit, A2_fit, sigma2_fit), label="1332.49 keV", color="red", linewidth=2)
    # plt.plot(xdata, gauss_calib2(xdata, a_fit, b_fit, A2_fit, sigma1_fit), label="1332.49 keV", color="red", linewidth=2)
    plt.plot(xdata, gauss_calib2(xdata, a_fit, b_fit, A1_fit, sigma1_fit), label="1332.49 keV", color="red", linewidth=2)
    xf = np.linspace(min(xdata), max(xdata), 2000)
    yf = two_gauss_calib3(xf, *popt)

    plt.plot(xf, yf, label="Fit", linewidth=2)

    plt.xlabel("Channel")
    plt.ylabel("Counts")
    # plt.yscale("log")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.close()
print("\n=== Summary of Calibration Parameters ===")
np.save("calibParams.npy", params)
