import ROOT
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import random

# -------- Known gamma energies (keV) ----------
E1 = 1.17323
E2 = 1.33249
# E1 = (E1+E2)/2
# E1 = 0.963
# E2 = 1.117
# E1 = 1.3
# E2 = 1.5
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
show = False

datapath = "../data/paperBeamtime"

if(useCutoff):
    f = ROOT.TFile.Open(f"{datapath}/notarget/output/Histograms.root")

cutoff = [4000, 3500, 4500, 1640, 3500, 1700, 1410, 1630, 5000, 1200, 1150, 4000, 900, 700, 500, 800, 4000, 4000, 4000, 500, 620, 1200, 0, 500, 600, 600, 5000, 5500, 5000, 5000, 5000, 5000]


highgain = [0,1,2,4,8, 29, 30, 31]

highcut = []
params = []
cutoffch = [3, 5,6,7,9,10,12]

channel7 = 0

channels = [7969.63, 6810.33, 8544.35, 8366.75, 8248.32, 8399.27, 7839.18, 8067.3, 9352.22, 8115.49, 7853.99, 3972.46, 4136.24, 3851.16, 4036.96, 4328.33, 3523.18, 3633.97, 4198.85, 3682.06, 4013.73, 4020.97, 3362.79, 3471.02, 4460.77, 4230.29, 5721, 9568.39, 4908.55]

channels_Err = [1176.55, 1032.72, 1230.33, 1199.13, 1169.02, 1202.06, 1155.05, 1157.6, 1213.28, 1140.51, 1130.95, 729.974, 767.702, 695.278, 754.469, 791.3, 691.12, 732.849, 868.613, 730.777, 775.843, 783.083, 730.004, 775.425, 1025.96, 978.059, 1431.56, 2485.46, 2373.64]

# channels = [8260.4, 7005.78, 8775.9, 8635.58, 8498.3, 8627.36, 8096.46, 8361.34, 9589.23, 8386.07, 8059.63, 4085.59, 4268.36, 3952.39, 4142.66, 4433.85, 3616.99, 3750.23, 4325.86, 3819.41, 4151, 4197.05, 3512.9, 3639.87, 4802.1, 4735.61, 7268.79, 9108.95, 4908.55]

# channels_Err = [1211.02, 1046.36, 1244.99, 1228.21, 1197.43, 1190.42, 1150.61, 1155.63, 1207.02, 1145.2, 1135.46, 732.577, 769.438, 705.058, 749.772, 821.845, 701.081, 732.514, 838.599, 738.325, 788.51, 812.071, 787.064, 829.456, 1124.78, 1145.54, 2143.24, 1811.97, 2373.64]

# energies = [5.93944, 6.04039, 6.15737, 6.27501, 6.40478, 6.54775, 6.70606, 6.87805, 7.06705, 7.28057, 7.51988, 5.16491, 5.29576, 5.4461, 5.61065, 5.79852, 6.01081, 6.25395, 6.53684, 6.87087, 7.28225, 7.79209, 8.45617, 9.3807, 10.6148, 13.246, 13.2788, 20.4859, 10.7988]

# energies_Err = [0.450472,0.449976,0.450388,0.446538,0.446791,0.445129,0.445779,0.441615,0.441219,0.444092,0.441296,0.360322,0.36038,0.359645,0.362348,0.36269,0.367195,0.372248,0.378452,0.389481,0.411361,0.454613,0.527966,0.676651,1.54338,2.72723,2.30265,2.6112,5.83477]

energies = [5.71815, 5.82682, 5.92648, 6.02858, 6.1768, 6.31649, 6.47338, 6.64359, 6.81676, 7.02119, 7.25093, 4.96042, 5.08051, 5.23889, 5.37492, 5.54862, 5.73106, 5.94857, 6.19156, 6.47118, 6.81441, 7.22038, 7.71979, 8.36851, 9.27216, 10.6161, 13.0496, 20.6016, 10.7988]


energies_Err = [0.444011, 0.436206, 0.434967, 0.433878, 0.428874, 0.430721, 0.43061, 0.426988, 0.418238, 0.428288, 0.435819, 0.346543, 0.349299, 0.351787, 0.34675, 0.345176, 0.348694, 0.346973, 0.341384, 0.331524, 0.339027, 0.336028, 0.345539, 0.366425, 0.384159, 0.468517, 0.704609, 2.72207, 5.83477]


channelWidth = 3
kB = 12.68*0.001
for idx, value in enumerate(channels):
    if(idx>10): channelWidth = 2
    quenched = energies[idx]/(1 + kB*energies[idx]/channelWidth)
    a_fit = quenched/value
    print(f"E {energies[idx]}, Quenched {quenched}")
    params.append([a_fit, 0,0,0])
    
for channel in range(29,32):
    # if channel == 21: channel-=1
    f = ROOT.TFile.Open(f"{datapath}/notarget/output/coincHistogram{channel}.root")
    if(channel in cutoffch ):
        useCutoff = True
        print("usecutoff") 
    else:
        useCutoff = False

    h = f.Get(f"h_coincCharge_{channel}")

    # channel += +1
    xdata, ydata = hist_to_numpy(h)
    a_guess = 1/500
    xdata = xdata/4
    if channel not in highgain:
        xdata = xdata/4
        a_guess = 1/50
    b_guess = 0.0005
    A1_guess = 30
    A2_guess = A1_guess
    sigma_guess = 700

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
    
    prefitchannels = [8200, 7110, 8900, 8750, 8600, 8750, 8250, 8521, 11830, 8550, 8240, 4230, 4300, 4000, 4250, 4500, 3668, 3780, 4470, 3871, 4200, 4260, 3540, 3700, 4870, 4900, 7725, 8381, 4700, 4472, 5500, 5400]
    
    popt2, pcov2 = curve_fit(gauss_calib1, xdata, ydata, p0=p02, bounds=(lower2, upper2), maxfev=1000000)

    a_fit2, b_fit2, A1_fit2, sigma1_fit2 = popt2[0], popt2[1], popt2[2], popt2[3]

    c3 = ((E1+E2)/2 - b_fit2) / a_fit2

    print(f"länge: {len(xdata)}" )
    # popt, pcov = curve_fit(two_gauss_calib2, xdata, ydata, p0=p0, bounds=(lower, upper), maxfev=1000000)

    print(xdata[int(cutoff[channel]/160)])
    diff = xdata[1]-xdata[0]
    if useCutoff:
        if channel in highcut:
            print("hightcut")
            popt, pcov = curve_fit(two_gauss_calib3, xdata[:int(cutoff[channel]/diff)], ydata[:int(cutoff[channel]/diff)], p0=p0, bounds=(lower, upper), maxfev=1000000)
        else: 
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
    
    print(f"a_fit: {a_fit}")
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
    energy = (a_fit*prefitchannels[channel]+b_fit)
    print(f"Energy {energy}")
    kB = 12.68*0.001
    channelWidth = 3
    if(channel>10): channelWidth = 2
    quenched = 1/(1/energy - kB/channelWidth)
    print(f"BeamPo: {quenched:.2f} MeV")


    if channel == 7:
        channel7 = (1 - b_fit) / a_fit
    if channel == 8:
        channel8 = (1 - b_fit) / a_fit
        print(f"Quotient: {channel8/channel7} ")    

    # ---------- Plot ----------
    plt.figure(figsize=(10,6))
    # plt.step(xdata, two_gauss_calib(xdata, *p0), where="mid", label="Fit", color="red", linewidth=2)
    plt.step(xdata, ydata, where="mid", label="Data")
    # plt.plot(xdata, gauss_calib1(xdata, a_fit, b_fit, A1_fit, sigma1_fit), label="1173.23 keV", color="green", linewidth=2)
    plt.plot(xdata, gauss_calib1(xdata, a_fit, b_fit, A1_fit, sigma1_fit), label="1173.23 keV", color="green", linewidth=2)
    # plt.plot(xdata, gauss_calib2(xdata, a_fit, b_fit, A2_fit, sigma2_fit), label="1332.49 keV", color="red", linewidth=2)
    # plt.plot(xdata, gauss_calib2(xdata, a_fit, b_fit, A2_fit, sigma1_fit), label="1332.49 keV", color="red", linewidth=2)
    plt.plot(xdata, gauss_calib2(xdata, a_fit, b_fit, A1_fit, sigma1_fit), label="1332.49 keV", color="red", linewidth=2)
    yf = two_gauss_calib3(xdata, *popt)

    plt.plot(xdata, yf, label="Fit", linewidth=2)

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

print("Individual")
for idx, value in enumerate(params):
    print(f"ch {idx}: {(1.2-value[2])/value[0]}")
np.save(f"{datapath}/detector/calibParams.npy", params)
