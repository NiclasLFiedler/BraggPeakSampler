import uproot
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import scienceplots
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

plt.style.use(['science','notebook','grid']) 
# List of ROOT files to read

alpha  = 2.585e-2
alpha_o  = 6.826e-4
p = 1.738
p_o = 4.95e-3


def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

def GetRange(E):
    return alpha*E**p

def GetSigmaR(E, sigmaE):
    # print()
    # print(E)
    # print(sigmaE)
    # print(sigmaE*alpha*p*E**(p-1))
    # print(alpha_o*E**p)
    # print(alpha*E**p*np.log(E)*p_o)
    return np.sqrt((sigmaE*alpha*p*E**(p-1))**2)

def GetConv(E, sigma):
    return sigma*alpha*p*E**(p-1)

colors = [
    "#1f77b4",  # muted blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # yellow-green
]


root_files = []
pmods = [100,200,300,400,500,600,700,800]
thicknesses = [50,100,150,200]

pmods = [200, 300, 400, 500]
#thicknesses = [200]

root_files.append(f"../data/modulation/output/{0}mu_{0}mmMeans.root")

combination = []
combination.append([0,0,0])

for index, pmod in enumerate(pmods):
    for thickness in thicknesses:
        combination.append([pmod, thickness, index+1])
        
notargetR0 = 0
notargetR0sigma = 0

# Create a figure and axis for plotting
fig, ax = plt.subplots(figsize=(30, 13))
range0 = 0
sigmaE0 = 0
sigmar00 = 0
for comb in combination:
    file = ROOT.TFile(f"../data/modulation/output/{comb[0]}mu_{comb[1]}mmMeans.root")
    
    bin_centers = []
    bin_values = []
    
    # Access the histogram (replace 'yourhist' with your histogram's name)
    hist = file.Get("h_total_edep")
    # Loop through the bins of the histogram (assuming it's a 1D histogram)
    for bin in range(1, hist.GetNbinsX() + 1):
        bin_centers.append(hist.GetBinCenter(bin))
        bin_values.append(hist.GetBinContent(bin))
    bin_centers = np.array(bin_centers)
    bin_values = np.array(bin_values)
    #bin_values = np.where(bin_values == 0, 0.1, bin_values)

    mask = bin_centers > 175
    bin_centers = bin_centers[mask]
    bin_values = bin_values[mask]

    peaks, _ = find_peaks(bin_values)
    
    # Find the index of the largest peak
    peak_index = peaks[np.argmax(bin_values[peaks])]
    
    # Define the region around the largest peak for Gaussian fitting (adjust the range as needed)
    fit_range = (bin_centers > bin_centers[peak_index] - 20) & (bin_centers < bin_centers[peak_index] + 20)
    
    # Fit a Gaussian to the data in the selected range
    if comb[2] != 0:
        popt, _ = curve_fit(gaussian, bin_centers[fit_range], bin_values[fit_range], 
                        p0=[bin_values[peak_index], bin_centers[peak_index], 2], maxfev=100000)
    else:
        popt = [1, 220, 0]
       
    R0 = GetRange(popt[1])
    sigmaR0 = GetSigmaR(popt[1], popt[2])
    
    #sigmaR0 = R0 - GetRange(popt[1]-popt[2])

    if notargetR0 == 0:
        notargetR0 = R0
        sigmar00 = sigmaR0
        t = 0
        sigmat = 0
        Pmod = 0
        sigmaE0 = popt[2]
    else:
        t = notargetR0-R0
        #sigmat = GetConv(popt[1], np.sqrt(popt[2]**2-sigmaE0**2))
        sigmat = np.sqrt(sigmaR0**2-sigmar00**2)
        Pmod = sigmat**2/t*1e3
        #print(f"{comb[0]} um, {comb[1]} mm, t {t}, sigma {sigmat}, pmod {Pmod}")

    texText = f"{comb[0]} & {comb[1]} & {popt[1]:.2f}" " $\\pm$ " f"{popt[2]:.2f} & {R0:.2f}" " $\\pm$ " f"{sigmaR0:.2f} & {Pmod:.2f} \\\\"
    labeltext = f"{comb[0]} um, {comb[1]} mm, E={popt[1]:.2f} +- {popt[2]:.2f} MeV, R0 = {R0:.2f} +- {sigmaR0:.2f} mm, Pmod={Pmod:.2f} um"
    print(texText)
    ax.step(bin_centers, bin_values, where='mid', label=labeltext, linewidth=1, color=colors[comb[2]])
    centeres = bin_centers[fit_range]
    gplot = gaussian(centeres, *popt)
    mask = gplot > 1
    
    ax.plot(centeres[mask], gplot[mask], 'r--')

ax.set_yscale('log')
#ax.set_xlim(175, 225)
# Add labels and title
ax.set_xlabel('Energy / MeV')
ax.set_ylabel('Counts')
ax.set_title('Modulation power analysis from total energy distribution')

# Add a legend
#ax.legend(fontsize=12)
#fig.tight_layout()
# Show the plot
plt.show()
