import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from numpy import log

# --- Define the exponential decay function ---
def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

def fit_uncertainty(x, A, tau, cov):
    dA2 = cov[0, 0]
    dtau2 = cov[1, 1]
    covAtau = cov[0, 1]

    df_dA = np.exp(-x / tau)
    df_dtau = A * x * np.exp(-x / tau) / tau**2

    sigma2 = (df_dA**2) * dA2 + (df_dtau**2) * dtau2 + 2 * df_dA * df_dtau * covAtau
    return np.sqrt(sigma2)

# --- Load histogram from ROOT file ---
dir = ["ej200_50x50", "bc408_50x50", "bc420_50x50", "bc408_100x100", "bc420_100x100"]
scint = ["EJ200 50x50x10 mm³", "BC408 50x50x3 mm³", "BC420 50x50x5 mm³", "BC408 100x100x3 mm³", "BC420 100x100x5 mm³"]
sipmPos = ["o", "c", "u"]
sourcePos = ["Center", "Top Right", "Bottom Right", "Top Left", "Bottom Left"]
sipmPosname = ["Left", "Center", "Right"]
hist_name = "h_ampl0"

sel_dir = 0

scintData = []

for sel_dir in range(0, 5):
    data = []
    dataPos = []
    for sipmindex in range(0,3):
        for sourceindex in range(1,6):
            file_path = f"data/{dir[sel_dir]}/{dir[sel_dir]}_{sipmPos[sipmindex]}_{sourceindex}.root"

            with uproot.open(file_path) as file:
                hist = file[hist_name]
                bin_edges = hist.axis().edges()
                bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                counts = hist.values()

            # --- Define fit range (tail only, e.g. last 20% of nonzero bins) ---
            #nonzero_indices = np.nonzero(counts)[0]
            #tail_start = int(nonzero_indices[-1] - 0.5 * len(nonzero_indices))
            tail_start = int(len(bin_centers) * 0.45)
            nonzero_indices = np.nonzero(counts)[0]
            print("Nonzero indices:", len(nonzero_indices))

            tail_start = int(nonzero_indices[-1] - 0.8 * len(nonzero_indices))

            total_counts = np.sum(counts)
            cum_counts = np.cumsum(counts)
            tail_threshold = 0.99  # Use last 15% of cumulative counts

            tail_start_index = np.searchsorted(cum_counts, tail_threshold * total_counts)
            x_tail = bin_centers[tail_start_index:]
            y_tail = counts[tail_start_index:]

            sigma_y = np.sqrt(np.maximum(y_tail, 1))        
            popt, pcov = curve_fit(exp_decay, x_tail, y_tail, p0=(1e7, 1000), bounds=([0, 10], [np.inf, 1e6]), sigma=sigma_y, 
            absolute_sigma=True, maxfev=1000000)        

            residuals = y_tail - exp_decay(x_tail, *popt)
            chi2 = np.sum((residuals / sigma_y)**2)
            ndf = len(x_tail) - len(popt)
            print(f"Chi²/ndf = {chi2:.2f}/{ndf} = {chi2/ndf:.2f}")

            perr = np.sqrt(np.diag(pcov))

            def exp_decay_upper(x, A, tau, dA, dtau):
                return (A + dA) * np.exp(-x / (tau - dtau))

            def exp_decay_lower(x, A, tau, dA, dtau):
                return (A - dA) * np.exp(-x / (tau + dtau))

            mask = exp_decay(x_tail, *popt) > 0.1

            y_upper = exp_decay_upper(x_tail[mask], *popt, perr[0], perr[1])
            y_lower = exp_decay_lower(x_tail[mask], *popt, perr[0], perr[1])

            plt.figure(figsize=(10, 6))
            plt.plot(bin_centers, counts, drawstyle='steps-mid', label="Data")
            plt.plot(x_tail[mask], exp_decay(x_tail[mask], *popt), 'r--', label=f"Fit: A*exp(-x/τ)\nA={popt[0]:.2e}, τ={popt[1]:.2f}")
            plt.fill_between(x_tail[mask], y_lower, y_upper, color='red', alpha=0.3, label="Fit ± 1σ")

            plt.yscale("log")
            # plt.xscale("log")
            plt.xlabel("Amplitude / ADC Channel")
            plt.ylabel("Counts")
            plt.title(f"$^{{90}}$Sr Source at {sourcePos[sourceindex-1]} Position with {scint[sel_dir]} Scintillator,\n Read Out by Hamamatsu SiPM S14160-1350 at {sipmPosname[sipmindex]} Position")

            # Calculate the intersection of the fit with y = 1
            A, tau = popt
            dA, dtau = perr

            # Intersection point (x value where y = 1)
            x_intersect = -tau * log(1 / A)

            # Propagate uncertainty
            dx_dA = tau / A  # Partial derivative of x with respect to A
            dx_dtau = -log(1 / A)  # Partial derivative of x with respect to tau

            x_intersect_std = np.sqrt((dx_dA * dA)**2 + (dx_dtau * dtau)**2)

            print(f"Intersection with y=1: x = {x_intersect:.2f} ± {x_intersect_std:.2f}")

            # Optionally, you can add this information to the plot
            plt.axvline(x_intersect, color='blue', linestyle='--', label=f"Intersection: x={x_intersect:.2f} ± {x_intersect_std:.2f}")

            plt.legend()
            plt.grid(True, which="both", ls="--")
            plt.tight_layout()

            plt.savefig(f"data/{dir[sel_dir]}/pdf/Sr90_fit_{dir[sel_dir]}_{sipmPos[sipmindex]}_{sourceindex}.pdf", format="pdf")
            #plt.show()
            plt.close()
            if sel_dir == 0:
                if sipmindex == 2:
                    if sourceindex == 5:
                        x_intersect = x_intersect * 4
                        x_intersect_std = x_intersect_std * 4

            if sel_dir == 2:
                if sipmindex == 1:
                    if sourceindex == 4:
                        x_intersect = x_intersect / 4
                        x_intersect_std = x_intersect_std / 4

            if sel_dir == 3:
                if sipmindex == 0:
                    if sourceindex == 3:
                        x_intersect = x_intersect * 4
                        x_intersect_std = x_intersect_std * 4
                if sipmindex == 2:
                    if sourceindex == 5:
                        x_intersect = x_intersect * 4
                        x_intersect_std = x_intersect_std * 4

            if sel_dir == 4:
                if sipmindex == 0:
                    if sourceindex == 3:
                        x_intersect = x_intersect * 4
                        x_intersect_std = x_intersect_std * 4
                if sipmindex == 2:
                    if sourceindex == 5:
                        x_intersect = x_intersect * 4
                        x_intersect_std = x_intersect_std * 4

            dataPos.append([x_intersect, x_intersect_std])
        data.append(dataPos)
        dataPos = []

    normalized_intersections = []

    x_ref = data[1][0][0]
    x_ref_std = data[1][0][1]

    for i, dataPos in enumerate(data):
        for j, (x_intersect, x_intersect_std) in enumerate(dataPos):
            # Normalize the intersection value
            x_norm = x_intersect / x_ref

            # Propagate the uncertainty
            x_norm_std = x_norm * np.sqrt((x_intersect_std / x_intersect)**2 + (x_ref_std / x_ref)**2)

            # Store the normalized value and its std dev
            normalized_intersections.append({
                "sipm_index": i + 1,
                "source_index": j + 1,
                "x_norm": x_norm,
                "x_norm_std": x_norm_std,
                "xintersect": x_intersect,
                "x_intersect_std": x_intersect_std
            })

    # Print the normalized intersections
    for result in normalized_intersections:
        print(f"SiPM {result['sipm_index']}, Source {result['source_index']}: channel_norm = {result['x_norm']:.2f} ± {result['x_norm_std']:.3f}, channel = {result['xintersect']:.2f} ± {result['x_intersect_std']:.2f}")


    plt.figure(figsize=(10, 6))
    plt.xticks([1,2,3,4,5], ["Bottom Left", "Top Left", "Center", "Top Right", "Bottom Right"])

    leftSipmvalue = [[normalized_intersections[i]['x_norm'], normalized_intersections[i]['x_norm_std']] for i in range(10, 15)]

    centerSipmvalue = [[normalized_intersections[i]['x_norm'], normalized_intersections[i]['x_norm_std']] for i in range(5, 10)]

    rightSipmvalue = [[normalized_intersections[i]['x_norm'], normalized_intersections[i]['x_norm_std']] for i in range(0, 5)]

    absleft = [[normalized_intersections[i]['xintersect'], normalized_intersections[i]['x_intersect_std']] for i in range(10, 15)]
    abscenter = [[normalized_intersections[i]['xintersect'], normalized_intersections[i]['x_intersect_std']] for i in range(5, 10)]
    absright = [[normalized_intersections[i]['xintersect'], normalized_intersections[i]['x_intersect_std']] for i in range(0, 5)]

    leftSipmvalue = [leftSipmvalue[4], leftSipmvalue[3], leftSipmvalue[0],leftSipmvalue[1], leftSipmvalue[2]] 
    centerSipmvalue = [centerSipmvalue[4], centerSipmvalue[3], centerSipmvalue[0],centerSipmvalue[1], centerSipmvalue[2]]
    rightSipmvalue = [rightSipmvalue[4], rightSipmvalue[3], rightSipmvalue[0],rightSipmvalue[1], rightSipmvalue[2]]

    absleft = [absleft[4], absleft[3], absleft[0],absleft[1], absleft[2]] 
    abscenter = [abscenter[4], abscenter[3], abscenter[0],abscenter[1], abscenter[2]]
    absright = [absright[4], absright[3], absright[0],absright[1], absright[2]]


    plt.errorbar([1,2,3,4,5], [val[0] for val in leftSipmvalue], yerr=[val[1] for val in leftSipmvalue], fmt='o', capsize=4, label='Left SiPM Position', color='blue')
    plt.plot([1,2,3,4,5], [val[0] for val in leftSipmvalue], linestyle='--', linewidth=1,color='blue')

    plt.errorbar([1,2,3,4,5], [val[0] for val in centerSipmvalue], yerr=[val[1] for val in centerSipmvalue], fmt='o', capsize=4, label='Center SiPM Position', color='red')
    plt.plot([1,2,3,4,5], [val[0] for val in centerSipmvalue], linestyle='--', linewidth=1,color='red')

    plt.errorbar([1,2,3,4,5], [val[0] for val in rightSipmvalue], yerr=[val[1] for val in rightSipmvalue], fmt='o', capsize=4, label='Right SiPM Position', color='orange')
    plt.plot([1,2,3,4,5], [val[0] for val in rightSipmvalue], linestyle='--', linewidth=1,color='orange')

    plt.xlabel("Source Position")
    plt.ylabel("Normalized Intersection / ADC Channel")
    plt.title(f"Normalized y-Intersection of $^{{90}}$Sr Source \n for {scint[sel_dir]} Scintillator, Read Out by Hamamatsu SiPM S14160-1350")
    plt.grid()
    plt.legend()

    plt.savefig(f"data/{dir[sel_dir]}/pdf/Sr90_fit_{dir[sel_dir]}.pdf", format="pdf")

    plt.yscale("log")
    plt.savefig(f"data/{dir[sel_dir]}/pdf/Sr90_fit_{dir[sel_dir]}_log.pdf", format="pdf")

    scintData.append([absleft, abscenter, absright])
    plt.close()

plt.figure(figsize=(15, 9))
plt.xticks([1,2,3,4,5], ["Bottom Left", "Top Left", "Center", "Top Right", "Bottom Right"])
colors = ["green", "red", "blue", "orangered", "cyan"]

spos = 1

for index, scintD in enumerate(scintData):
    plt.errorbar([1,2,3,4,5], [val[0]/scintData[0][spos][2][0] for val in scintD[spos]], yerr=[val[1]/scintData[0][spos][2][0] for val in scintD[spos]], fmt='o', capsize=4, label=f'{scint[index]}', color=colors[index])
    plt.plot([1,2,3,4,5], [val[0]/scintData[0][spos][2][0] for val in scintD[spos]], linestyle='--', linewidth=1,color=colors[index])
plt.xlabel("Source Position")
plt.ylabel("EJ200 Center Normed Intersection / ADC Channel")
plt.title(f"y-Intersection of $^{{90}}$Sr Source \n for Different Scintillator Types, Read Out by Hamamatsu SiPM S14160-1350")
plt.grid()  
plt.legend()
plt.savefig(f"data/All_Sr90_fit_left.pdf", format="pdf")
plt.yscale("log")
plt.savefig(f"data/All_Sr90_fit_log_left.pdf", format="pdf")
plt.show()
plt.close()

plt.figure(figsize=(15, 9))
plt.xticks([1,2,3,4,5], ["Bottom Left", "Top Left", "Center", "Top Right", "Bottom Right"])
colors = ["green", "red", "blue", "orangered", "cyan"]
for index, scintD in enumerate(scintData):
    if index == 0:
        continue
    plt.errorbar([1,2,3,4,5], [(scintData[0][1][idx][0]-val[0])/scintData[0][1][2][0] for idx, val in enumerate(scintD[1])], yerr=[val[1]/scintData[0][1][2][0] for val in scintD[1]], fmt='o', capsize=4, label=f'{scint[index]}', color=colors[index])
    plt.plot([1,2,3,4,5], [(scintData[0][1][idx][0]-val[0])/scintData[0][1][2][0] for idx, val in enumerate(scintD[1])], linestyle='--', linewidth=1,color=colors[index])
plt.xlabel("Source Position")
plt.ylabel("Delta EJ200 Center Normed Intersection / ADC Channel")
plt.title(f"Delta y-Intersection of $^{{90}}$Sr Source \n for Different Scintillator Types, Read Out by Hamamatsu SiPM S14160-1350")
plt.grid()  
plt.legend()
plt.savefig(f"data/Diff_Sr90_fit.pdf", format="pdf")
plt.show()
plt.close()