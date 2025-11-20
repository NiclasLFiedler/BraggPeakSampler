import matplotlib.pyplot as plt
import numpy as np

ped = 85.4
pe = 110.53
#na22 = 0.059409
na22 = 0.511

mean_SPE = 110.342
sigma_SPE= 10.0004
mean_ped= 85.5295
sigma_ped= 0.48595
PDE = 0.2316
na22 = 0.511

def normalize(data):
    # Check if the first element is zero to avoid division by zero
    if data[0] == 0:
        raise ValueError("The first element of the list is zero, cannot normalize.")
    return [value / data[0] for value in data]

def convert_to_secondary_units(y):
    return y*50.65/11.73  # Example: converting y to another unit system

def GetPE(channel):
    return abs(((channel-ped)/(pe-ped))/na22)

def GetPH(channel):
    return abs(((channel-ped)/(pe-ped))/na22/0.2316)

def GetPEDev(channel):
    return abs(((channel)/(pe-ped))/na22)

def GetPHDev(channel):
    return abs(((channel)/(pe-ped))/na22/0.2316)

def mean_and_std_of_mean(values, std_devs):
    """
    Calculate the mean and the standard deviation of the mean
    (not weighted).
    
    Parameters:
        values (array-like): Measured values.
        std_devs (array-like): Standard deviations of each value.
        
    Returns:
        (float, float): (mean, standard deviation of the mean)
    """
    values = np.asarray(values)
    std_devs = np.asarray(std_devs)
    
    mean = np.mean(values)
    # print(np.sum(std_devs**2))
    test = 0
    for i in std_devs:
        test = test + i**2

    std_of_mean = np.sqrt(np.sum(std_devs**2)) / len(values)
    
    return mean, std_of_mean

def main():        
    showall = False
    energyNa22 = 0.511
    flatparams = np.load('data/2504/flat/output/fit_params.npy', allow_pickle=True)
    lateralparams = np.load('data/2504/lateral/output/fit_params.npy', allow_pickle=True)
    flatLightyield = []
    lateralLightyield = []
    mean = []
    stdDev = []
    meanErr = []
    stdDevErr = []
    for popt, perr in flatparams:
        mean.append(popt[1]/energyNa22)
        stdDev.append(popt[2]/energyNa22)
        meanErr.append(perr[1]/energyNa22)
        stdDevErr.append(perr[2]/energyNa22)

    flatLightyield.append(mean)
    flatLightyield.append(stdDev)
    flatLightyield.append(meanErr)
    flatLightyield.append(stdDevErr)

    mean = []
    stdDev = []
    meanErr = []
    stdDevErr = []
    for popt, perr in lateralparams:
        mean.append(popt[1]/energyNa22)
        stdDev.append(popt[2]/energyNa22)
        meanErr.append(perr[1]/energyNa22)
        stdDevErr.append(perr[2]/energyNa22)

    lateralLightyield.append(mean)
    lateralLightyield.append(stdDev)
    lateralLightyield.append(meanErr)
    lateralLightyield.append(stdDevErr)

    flatLightyield = np.array(flatLightyield)
    lateralLightyield = np.array(lateralLightyield)

    flat_means = np.array(flatLightyield[0])
    flat_indices = np.arange(len(flat_means))
    thickflat_sorted_indices = np.argsort(flat_means[:15])[::-1][:11]
    thinnflat_sorted_indices = np.argsort(flat_means[15:])[::-1] + 15
    thinnflat_sorted_indices= np.concatenate((thinnflat_sorted_indices[3:], thinnflat_sorted_indices[:3]))


    thickFlatMean = flatLightyield[0][thickflat_sorted_indices]
    thickFlatstddev = flatLightyield[1][thickflat_sorted_indices]
    thickFlatErr = flatLightyield[2][thickflat_sorted_indices]
    thickFlatstddevErr = flatLightyield[3][thickflat_sorted_indices]
    
    thinnFlatMean = flatLightyield[0][thinnflat_sorted_indices]
    thinnFlatstddev = flatLightyield[1][thinnflat_sorted_indices]
    thinnFlatErr = flatLightyield[2][thinnflat_sorted_indices]
    thinnFlatstddeverr = flatLightyield[3][thinnflat_sorted_indices]
    
    thinnFlatMean = flatLightyield[0][thinnflat_sorted_indices]
    thinnFlatErr = flatLightyield[2][thinnflat_sorted_indices]
    
    lateral_means = np.array(lateralLightyield[0])    
    lateral_indices = np.arange(len(lateral_means))
    thicklateral_sorted_indices = np.argsort(lateral_means[:15])[::-1][:11]
    thinnlateral_sorted_indices = np.argsort(lateral_means[15:])[::-1] + 15
    thinnlateral_sorted_indices= np.concatenate((thinnlateral_sorted_indices[3:], thinnlateral_sorted_indices[:3]))

    thicklateralMean = lateralLightyield[0][thicklateral_sorted_indices]
    thicklateralstddev = lateralLightyield[1][thicklateral_sorted_indices]
    thicklateralErr = lateralLightyield[2][thicklateral_sorted_indices]
    thicklateralstddeverr = lateralLightyield[3][thicklateral_sorted_indices]
    
    thinnlateralMean = lateralLightyield[0][thinnlateral_sorted_indices]
    thinnlateralstddev = lateralLightyield[1][thinnlateral_sorted_indices]
    thinnlateralErr = lateralLightyield[2][thinnlateral_sorted_indices]
    thinnlateralstddeverr = lateralLightyield[3][thinnlateral_sorted_indices]

    first15_flat_mean = np.mean(flat_means[:15])
    last21_flat_mean = np.mean(flat_means[-21:])
    first15_flat_std = np.std(flat_means[:15], ddof=1)
    last21_flat_std = np.std(flat_means[-21:], ddof=1)

    first15_lateral_mean = np.mean(lateral_means[:15])
    last21_lateral_mean = np.mean(lateral_means[-21:])
    first15_lateral_std = np.std(lateral_means[:15], ddof=1)
    last21_lateral_std = np.std(lateral_means[-21:], ddof=1)
    
    
    thickflataverage = np.mean(thickFlatMean)
    thinnflataverage = np.mean(thinnFlatMean)
    thickflatstd = np.std(thickFlatMean, ddof=1)
    thinnflatstd = np.std(thinnFlatMean, ddof=1)

    thicklateralaverage = np.mean(thicklateralMean)
    thinnlateralaverage = np.mean(thinnlateralMean)
    thicklateralstd = np.std(thicklateralMean, ddof=1)
    thinnlateralstd = np.std(thinnlateralMean, ddof=1)

    plt.rcParams.update({'font.size': 28})
    plt.figure(figsize=(18, 15))
    
    if(showall):
        plt.errorbar(range(len(flatLightyield[0])), flatLightyield[0], yerr=flatLightyield[1], fmt='none',linewidth = 3, ecolor='#FFD580',capsize=5, label='Std. dev.')

        plt.errorbar(range(len(lateralLightyield[0])), lateralLightyield[0], yerr=lateralLightyield[2], fmt='o', color = "green", capsize=5, label='Mean ± error')
        plt.errorbar(range(len(lateralLightyield[0])), lateralLightyield[0], yerr=lateralLightyield[1], fmt='none',linewidth = 3, ecolor='#90EE90',capsize=5,  label='Std. dev.')

        plt.text(-1, first15_flat_mean-3, f'{first15_flat_mean:.2f} \n' f'$\\pm${first15_flat_std:.2f}', color='dodgerblue', fontsize=16, va='bottom', ha='center')
        plt.text(36, last21_flat_mean-3, f'{last21_flat_mean:.2f}\n' f'$\\pm${last21_flat_std:.2f}', color='dodgerblue', fontsize=16, va='bottom', ha='center')

        plt.text(-1, first15_lateral_mean-3, f'{first15_lateral_mean:.2f}\n' f'$\\pm${first15_lateral_std:.2f}', color='crimson', fontsize=16,va='bottom', ha='center')
        plt.text(36, last21_lateral_mean-3, f'{last21_lateral_mean:.2f}\n' f'$\\pm${last21_lateral_std:.2f}', color='crimson', fontsize=16, va='bottom', ha='center')

        plt.plot(range(15), [first15_flat_mean] * 15, color='dodgerblue', linestyle='--')
        plt.plot(range(15, 36), [last21_flat_mean] * 21, color='dodgerblue', linestyle='--')

        plt.plot(range(15), [first15_lateral_mean] * 15, color='crimson', linestyle='--') 
        plt.plot(range(15, 36), [last21_lateral_mean] * 21, color='crimson', linestyle='--')
        
        plt.xlabel('Crystal index')
        plt.ylabel('Light yield / ph/MeV')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    else:
        plt.errorbar(range(len(thickFlatMean)), thickFlatMean, yerr=thickFlatErr, fmt='o', color = "orange", markersize=8, capsize=8, label='Flat: Mean ± error')
        
        plt.errorbar(range(len(thickFlatMean)), thickFlatMean, yerr=thickFlatstddev, fmt='none',linewidth = 6, ecolor='#FFD580',capsize=8, label='Flat: Std. dev.')

        plt.errorbar(range(len(thickFlatMean), len(thickFlatMean)+len(thinnFlatMean)), thinnFlatMean, yerr=thinnFlatErr, fmt='o', color = "orange", markersize=8,  capsize=8)
        
        plt.errorbar(range(len(thickFlatMean), len(thickFlatMean)+len(thinnFlatMean)), thinnFlatMean, yerr=thinnFlatstddev, fmt='none',linewidth = 6, ecolor='#FFD580',capsize=8)


        plt.errorbar(range(len(thicklateralMean)), thicklateralMean, yerr=thicklateralErr, fmt='o', color = "green", markersize=8, capsize=8, label='Vertical: Mean ± error')
        
        plt.errorbar(range(len(thicklateralMean)), thicklateralMean, yerr=thicklateralstddev, fmt='none',linewidth = 6, ecolor='#90EE90',capsize=8, label='Vertical: Std. dev.')

        plt.errorbar(range(len(thicklateralMean), len(thicklateralMean)+len(thinnlateralMean)), thinnlateralMean, yerr=thinnlateralErr, fmt='o', color = "green", markersize=8, capsize=8)
        
        plt.errorbar(range(len(thicklateralMean), len(thicklateralMean)+len(thinnlateralMean)), thinnlateralMean, yerr=thinnlateralstddev, fmt='none',linewidth = 6, ecolor='#90EE90',capsize=8)

        plt.plot(range(len(thickFlatMean)), [thickflataverage] * len(thickFlatMean), color='dodgerblue', linestyle='--', linewidth=6, label=f'Flat average: 3 mm: {thickflataverage:.2f}' fr'$\pm$' f'{thickflatstd:.2f} & 2 mm: ' f'{thinnflataverage:.2f}' fr'$\pm$' f'{thinnflatstd:.2f}')
        plt.plot(range(len(thickFlatMean), len(thickFlatMean)+len(thinnFlatMean)), [thinnflataverage] * len(thinnFlatMean), linewidth=6, color='dodgerblue', linestyle='--')

        plt.plot(range(len(thicklateralMean)), [thicklateralaverage] * len(thicklateralMean), color='crimson', linewidth=6,linestyle='--', label=f'Vertical average: 3 mm: {thicklateralaverage:.2f}' fr'$\pm$' f'{thicklateralstd:.2f} & 2 mm: ' f'{thinnlateralaverage:.2f}' fr'$\pm$' f'{thinnlateralstd:.2f}') 
        plt.plot(range(len(thicklateralMean), len(thicklateralMean)+len(thinnlateralMean)), [thinnlateralaverage] * len(thinnlateralMean), color='crimson', linewidth=6, linestyle='--')
        
        plt.xlabel('Crystal index')
        plt.ylabel('Light yield / ph/MeV')
        plt.ylim(0)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("lightyield_comparison.pdf", format="pdf", dpi=300)
        plt.show()
if __name__ == "__main__":
    main()
