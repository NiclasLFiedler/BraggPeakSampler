import matplotlib.pyplot as plt

def normalize(data):
    # Check if the first element is zero to avoid division by zero
    if data[0] == 0:
        raise ValueError("The first element of the list is zero, cannot normalize.")
    return [value / data[0] for value in data]

def convert_to_secondary_units(y):
    return y*50.65/11.73  # Example: converting y to another unit system

def main():
    LY = [354.7, 337.1, 335.4, 403.1, 437.1, 483.9, 431.5, 411.1, 420.4, 435.7, 449.7, 463.9, 405.6, 422.8, 390.5]
    LY_err = [118.3, 117.2, 115.5, 123.8, 129.2, 135.8, 127.1, 125.5, 124.2, 126.7, 128.8, 137.4, 124.3, 127.9, 122.3]
    #LY = [LY_old[0],LY_old[1],LY_old[2],LY_old[3],LY_old[4],LY_old[5],LY_old[6],LY_old[8],LY_old[9],LY_old[10],LY_old[11],LY_old[14],LY_old[12],LY_old[13],LY_old[7]]
    muon = [970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8]
    co60 = [405, 410.2, 415.6, 447.3, 441.77, 438.24, 416.8, 48.9, 408.8, 422.34, 461.6, 51.59, 51.39, 49.26, 36.42]
    old_config = range(15)
    new_config = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 14, 12, 13, 7]
    # Normalize muon and co60 data
    normalized_muon = normalize(muon)
    normalized_co60 = normalize(co60)

    ped = 149.4
    pe = 175.85
    cs137 = 0.6617
    LY_pe = []
    LY_pe_err = []
    LY_ph = []
    LY_ph_err = []
    ly_average = 0
    ly_average_err = 0
    new = 0
    for index in new_config:
        LY_pe.append(((LY[index]-ped)/(pe-ped))/cs137)
        LY_pe_err.append((LY_err[index]/(pe-ped))/cs137)
        LY_ph.append((((LY[index]-ped)/(pe-ped))/cs137)/0.2316)
        LY_ph_err.append(((LY_err[index]/(pe-ped))/cs137)/0.2316)
        print(f"CHold {index} CH new {new} Light yield: {LY[index]} ± {LY_err[index]:.2f} LY_pe: {LY_pe[-1]:.2f} ± {LY_pe_err[-1]:.2f}  LY_ph: {LY_ph[-1]:.2f} ± {LY_ph_err[-1]:.2f}")
        ly_average = ly_average + LY_ph[-1]
        ly_average_err = ly_average_err + LY_ph_err[-1]
        new = new + 1
    print(f"\nAverage Light yield: {ly_average/15}")
    print(f"Average Light yield Error: {ly_average_err/15}")
    #normalized_LY_pe = normalize(LY_pe)
    #normalized_LY_ph = normalize(LY_ph)
    # Plotting
    
    # Plotting
    
    # Create the primary plot
    plt.rcParams.update({'font.size': 14})
    fig, ax1 = plt.subplots(figsize=(10, 6))
    color = 'tab:blue'
    ax1.set_xlabel('Crystal Index')
    ax1.set_ylabel('Light yield  [ph/MeV]')
    ax1.errorbar(range(15), LY_ph, yerr=LY_ph_err, fmt='o', capsize=5, elinewidth=1, capthick=1, color="black", label='Light yield data points')
    ax1.tick_params(axis='y')
    ax1.grid()

    ax1.set_ylim(0,120)

    twin1 = ax1.twinx()
    twin1.set_ylim(0,120*11.73/50.65)
    twin1.set_ylabel('Light yield  [p.e./MeV]')

    #twin2 = ax1.twinx()
    #twin2.set_ylim(0,120*354.7/50.65)
    #twin2.set_ylabel('Cs-137 Peak [ADC]')

    # Create the first secondary y-axis
    # ax2 = ax1.secondary_yaxis("right", functions=(lambda x: x * 11.73/50.65, lambda x: x * 50.65/11.73))
    # ax2.set_ylabel('Lightyield  [p.e./MeV]')
    
    # ax3 = ax1.secondary_yaxis("right", functions=(lambda x: x * 354.7/50.65, lambda x: x * 50.65/354.7))
    # ax3.set_ylabel('Lightyield  [ADC/MeV]')
    #twin2.spines.right.set_position(("axes", 1.1))
    
    fig.tight_layout()
    ax1.legend()
    plt.savefig("lightyield_plot.pdf", format="pdf")
    plt.show()
    
    # Plot normalized_muon
    #plt.plot(normalized_muon, label='Normalized Muon', marker='o', linestyle='-', color='blue')

    # Plot normalized_co60
    #plt.plot(normalized_co60, label='Normalized Co60', marker='s', linestyle='--', color='green')

    # Plot LY_pe
    #plt.plot(normalized_LY_pe, label='LY_pe', marker='^', linestyle='-.', color='red')

    # Plot LY_ph
    #plt.plot(normalized_LY_ph, label='LY_ph', marker='d', linestyle=':', color='purple')

if __name__ == "__main__":
    main()
