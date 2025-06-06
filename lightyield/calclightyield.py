import matplotlib.pyplot as plt
import numpy as np

def normalize(data):
    # Check if the first element is zero to avoid division by zero
    if data[0] == 0:
        raise ValueError("The first element of the list is zero, cannot normalize.")
    return [value / data[0] for value in data]

def convert_to_secondary_units(y):
    return y*50.65/11.73  # Example: converting y to another unit system

ped = 85.4
pe = 110.53
# na22 = 0.059409
na22 = 0.511    
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
    print(np.sum(std_devs**2))
    test = 0
    for i in std_devs:
        test = test + i**2
    print(test)
    print(np.sum(std_devs**2))
    print(np.sqrt(np.sum(std_devs**2)))
    print(len(values))
    std_of_mean = np.sqrt(np.sum(std_devs**2)) / len(values)
    
    return mean, std_of_mean

def main():    
    #old flat
    LY_old = [np.float64(442.18853875881183), np.float64(412.7240577897924), np.float64(438.36596726125043), np.float64(445.639697109856), np.float64(439.1334908732424), np.float64(439.1418069864141), np.float64(444.9552486293444), np.float64(431.3836981207205), np.float64(416.2072835187019), np.float64(430.3227828625882), np.float64(412.61192965820135), np.float64(450.87171511477146), np.float64(445.3323072324525), np.float64(433.221253619175), np.float64(444.1094086965365), np.float64(305.5188579947859), np.float64(289.5574462218171), np.float64(303.58029062085296), np.float64(293.0680543882572), np.float64(346.1524491452491), np.float64(353.8045692495503), np.float64(354.66877605055777), np.float64(343.7406251010896), np.float64(285.9187648093916), np.float64(344.59569533221736), np.float64(330.34242904149585), np.float64(301.83627459192377), np.float64(304.20053759205473), np.float64(382.14857044549683), np.float64(344.13649949009846), np.float64(391.4047214439339), np.float64(404.5485256861425), np.float64(414.6991612973868), np.float64(393.7214815345912), np.float64(399.27750939720016), np.float64(398.9512033989302)]
    LY_err_old = [np.float64(132.1816042675276), np.float64(129.73259992086594), np.float64(131.11876731678927), np.float64(132.43579645862093), np.float64(130.59866715814664), np.float64(132.15449610782647), np.float64(132.6555697697036), np.float64(131.0620965528852), np.float64(128.74867957532263), np.float64(129.7249706476238), np.float64(129.3683416038428), np.float64(130.61875915017876), np.float64(131.97778685142865), np.float64(131.0589564180455), np.float64(133.69902525458764), np.float64(125.91175444785588), np.float64(125.66714018538698), np.float64(124.86467478981392), np.float64(126.90801213861991), np.float64(128.86630358015583), np.float64(129.42676510647516), np.float64(132.7940884625936), np.float64(127.87470686493678), np.float64(121.6867900395244), np.float64(126.62112363798262), np.float64(128.84257735142978), np.float64(123.79852102350259), np.float64(123.30387482049353), np.float64(133.12130875474531), np.float64(128.01592205417637), np.float64(134.07027727293803), np.float64(135.77631854384015), np.float64(136.30171604988905), np.float64(135.71645714433086), np.float64(137.49111669028878), np.float64(134.99554275900556)]
    
    #flat
    LY = [580.8736157052743, 558.5195832257807, 571.796531005362, 576.9528373856286, 581.7088213974404, 586.2622705034566, 573.5793128844257, 565.3303270928747, 554.7419921907729, 564.7257886906036, 538.7889547317318, 588.0257534263993, 584.7716124489042, 573.7353661896567, 585.5823262525396, 417.0285768032109, 402.955441148207, 425.66970208316394, 407.06210026235993, 478.7078721131375, 480.24524947929166, 498.7221735173861, 479.11535774632296, 412.3668762371207, 467.5596796449048, 470.55636778080645, 414.2348792379325, 421.9478319974132, 530.3742362501756, 462.0656808184509, 537.5481311803621, 522.608248638223, 526.8485577417683, 522.1902140873837, 528.1508335227278, 569.2556589126506]
    LY_err = [152.41627883117627, 146.19071705929454, 145.9942042895179, 151.95897799607252, 151.1275107286317, 152.59661217075373, 145.64515011341533, 146.35429211666116, 143.12000232187242, 145.32571479615112, 143.74734084181577, 148.4629563655851, 145.95933519107834, 143.59395337382034, 148.32415844458603, 137.18618862882232, 133.84921507450196, 137.9784795402983, 138.0449497850902, 144.4276871002791, 141.60673651274504, 142.54986144632, 142.3778908671625, 133.15368102894138, 139.92401652362534, 143.3781761502344, 136.39231347801788, 141.45312621528737, 148.144104072922, 140.5266152100586, 150.71895214029283, 147.03302339313316, 146.19217594025977, 150.61858583846637, 150.19702944041194, 156.5074285469471]

    #lateral
    LY_lat = [457.3293200657901, 473.1145670080971, 460.2227469560361, 451.4789680411492, 479.0370633469254, 484.41893790520476, 457.58561819358823, 452.27247185213486, 467.0652346924665, 486.1498828479076, 457.9884315099127, 479.81885244448364, 485.6655351647721, 495.98592633091937, 467.39202022687977, 313.08219073193203, 313.3490354833067, 301.2406239993116, 301.9498689381933, 352.4771095133145, 375.95712819446646, 396.59781031823076, 381.18803950738703, 325.10196567755065, 380.48697262401015, 347.46567514220203, 316.43229799388115, 312.16901475677275, 417.9609414962273, 360.0856903237679, 430.4278511659892, 420.0990384588279, 412.9547851527781, 391.5824987128114, 397.28296194386303, 436.51222800479803, 309.2099430572384, 262.9420838345245]
    LY_lat_err = [135.15219777138492, 137.9848543055064, 124.73222735843882, 123.77041412695226, 131.36498055446114, 133.29327635051513, 133.29125319115678, 126.38899739123876, 125.16996504359446, 130.65907950528802, 125.17013060330986, 131.78206645838952, 128.72366486072883, 133.36027504155135, 127.38872636454161, 103.5534501088127, 102.8606319617789, 98.87466583231763, 96.62153505146449, 113.16439535989089, 120.34496590629885, 127.2955910933667, 121.72637157473348, 107.71730027027988, 123.33702875892618, 113.69521031763921, 106.64027555610934, 94.92042698482652, 126.93579870849513, 125.62193365805214, 130.62142571901975, 128.79891881010582, 126.00432630801515, 125.09478683454286, 127.32259101682668, 132.75867345516508, 98.41133605624921, 57.16752058574868]
    
    LY_ej = [3065.8680562948466, 2111.5204682089284, 2195.380253588704, 947.3111612280909, 973.588175501902]
    LY_ej_err = [445.2448683651481, 333.403234965277, 367.7696271726769, 212.18734149175012, 241.85514481819186]
    
    LY_pe = []
    LY_pe_err = []
    LY_ph = []
    LY_ph_err = []
    
    LY_pe_lat = []
    LY_pe_err_lat = []
    LY_ph_lat = []
    LY_ph_err_lat = []
    
    ly_average = 0
    ly_average_err = 0
    ly_average_lat = 0
    ly_average_err_lat = 0

    # for index in range(5):
    #     LY_pe.append(GetPE(LY_ej[index]))
    #     LY_pe_err.append(GetPEDev(LY_ej_err[index]))
    #     LY_ph.append(GetPH(LY_ej[index]))
    #     LY_ph_err.append(GetPHDev(LY_ej_err[index]))
        
    #     print(f"Flat: CH {index} Light yield: {LY[index]} ± {LY_err[index]:.2f} LY_pe: {LY_pe[-1]:.2f} ± {LY_pe_err[-1]:.2f}  LY_ph: {LY_ph[-1]:.2f} ± {LY_ph_err[-1]:.2f}")

    # exit()

    for index in range(36):
        LY_pe.append(GetPE(LY[index]))
        LY_pe_err.append(GetPEDev(LY_err[index]))
        LY_ph.append(GetPH(LY[index]))
        LY_ph_err.append(GetPHDev(LY_err[index]))
        
        print(f"Flat: CH {index} Light yield: {LY[index]} ± {LY_err[index]:.2f} LY_pe: {LY_pe[-1]:.2f} ± {LY_pe_err[-1]:.2f}  LY_ph: {LY_ph[-1]:.2f} ± {LY_ph_err[-1]:.2f}")
        ly_average = ly_average + LY_ph[-1]
        ly_average_err = ly_average_err + LY_ph_err[-1]

    for index in range(38):
        LY_pe_lat.append(GetPE(LY_lat[index]))
        LY_ph_lat.append(GetPH(LY_lat[index]))
        LY_pe_err_lat.append(GetPEDev(LY_lat_err[index]))
        LY_ph_err_lat.append(GetPHDev(LY_lat_err[index]))
        
        print(f"Lateral: CH {index} Light yield: {LY_lat[index]} ± {LY_lat_err[index]:.2f} LY_pe: {LY_pe_lat[-1]:.2f} ± {LY_pe_err_lat[-1]:.2f}  LY_ph: {LY_ph_lat[-1]:.2f} ± {LY_ph_err_lat[-1]:.2f}")
        ly_average_lat = ly_average_lat + LY_ph_lat[-1]
        ly_average_err_lat = ly_average_err_lat + LY_ph_err_lat[-1]


    print(f"\nAverage flat light yield: {ly_average/36}")
    print(f"Average flat light yield Error: {ly_average_err/36}")
    print()
    print(f"\nAverage lateral light yield: {ly_average_lat/38}")
    print(f"Average lateral light yield Error: {ly_average_err_lat/38}")
    # Create the primary plot
    plt.rcParams.update({'font.size': 14})
    fig, ax1 = plt.subplots(figsize=(15, 9))
    ax1.set_title('Measured light yield of PbWO4 crystals 3mm and 2mm thick using a Na22 source at 20°C')
    ax1.set_xlabel('Crystal Index')
    ax1.set_ylabel('Light yield / ph/MeV')

    xaxis = [i for i in range(1)] + [i for i in range(2, 23)] + [i for i in range(24, 38)]

    ax1.errorbar(xaxis, LY_ph, yerr=LY_ph_err, fmt='o', capsize=5, elinewidth=2, capthick=2, markersize=8, color="orange", label='Flat position')
    ax1.errorbar(xaxis, LY_ph_lat[:36], yerr=LY_ph_err_lat[:36], fmt='o', capsize=5, elinewidth=2, capthick=2, markersize=8, color="green", label='Vertical position')
    ax1.errorbar([1,23], [LY_ph_lat[36], LY_ph_lat[37]], yerr=[LY_ph_err_lat[36], LY_ph_err_lat[37]], fmt='o', capsize=5, elinewidth=2, capthick=2, markersize=8, color="black", label='Vertical position & window cut-out')

    ax1.grid()
    ax1.tick_params(axis='y')
    print()
    print()
    print(LY_ph_err[:15])
    print()
    print()
    average_first_15, average_first_15err = mean_and_std_of_mean(LY_ph[:15], LY_ph_err[:15])
    print()
    print()
    average_last_21, average_last_21err = mean_and_std_of_mean(LY_ph[15:36], LY_ph_err[15:36])
    average_first_15_lat, average_first_15_laterr = mean_and_std_of_mean(LY_ph_lat[:15], LY_ph_err_lat[:15])
    average_last_21_lat, average_last_21_latererr = mean_and_std_of_mean(LY_ph_lat[15:36], LY_ph_err_lat[15:36])
    
    ax1.plot(range(16), [average_first_15] * 16, color='dodgerblue', linestyle='--')
    ax1.plot(range(16, 38), [average_last_21] * 22, color='dodgerblue', linestyle='--')

    ax1.text(-1, average_first_15-3, f'{average_first_15:.2f} \n' f'$\\pm${average_first_15err:.2f}', color='dodgerblue', fontsize=16, va='bottom', ha='center')
    ax1.text(38, average_last_21-3, f'{average_last_21:.2f}\n' f'$\\pm${average_last_21err:.2f}', color='dodgerblue', fontsize=16, va='bottom', ha='center')

    ax1.plot(range(16), [average_first_15_lat] * 16, color='crimson', linestyle='--')
    ax1.plot(range(16, 38), [average_last_21_lat] * 22, color='crimson', linestyle='--')

    ax1.text(-1, average_first_15_lat-3, f'{average_first_15_lat:.2f}\n' f'$\\pm${average_first_15_laterr:.2f}', color='crimson', fontsize=16,va='bottom', ha='center')
    ax1.text(38, average_last_21_lat-3, f'{average_last_21_lat:.2f}\n' f'$\\pm${average_last_21_latererr:.2f}', color='crimson', fontsize=16, va='bottom', ha='center')

    axX = 240
    axY = 0
    ax1.set_ylim(axY, axX)

    twin1 = ax1.twinx()
    twin1.set_ylim(axY * 11.73 / 50.65, axX * 11.73 / 50.65)
    twin1.set_ylabel('Light yield / p.e./MeV')

    # Define custom x-tick labels
    xtick_labels = [f'{0}'] + [f'window {0}'] + [f'{i+1}' for i in range(0,21)] +[f'window {21}'] +[f'{i+1}' for i in range(21,35)]
    print(xtick_labels)
    ax1.set_xticks(range(38))  # Set x-ticks for all 38 entries
    ax1.set_xticklabels(xtick_labels, rotation=90, ha='right')  # Rotate labels for better readability

    fig.tight_layout()
    ax1.legend()
    plt.savefig("lightyield_plot.pdf", format="pdf")
    plt.show()
    
if __name__ == "__main__":
    main()
