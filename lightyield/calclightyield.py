import matplotlib.pyplot as plt
import numpy as np

def normalize(data):
    # Check if the first element is zero to avoid division by zero
    if data[0] == 0:
        raise ValueError("The first element of the list is zero, cannot normalize.")
    return [value / data[0] for value in data]

def convert_to_secondary_units(y):
    return y*50.65/11.73  # Example: converting y to another unit system

def main():    
    LY = [np.float64(442.18853875881183), np.float64(412.7240577897924), np.float64(438.36596726125043), np.float64(445.639697109856), np.float64(439.1334908732424), np.float64(439.1418069864141), np.float64(444.9552486293444), np.float64(431.3836981207205), np.float64(416.2072835187019), np.float64(430.3227828625882), np.float64(412.61192965820135), np.float64(450.87171511477146), np.float64(445.3323072324525), np.float64(433.221253619175), np.float64(444.1094086965365), np.float64(305.5188579947859), np.float64(289.5574462218171), np.float64(303.58029062085296), np.float64(293.0680543882572), np.float64(346.1524491452491), np.float64(353.8045692495503), np.float64(354.66877605055777), np.float64(343.7406251010896), np.float64(285.9187648093916), np.float64(344.59569533221736), np.float64(330.34242904149585), np.float64(301.83627459192377), np.float64(304.20053759205473), np.float64(382.14857044549683), np.float64(344.13649949009846), np.float64(391.4047214439339), np.float64(404.5485256861425), np.float64(414.6991612973868), np.float64(393.7214815345912), np.float64(399.27750939720016), np.float64(398.9512033989302)]

    LY_err = [np.float64(132.1816042675276), np.float64(129.73259992086594), np.float64(131.11876731678927), np.float64(132.43579645862093), np.float64(130.59866715814664), np.float64(132.15449610782647), np.float64(132.6555697697036), np.float64(131.0620965528852), np.float64(128.74867957532263), np.float64(129.7249706476238), np.float64(129.3683416038428), np.float64(130.61875915017876), np.float64(131.97778685142865), np.float64(131.0589564180455), np.float64(133.69902525458764), np.float64(125.91175444785588), np.float64(125.66714018538698), np.float64(124.86467478981392), np.float64(126.90801213861991), np.float64(128.86630358015583), np.float64(129.42676510647516), np.float64(132.7940884625936), np.float64(127.87470686493678), np.float64(121.6867900395244), np.float64(126.62112363798262), np.float64(128.84257735142978), np.float64(123.79852102350259), np.float64(123.30387482049353), np.float64(133.12130875474531), np.float64(128.01592205417637), np.float64(134.07027727293803), np.float64(135.77631854384015), np.float64(136.30171604988905), np.float64(135.71645714433086), np.float64(137.49111669028878), np.float64(134.99554275900556)]

    ped = 85.4
    pe = 110.53
    na22 = 0.511
    LY_pe = []
    LY_pe_err = []
    LY_ph = []
    LY_ph_err = []
    ly_average = 0
    ly_average_err = 0
    new = 0
    for index in range(36):
        LY_pe.append(((LY[index]-ped)/(pe-ped))/na22)
        LY_pe_err.append((LY_err[index]/(pe-ped))/na22)
        LY_ph.append((((LY[index]-ped)/(pe-ped))/na22)/0.2316)
        LY_ph_err.append(((LY_err[index]/(pe-ped))/na22)/0.2316)
        print(f"CHold {index} CH new {new} Light yield: {LY[index]} ± {LY_err[index]:.2f} LY_pe: {LY_pe[-1]:.2f} ± {LY_pe_err[-1]:.2f}  LY_ph: {LY_ph[-1]:.2f} ± {LY_ph_err[-1]:.2f}")
        ly_average = ly_average + LY_ph[-1]
        ly_average_err = ly_average_err + LY_ph_err[-1]
        new = new + 1
    print(f"\nAverage Light yield: {ly_average/36}")
    print(f"Average Light yield Error: {ly_average_err/36}")

    # Create the primary plot
    plt.rcParams.update({'font.size': 14})
    fig, ax1 = plt.subplots(figsize=(10, 6))
    color = 'tab:blue'
    ax1.set_xlabel('Crystal Index')
    ax1.set_ylabel('Light yield  [ph/MeV]')
    ax1.errorbar(range(36), LY_ph, yerr=LY_ph_err, fmt='o', capsize=5, elinewidth=1, capthick=1, color="black", label='Light yield data points')
    ax1.tick_params(axis='y')
    ax1.grid()

    ax1.set_ylim(0,180)

    twin1 = ax1.twinx()
    twin1.set_ylim(0, 180*11.73/50.65)
    twin1.set_ylabel('Light yield  [p.e./MeV]')

    fig.tight_layout()
    ax1.legend()
    plt.savefig("lightyield_plot.pdf", format="pdf")
    plt.show()
    
if __name__ == "__main__":
    main()
