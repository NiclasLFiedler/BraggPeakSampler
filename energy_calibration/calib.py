ch_na = [187.95, 182.31, 202.14, 185.35, 218.14, 175.91, 178.96, 187.66, 168.23, 173.47, 190.77, 116.02, 121.70, 114.04, 90.49]
ch_muon = [652.10, 633.39, 733.54, 667.98, 741.10, 712.55, 696.39, 715.85, 697.87, 725.30, 712.13, 387.25, 408.65, 404.82, 392.22]
e_na = 1.2745
e_muon = 4.95

def get_linear_interpolation(ch):
    calib = []
    calib.append(abs((e_na-e_muon)/(ch_na[ch]-ch_muon[ch])))
    calib.append(e_na-abs((e_na-e_muon)/(ch_na[ch]-ch_muon[ch]))*ch_na[ch])
    print(f"{calib[1]},", end="")

for i in range(15):
    get_linear_interpolation(i)