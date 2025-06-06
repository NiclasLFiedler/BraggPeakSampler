import numpy as np
import matplotlib.pyplot as plt

mat = {
    "pwo": {
        "name": "Lead Tungstate",
        "formula": "PbWO4",
        "density": 8.28,         # g/cm^3
        "I": 600.7e-6,                # eV
        "Z/A": 188/455.04,             # g/mol (approx)
    },
    "h2o": {
        "name": "Water",
        "formula": "H2O",
        "density": 1.0,
        "I": 75e-6,
        "Z/A": 10/18.01528,            
    },
    "pmma": {
        "name": "PMMA",
        "formula": "C5H8O2",
        "density": 1.19,
        "I": 74e-6,
        "Z/A": 54/100.18,             
    },
    "air": {
        "name": "AIR",
        "formula": "AIR",
        "density": 1.205e-3,
        "I": 85.7e-6,                # eV
        "Z/A": 0.49919,             
    },
    "lung_def": {
        "name": "Lung",
        "formula": "Lung",
        "density": 1.05,              # g/cm^3
        "I": 75.3e-6,                # eV
        "Z/A": 0.5496,             
    },
}



comp = {
    "air": {
        "C": [0.000124, 0.4995],
        "N": [0.755267, 0.4998],
        "O": [0.231781, 0.5000],
        "Ar": [0.012827, 0.4506]
    },
    "lung": {
        "H": [0.101278, 0.9923],
        "C": [0.102310, 0.4995],
        "N": [0.02865, 0.4998],
        "O": [0.757072, 0.5000],
        "Na": [0.001840, 0.4785],
        "Mg": [0.000730, 0.4937],
        "P": [0.0008, 0.4843],
        "S": [0.002250, 0.4990],
        "Cl": [0.002660, 0.4795],
        "K": [0.001940, 0.4860],
        "Ca": [0.000090, 0.4990],
        "Fe": [0.000370, 0.4656],
        "Zn": [0.000010, 0.4589]
    }
}

def SP(mat):
    return mat["density"] * (mat["Z/A"]) * np.log(2 * 0.511 / mat["I"])

def mSP(mat):
    return SP(mat) * 1/mat["density"]

def SPRatio(mat1, mat2):
    return SP(mat1) / SP(mat2)

def mSPRatio(mat1, mat2):
    return mSP(mat1) / mSP(mat2)

def SPRatioComb1(mat1, mat2, mat3, w, rho1):
    return SPComp(mat1, mat2, w, rho1) / mSP(mat3)

def mSPRatioComb1(mat1, mat2, mat3, w):
    return mSPComp(mat1, mat2, w) / mSP(mat3)

def mSPComp(mat1, mat2, w):
    return mSP(mat1) * w + mSP(mat2) * (1 - w)

def SPComp(mat1, mat2, w, rho):
    return rho*mSPComp(mat1, mat2, w)

def Pmod(lungdef, air, material, water, d):
    rhoL = 0.26
    mSL_SA = mSPRatioComb1(lungdef, air, air, rhoL/lungdef["density"])
    mSM_SA = mSPRatio(material, air)
    SM_SW = SPRatio(material, water)
    SM_SL = SPRatioComb1(lungdef, air, air, rhoL/lungdef["density"], rhoL)**(-1)
    
    return d*(mSL_SA-1)*(mSM_SA-mSL_SA)/(mSM_SA-1)**2*SM_SW*SM_SL

def dConstant(lungdef, air, material, water, Pmod):
    rhoL = 0.26
    
    mSL_SA = mSPRatioComb1(lungdef, air, air, rhoL/lungdef["density"])

    mSM_SA = mSPRatio(material, air)
    SM_SW = SPRatio(material, water)
    SM_SL = SPRatioComb1(lungdef, air, material, rhoL/lungdef["density"], rhoL)**(-1)

    # print("mSAir: ", mSP(air))
    # print("mSwater: ", mSP(water))
    # print("mSmaterial: ", mSP(material))

    # print("mSL_SA:", mSL_SA)
    # print("mSM_SA:", mSM_SA)
    # print("SM_SW:", SM_SW)
    # print("SM_SL:", SM_SL)  
    # print("Pmod:", Pmod)
    # print("(mSL_SA-1):", mSL_SA-1)
    # print("(mSM_SA-mSL_SA):", mSM_SA-mSL_SA)
    # print("(mSM_SA-1):", mSM_SA-1)
    print("w_m: ", (mSL_SA-1)/(mSM_SA-1))
    print("d:", Pmod/((mSL_SA-1)*(mSM_SA-mSL_SA)/(mSM_SA-1)**2*SM_SW*SM_SL))

    return ((mSL_SA-1)*(mSM_SA-mSL_SA)/(mSM_SA-1)**2*SM_SW*SM_SL)

def sumZA(comp):
    ZA = 0
    for c in comp:
        ZA += comp[c][1] * comp[c][0]
    return ZA

print("pwo/water", mSPRatio(mat["lung_def"], mat["air"]))
print("air/water", mSPRatio(mat["air"], mat["h2o"]))

print("lung/water", SPRatioComb1(mat["lung_def"], mat["air"], mat["h2o"], 0.26/mat["lung_def"]["density"], 0.26))

print(dConstant(mat["lung_def"], mat["air"], mat["h2o"], mat["h2o"], 300))

print((1-0.2476)*1.05)