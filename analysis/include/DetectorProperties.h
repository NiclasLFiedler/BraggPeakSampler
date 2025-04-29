#ifndef DETECTORPROPERTIES_H
#define DETECTORPROPERTIES_H

#include <vector>
#include <unordered_map>
#include <string>

class Calibration;  // Tells the compiler that Calibration exists

struct MaterialProperties {
    double alpha;
    double p;
    double kB;
};

struct DetectorProperties {

    std::string scintillator;
    int nLayers;
    int nSecLayers;
    double energy;
    std::vector<double> crystalSize;
    std::vector<double> TeflonThickness;
    std::vector<double> LayerSizeZ;
    std::vector<double> AluThickness;
    double absorberSizeZ;
    double gapSizeZ;
    double SeclayerSizeZ;

    const char* Target;
    Calibration* calibration=nullptr;
    
    bool SimulationStatus;
    bool NormStatus;
    bool CalibStatus;
    bool ReversedStatus;

    std::unordered_map<std::string, MaterialProperties> materials;

    double alpha_h2o = 2.585e-2;
    double p_h2o = 1.738;
    double kB_h2o = 1.738; //https://link.springer.com/article/10.1140/epjc/s10052-023-11242-2

    double alpha_pbwo4 = 7.275e-3;
    double p_pbwo4 = 1.690;
    double kB_pbwo4 = 0.01268; //https://arxiv.org/pdf/0911.3041

    double alpha_dsb = 1.030e-2;
    double p_dsb = 1.713;
    double kB_dsb = 0.01268; //https://arxiv.org/pdf/0911.3041

    double alpha_ej256 = 2.393e-2;
    double p_ej256 = 1.742;
    double kB_ej256 = 0.154; //cite pöschl
    
    double alpha_ej254 = 2.495e-2;
    double p_ej254 = 1.743;
    double kB_ej254 = 0.154;

    double alpha_ej212 = 2.483e-2;
    double p_ej212 = 1.743;
    double kB_ej212 = 0.154;

    double alpha_teflon = 1.459e-2;
    double p_teflon = 1.733;
    double kB_teflon = 0;

    double alpha_alu = 1.326-2;
    double p_alu = 1.723;
    double kB_alu = 0;

    double alpha_pmma = 2.173e-2;
    double p_pmma = 1.742;
    double kB_pmma = 0;

    double alpha_air = 2.456e1;
    double p_air = 1.736;
    double kB_air = 0;

    DetectorProperties(){
        materials = {
            {"h2o", {alpha_h2o, p_h2o, kB_h2o}},
            {"pbwo4", {alpha_pbwo4, p_pbwo4, kB_pbwo4}},
            {"dsb", {alpha_dsb, p_dsb, kB_dsb}},
            {"ej256", {alpha_ej256, p_ej256, kB_ej256}},
            {"ej254", {alpha_ej254, p_ej254, kB_ej254}},
            {"ej212", {alpha_ej212, p_ej212, kB_ej212}},
            {"teflon", {alpha_teflon, p_teflon, kB_teflon}},
            {"alu", {alpha_alu, p_alu, kB_alu}},
            {"pmma", {alpha_pmma, p_pmma, kB_pmma}}
        };

    }
    ~DetectorProperties() = default;

    void SetScintillator(std::string fscintillator){
        scintillator = fscintillator;
        return;
    }

    void SetBeamEnergy(double fenergy){
        energy = fenergy;
        return;
    }

    void SetNLayers(int fnLayers){
        nLayers = fnLayers;
        return;
    }

    void SetScintillatorDimensions(std::vector<double> fcrystalSize){
        crystalSize = fcrystalSize;
        return;
    }

    void SetGapSizeZ(double fgapSizeZ){
        gapSizeZ = fgapSizeZ;
        return;
    }

    void SetTeflonThickness(std::vector<double> fTeflonThickness){
        if (fTeflonThickness.size() < nLayers) {
            fTeflonThickness.resize(nLayers, 0);
        }
        TeflonThickness = fTeflonThickness;
        return;
    }

    void SetAluThickness(std::vector<double> fAluThickness){
        if (fAluThickness.size() < nLayers) {
            fAluThickness.resize(nLayers, 0);
        }
        AluThickness = fAluThickness;
        return;
    }

    void SetAbsorberSize(double fabsorberSizeZ){
        absorberSizeZ = fabsorberSizeZ;
        return;
    }

    void SetNSecondaryLayers(int fnSecLayers){
        nSecLayers = fnSecLayers;
        return;
    }

    void SetSecondaryLayerSizeZ(double fSeclayerSizeZ){
        SeclayerSizeZ = fSeclayerSizeZ;
        return;
    }

    void SetTarget(const char* fTarget){
        Target = fTarget;
        return;
    }

    void SetCalibration(Calibration* fcalibration){
        calibration = fcalibration;
        return;
    }

    void SetNormStatus(bool fNormStatus){
        NormStatus = fNormStatus;
        return;
    }

    void SetCalibrationStatus(bool fCalibStatus){
        CalibStatus = fCalibStatus;
        return;
    }

    void SetSimulationStatus(bool fSimulationStatus){
        SimulationStatus = fSimulationStatus;
        return;
    }

    void SetReversedStatus(bool fReversedStatus){
        ReversedStatus= fReversedStatus;
        return;
    }

    std::string GetScintillator(){return scintillator; }

    double GetBeamEnergy(){return energy; }

    int GetNLayers(){return nLayers; }

    std::vector<double> GetScintillatorDimensions(){return crystalSize; }

    double GetLayerSizeX(){return crystalSize.at(0); }
    double GetLayerSizeY(){return crystalSize.at(1); }
    double GetLayerSizeZ(){return crystalSize.at(2); }

    double GetLayerSizeZ(int channel){
        if(channel < 0){
            return 0;
        }
        else if(channel < GetNSecondaryLayers()){
            return GetSecondaryLayerSizeZ();
        }
        else{
            return crystalSize.at(2);
        }
    }

    double GetGapSizeZ(){return gapSizeZ; }

    std::vector<double> GetTeflonThickness(){return TeflonThickness; }
    double GetTeflonThickness(int layer){return TeflonThickness.at(layer)*1/10000; }

    std::vector<double> GetAluThickness(){return AluThickness; }
    double GetAluThickness(int layer){return AluThickness.at(layer)*1/10000; }

    double GetAbsorberSize(){return absorberSizeZ; }

    int GetNSecondaryLayers(){return nSecLayers; }

    double GetSecondaryLayerSizeZ(){return SeclayerSizeZ; }

    const char*  GetTarget(){return Target; }

    Calibration* GetCalibration(){return calibration; }

    bool GetNormStatus(){return NormStatus; }

    bool GetCalibrationStatus(){return CalibStatus; }

    bool GetSimulationStatus(){return SimulationStatus; }

    bool GetReversedStatus(){return ReversedStatus; }

    double GetBirksConstant(){return materials[GetScintillator()].kB; }

    double GetP(){
        return materials[GetScintillator()].p; 
    }

    double GetP(std::string mat){
        return materials[mat].p; 
    }

    double GetAlpha(){
        return materials[GetScintillator()].alpha; 
    }

    double GetAlpha(std::string mat){
        return materials[mat].alpha; 
    }
};

#endif
