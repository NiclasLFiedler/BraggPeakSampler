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
    double R;
};

struct DetectorProperties {

    std::string scintillator;
    int nLayers;
    int nSecLayers;
    double energy;
    std::vector<double> crystalSize;
    std::vector<double> LayerSizeZ;
    double TeflonThickness;
    double AluThickness;
    std::string absorberType;
    double absorberSizeZ;
    double gapSizeZ;
    double SeclayerSizeZ;

    const char* Target;
    Calibration* calibration=nullptr;
    
    bool SimulationStatus;
    bool NormStatus;
    bool CalibStatus;
    bool ReversedStatus;
    bool AbsorberStatus;

    std::unordered_map<std::string, MaterialProperties> materials;
    
    DetectorProperties(){}
    ~DetectorProperties() = default;

    void Process(){
        double alpha_h2o = 2.585e-2;
        double p_h2o = 1.738;
        double kB_h2o = 1.738; //https://link.springer.com/article/10.1140/epjc/s10052-023-11242-2
        double R_h2o = CalcRange(alpha_h2o, p_h2o, energy);
        
        double alpha_pbwo4 = 7.275e-3;
        double p_pbwo4 = 1.690;
        double kB_pbwo4 = 0.01268; //https://arxiv.org/pdf/0911.3041
        double R_pbwo4 = CalcRange(alpha_pbwo4, p_pbwo4, energy);
        
        double alpha_dsb = 1.030e-2;
        double p_dsb = 1.713;
        double kB_dsb = 0.01268; //https://arxiv.org/pdf/0911.3041
        double R_dsb = CalcRange(alpha_dsb, p_dsb, energy);
    
        double alpha_ej256 = 2.393e-2;
        double p_ej256 = 1.742;
        double kB_ej256 = 0.154; //cite pöschl
        double R_ej256 = CalcRange(alpha_ej256, p_ej256, energy);
        
        double alpha_ej254 = 2.495e-2;
        double p_ej254 = 1.743;
        double kB_ej254 = 0.154;
        double R_ej254 = CalcRange(alpha_ej254, p_ej254, energy);
        
        double alpha_ej212 = 2.483e-2;
        double p_ej212 = 1.743;
        double kB_ej212 = 0.154;
        double R_ej212 = CalcRange(alpha_ej212, p_ej212, energy);
        
        double alpha_ej200 = 2.483e-2;
        double p_ej200 = 1.743;
        double kB_ej200 = 0.154;
        double R_ej200 = CalcRange(alpha_ej200, p_ej200, energy);
        
        double alpha_teflon = 1.459e-2;
        double p_teflon = 1.733;
        double kB_teflon = 0;
        double R_teflon = CalcRange(alpha_teflon, p_teflon, energy);
        
        double alpha_alu = 1.326e-2;
        double p_alu = 1.723;
        double kB_alu = 0;
        double R_alu = CalcRange(alpha_alu, p_alu, energy);
        
        double alpha_pmma = 2.173e-2;
        double p_pmma = 1.742;
        double kB_pmma = 0;
        double R_pmma = CalcRange(alpha_pmma, p_pmma, energy);
        
        double alpha_air = 2.456e1;
        double p_air = 1.736;
        double kB_air = 0;
        double R_air = CalcRange(alpha_air, p_air, energy);

        materials = {
            {"h2o", {alpha_h2o, p_h2o, kB_h2o, R_h2o}},
            {"pbwo4", {alpha_pbwo4, p_pbwo4, kB_pbwo4, R_pbwo4}},
            {"dsb", {alpha_dsb, p_dsb, kB_dsb, R_dsb}},
            {"ej256", {alpha_ej256, p_ej256, kB_ej256, R_ej256}},
            {"ej254", {alpha_ej254, p_ej254, kB_ej254, R_ej254}},
            {"ej212", {alpha_ej212, p_ej212, kB_ej212, R_ej212}},
            {"ej200", {alpha_ej200, p_ej200, kB_ej200, R_ej200}},
            {"teflon", {alpha_teflon, p_teflon, kB_teflon, R_teflon}},
            {"alu", {alpha_alu, p_alu, kB_alu, R_alu}},
            {"pmma", {alpha_pmma, p_pmma, kB_pmma, R_pmma}},
            {"air", {alpha_air, p_air, kB_air, R_air}}
        };

    }
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

    void SetTeflonThickness(double fTeflonThickness){
        TeflonThickness = fTeflonThickness;
        return;
    }

    void SetAluThickness(double fAluThickness){
        AluThickness = fAluThickness;
        return;
    }

    void SetAbsorberType(std::string fabsorberType){
        absorberType = fabsorberType;
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

    void SetAbsorberStatus(bool fAbsorberStatus){
        AbsorberStatus= fAbsorberStatus;
        return;
    }

    std::string GetScintillator(){return scintillator; }

    std::string GetAbsorberType(){return absorberType; }

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

    double GetTeflonThickness(){return TeflonThickness; }

    double GetAluThickness(){return AluThickness; }

    double GetAbsorberSize(){return absorberSizeZ; }

    int GetNSecondaryLayers(){return nSecLayers; }

    double GetSecondaryLayerSizeZ(){return SeclayerSizeZ; }

    const char*  GetTarget(){return Target; }

    Calibration* GetCalibration(){return calibration; }

    bool GetNormStatus(){return NormStatus; }

    bool GetCalibrationStatus(){return CalibStatus; }

    bool GetSimulationStatus(){return SimulationStatus; }

    bool GetReversedStatus(){return ReversedStatus; }

    bool GetAbsorberStatus(){return AbsorberStatus; }

    double GetBirksConstant(){return materials[GetScintillator()].kB; }

    double GetR(){
        return materials[GetScintillator()].R; 
    }

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

    double CalcRange(double alpha, double p, double energy) {
        double range = (alpha * std::pow(energy, p));
        return range;
    }

    MaterialProperties GetMaterial(std::string mat){
        return materials[mat];
    }
};

#endif
