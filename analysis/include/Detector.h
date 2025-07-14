#ifndef DETECTOR_H
#define DETECTOR_H

#include <Rtypes.h> // For Double_t
#include "Calibration.h"
#include "DetectorProperties.h"


struct EnergyDose {
    double dose = 0.0;
    double stddev = 0.0;
};

struct Position {
    double depth = 0.0;
    double stddev = 0.0;
};

struct Crystal {
    Position pos = {0,0};
    EnergyDose dose = {0,0};
    double size = 0;
};

class Detector {
public:
    Detector(DetectorProperties * detectorProperties);

    ~Detector() = default;

    TH1D* EnergyHist(const int channel);
    TH1D* CoincEnergyHist(const int channel);
    TH1D* StoppedEnergyHist(const int channel);
    TH1D* PhotonHist(const int channel);
    TH1D* CoincPhotonHist(const int channel);
    TH1D* TotalEnergyHist();
    TH1D* EntryHist(const int channel);
    TH1D* ExitHist(const int channel);
    TH1D* AngleHist(const int channel);
    TGraphErrors* MeansGraph();
    double CalcWET(double depth);
    double CalcWET(double z0, double delta, MaterialProperties mat);
    void Process();
    void CalcDose(int layer);
    void CalcPosition(int layer);
    void CalcWETConv();
    void FillMeansGraph();
    void Info();
    void CalcAbsorber();
    int nLayers = 0;
    int nSecLayers= 0;

    double WETConv = 0;
    double WETConvTeflon = 0; 
    double WETConvAlu = 0;
    double WETConvPMMA = 0;
   
    double Rw = 0;
    double alphaw = 0;
    double pw = 0;
    
    double Rm = 0;
    double alpham = 0;
    double pm = 0;
   
    double absorberConv = 0;
    double energy = 0;
    double layerSizeX= 0;
    double layerSizeY= 0;
    double layerSizeZ= 0;
    double gapSizeZ= 0;
    double teflonSizeZ= 0;
    double aluSizeZ= 0;
    double absorberSizeZ= 0;
    double SeclayerSizeZ= 0;
    std::string scintillator;
    const char* target;
    Calibration* calib=nullptr;

    DetectorProperties* detectorProperties;

    double alpha_m;
    double p_m;
    
    TH1D *h_total_edep;
    TGraphErrors *g_means;
    std::vector<TH1D*> h_edep_all;
    std::vector<TH1D*> h_edep_coinc;
    std::vector<TH1D*> h_stopped;
    std::vector<TH1D*> h_photons;
    std::vector<TH1D*> h_photons_coinc;
    std::vector<TH1D*> h_entry;
    std::vector<TH1D*> h_exit;
    std::vector<TH1D*> h_angle;

    std::vector<Crystal> crystals;

private:
    void PrintVector(const std::vector<double>& vec);
};

#endif