#ifndef PARTICLE_H
#define PARTICLE_H

#include <Rtypes.h> 
#include "Calibration.h"

class Particle {
public:
    int coincidence_time=0;
    int coincidence_layer = 0;
    int nbDetectors=0;
    Double_t total_edep = 0; 
    Double_t total_edep_err = 0; 
    Double_t coinc_layer = 0;
    
    
    std::vector<TraceProperties> traces = {};
    std::vector<Double_t> QuenchedDeposition = {};
    std::vector<Double_t> Deposition = {};
    std::vector<Double_t> DepositionStdDev = {};

    bool missingChannel = false;
    bool pileupStatus = false;
    bool ampOffsetStatus = false;
    
    Calibration* calibration=nullptr;
    //sim
    std::vector<std::vector<double>> dE={};
    std::vector<Double_t> NPhotons={};
    std::vector<std::vector<Double_t>> EntryPos = {}; 

    Particle(int nbDetectors, int coincidence_time, int coincidence_layer, Calibration* fcalibration);

    ~Particle() = default;
    void SumEDep();
    void InsertInitial(TraceProperties trace);
    void Insert(TraceProperties trace);
    void Coincidence(TraceProperties trace);
    bool Coincidence(int layer);
    bool CoincidencePhotons(int layer);
    void Test();
    void TestCoincLayer();
    void TestBuffer();
    void Clear();
    void SetdE(const int channel, double dEdep);
    void SetEDep(const int channel, double dEdep);
    void SetNPhotons(const int channel, double fNPhotons);
    double GetNPhotons(const int channel);
    void SetHitPosition(const int channel, double x, double y, double z);
    void CalculateEDep(const int channel);
    void ProcessEDep();

    std::vector<double> GetHitPosition(const int channel);
    std::vector<Double_t> GetTrace(int channel);
    Double_t GetEDep(int channel);
    Double_t GetTimeNS(int channel);
    Double_t GetTimePS(int channel);
    Double_t GetAmplitude(int channel);
    Double_t GetEDepDeviation(int channel);
};

#endif