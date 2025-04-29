#ifndef PARTICLE_H
#define PARTICLE_H

#include <Rtypes.h> 

class Particle {
public:
    Double_t coincidence_time = 20;
    Double_t discard_index = 50;
    
    Double_t total_edep = 0; 
    Double_t total_edep_err = 0; 
    Double_t coinc_layer = 0;
    std::vector<Double_t> amplitude = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> edep = std::vector<Double_t>(15, 0.0); 
    std::vector<Double_t> edep_err = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> time = std::vector<Double_t>(15, 0.0); //ns

    bool det_prepeak_step = false; // flag for prepeak step detected
    bool det_pileup = false; // flag for pileup detection
    bool missing_buffer = false; // flag for missing buffer entry

    // Member functions
    Double_t GetEDep(int channel);
    void sum_edep();
    void coincidence(TraceProperties trace, int channel);
    void TestCoincLayer();
    void TestBuffer();
    void Clear();
};

#endif