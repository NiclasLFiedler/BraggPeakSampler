#include "Particle.h"
#include "TraceProperties.h"

void Particle::Clear() {
    total_edep = 0; 
    total_edep_err = 0; 
    coinc_layer = 0;
    amplitude = std::vector<Double_t>(15, 0.0);
    edep = std::vector<Double_t>(15, 0.0); 
    edep_err = std::vector<Double_t>(15, 0.0);
    time = std::vector<Double_t>(15, 0.0); //ns
    
    det_prepeak_step = false; // flag for prepeak step detected
    det_pileup = false; // flag for pileup detection
    missing_buffer = false; // flag for missing buffer entry
}

void Particle::sum_edep(){
    total_edep = 0;
    for(int i=0; i<edep.size(); i++){
        total_edep += GetEDep(i);
    }
}

Double_t Particle::GetEDep(int channel){
    // Double_t kB = 0.008694; //0.0333333
    return edep.at(channel); // /(1-kB*edep.at(channel)/5);
}

void Particle::coincidence(TraceProperties trace, int channel){
    if (time.at(0) != 0) {
        // cout << std::abs(time.at(0) - trace.time_ns) << endl;
        if(std::abs(time.at(0) - trace.time_ns) <= coincidence_time){
            time.at(channel) = trace.time_ns;
            edep.at(channel) = trace.energy;
            edep_err.at(channel) = trace.energy_err;
            amplitude.at(channel) = trace.amp;
            det_pileup = trace.det_pileup;
            if(trace.amp_idx>discard_index){
                det_prepeak_step = true;
            }
        }
        return;
    }
    return;
}

void Particle::TestCoincLayer(){
    for(int channel = edep.size()-1; channel>=0; channel--){
        if(edep.at(channel) != 0){
            coinc_layer = channel;
            return;
        }
    }
}
void Particle::TestBuffer(){
    bool emptyChannel = false;
    for(int channel = 0; channel<edep.size(); channel++){
        if(edep.at(channel) != 0){
            if(emptyChannel == false){
                continue;
            }
            else{
                missing_buffer = true;
            }
        }
        else{
            emptyChannel = true;
        }
    }
}