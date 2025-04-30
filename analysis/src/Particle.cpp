#include "../include/Particle.h"
#include "../include/TraceProperties.h"

Particle::Particle(int fnbDetectors, int fcoincidence_time, int fcoincidence_layer, Calibration* fcalibration)
    : nbDetectors(fnbDetectors), coincidence_time(fcoincidence_time), coincidence_layer(fcoincidence_layer), calibration(fcalibration){
    
    Deposition = std::vector<Double_t>(nbDetectors, 0.0);
    DepositionStdDev = std::vector<Double_t>(nbDetectors, 0.0);
    traces = std::vector<TraceProperties>(nbDetectors);
    QuenchedDeposition = std::vector<Double_t>(nbDetectors, 0.0);
    dE = std::vector<std::vector<double>>(nbDetectors);
    NPhotons = std::vector<Double_t>(nbDetectors, 0.0);
    EntryPos = std::vector<std::vector<Double_t>>(nbDetectors, std::vector<Double_t>(3,0.0)); 

}

void Particle::Clear() {
    total_edep = 0; 
    total_edep_err = 0;
    coinc_layer = 0;
    pileupStatus = false;
    ampOffsetStatus = false;
    missingChannel = false;

    Deposition = std::vector<Double_t>(nbDetectors, 0.0);
    DepositionStdDev = std::vector<Double_t>(nbDetectors, 0.0);
    traces = std::vector<TraceProperties>(nbDetectors);
    QuenchedDeposition = std::vector<Double_t>(nbDetectors, 0.0);
    dE = std::vector<std::vector<double>>(nbDetectors);
    NPhotons = std::vector<Double_t>(nbDetectors, 0.0); 
    EntryPos = std::vector<std::vector<Double_t>>(nbDetectors, std::vector<Double_t>(3,0.0)); 
}

void Particle::SumEDep(){
    total_edep = 0;
    for(int i=0; i<traces.size(); i++){
        total_edep += GetEDep(i);
    }
}

void Particle::InsertInitial(TraceProperties trace){
    if(trace.channel != 0){
        std::cout << "Initial trace not channel 0" << std::endl;
    }
    Insert(trace);
    return;
}

void Particle::Insert(TraceProperties trace){
    int fchannel = trace.channel;
    traces.at(fchannel) = trace;
    QuenchedDeposition.at(fchannel) = calibration->GetQuenchedEnergy(fchannel, trace.amp);
    Deposition.at(fchannel) = calibration->GetEnergy(fchannel, QuenchedDeposition.at(fchannel));
    DepositionStdDev.at(fchannel) = calibration->GetQuenchedStdDev(fchannel, QuenchedDeposition.at(fchannel));
    return;
}

Double_t Particle::GetEDep(int channel){
    return Deposition.at(channel);
}

void Particle::SetEDep(const int channel, double dEdep){
    Deposition.at(channel) = dEdep;
    return;
}

Double_t Particle::GetTimeNS(int channel){
    return traces.at(channel).time_ns;
}

Double_t Particle::GetTimePS(int channel){
    return traces.at(channel).time_ps;
}

Double_t Particle::GetAmplitude(int channel){
    return traces.at(channel).amp;
}

std::vector<Double_t> Particle::GetTrace(int channel){
    return traces.at(channel).trace;
}

Double_t Particle::GetEDepDeviation(int channel){
    return DepositionStdDev.at(channel);
}

void Particle::Coincidence(TraceProperties trace){
    
    if(GetTimePS(trace.channel) != 0){
        std::cout << "Already filled channel " << trace.channel << ", duplicate events" << std::endl;
        return;
    }
    if(std::abs(GetTimePS(0) - trace.time_ps) <= coincidence_time){
        Insert(trace);
    }
    else{
        std::cout << "No coincidence found" << std::endl;
    }
    return;
}

void Particle::Test(){
    for(auto trace : traces){
        if(trace.pileupStatus) pileupStatus = true;
        if(trace.ampOffsetStatus) ampOffsetStatus = true;
    }
    TestBuffer();
    TestCoincLayer();
    return;
}

void Particle::TestCoincLayer(){
    for(int channel = traces.size()-1; channel>=0; channel--){
        if(GetEDep(channel) != 0){
            coinc_layer = channel;
            return;
        }
    }
}
void Particle::TestBuffer(){
    bool emptyChannel = false;
    for(int channel = 0; channel<traces.size(); channel++){
        if(GetEDep(channel) != 0){
            if(emptyChannel == false){
                continue;
            }
            else{
                missingChannel = true;
            }
        }
        else{
            emptyChannel = true;
        }
    }
}

// simulation specific
void Particle::SetdE(const int channel, double dEdep){
    dE.at(channel).push_back(dEdep);
    return;
}

void Particle::SetHitPosition(const int channel, double x, double y, double z){
    EntryPos.at(channel).at(0) = x;
    EntryPos.at(channel).at(1) = y;
    EntryPos.at(channel).at(2) = z;
    return;
}

std::vector<double> Particle::GetHitPosition(const int channel){
    return EntryPos.at(channel);
}

void Particle::SetNPhotons(const int channel, double fNPhotons){
    NPhotons.at(channel) = fNPhotons;
    return;
}

double Particle::GetNPhotons(const int channel){
    return NPhotons.at(channel);
}

void Particle::CalculateEDep(const int channel){
    double fedep = 0;
    for(const double& dE_i :dE.at(channel)){
        fedep += dE_i;
    }
    if(fedep != 0){
        SetEDep(channel, fedep);
    }
    return;
}

void Particle::ProcessEDep(){
    for(int ch=0; ch<traces.size(); ch++){            
        CalculateEDep(ch);
    }
    SumEDep();
    return;
}

bool Particle::Coincidence(int layer){
    if(layer <= coincidence_layer){
        layer = coincidence_layer;
    }
    for(int i = 0; i < layer; i++){
        if(GetEDep(i) == 0){
            return false;
        }
    }
    return true;
}

bool Particle::CoincidencePhotons(int layer){
    if(layer <= coincidence_layer){
        layer = coincidence_layer;
    }
    for(int i = 0; i < layer; i++){
        if(NPhotons.at(i) == 0){
            return false;
        }
    }
    return true;
}