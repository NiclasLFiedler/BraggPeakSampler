#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"
#include "src/TraceProperties.cpp"
#include "src/Particle.cpp"
#include "src/Plotter.cpp"
#include "src/Detector.cpp"
#include "include/EnergyCh.h"
#include "include/Calibration.h"
#include "include/DetectorProperties.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

void HistTraces(){
    ROOT::EnableImplicitMT();

    std::ifstream configFile("config.json");
    if (!configFile) {
        std::cerr << "Error opening config file!" << std::endl;
        return;
    }

    json allConfigs;
    configFile >> allConfigs;

    int fileSelect                      = allConfigs["targetSelect"];
    
    const auto& config = allConfigs["detectors"][int(allConfigs["detectorSelect"])];
    
    int datasetSelect                   = config["datasetSelect"];
    std::string detectortype            = config["detectorType"];
    std::string absorberType            = config["absorberType"];
    double beamEnergy                   = config["beamEnergy"];
    int nLayers                         = config["nLayers"];
    std::vector<double> crystalSize     = config["crystalSize"];
    double gapSizeZ                     = config["gapSizeZ"];
    bool secondaryLayerStatus           = config["secondaryLayerStatus"];
    int nSecLayers                      = config["nSecondaryLayers"];
    double secLayerSizeZ                = config["secLayerSizeZ"];
    bool absorberStatus                 = config["absorberStatus"];
    double absorberSize                 = config["absorberSize"];
    bool reversedStatus                 = config["reversedStatus"];
    bool normStatus                     = config["normStatus"];
    bool simulationStatus               = config["simulationStatus"];
    std::vector<double> teflonThickness = config["teflonThickness"];
    std::vector<double> aluThickness    = config["aluThickness"];
    int coincidenceTime                 = config["coincidenceTime"];
    int coincidenceLayer                = config["coincidenceLayer"];
    int discard_index                   = config["discardIndex"];
    int pmod                            = allConfigs["pmod"];
    int heteroThickness                 = allConfigs["heteroThickness"];


    const char* datasets[2] = {"MIT_05_2024", "simulation"};
    const char* in_data[3] = {"notarget", "homotarget", "heterotarget"};
    const char* target_data[3] = {"without a target", "with the homogeneous target", "with the heterogeneous target"};

    Char_t dataset[200];
    Char_t filename[100];
    sprintf(dataset, "%s", datasets[datasetSelect]);
    sprintf(filename, "%s", in_data[fileSelect]);

    
    Char_t histdesc[100];
    Char_t histname[100];
    Char_t title[100];
    Char_t file[100];
    Char_t in_path[200];
    Char_t out_path[200];
    
    Int_t eventCounter = 0; //store tree values measurement
    Int_t channel = 0;
    uint32_t timestamp_ns = 0;
    Long64_t timestamp_ps = 0;
    std::vector<std::vector<Double_t>> *trace = 0; //old
    //std::vector<Double_t> *trace = 0; new
    
    
    sprintf(in_path, "../data/%s/%s/input/", dataset, filename);
    sprintf(out_path, "../data/%s/%s/output/", dataset, filename);
    sprintf(file, "%s%s.root", in_path, filename);
    cout << "In path: " << file << endl;
    
    TFile *input = new TFile(file, "READ");
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }
    TTree *datatree;
    datatree = (TTree*)input->Get("vtree");
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }
    
    datatree->SetBranchAddress("BufferCounter", &eventCounter);
    datatree->SetBranchAddress("Channel", &channel);
    datatree->SetBranchAddress("TimeTag", &timestamp_ns);
    datatree->SetBranchAddress("TimeStampPico", &timestamp_ps);
    datatree->SetBranchAddress("Trace", &trace);
    
    energy_ch muon;
    muon.CH = {970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8};
    muon.o_CH = {100.52, 90.4, 100.65, 98.31, 107.62, 117.89, 98.19, 23.22, 99.15, 96.28, 116.67, 17.89, 19.18, 15.76, 9.33}; //old
    muon.E = 4.95;
    muon.o_E = 0.3433;
    
    energy_ch co60;
    co60.CH = {402.56, 401.75, 421.52, 444.63, 437.76, 428.52, 415.91, 49.03, 403.4, 417.63, 457.43, 52.08, 52.11, 47.11, 35.55};
    co60.o_CH = {97.79, 97.82, 101.99, 99.31, 115.65, 102.12, 90.23, 11.05, 96.9, 97.69, 111.66, 20.22, 23.25, 19.32, 12.67};
    co60.E = 1.25275;
    co60.o_E = 0.068;
    
    
    Double_t entries = datatree->GetEntries();
    // entries = 100000;
    Int_t prevEvent = -1;
    int missing_buffer_counter = 0;
    int pileup_counter = 0;
    int prepeak_step_counter = 0;
    int coinc_layer_counter = 0;
    
    TraceProperties trace_props;
    DetectorProperties* detectorProperties = new DetectorProperties();
    
    detectorProperties->SetScintillator(detectortype);
    detectorProperties->SetBeamEnergy(beamEnergy);
    detectorProperties->SetNLayers(nLayers);
    detectorProperties->SetScintillatorDimensions(crystalSize);
    detectorProperties->SetGapSizeZ(gapSizeZ);
    detectorProperties->SetTeflonThickness(teflonThickness);
    detectorProperties->SetAluThickness(aluThickness);
    detectorProperties->SetAbsorberType(absorberType);
    detectorProperties->SetAbsorberSize(absorberSize);
    detectorProperties->SetNSecondaryLayers(nSecLayers);
    detectorProperties->SetSecondaryLayerSizeZ(secLayerSizeZ);
    detectorProperties->SetTarget(target_data[fileSelect]);
    detectorProperties->SetNormStatus(normStatus);
    detectorProperties->SetReversedStatus(reversedStatus);
    detectorProperties->SetSimulationStatus(simulationStatus);
    detectorProperties->SetCalibrationStatus(!simulationStatus);
    detectorProperties->SetAbsorberStatus(absorberStatus);
    
    detectorProperties->Process();

    Calibration* calib = new Calibration(detectorProperties);
    calib->energy_extrapolation(co60, muon);
    
    detectorProperties->SetCalibration(calib);
    
    
    Detector* detector = new Detector(detectorProperties);
    detector->Info();
    
    std::vector<Particle> particles;
    std::vector<Particle> ScintParticles;
    
    Particle proton(nLayers, coincidenceTime, coincidenceLayer, calib);
    std::unordered_map<double, TraceProperties> initialEvents;
    std::multimap<double, TraceProperties> postEvents;
    bool bsaveTrace = false; 
    Plotter plotter(coincidenceTime);

    cout << "Dataset: " << dataset << endl;
    cout << "Input file: " << filename << endl;
    
    cout << "Measurement: Getting raw data." << endl;
    for (double e = 0; e<entries; e++){
    	datatree->GetEntry(e);
        channel--;
        trace_props.Clear();    
        trace_props.SetParameters(trace->at(0), channel, static_cast<double>(timestamp_ns)*2, static_cast<double>(timestamp_ps), discard_index, bsaveTrace);
        detector->EnergyHist(channel)->Fill(calib->GetQuenchedEnergy(channel, trace_props.amp)); //todo set birks to 0 if no quenching // sometimes double accounted for
        if(trace_props.channel == 0){
            initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            postEvents.insert({trace_props.time_ps, trace_props});
        }
    }
    cout << "Measurement: Processing raw data." << endl;
    for(const auto& [initialTime, initialTrace] : initialEvents) {
        proton.Clear();
        proton.InsertInitial(initialTrace);
        auto lowerTraces = postEvents.lower_bound(initialTime - coincidenceTime);
        auto upperTraces = postEvents.upper_bound(initialTime + coincidenceTime);            
        for (auto coincTrace = lowerTraces; coincTrace != upperTraces; ++coincTrace) {                                
            proton.Coincidence(coincTrace->second);
        }            
        particles.push_back(proton);
    }
    for(auto coinProton : particles){
        coinProton.Test();
        if(coinProton.missingChannel == true){
            missing_buffer_counter++;
            continue;
        }
        if(coinProton.pileupStatus == true){
            pileup_counter++;
            continue;
        }
        if(coinProton.ampOffsetStatus == true){
            prepeak_step_counter++;
            continue;
        }
        if(coinProton.coinc_layer < coincidenceLayer){
            coinc_layer_counter++;
            continue;
        }
        
        coinProton.SumEDep();
        detector->TotalEnergyHist()->Fill(coinProton.total_edep);
        for(int i=0; i<coinProton.traces.size(); i++){
            if(coinProton.GetEDep(i)!=0){
                detector->CoincEnergyHist(i)->Fill(coinProton.GetEDep(i));
            }
        }
        for (int i = coinProton.traces.size() - 1; i >= 0; --i) {
            if(coinProton.GetEDep(i)!=0){
                detector->StoppedEnergyHist(i)->Fill(coinProton.GetEDep(i));
                break;
            }
        }
    }
}