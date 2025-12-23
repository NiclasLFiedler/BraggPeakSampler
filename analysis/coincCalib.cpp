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

void coincCalib(){
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
    double teflonThickness              = config["teflonThickness"];
    double aluThickness                 = config["aluThickness"];
    int coincidenceTime                 = config["coincidenceTime"];
    int coincidenceLayer                = config["coincidenceLayer"];
    int discard_index                   = config["discardIndex"];
    int pmod                            = allConfigs["pmod"];
    int heteroThickness                 = allConfigs["heteroThickness"];

    bool bScintSim = false;
    bool bSavepdf = true;
    bool bPhotons = false;

    const char* datasets[3] = {"MIT_05_2024", "simulation", "paperBeamtime"};
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
    Char_t outFile_path[200];
    u_int32_t eventCounter = 0; //store tree values measurement
    Int_t channel = 0;
    Int_t board = 0;
    Int_t charge = 0;
    uint64_t timestamp_ns = 0;
    // uint64_t timestamp_ps = 0;
    Long64_t timestamp_ps = 0;

    // std::vector<std::vector<Double_t>> *trace = 0; //old
    std::vector<u_int16_t> *trace = 0;

    Int_t event = 0;
    int eventPhotons = 0;
    Int_t NDet = 0;
    int NDetPhotons = 0;
    int NPhotons = 0;
    Double_t EDep = 0;
    Double_t EntryPosX = 0;
    Double_t EntryPosY = 0;
    Double_t EntryPosZ = 0;
    Double_t ExitPosX = 0;
    Double_t ExitPosY = 0;
    Double_t ExitPosZ = 0;
    Double_t phi = 0;
    Int_t TrackID;
    sprintf(in_path, "../data/%s/%s/input/", dataset, filename);
    sprintf(out_path, "../data/%s/%s/output/", dataset, filename);

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
    detectorProperties->SetCalibrationStatus(false);
    detectorProperties->SetAbsorberStatus(absorberStatus);

    detectorProperties->Process();

    Calibration* calib = new Calibration(detectorProperties);
    // calib->energy_extrapolation(co60, muon);

    detectorProperties->SetCalibration(calib);


    Detector* detector = new Detector(detectorProperties);
    detector->Info();

    std::vector<Particle> particles;

    Particle proton(nLayers, coincidenceTime, coincidenceLayer, calib);
    std::unordered_map<double, TraceProperties> initialEvents;
    std::multimap<double, TraceProperties> postEvents;
    bool bsaveTrace = false;
    Plotter plotter(coincidenceTime);

    cout << "Dataset: " << dataset << endl;
    cout << "Input file: " << filename << endl;

    cout << "Measurement: Getting raw data." << endl;

    double ylim = 65535;
    std::vector<TH1D*> h_coincCharge;
    std::vector<TH1D*> h_timediff;
    for (int i = 0; i < nLayers; i++){
        sprintf(histname, "h_coincCharge_%d", i);
        sprintf(histdesc, "h_coincCharge_%d", i);
        TH1D* h_coincCharge_i = new TH1D(histname, histdesc, 65536/256, 0, ylim);
        h_coincCharge.push_back(h_coincCharge_i);
    }

    for(int i = 0; i < nLayers/2; i++){
        sprintf(histname, "h_timediff_%d%d", i*2, i*2+1);
        sprintf(histdesc, "h_timediff_%d%d", i*2, i*2+1);
        TH1D* h_timediff_i = new TH1D(histname, histdesc, 800, -200, 200);
        h_timediff.push_back(h_timediff_i);
    }

    bool useCharge = true;
    bool useCoinc = false;
    Double_t entries;

    TFile *input;
    TTree *datatree;
    TTree *datatree2;
    TFile *input2;

    int coincChannel = 28;
    bool second = 1;
    particles.clear();
    initialEvents.clear();
    postEvents.clear();

    //0 co60_%d%d.root 50ns 5.99
    //1 co60_%d%d.root nocinc 5.76
    //2 co60_%d%d.root nocoinc 6.04
    //3 ch%d%d_co60_2.root 2ns 6.22 1640
    //4 co60_%d%d.root 4ns 6.3 
    //5 ch%d%d_co60.root 25ns 6.4  1700
    //6 %sch%d%d_co60  1ns 6.52 cutoff 1410
    //7 ch%d%d_co60_2.root 25ns 6.93  1630
    //8 co60_%d%d.root 2ns 7.02               
    //9 sch%d%d_co60.root 1ns 7.16 
    //10 co60_%d%d.root 25ns 6.7             ch%d%d_co60_2.root 1ns 7.4 cutoff 1150
    //11                                     ch%d%d_co6.root 3ns 5.1
    //12                                     co60ch%d%d.root 4ns 5.25 cutoff 900
    //13                                     co60ch%d%d.root 20ns 5.38 
    //14                                     co60ch%d%d.root 1ns 5.56
    //15                                     co60ch%d%d.root 2ns 5.82 
    //16                                     co60ch%d%d.root 1ns 5.95 
    //17                                     co60ch%d%d.root 1ns 6.14
    //18                                     co60ch%d%d.root 10ns 6.47  
    //19                                     co60ch%d%d.root 40ns 6.84  
    //20                                     co60ch%d%d.root 30ns 7.22 
    //21                                     co60ch%d%d.root 40ns 7.73
    //22                                     sch%d%d_co60.root 2ns 8.24
    //23                                     ch%d%d_co60.root 10ns 9.09
    //24                                     %sch%d%d_co60.root 14ns 10.12
    //25                                     %sch%d%d_co60.root 2ns 12.24
    //26                                     ch%d%d_co60.root 12ns 17.32
    //27                                     ch%d%d_co60.root 12ns 19.24

    //28                                     sch%d%d_co60 15ns 12.09

    //29                                     %sco60_%d%d.root 2ns 2.9

    //30                                     %sco60_%d%d.root 2ns 3.32
    //31                                     %sco60_%d%d.root 2ns 3.2

    // sprintf(file, "%sch%d%d_co60.root", in_path, coincChannel, coincChannel+1);
    sprintf(file, "%sco60_%d%d.root", in_path, coincChannel, coincChannel+1);
    // sprintf(file, "%sco60ch%d%d.root", in_path, coincChannel, coincChannel+1);
    // sprintf(file, "%sbeamtimebackground.root", in_path); 
    cout << "In path: " << file << endl;
    input = new TFile(file, "READ");
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }
    datatree = (TTree*)input->Get("RawData");
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }
    datatree->SetBranchAddress("EventCounter", &eventCounter);
    datatree->SetBranchAddress("Channel", &channel);
    datatree->SetBranchAddress("Board", &board);
    datatree->SetBranchAddress("Charge", &charge);
    datatree->SetBranchAddress("TimeStamp", &timestamp_ps);
    datatree->SetBranchAddress("Trace", &trace);
    entries = datatree->GetEntries();
    for (double e = 0; e<entries; e++){
        trace_props.Clear();
        datatree->GetEntry(e);
        if(charge < 0) charge += 65536;
        if (board == 1) {
            channel += 16;
            timestamp_ps += 32*1000;
        };
        trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000, static_cast<double>(timestamp_ps), discard_index, bsaveTrace, true);

        if (!useCoinc){
            if(channel == coincChannel){
                h_coincCharge.at(coincChannel)->Fill(trace_props.charge);
            }
            else if(channel == coincChannel+1)
            {
                h_coincCharge.at(coincChannel+1)->Fill(trace_props.charge);
            }
        }

        if(trace_props.channel == coincChannel){
            initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            postEvents.insert({trace_props.time_ps, trace_props});
        }
    }
    cout << "Measurement: Processing raw data." << endl;
    for(const auto& [initialTime, initialTrace] : initialEvents) {
        proton.Clear();
        proton.Insert(initialTrace);
        auto lowerTraces = postEvents.lower_bound(initialTime - coincidenceTime*1000);
        auto upperTraces = postEvents.upper_bound(initialTime + coincidenceTime*1000);
        for (auto coincTrace = lowerTraces; coincTrace != upperTraces; ++coincTrace) {
            proton.Coincidence(coincTrace->second, coincChannel);
            if(coincTrace->second.channel == coincChannel+1) h_timediff.at(coincChannel/2)->Fill(static_cast<double>(coincTrace->first-initialTime)/1000);
        }
        
        particles.push_back(proton);
    }
    if (useCoinc){
        for(auto coinProton : particles){
            coinProton.Test();
            if(coinProton.GetCharge(coincChannel+1) != 0){
                if(std::abs(coinProton.GetCharge(coincChannel)-coinProton.GetCharge(coincChannel+1)) > 0){
                    h_coincCharge.at(coincChannel)->Fill(coinProton.GetCharge(coincChannel));
                    h_coincCharge.at(coincChannel+1)->Fill(coinProton.GetCharge(coincChannel+1));
                }
            }
        }
    }
    
    if(second) coincChannel++;
    sprintf(outFile_path, "../data/%s/%s/output/coincHistogram%d.root", dataset, filename, coincChannel);
    TFile* outputFile = new TFile(outFile_path, "RECREATE");
    // if(coincChannel == 21){
    //     h_coincCharge.at(coincChannel) = h_coincCharge.at(coincChannel+1)    
    // }
    h_coincCharge.at(coincChannel)->Write();
    h_timediff.at(int(coincChannel/2))->Write();
    outputFile->Close();
    particles.clear();
}
