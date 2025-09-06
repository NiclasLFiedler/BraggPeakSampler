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

void sourceCalib(){
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

    const char* datasets[3] = {"MIT_05_2024", "simulation", "test"};
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
    sprintf(file, "%s%s.root", in_path, filename);
    cout << "In path: " << file << endl;
    
    TFile *input = new TFile(file, "READ");
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }
    TTree *datatree;
    TTree *datatree2;
    TTree *photontree;
    
    datatree = (TTree*)input->Get("RawData");

    Char_t backgroundpath[200];
    sprintf(backgroundpath, "%snotargetBack.root", in_path);
    TFile *input2 = new TFile(backgroundpath, "READ");
    datatree2 = (TTree*)input2->Get("RawData");

    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

    if(!simulationStatus){
        datatree->SetBranchAddress("EventCounter", &eventCounter);
        datatree->SetBranchAddress("Channel", &channel);
        datatree->SetBranchAddress("Board", &board);
        datatree->SetBranchAddress("Charge", &charge);
        datatree->SetBranchAddress("TimeStamp", &timestamp_ps);
        datatree->SetBranchAddress("Trace", &trace);

        datatree2->SetBranchAddress("EventCounter", &eventCounter);
        datatree2->SetBranchAddress("Channel", &channel);
        datatree2->SetBranchAddress("Board", &board);
        datatree2->SetBranchAddress("Charge", &charge);
        datatree2->SetBranchAddress("TimeStamp", &timestamp_ps);
        datatree2->SetBranchAddress("Trace", &trace);
    }
    else{
        photontree = (TTree*)input->Get("ptree");
        
        datatree->SetBranchAddress("event", &event);
        datatree->SetBranchAddress("track", &TrackID);
	    datatree->SetBranchAddress("NDet", &NDet);
	    datatree->SetBranchAddress("EDep", &EDep);
        
        if(bScintSim){
            if (photontree) {
                photontree->SetBranchAddress("eventPhotons", &eventPhotons);
	            photontree->SetBranchAddress("NDetPhotons", &NDetPhotons);
                photontree->SetBranchAddress("NPhotons", &NPhotons);
                photontree->SetBranchAddress("EntryPosX", &EntryPosX);
                photontree->SetBranchAddress("EntryPosY", &EntryPosY);
                photontree->SetBranchAddress("EntryPosZ", &EntryPosZ);
                photontree->SetBranchAddress("ExitPosX", &ExitPosX);
                photontree->SetBranchAddress("ExitPosY", &ExitPosY);
                photontree->SetBranchAddress("ExitPosZ", &ExitPosZ);
                photontree->SetBranchAddress("AngleZ", &phi);
            }
            else{
                cout << "Error: Failed to retrieve tree 'ptree' from file!" << endl;
            }
        }   
    }
    
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
    detectorProperties->SetCalibrationStatus(false);
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


    TH1D* h_amp_0 = new TH1D("h_amp_0", "h_amp_0", 1000, 0, 5000);
    TH1D* h_amp_1 = new TH1D("h_amp_1", "h_amp_1", 1000, 0, 5000);
    int coincChannel = 10;
    bool useCharge = true;
    for (double e = 0; e<entries; e++){
        trace_props.Clear();    
        datatree->GetEntry(e);
        if (board == 1) {
            channel += 16;
            timestamp_ps += 32*1000;
        };
        
        trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000, static_cast<double>(timestamp_ps), discard_index, bsaveTrace);
        detector->EnergyHist(channel)->Fill(trace_props.amp);

        if(channel == coincChannel){
            if(useCharge) h_amp_0->Fill(trace_props.charge);
            else h_amp_0->Fill(trace_props.amp);
        }
        else if(channel == coincChannel+1)
        {
            if(useCharge) h_amp_1->Fill(trace_props.charge);
            else h_amp_1->Fill(trace_props.amp);
        }

        if(trace_props.channel == coincChannel){
            initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            postEvents.insert({trace_props.time_ps, trace_props});
        }
    }

    TH1D* h_amp_2 = new TH1D("h_amp_2", "h_amp_2", 1000, 0, 6000);
    TH1D* h_amp_3 = new TH1D("h_amp_3", "h_amp_3", 1000, 0, 6000);

    for (double e = 0; e<entries; e++){
        trace_props.Clear();    
        datatree2->GetEntry(e);
        if (board == 1) {
            channel += 16;
            timestamp_ps += 32*1000;
        };
        
        trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000, static_cast<double>(timestamp_ps), discard_index, bsaveTrace);
        //detector->EnergyHist(channel)->Fill(trace_props.amp);

        if(channel == coincChannel){
            h_amp_0->Fill(trace_props.amp,-1);
        }
        else if(channel == coincChannel+1)
        {
            h_amp_1->Fill(trace_props.amp,-1);
        }

        if(trace_props.channel == coincChannel){
            // initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            // postEvents.insert({trace_props.time_ps, trace_props});
        }
    }


    cout << "Measurement: Processing raw data." << endl;
    TH1D* h_timediff = new TH1D("h_timediff", "h_timediff", 800, -10000, 10000);
    for(const auto& [initialTime, initialTrace] : initialEvents) {
        proton.Clear();
        proton.InsertInitial(initialTrace);
        auto lowerTraces = postEvents.lower_bound(initialTime - coincidenceTime*1000);
        auto upperTraces = postEvents.upper_bound(initialTime + coincidenceTime*1000);            
        for (auto coincTrace = lowerTraces; coincTrace != upperTraces; ++coincTrace) {                                
            proton.Coincidence(coincTrace->second, coincChannel);
            if(coincTrace->second.channel == coincChannel+1) h_timediff->Fill(static_cast<double>(coincTrace->first-initialTime)/1000);
        }            
        particles.push_back(proton);
    }


    for(auto coinProton : particles){
        coinProton.Test();
        if(coinProton.GetCharge(coincChannel+1) != 0){
            if(useCharge){
                // h_amp_0->Fill(coinProton.GetCharge(coincChannel));
                // h_amp_1->Fill(coinProton.GetCharge(coincChannel+1));    
            }
            else{
                // h_amp_0->Fill(coinProton.GetAmplitude(coincChannel));
                // h_amp_1->Fill(coinProton.GetAmplitude(coincChannel+1));
            }
        }

    }

    // TF1 *gaus = new TF1("gaus", "gaus");
    // h_amp_0->Fit(gaus, "Q"); // "Q" = quiet, remove for verbose

    // double mean  = gaus->GetParameter(1);
    // double sigma = gaus->GetParameter(2);
    // double ampl  = gaus->GetParameter(0);


    // std::cout << "Fit results:\n";
    // std::cout << "  Amplitude = " << ampl  << "\n";
    // std::cout << "  Mean      = " << mean  << "\n";
    // std::cout << "  Sigma     = " << sigma << "\n";

    // h_amp_1->Fit(gaus, "Q"); // "Q" = quiet, remove for verbose
    
    // mean  = gaus->GetParameter(1);
    // sigma = gaus->GetParameter(2);
    // ampl  = gaus->GetParameter(0);

    // std::cout << "Fit results:\n";
    // std::cout << "  Amplitude = " << ampl  << "\n";
    // std::cout << "  Mean      = " << mean  << "\n";
    // std::cout << "  Sigma     = " << sigma << "\n";

    cout << "Processing Detector" << endl;
    detector->Process();

   //------------------------Plots-----------------------------//
    sprintf(title, "Bragg Sampler %s Analysis", filename);
    TCanvas *c1 = new TCanvas("c1", title, 10, 10, 1900, 1000);
    c1->Divide(3,2); // two pads side by side

    // Draw first amplitude histogram
    c1->cd(1);
    h_amp_0->SetLineColor(kBlue);
    h_amp_0->SetTitle("Amplitude Spectrum Channel 0");
    h_amp_0->GetXaxis()->SetTitle("Amplitude");
    h_amp_0->GetYaxis()->SetTitle("Counts");
    h_amp_0->Draw();

    // Draw second amplitude histogram
    c1->cd(2);
    h_amp_1->SetLineColor(kRed);
    h_amp_1->SetTitle("Amplitude Spectrum Channel 1");
    h_amp_1->GetXaxis()->SetTitle("Amplitude");
    h_amp_1->GetYaxis()->SetTitle("Counts");
    h_amp_1->Draw();
    
    c1->cd(3);
    h_timediff->SetLineColor(kRed);
    h_timediff->SetTitle("Time Spectrum Channel 1");
    h_timediff->GetXaxis()->SetTitle("Time");
    h_timediff->GetYaxis()->SetTitle("Counts");
    h_timediff->Draw();

        // Draw first amplitude histogram
    c1->cd(4);
    h_amp_2->SetLineColor(kBlue);
    h_amp_2->SetTitle("Amplitude Background Spectrum Channel 0");
    h_amp_2->GetXaxis()->SetTitle("Amplitude");
    h_amp_2->GetYaxis()->SetTitle("Counts");
    h_amp_2->Draw();

    // Draw second amplitude histogram
    c1->cd(5);
    h_amp_3->SetLineColor(kRed);
    h_amp_3->SetTitle("Amplitude Background Spectrum Channel 1");
    h_amp_3->GetXaxis()->SetTitle("Amplitude");
    h_amp_3->GetYaxis()->SetTitle("Counts");
    h_amp_3->Draw();
}
