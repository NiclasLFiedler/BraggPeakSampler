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

void idealAnalysis(){
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

    bool bScintSim = true;
    bool bSavepdf = false;
    bool bPhotons = false;


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
    datatree = (TTree*)input->Get("vtree");
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

    datatree->SetBranchAddress("event", &event);
	datatree->SetBranchAddress("NDet", &NDet);
	datatree->SetBranchAddress("EDep", &EDep);
    
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
    cout << "Entries " << entries << endl; 
    
    TStopwatch timer;
    timer.Start();
    for(int64_t e = 0; e<entries; e++){
	    datatree->GetEntry(e);
        if(prevEvent != event){
            proton.ProcessSimulationEDep();
            for(int ch = 0; ch<nLayers; ch++){
                double FillEnergy = proton.GetEDep(ch);
                if(FillEnergy > 0.0){
                    detector->EnergyHist(ch)->Fill(FillEnergy);
                    detector->CoincEnergyHist(ch)->Fill(FillEnergy);
                }
                else{
                    break;
                }
            }
            detector->TotalEnergyHist()->Fill(proton.total_edep);
            proton.Clear();
        }
        proton.SetEDep(NDet, EDep);
        prevEvent = event; 
    }
    timer.Stop();
    timer.Print();
    
    cout << "Data acquisition & Histograms finished" << endl;
  
    cout << "Processing Detector" << endl;
    detector->Process();

    for(int i = 0;  i < nLayers; i++){
        if(detector->crystals.at(i).dose.dose > 0){
        std::cout << "Channel: " << i << " Total Energy Dose: " << detector->crystals.at(i).dose.dose << " +- " << detector->crystals.at(i).dose.stddev << " ~" << detector->crystals.at(i).dose.stddev/detector->crystals.at(i).dose.dose*100 << "%" <<std::endl;
        }
        else{
            std::cout << "Channel: " << i << " Total Energy Dose: " << detector->crystals.at(i).dose.dose << " +- " << detector->crystals.at(i).dose.stddev << std::endl;
        }
    }
    //------------------------End of Analysis-------------------//

    //------------------------Plots-----------------------------//
    sprintf(title, "Bragg Sampler %s Analysis", filename);
    TCanvas *c1 = new TCanvas("c1", title, 10, 10, 1900, 1000);
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    int coloums = 4;
    int rows = (nLayers+2)/coloums;
    int modCanvas = (nLayers+2)%coloums;
    if(modCanvas != 0) rows++;
    c1->Divide(coloums, rows);

    cout << setfill(' ');
    cout << setw(10) << "Channel" << setw(2) << "|" << setw(15) << "Unfiltered" << setw(2) << "|"  << setw(25) << "Coinc. /w channel 1" << setw(2) << "|" << setw(18) << " Stopped particle" << endl;
    cout << setfill('-') << setw(80) << "-" << endl;
    cout << setfill(' ');

    // for(int i = 0; i<nLayers; i++){
    //     c1->cd(i+1);
    //     plotter.Histogram1D(detector->EnergyHist(i), "EDep [MeV]", "Counts");
    //     detector->CoincEnergyHist(i)->SetLineColor(kRed+1);
    //     detector->CoincEnergyHist(i)->Draw("HIST SAME");
    //     detector->StoppedEnergyHist(i)->SetLineColor(kOrange+1);
    //     detector->StoppedEnergyHist(i)->Draw("HIST SAME");
    //     plotter.Legend(detector->EnergyHist(i));

    //     cout << setw(10) << i << setw(2) << "|" << setw(15) << detector->EnergyHist(i)->GetEntries() << setw(2) << "|" << setw(25) << detector->CoincEnergyHist(i)->GetEntries() << setw(2) << "|" << setw(18) <<  detector->StoppedEnergyHist(i)->GetEntries() << endl;
    // }
    cout << setfill('-') << setw(80) << "-" << endl;
    cout << setfill(' ');
    cout << "Detected missing buffer entries: " << missing_buffer_counter << endl;
    cout << "Detected pileups: " << pileup_counter << endl;
    cout << "Detected prepeak steps: " << prepeak_step_counter << endl;
    cout << "Number of lower coincidence particles: " << coinc_layer_counter << endl;
    c1->cd(nLayers+1);    

    sprintf(histdesc, "Norm. energy depth dose distribution %s", target_data[fileSelect]);
    plotter.GraphError(detector->MeansGraph(), "Depth [cm]", "Norm. Energy Dose [MeV]", histdesc);
    plotter.Legend(detector->MeansGraph());

    c1->cd(nLayers+2);
    detector->TotalEnergyHist()->SetLineColor(kGreen+1);
    plotter.Histogram1D(detector->TotalEnergyHist(), "EDep [MeV]", "Counts");
    
    //outputfile
    sprintf(file, "%splots.root", out_path);
    TFile *	hfile = new TFile(file,"RECREATE");

    c1->Write("ALL");
    //c1->Close();
    TCanvas *c2 = new TCanvas("c2", title, 4000, 2000);
    if(bSavepdf){
        for(int i = 0; i < nLayers; i++){
            c2->SetFillColor(0);
            c2->SetGrid();
	        c2->SetBorderMode(0);
	        c2->SetBorderSize(2);
	        c2->SetFrameBorderMode(0);

            plotter.Histogram1D(detector->CoincEnergyHist(i), "EDep [MeV]", "Counts");
            gStyle->SetTitleFontSize(0.07); // Title font size for the canvas
            plotter.Legend(detector->CoincEnergyHist(i));
            sprintf(file, "%spdf/h_edep_%i.pdf", out_path, i);
            c2->Update();
            c2->Print(file);
            c2->Clear();
        }

        c2->SetFillColor(0);
        c2->SetGrid();
	    c2->SetBorderMode(0);
	    c2->SetBorderSize(2);
	    c2->SetFrameBorderMode(0);

        gStyle->SetTitleFontSize(0.07); // Title font size for the canvas
        sprintf(histdesc, "Norm. energy depth dose distribution %s", target_data[fileSelect]);
        plotter.GraphError(detector->MeansGraph(), "Depth [cm]", "Norm. Energy Dose [MeV]", histdesc);
        plotter.Legend(detector->MeansGraph());
        c2->Update();
        sprintf(file, "%spdf/g_means.pdf", out_path);
        c2->Print(file);

        plotter.Histogram1D(detector->TotalEnergyHist(), "EDep [MeV]", "Counts");
        gStyle->SetTitleFontSize(0.07); // Title font size for the canvas
        plotter.Legend(detector->TotalEnergyHist());
        sprintf(file, "%spdf/h_total_edep.pdf", out_path);
        c2->Print(file);

        for(int i = 0; i < nLayers; i++){
            detector->CoincEnergyHist(i)->Write();
        }
        hfile->Close();
    }
    //store hist means for bortfield fit
    sprintf(file, "%s%sMeans.root", out_path, in_data[fileSelect]);
    std::cout << "Outpath Means file: " << file << std::endl;
    hfile = new TFile(file, "RECREATE");
    
    TTree *meantree = new TTree("meantree", "Tree storing histogram means");
    
    Double_t layer;
    Double_t layererr;
    Double_t mean;
    Double_t error;
    
    meantree->Branch("x", &layer);
    meantree->Branch("x_sigma", &layererr);
    meantree->Branch("mean", &mean);
    meantree->Branch("error", &error);

    for (int ch = 0; ch < nLayers; ++ch) {
        layer = detector->crystals.at(ch).pos.depth;
        layererr = detector->crystals.at(ch).pos.stddev;
        mean = detector->crystals.at(ch).dose.dose;
        error = detector->crystals.at(ch).dose.stddev;
        meantree->Fill();        
    }
    meantree->Write(); 


    hfile->Close();
    
}