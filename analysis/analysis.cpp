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

void analysis(){
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
    bool bSavepdf = false;
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
    
    u_int32_t eventCounter = 0; //store tree values measurement
    Int_t channel = 0;
    Int_t board = 0;
    Int_t charge = 0;
    uint64_t timestamp_ns = 0;
    Long64_t timestamp_ps = 0;

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
    // sprintf(file, "%s%s.root", in_path, "beamtimebackground");
    cout << "In path: " << file << endl;
    
    TFile *input = new TFile(file, "READ");
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }
    TTree *datatree;
    TTree *photontree;

    datatree = (TTree*)input->Get("RawData");
    // datatree = (TTree*)input->Get("vtree");
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
    
    energy_ch co60_1;
    co60_1.CH = {6410.23, 6008.26, 7291.29, 7270.89, 6354.13, 6787.23, 6553.78, 6340.23, 7951.94, 6043.47, 5709.45, 5186.46, 5187.57, 4614.36, 4608.14, 4652.09, 4479.48, 4501.15, 4317.02, 4119.14, 4686.39, 4788.02, 5525.51, 4077.01, 5684.99, 5131.04, 6410.33, 6605.79, 5964.62, 7388.72, 7951.13, 9030.46};
    for (int i = 0; i<32; i++){
        co60_1.CH.at(i)=co60_1.CH.at(i)/4;
        co60_1.E.at(i) = 1.17323;
        co60_1.o_CH.at(i) = 0;
        co60_1.o_E.at(i) = 0;
    }
    
    energy_ch co60_2;
    co60_2.CH = {7280.39, 6823.85, 8281.04, 8257.88, 7216.67, 7708.56, 7443.42, 7200.89, 9031.37, 6863.84, 6484.48, 5890.50, 5891.76, 5240.74, 5233.67, 5283.58, 5087.55, 5112.15, 4903.04, 4678.29, 5322.54, 5437.97, 6282.01, 4635.20, 6456.70, 5833.00, 7280.50, 7502.49, 6774.28, 8391.70, 8079.49, 9176.24};
    for (int i = 0; i<32; i++){
        co60_2.CH.at(i)=co60_2.CH.at(i)/4;
        co60_2.E.at(i) = 1.33249;
        co60_2.o_CH.at(i) = 0;
        co60_2.o_E.at(i) = 0;
    }

    Double_t entries = datatree->GetEntries();
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
    calib->energy_extrapolation(co60_1, co60_2);
    detectorProperties->SetCalibration(calib);
    
    
    Detector* detector = new Detector(detectorProperties);
    detector->Info();
    
    std::vector<Particle> particles;
    std::vector<Particle> ScintParticles;
    
    Particle proton(nLayers, coincidenceTime, coincidenceLayer, calib);
    std::map<Long64_t, TraceProperties> initialEvents;
    std::map<Long64_t, TraceProperties> secondaryEvents;
    std::multimap<Long64_t, TraceProperties> postEvents;
    bool bsaveTrace = false; 
    Plotter plotter(coincidenceTime);

    //tests
    TH1D* h_timediff = new TH1D("h_timediff", "h_timediff", 2000, 0, 1000);
    TH1D* h_inittime = new TH1D("h_inittime", "h_inittime", 200000, 0, 1e11);

    cout << "Dataset: " << dataset << endl;
    cout << "Input file: " << filename << endl;
    bool bTraceDisable = true;

    sprintf(file, "%scoincidence_output.root", out_path);
    TFile* f = new TFile(file, "RECREATE");


    TTree* t = new TTree("events", "Coincident events");

    double Q[32] = {0.0};
    t->Branch("Charge", Q, "Charge[32]/D");

    if(!simulationStatus){
        cout << "Measurement: Getting raw data." << endl;
        for (double e = 0; e<entries; e++){
            datatree->GetEntry(e);
            if(charge < 0) charge = 0; // charge += 65536;
            // printf("Timestamp: %llu ps, Timestamp<static> %f ps, \n", timestamp_ps, static_cast<double>(timestamp_ps));
            if (board == 1) {
                channel += 16;
                timestamp_ps += 32*1000;
            };
            trace_props.Clear();    
            trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000,  timestamp_ps, discard_index, bsaveTrace, bTraceDisable);
            
            if(trace_props.charge > 0) {
                // detector->EnergyHist(channel)->Fill(trace_props.charge);
                detector->EnergyHist(channel)->Fill(calib->GetQuenchedEnergy(channel, trace_props.charge));
            }
            else{
                // std::cout << "Energy zero or negative!" << " WARNING " << trace_props.charge << std::endl;
            }
            if(trace_props.channel == 0){
                initialEvents.insert({trace_props.time_ps, trace_props});
            }
            else{
                postEvents.insert({trace_props.time_ps, trace_props});
            } 

            if(channel == 16 && board == 1){
                secondaryEvents.insert({timestamp_ps, trace_props});
            }
        }
        cout << "Measurement: Processing raw data." << endl;
        
        Long64_t time1 = 0;
        Long64_t time2 = 0;
        Long64_t timeDiff = 0;
        int iMod = 0;
        // for(const auto& [secondaryTime, secondaryTrace] : secondaryEvents) {
        //     //h_inittime->Fill(secondaryTime/1000);
        //     if(time1 == 0){
        //         time1 = secondaryTime;
        //     }
        //     else{
        //         time2 = secondaryTime;
        //         timeDiff = (time2-time1)/1000;
        //         printf("Time1 %lld Tim2: %lld ps Timediff: %lld ns \n", time1, time2, timeDiff);
        //         h_timediff->Fill(timeDiff);
        //         time1 = time2;
        //     }
        // }

        for(const auto& [initialTime, initialTrace] : initialEvents) {
            // h_inittime->Fill(initialTime/1000);
            // if(time1 == 0){
            //     time1 = initialTime;
            // }
            // else{
            //     time2 = initialTime;
            //     timeDiff = (time2-time1)/1000;
            //     h_timediff->Fill(timeDiff);
            //     time1 = time2;
            // }
            // printf("Time1 %lld Tim2: %lld ps Timediff: %lld ns \n", time1, time2, timeDiff);

            proton.Clear();
            proton.Insert(initialTrace);
            // if((timeDiff)<1200){
            //     // cout << "Short time difference: " << timeDiff/1000 << " us, Event: " << iMod << endl;
            //     detector->StoppedEnergyHist(0)->Fill(proton.GetEDep(0));
            // }
            auto lowerTraces = postEvents.lower_bound(initialTime - coincidenceTime*1000);
            auto upperTraces = postEvents.upper_bound(initialTime + coincidenceTime*1000);            
            for (auto coincTrace = lowerTraces; coincTrace != upperTraces; ++coincTrace) {                                
                proton.Coincidence(coincTrace->second);
            }            
            particles.push_back(proton);
        }
        for(auto coinProton : particles){
            std::fill(Q, Q + 32, 0.0);
            coinProton.Test();
            if(coinProton.missingChannel == true){
                missing_buffer_counter++;
                continue;
            }
            if(coinProton.coinc_layer < coincidenceLayer){
                coinc_layer_counter++;
                continue;
            }
            
            coinProton.SumEDep();
            detector->TotalEnergyHist()->Fill(coinProton.total_edep);
            for(int i=0; i<coinProton.traces.size(); i++){
                if(coinProton.GetEDep(i)>0){
                    Q[i] = coinProton.GetEDep(i);
                    detector->CoincEnergyHist(i)->Fill(coinProton.GetEDep(i));
                }
            }
            t->Fill();
            for (int i = coinProton.traces.size() - 1; i >= 0; --i) {
                if(coinProton.GetEDep(i)!=0){
                    //detector->StoppedEnergyHist(i)->Fill(coinProton.GetEDep(i));
                    break;
                }
            }
        }
        t->Write();
        f->Close();
        cout << "Data acquisition & Histograms finished" << endl;
    }
    else{
        cout << "Simulation: Getting raw data from ROOT file." << endl;
        cout << "Entries " << entries << endl; 
        for(int64_t e = 0; e<entries; e++){
		    datatree->GetEntry(e);

            if(prevEvent != event){
                proton.ProcessEDep();
                for(int ch = 0; ch<nLayers; ch++){
                    if(proton.GetEDep(ch) > 0.0){
                        detector->EnergyHist(ch)->Fill(proton.GetEDep(ch));
                    }
                    if(proton.Coincidence(ch) && proton.GetEDep(ch) > 0.0){
                        detector->CoincEnergyHist(ch)->Fill(proton.GetEDep(ch));
                    }
                }
                detector->TotalEnergyHist()->Fill(proton.total_edep);
                proton.Clear();
            }

            if(TrackID == 1){
                proton.SetdE(NDet, EDep);
            }
            prevEvent = event; 
        }
        cout << "Data acquisition & Histograms finished" << endl;
    }
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
    TCanvas *c1 = new TCanvas("c1", title, 10, 10, 800, 500);
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

    for(int i = 0; i<nLayers; i++){
        c1->cd(i+1);
        // plotter.Histogram1D(detector->EnergyHist(i), "EDep [MeV]", "Counts");
        // detector->CoincEnergyHist(i)->SetLineColor(kRed+1);
        // detector->CoincEnergyHist(i)->Draw("HIST SAME");

        plotter.Histogram1D(detector->CoincEnergyHist(i), "EDep [MeV]", "Counts");

        // detector->StoppedEnergyHist(i)->SetLineColor(kOrange+1);
        // detector->StoppedEnergyHist(i)->Draw("HIST SAME");
        plotter.Legend(detector->CoincEnergyHist(i));

        cout << setw(10) << i << setw(2) << "|" << setw(15) << detector->EnergyHist(i)->GetEntries() << setw(2) << "|" << setw(25) << detector->CoincEnergyHist(i)->GetEntries() << setw(2) << "|" << setw(18) <<  detector->StoppedEnergyHist(i)->GetEntries() << endl;
    }
    cout << setfill('-') << setw(80) << "-" << endl;
    cout << setfill(' ');
    cout << "Detected missing buffer entries: " << missing_buffer_counter << endl;
    cout << "Detected pileups: " << pileup_counter << endl;
    cout << "Detected prepeak steps: " << prepeak_step_counter << endl;
    cout << "Number of lower coincidence particles: " << coinc_layer_counter << endl;
    c1->cd(nLayers+1);    

    for(int i = 0; i<nLayers; i++){
        std::cout << "Channel " << i << " Mean " <<  detector->CoincEnergyHist(i)->GetMean() << " Error " << detector->CoincEnergyHist(i)->GetStdDev() << std::endl;
    }
    sprintf(histdesc, "Norm. energy depth dose distribution %s", target_data[fileSelect]);
    plotter.GraphError(detector->MeansGraph(), "Depth [cm]", "Norm. Energy Dose [MeV]", histdesc);
    plotter.Legend(detector->MeansGraph());

    c1->cd(nLayers+2);
    detector->TotalEnergyHist()->SetLineColor(kGreen+1);
    plotter.Histogram1D(detector->TotalEnergyHist(), "EDep [MeV]", "Counts");
    
    c1->cd(nLayers+3);
    plotter.Histogram1D(h_timediff, "time / ns", "Counts");

    c1->cd(nLayers+4);
    plotter.Histogram1D(h_inittime, "time / ns", "Counts");

    //outputfile
    sprintf(file, "%splots.root", out_path);
    TFile *	hfile = new TFile(file,"RECREATE");

    c1->Write("ALL");
    //c1->Close();
    TCanvas *c2 = new TCanvas("c2", title, 800, 500);
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

    if(bPhotons){
        TCanvas* c3 = new TCanvas("c3", title, 10, 10, 1900, 1000);
        c3->SetFillColor(0);
	    c3->SetGrid();
	    c3->SetBorderMode(0);
	    c3->SetBorderSize(2);
	    c3->SetFrameBorderMode(0);
        c3->Divide(coloums, rows);
        
        for(int i = 0; i<nLayers; i++){
            c3->cd(i+1);
            plotter.Histogram1D(detector->PhotonHist(i), "# Photons", "Counts");
            detector->CoincPhotonHist(i)->SetLineColor(kRed+1);
            detector->CoincPhotonHist(i)->Draw("HIST SAME");
        }
    }

    hfile->Close();
    
    // Char_t heteroPath[200];
    // sprintf(heteroPath, "../data/modulation/output/%ium_%immMeans.root", pmod, heteroThickness);
    // std::cout << "Out Path: " << heteroPath << std::endl;
    // TFile* heteroFile = new TFile(heteroPath, "RECREATE");
    // detector->TotalEnergyHist()->Write();
    // heteroFile->Close();
    
    // c2->Close();   
}
