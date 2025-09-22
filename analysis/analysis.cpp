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


    const char* datasets[3] = {"MIT_05_2024", "simulation", "beamtime"};
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
    
    energy_ch base;
    base.CH = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    base.o_CH = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    base.E = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    base.o_E = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    energy_ch na22;
    na22.CH = {1572.17, 1540.67, 1394.92, 49.7185, 1853.25, 2078.32, 1007.81, 324.048, 500, 1374.37, 726.037, 900, 514.404, 491.578, 1086.86, 921.839, 500, 500, 1041.1, 217.159, 500, 1008.07, 742.082, 756.624, 500, 353.818, 500, 881.16, 500, 500, 1126.8, 203.876};
    // na22.o_CH = {785.518, 623.046, 681.057, 361.922, 586.363, 534.047, 811.318, 1032.56, 1000, 648.283, 618.388, 1000, 334.118, 874.19, 678.237, 459.7, 1000, 1000, 517.583, 337.498, 1000, 527.741, 498.61, 484.256, 1000, 422.4, 1000, 510.73, 1000, 1000, 604.186, 317.11};

    // na22.E = 1.275;
    // na22.o_E = 0.01;
    
    energy_ch beam;
    beam.CH = {5419.79, 5673.85, 6221.08, 5718.44, 6121.69, 6272.27, 5510.53, 24060.1, 5106.76, 6495.73, 5359.34, 1369.69, 1742.07, 1350.71, 4038.64, 3880.83, 1902.56, 1303.45, 3584.7, 2008.28, 3627.02, 3616.78, 2946.28, 3312.72, 3856.61, 6518.53, 7242.99, 5724.06, 10907.4, 5587.3, 5645.68, 11026.2};

    beam.o_CH = {1278.59, 1066.01, 1145.18, 1114.99, 1179.06, 1134.02, 1123.82, 4018.15, 1124.86, 1167.49, 1184.63, 649.625, 702.544, 659.381, 653.853, 1328.2, 799.764, 760.142, 958.492, 756.957, 1018.61, 966.635, 941.913, 934.82, 2004, 1385.32, 1179.33, 1981.61, 2367.1, 1445.27, 1406.13, 2694.97};

    beam.E = {5.7038, 5.79716, 5.94443, 6.03835, 6.16728, 6.3360, 6.45465, 6.6303, 6.83547, 7.03770, 7.2600, 4.93424, 5.07636, 5.20989, 5.36086, 5.57017, 5.73576, 5.96942, 6.21447, 6.52709, 6.89378, 7.29601, 7.88811, 8.56514, 9.53018, 10.9653, 13.5613, 20.5618, 11.2318, 5.7038, 5.7038, 5.7038};

    beam.o_E = {0.14336, 0.149442, 0.14427, 0.143871, 0.154516, 0.143849, 0.14975, 0.149520, 0.141375, 0.146101, 0.13765, 0.13335, 0.139438, 0.133682, 0.130915, 0.112131, 0.116986, 0.112732, 0.118858, 0.118496, 0.120264, 0.126596, 0.143783, 0.164512, 0.229532, 0.378573, 0.843234, 3.08753, 6.45036, 0.14336, 0.14336, 0.14336};
    //4 - 93.01 + 19.97; 1779.12+690.66
    //5 - 78.576 + 9.49; 2034.19+675.6

    //10 -  85.4175 + 19.5839- 1137.05 + 696.43
    //11 -12.6472 +12.0657-
    //12 - 28.2955 + 9.40731 +179.523 + 399.546
    //13 -52.445 +13.3292 + 370.058 + 340.344
    //14  - 73.2931 + 15.3272 - 1277.45 +613.614
    //15 - 69.6274 + 16.6303 - 903.729 + 604.97
    //16 - 98.5465 + 23.1831 - 2342.06 +991.668
    //17 - 12.4756 + 4.49776 - 251.097 + 383.081
    //18 - 74.8514 + 20.607 - 1047.09+596.422
    //19 - 45.3675 + 17.5071 - 324.869 - 361.442
    //20 - 71.723 + 20.4579 - 1022.3 + 576.937
    //21 - 59.3453 + 33.7472 - 883.589 + 656.023
    //22 - 4.50820e+01 + 1.44533e+01 - 691.137 + 469.948
    //23 - 22.4054 + 4.1494 - 724.823 + 479.132
    //24 - 5.15129e+01+1.87916e+01 -630.284 + 545.392
    //25 - 26.7939 + 13.1999 - 310.344 + 438.827
    //26 - 19.8751 + 5.25402 +  630.666 + 523.842
    //27 - 55.379 + 22.6026 - 672.635 + 541.671
    //28 - 96.4732 + 22.9196 - 2577.5 + 824.663
    //29 - 61.937 + 21.9025 -  989.291 + 605.193
    //30 - 85.16 + 21.75 - 1574.37 + 590.742
    //31 - 102.57 + 17.25 - 2.97269e+03 - 1.03523e+03



    //new na22
    // 0 1599.14 + 882.484
    // 1 1496.94 + 647.818
    // 2
    // 3 

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
    detectorProperties->SetCalibrationStatus(true);
    detectorProperties->SetAbsorberStatus(absorberStatus);
    
    detectorProperties->Process();

    Calibration* calib = new Calibration(detectorProperties);
    calib->energy_extrapolation(base, beam);
    for(int i=0; i<nLayers; i++){
        printf("Layer %i: Channel Position %f Current: %f \n", i, calib->GetChannel(i, 1.275), na22.CH[i]);    
    }
    detectorProperties->SetCalibration(calib);
    
    
    Detector* detector = new Detector(detectorProperties);
    detector->Info();
    
    std::vector<Particle> particles;
    std::vector<Particle> ScintParticles;
    
    Particle proton(nLayers, coincidenceTime, coincidenceLayer, calib);
    std::map<double, TraceProperties> initialEvents;
    std::map<Long64_t, TraceProperties> secondaryEvents;
    std::multimap<double, TraceProperties> postEvents;
    bool bsaveTrace = false; 
    Plotter plotter(coincidenceTime);

    //tests
    TH1D* h_timediff = new TH1D("h_timediff", "h_timediff", 2000, 0, 1000);
    TH1D* h_inittime = new TH1D("h_inittime", "h_inittime", 200000, 0, 1e11);

    cout << "Dataset: " << dataset << endl;
    cout << "Input file: " << filename << endl;
    bool bTraceDisable = true;
    if(!simulationStatus){
        cout << "Measurement: Getting raw data." << endl;
        for (double e = 0; e<entries; e++){
            datatree->GetEntry(e);
            // printf("Timestamp: %llu ps, Timestamp<static> %f ps, \n", timestamp_ps, static_cast<double>(timestamp_ps));
            if (board == 1) {
                channel += 16;
                timestamp_ps += 32*1000;
            };
            trace_props.Clear();    
            trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000,  static_cast<double>(timestamp_ps), discard_index, bsaveTrace, bTraceDisable);
            
            if(trace_props.charge > 0) {
                detector->EnergyHist(channel)->Fill(calib->GetQuenchedEnergy(channel, trace_props.charge)); //todo set birks to 0 if no quenching // sometimes double accounted for
            }
            else{
                detector->EnergyHist(channel)->Fill(-1); //todo set birks to 0 if no quenching // sometimes double accounted for
            }
            if(trace_props.channel == 0){
                initialEvents.insert({trace_props.time_ps, trace_props});
            }
            else{
                postEvents.insert({trace_props.time_ps, trace_props});
            } 

            // if(channel == 0 && board == 0){
            //     secondaryEvents.insert({timestamp_ps, trace_props});
            // }
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
        //         timeDiff = (time2-time1);
        //         printf("Time1 %lld Tim2: %lld Timediff: %lld ns \n", time1, time2, timeDiff);
        //         h_timediff->Fill(timeDiff);
        //         time1 = time2;
        //     }
        // }

        for(const auto& [initialTime, initialTrace] : initialEvents) {
            //h_inittime->Fill(initialTime/1000);
            if(time1 == 0){
                time1 = initialTime;
            }
            else{
                time2 = initialTime;
                timeDiff = (time2-time1)/1000;
                //h_timediff->Fill((time2-time1)/1000);
                time1 = time2;
            }
            proton.Clear();
            proton.Insert(initialTrace);
            if((timeDiff)<1200){
                // cout << "Short time difference: " << timeDiff/1000 << " us, Event: " << iMod << endl;
                detector->StoppedEnergyHist(0)->Fill(proton.GetEDep(0));
            }
            auto lowerTraces = postEvents.lower_bound(initialTime - coincidenceTime*1000);
            auto upperTraces = postEvents.upper_bound(initialTime + coincidenceTime*1000);            
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
                //continue;
            }
            if(coinProton.ampOffsetStatus == true){
                prepeak_step_counter++;
                //continue;
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
                    //detector->StoppedEnergyHist(i)->Fill(coinProton.GetEDep(i));
                    break;
                }
            }
        }
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
        if(bScintSim){
            photontree->GetEntry(0);
            prevEvent = eventPhotons;
            proton.Clear();
            int64_t photonEntries = photontree->GetEntries();
            cout << "Photon Entries " << photonEntries << endl; 
            for(int64_t e = 0; e<photonEntries; e++){
		        photontree->GetEntry(e);
                detector->EntryHist(NDetPhotons)->Fill(EntryPosX);
                if( EntryPosX > 10){
                    // std::cout << "EntryPosX: " << EntryPosX << " NDetPhotons: " << NDetPhotons << " Event: " << eventPhotons << std::endl;
                }
                // std::cout << "Photon Event: " << eventPhotons << " NDetPhotons: " << NDetPhotons << " EntryPosX: " << EntryPosX << std::endl;
                detector->ExitHist(NDetPhotons)->Fill(ExitPosX);
                detector->AngleHist(NDetPhotons)->Fill(phi);
            }    
            
            
            
            
            
                // if(prevEvent != eventPhotons){
                    // ScintParticles.push_back(proton);
                    // proton.Clear();
                // }
                // proton.SetNPhotons(NDetPhotons, NPhotons);
                // prevEvent = event; 
            // }
            // for(Particle p : ScintParticles){
                // for(int ch = 0; ch<nLayers; ch++){
                    // if(p.GetNPhotons(ch) > 0){
                        // detector->PhotonHist(ch)->Fill(p.GetNPhotons(ch));
                        // if(p.CoincidencePhotons(ch)){
                            // detector->CoincPhotonHist(ch)->Fill(p.GetNPhotons(ch));
                        // }
                    // }
                // }
            // }
        }
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
        plotter.Histogram1D(detector->EnergyHist(i), "EDep [MeV]", "Counts");
        // detector->CoincEnergyHist(i)->SetLineColor(kRed+1);
        // detector->CoincEnergyHist(i)->Draw("HIST SAME");
        // detector->StoppedEnergyHist(i)->SetLineColor(kOrange+1);
        // detector->StoppedEnergyHist(i)->Draw("HIST SAME");
        plotter.Legend(detector->EnergyHist(i));

        cout << setw(10) << i << setw(2) << "|" << setw(15) << detector->EnergyHist(i)->GetEntries() << setw(2) << "|" << setw(25) << detector->CoincEnergyHist(i)->GetEntries() << setw(2) << "|" << setw(18) <<  detector->StoppedEnergyHist(i)->GetEntries() << endl;
    }
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
    
    Char_t heteroPath[200];
    sprintf(heteroPath, "../data/modulation/output/%ium_%immMeans.root", pmod, heteroThickness);
    std::cout << "Out Path: " << heteroPath << std::endl;
    TFile* heteroFile = new TFile(heteroPath, "RECREATE");
    detector->TotalEnergyHist()->Write();
    heteroFile->Close();
    
    c2->Close();   
}
