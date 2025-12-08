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

    const char* datasets[3] = {"MIT_05_2024", "simulation", "postCalib"};
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
    

    // background.root            na1617_300s.root (X)
    // background_1516.root       na1819_co23_340s.root (X)
    // co1617_na01_200s.root      na2021_co45_300s.root (X)
    // co1819_na23_200s.root  (X)    na2223_co67_300s.root  (X)
    // co2021_na45_200s.root  (X)    na23.root (X)
    // co2021_na45_200s_2nd.root (X) na23_200s.root (X)
    // co2223_na67_200s.root   (X)  na2425_co89_300s.root (X)
    // co2425_na89_200s.root   (?)   na2627_co1011_300s.root (X)
    // co2627_na1011_200s.root (X)  na2829_co1213_300s.root (X)
    // co2829_na1213_200s.root (X)  na3031_co1415_300s.root (X)
    // co3031_na1415_220s.root (X)   na3031_co1516_200s.root (X)
    // na01.root         (X)        na45.root (X)
    // na1011_300s.root  (X)        na67_200s.root (X)
    // na1213_320s.root  (X)         na89_300s.root (X)
    // na1415_co01_300s.root(X)


    // na22.CH = {1572.17, 1540.67, 1394.92, 49.7185, 1853.25, 2078.32, 1007.81, 324.048, 500, 1374.37, 726.037, 900, 514.404, 491.578, 1086.86, 921.839, 500, 500, 1041.1, 217.159, 500, 1008.07, 742.082, 756.624, 500, 353.818, 500, 881.16, 500, 500, 1126.8, 203.876};
    // na22.o_CH = {785.518, 623.046, 681.057, 361.922, 586.363, 534.047, 811.318, 1032.56, 1000, 648.283, 618.388, 1000, 334.118, 874.19, 678.237, 459.7, 1000, 1000, 517.583, 337.498, 1000, 527.741, 498.61, 484.256, 1000, 422.4, 1000, 510.73, 1000, 1000, 604.186, 317.11};

    //          na22                    c060
    //0   1572.17+-785.518              1718.21+-651.8
    //1   1540.67+-623.046              1417.37 +-611.446
    //2   1394.92+-681.057 (1438.3+-660.218)              1497.63 +-650.941 (coinc: 169.083+- 418.9)
    //3  49.7185+-361.922               96.751 +- 34.684
    //4  1853.25+-586.363  (1594.97+-664.613)             1799.13 +- 641.491(coinc?: 218.933+- 440.838)
    //5  2078.32+-534.047  (2013.39+-442.469)             2117.79 +- 625.909(coinc?: 187.256+-453.846)
    //6  1007.81+- 811.318 (1152.73+-1334.39)        ?1820.87+- 529.955
    //7  94.8475+-339.275 (coinc: 324.048 +- 1032.56)
    //8  88.355+- 267.49                 1141.68+- 1141.68
    //9  1374.37+-648.283   (1658.43 +-621.172)
    //10 726.037+- 618.388  (951.97 +-1145.0)             674.796+-901.862
    //11 103.091 +- 298.051 (900+1000)
    //12 35.938 +- 292.8 (514.404+-334.118)
    //13 93.93+- 273.44 (491.578+-874.19)
    //14 1448.6+-483.7(coinc?: 204.775+-392.306;  1086.86+-678.237)       1478.07 +- 529.0 (1124.81+-501.79)
    //15 1373.48+-491.068(coinc?: 188.322+-404.762; 921.839+-459.7)    1703.42 +- 624.532 (1126.0+-449.60)
    //16 51.33 +- 36.478
    //17 92.674+-331.529
    //18 1041.1+- 517.583 (1324.47+-356.078)
    //19 217.159+-337.498 (96.935+-364.724)
    //20 500+1000
    //21 ? 1008.07 +-  527.741
    //22 742.082+-498.61                             653.698+-464.584
    //23 756.624+-484.256                            688.189+-487.616
    //24 500+1000
    //25 353.818+-422.4
    //26 78.3982+-352.849 
    //27 881.16+-510.73 (785.04+-352.849)
    //28                                             2423.1+-1070.58
    //29                                             696.856+-545.819
    //30 1126.8+- 604.186 ()                         1109.41+-500.988
    //31 203.876+-317.11 (294.426+-582.268)          314.728 +-604.677


    //simulated energies 

    //0 5.7038 +-   0.14336  5419.79+-1278.59
    //1 5.79716 +-  0.149442 5673.85+-1066.01
    //2 5.94443 +-  0.14427  6221.08+-1145.18
    //3 6.03835 +-  0.143871 5718.44+-1114.99
    //4 6.16728 +-  0.154516 6121.69+-1179.06
    //5 6.3360  +-  0.143849 6272.27+-1134.02
    //6 6.45465 +-  0.14975  5510.53+-1123.82
    //7 6.6303 +-   0.149520 24060.1+-4018.15
    //8 6.83547 +-  0.141375 5106.76+-1124.86
    //9 7.03770 +-  0.146101 6495.73+-1167.49
    //10 7.2600 +-  0.13765  5359.34+-1184.63
    //11 4.93424 +- 0.13335  1369.69+-649.625
    //12 5.07636 +- 0.139438 1742.07+-702.544
    //13 5.20989 +- 0.133682 1350.71+-659.381
    //14 5.36086 +- 0.130915 4038.64+-653.853
    //15 5.57017 +- 0.112131 3880.83+-1328.2    
    //16 5.73576 +- 0.116986 1902.56+-799.764
    //17 5.96942 +- 0.112732 1303.45+-760.142
    //18 6.21447 +- 0.118858 3584.7+- 958.492
    //19 6.52709 +- 0.118496 2008.28+-756.957 
    //20 6.89378 +- 0.120264 3627.02+-1018.61
    //21 7.29601 +- 0.126596 3616.78+-966.635
    //22 7.88811 +- 0.143783 2946.28+-941.913
    //23 8.56514 +- 0.164512 3312.72+-934.82
    //24 9.53018 +- 0.229532 3856.61+-2004
    //25 10.9653 +- 0.378573 6518.53+-1385.32
    //26 13.5613 +- 0.843234 7242.99+-1179.33
    //27 20.5618 +- 3.08753  5724.06+-1981.61
    //28 11.2318 +- 6.45036  10907.4+-2367.1
    //29 5.7038 +-  0.14336  5587.3+- 1445.27
    //30 5.7038 +-  0.14336  5645.68+-1406.13
    //31 5.7038 +-  0.14336  11026.2+-2694.97


// 5.7038, 5.79716, 5.94443, 6.03835, 6.16728, 6.3360, 6.45465, 6.6303, 6.83547, 7.03770, 7.2600, 4.93424, 5.07636, 5.20989, 5.36086, 5.57017, 5.73576, 5.96942, 6.21447, 6.52709, 6.89378, 7.29601, 7.88811, 8.56514, 9.53018, 10.9653, 13.5613, 20.5618, 11.2318, 5.7038, 5.7038, 5.7038

// 0.14336, 0.149442, 0.14427, 0.143871, 0.154516, 0.143849, 0.14975, 0.149520, 0.141375, 0.146101, 0.13765, 0.13335, 0.139438, 0.133682, 0.130915, 0.112131, 0.116986, 0.112732, 0.118858, 0.118496, 0.120264, 0.126596, 0.143783, 0.164512, 0.229532, 0.378573, 0.843234, 3.08753, 6.45036, 0.14336, 0.14336, 0.14336

// 5419.79, 5673.85, 6221.08, 5718.44, 6121.69, 6272.27, 5510.53, 24060.1, 5106.76, 6495.73, 5359.34, 1369.69, 1742.07, 1350.71, 4038.64, 3880.83, 1902.56, 1303.45, 3584.7, 2008.28, 3627.02, 3616.78, 2946.28, 3312.72, 3856.61, 6518.53, 7242.99, 5724.06, 10907.4, 5587.3, 5645.68, 11026.2

// 1278.59, 1066.01, 1145.18, 1114.99, 1179.06, 1134.02, 1123.82, 4018.15, 1124.86, 1167.49, 1184.63, 649.625, 702.544, 659.381, 653.853, 1328.2, 799.764, 760.142, 958.492, 756.957, 1018.61, 966.635, 941.913, 934.82, 2004, 1385.32, 1179.33, 1981.61, 2367.1, 1445.27, 1406.13, 2694.97

    sprintf(in_path, "../data/%s/%s/input/", dataset, filename);
    sprintf(out_path, "../data/%s/%s/output/", dataset, filename);
    sprintf(file, "%sco60test.root", in_path);//, filename);
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
    sprintf(backgroundpath, "%sbackground.root", in_path);
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
    
    // energy_ch muon;
    // muon.CH = {970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8};
    // muon.o_CH = {100.52, 90.4, 100.65, 98.31, 107.62, 117.89, 98.19, 23.22, 99.15, 96.28, 116.67, 17.89, 19.18, 15.76, 9.33}; //old
    // muon.E = 4.95;
    // muon.o_E = 0.3433;
    
    // energy_ch co60;
    // co60.CH = {402.56, 401.75, 421.52, 444.63, 437.76, 428.52, 415.91, 49.03, 403.4, 417.63, 457.43, 52.08, 52.11, 47.11, 35.55};
    // co60.o_CH = {97.79, 97.82, 101.99, 99.31, 115.65, 102.12, 90.23, 11.05, 96.9, 97.69, 111.66, 20.22, 23.25, 19.32, 12.67};
    // co60.E = 1.25275;
    // co60.o_E = 0.068;
    
    
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
    // calib->energy_extrapolation(co60, muon);
    
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

    double ylim = 10000;
    TH1D* h_charge_0 = new TH1D("h_charge_0", "h_charge_0", 1000, 0, ylim);
    TH1D* h_charge_1 = new TH1D("h_charge_1", "h_charge_1", 1000, 0, ylim);
    TH1D* h_charge_2 = new TH1D("h_charge_2", "h_charge_2", 1000, 0, ylim);
    TH1D* h_charge_3 = new TH1D("h_charge_3", "h_charge_3", 1000, 0, ylim);
    TH1D* h_chargeDiff = new TH1D("h_chargeDiff", "h_chargeDiff", 1000, 0, 3000);
    //co2627_na1011_200s
    int coincChannel = 0;
    bool useCharge = true;
    bool useCoinc = true;

    for (double e = 0; e<entries; e++){
        trace_props.Clear();    
        datatree->GetEntry(e);
        if (board == 1) {
            channel += 16;
            timestamp_ps += 32*1000;
        };
        trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000, static_cast<double>(timestamp_ps), discard_index, bsaveTrace, true);
        // printf("test\n");

        if (!useCoinc){
            if(channel == coincChannel){
                if(useCharge) h_charge_0->Fill(trace_props.charge);
                else h_charge_0->Fill(trace_props.amp);
            }
            else if(channel == coincChannel+1)
            {
                if(useCharge) h_charge_1->Fill(trace_props.charge);
                else h_charge_1->Fill(trace_props.amp);
            }
        }
        if(trace_props.channel == coincChannel){
            initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            postEvents.insert({trace_props.time_ps, trace_props});
        }
    }


    for (double e = 0; e<entries; e++){
        trace_props.Clear();    
        datatree2->GetEntry(e);
        if (board == 1) {
            channel += 16;
            timestamp_ps += 32*1000;
        };
        
        trace_props.SetParameters(*trace, channel, charge, static_cast<double>(timestamp_ps)/1000, static_cast<double>(timestamp_ps), discard_index, bsaveTrace, true);
        //detector->EnergyHist(channel)->Fill(trace_props.amp);

        if (!useCoinc){
            if(channel == coincChannel){
                h_charge_2->Fill(trace_props.charge);
            }
            else if(channel == coincChannel+1)
            {
                h_charge_3->Fill(trace_props.charge);
            }
        }
        if(trace_props.channel == coincChannel){
            // initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            // postEvents.insert({trace_props.time_ps, trace_props});
        }
    }


    cout << "Measurement: Processing raw data." << endl;
    TH1D* h_timediff = new TH1D("h_timediff", "h_timediff", 800, -200, 200);
    for(const auto& [initialTime, initialTrace] : initialEvents) {
        proton.Clear();
        proton.Insert(initialTrace);
        auto lowerTraces = postEvents.lower_bound(initialTime - coincidenceTime*1000);
        auto upperTraces = postEvents.upper_bound(initialTime + coincidenceTime*1000);            
        for (auto coincTrace = lowerTraces; coincTrace != upperTraces; ++coincTrace) {                                
            proton.Coincidence(coincTrace->second, coincChannel);
            if(coincTrace->second.channel == coincChannel+1) h_timediff->Fill(static_cast<double>(coincTrace->first-initialTime)/1000);
        }            
        particles.push_back(proton);
    }

    if (useCoinc){
        for(auto coinProton : particles){
            coinProton.Test();
            if(coinProton.GetCharge(coincChannel+1) != 0){
                if(useCharge){
                    
                    if(std::abs(coinProton.GetCharge(coincChannel)-coinProton.GetCharge(coincChannel+1)) > 0){
                        h_charge_0->Fill(coinProton.GetCharge(coincChannel));
                        h_charge_1->Fill(coinProton.GetCharge(coincChannel+1)); 
                        h_chargeDiff->Fill(std::abs(coinProton.GetCharge(coincChannel)-coinProton.GetCharge(coincChannel+1)));
                    }    
                }
                else{
                    h_charge_0->Fill(coinProton.GetAmplitude(coincChannel));
                    h_charge_1->Fill(coinProton.GetAmplitude(coincChannel+1));
                }
            }
        }
    }

    if(!useCoinc){
        h_charge_0->Add(h_charge_2, -1);
        h_charge_1->Add(h_charge_3, -1);

        // h_charge_0->Add(h_charge_2, -0.1666);
        // h_charge_1->Add(h_charge_3, -0.1666);
        // h_charge_0->Add(h_charge_2, -0.333);
        // h_charge_1->Add(h_charge_3, -0.333);
        // h_charge_0->Add(h_charge_2, -0.366);
        // h_charge_1->Add(h_charge_3, -0.366);
        // h_charge_0->Add(h_charge_2, -0.5);
        // h_charge_1->Add(h_charge_3, -0.5);
        // h_charge_0->Add(h_charge_2, -0.5333);
        // h_charge_1->Add(h_charge_3, -0.5333);
    }

    // TF1 *gaus = new TF1("gaus", "gaus");
    // h_charge_0->Fit(gaus, "Q"); // "Q" = quiet, remove for verbose

    // double mean  = gaus->GetParameter(1);
    // double sigma = gaus->GetParameter(2);
    // double ampl  = gaus->GetParameter(0);


    // std::cout << "Fit results:\n";
    // std::cout << "  Charge = " << ampl  << "\n";
    // std::cout << "  Mean      = " << mean  << "\n";
    // std::cout << "  Sigma     = " << sigma << "\n";

    // h_charge_1->Fit(gaus, "Q"); // "Q" = quiet, remove for verbose
    
    // mean  = gaus->GetParameter(1);
    // sigma = gaus->GetParameter(2);
    // ampl  = gaus->GetParameter(0);

    // std::cout << "Fit results:\n";
    // std::cout << "  Charge = " << ampl  << "\n";
    // std::cout << "  Mean      = " << mean  << "\n";
    // std::cout << "  Sigma     = " << sigma << "\n";

    cout << "Processing Detector" << endl;
    detector->Process();

   //------------------------Plots-----------------------------//
    sprintf(title, "Bragg Sampler %s Analysis", filename);
    TCanvas *c1 = new TCanvas("c1", title, 10, 10, 1000, 700);
    c1->Divide(3,2); // two pads side by side

    // Draw first Charge histogram
    c1->cd(1);
    h_charge_0->SetLineColor(kBlue);
    h_charge_0->SetTitle("Charge Spectrum Channel 0");
    h_charge_0->GetXaxis()->SetTitle("Charge");
    h_charge_0->GetYaxis()->SetTitle("Counts");
    h_charge_0->Draw();

    // Draw second Charge histogram
    c1->cd(2);
    h_charge_1->SetLineColor(kRed);
    h_charge_1->SetTitle("Charge Spectrum Channel 1");
    h_charge_1->GetXaxis()->SetTitle("Charge");
    h_charge_1->GetYaxis()->SetTitle("Counts");
    h_charge_1->Draw();
    
    c1->cd(3);
    h_timediff->SetLineColor(kRed);
    h_timediff->SetTitle("Time Spectrum Channel 1");
    h_timediff->GetXaxis()->SetTitle("Time");
    h_timediff->GetYaxis()->SetTitle("Counts");
    h_timediff->Draw();

        // Draw first Charge histogram
    c1->cd(4);
    h_charge_2->SetLineColor(kBlue);
    h_charge_2->SetTitle("Charge Background Spectrum Channel 0");
    h_charge_2->GetXaxis()->SetTitle("Charge");
    h_charge_2->GetYaxis()->SetTitle("Counts");
    h_charge_2->Draw();

    // Draw second Charge histogram
    c1->cd(5);
    h_charge_3->SetLineColor(kRed);
    h_charge_3->SetTitle("Charge Background Spectrum Channel 1");
    h_charge_3->GetXaxis()->SetTitle("Charge");
    h_charge_3->GetYaxis()->SetTitle("Counts");
    h_charge_3->Draw();

    c1->cd(6);
    h_chargeDiff->SetLineColor(kOrange);
    h_chargeDiff->SetTitle("Charge Difference");
    h_chargeDiff->GetXaxis()->SetTitle("Charge");
    h_chargeDiff->GetYaxis()->SetTitle("Counts");
    h_chargeDiff->Draw();
}
