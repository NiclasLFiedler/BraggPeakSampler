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

double computeR2(TGraphErrors* graph, TF1* fit) {
    int n = graph->GetN();
    double *x = graph->GetX();
    double *y = graph->GetY();

    // Mean of y
    double y_mean = 0;
    for (int i = 0; i < n; ++i) y_mean += y[i];
    y_mean /= n;

    // Total sum of squares (SST) and sum of squares of residuals (SSR)
    double ss_tot = 0;
    double ss_res = 0;
    for (int i = 0; i < n; ++i) {
        double yi_fit = fit->Eval(x[i]);
        ss_res += TMath::Power(y[i] - yi_fit, 2);
        ss_tot += TMath::Power(y[i] - y_mean, 2);
    }

    double R2 = 1.0 - ss_res/ss_tot;
    return R2;
}


void photonAnalysis(){
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
    double teflonThickness = config["teflonThickness"];
    double aluThickness    = config["aluThickness"];
    int coincidenceTime                 = config["coincidenceTime"];
    int coincidenceLayer                = config["coincidenceLayer"];
    int discard_index                   = config["discardIndex"];
    int pmod                            = allConfigs["pmod"];
    int heteroThickness                 = allConfigs["heteroThickness"];

    bool bScintSim = true;
    bool bSavepdf = false;
    bool bPhotons = false;


    const char* datasets[3] = {"MIT_05_2024", "photon", "beamtime"};
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
    Char_t file2[100];
    Char_t in_path[200];
    Char_t out_path[200];

    
    Int_t event = 0;
    int eventPhotons = 0;
    Int_t NDet = 0;
    int NDetPhotons = 0;
    int NPhotons = 0;
    Double_t EDep = 0;
    Int_t photonScint = 0;
    Int_t photonSiPM = 0;
    Double_t phi = 0;
    Int_t TrackID;
    
    
    sprintf(in_path, "../data/%s/%s/input/", dataset, filename);
    sprintf(out_path, "../data/%s/%s/output/", dataset, filename);
    sprintf(file, "%s%s.root", in_path, filename);
    sprintf(file2, "%snoQuench.root", in_path);
    cout << "In path: " << file << endl;
    
    TFile *input = new TFile(file, "READ");
    TFile *input2 = new TFile(file2, "READ");

    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }
    TTree *datatree;
    TTree *datatree2;
    datatree = (TTree*)input->Get("RawData");
    datatree2 = (TTree*)input2->Get("RawData");
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

    datatree->SetBranchAddress("event", &event);
	datatree->SetBranchAddress("NDet", &NDet);
	datatree->SetBranchAddress("EDep", &EDep);
    datatree->SetBranchAddress("PhotonScint", &photonScint);
    datatree->SetBranchAddress("PhotonSiPM", &photonSiPM);
   
    datatree2->SetBranchAddress("event", &event);
	datatree2->SetBranchAddress("NDet", &NDet);
	datatree2->SetBranchAddress("EDep", &EDep);
    datatree2->SetBranchAddress("PhotonScint", &photonScint);
    datatree2->SetBranchAddress("PhotonSiPM", &photonSiPM);

    double resolution = 0.5;
    double edep_min = 5;
    double edep_max = 60;

    std::vector<TH1D*> h_scint;
    std::vector<TH1D*> h_sipm;

    for(int layer = 0; layer<(edep_max-edep_min)/resolution; layer++){
        sprintf(histdesc, "Created Photons edep %f", edep_min+layer*resolution);
        sprintf(histname, "h_scint_%i", layer);
        h_scint.push_back(new TH1D(histname, histdesc, 100, 0, 10000));
        sprintf(histdesc, "Measured Photons edep %f", edep_min+layer*resolution);
        sprintf(histname, "h_sipm_%i", layer);
        h_sipm.push_back(new TH1D(histname, histdesc, 100, 0, 2000));
    }

    std::vector<TH1D*> h_scint2;
    std::vector<TH1D*> h_sipm2;

    for(int layer = 0; layer<(edep_max-edep_min)/resolution; layer++){
        sprintf(histdesc, "Created Photons edep %f", edep_min+layer*resolution);
        sprintf(histname, "h_scint2_%i", layer);
        h_scint2.push_back(new TH1D(histname, histdesc, 100, 0, 10000));
        sprintf(histdesc, "Measured Photons edep %f", edep_min+layer*resolution);
        sprintf(histname, "h_sipm2_%i", layer);
        h_sipm2.push_back(new TH1D(histname, histdesc, 100, 0, 2000));
    }

    for(Long64_t i=0; i<datatree->GetEntries(); i++){
        datatree->GetEntry(i);
        if(i%1000 == 0){
            cout << "Processing event " << i << " / " << datatree->GetEntries() << "\r" << flush;
        }
        if(EDep > edep_min && EDep < edep_max){
            int index = (EDep-edep_min)/resolution;
            h_scint.at(index)->Fill(photonScint);
            h_sipm.at(index)->Fill(photonSiPM);
        }
        else{
            // std::cout << "Warning: EDep out of range: " << EDep << std::endl;
        }
    }

    for(Long64_t i=0; i<datatree2->GetEntries(); i++){
        datatree2->GetEntry(i);
        if(i%1000 == 0){
            cout << "Processing event " << i << " / " << datatree2->GetEntries() << "\r" << flush;
        }
        if(EDep > edep_min && EDep < edep_max){
            int index = (EDep-edep_min)/resolution;
            h_scint2.at(index)->Fill(photonScint);
            h_sipm2.at(index)->Fill(photonSiPM);
        }
        else{
            // std::cout << "Warning: EDep out of range: " << EDep << std::endl;
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Histograms", 1200, 800);
    c1->Divide(4, 10); // adjust grid depending on number of histograms

    for (size_t i = 0; i < h_scint.size(); ++i) {
        c1->cd(i+1);
        h_scint[i]->SetLineColor(kBlue);
        h_sipm[i]->SetLineColor(kRed);
        h_scint[i]->Draw();
        h_sipm[i]->Draw("SAME");

        gPad->BuildLegend();
    }
    
    std::vector<std::vector<double>> scint_means;
    std::vector<std::vector<double>> scint_means2;

    for(int layer = 0; layer<(edep_max-edep_min)/resolution; layer++){
        std::vector<double> means;
        if(h_scint.at(layer)->GetEntries() < 10){
            continue;
        }
        means.push_back(edep_min+layer*resolution);
        means.push_back(h_scint.at(layer)->GetMean());
        means.push_back(h_scint.at(layer)->GetStdDev());
        means.push_back(h_sipm.at(layer)->GetMean());
        means.push_back(h_sipm.at(layer)->GetStdDev());

        scint_means.push_back(means);
    }

    for(int layer = 0; layer<(edep_max-edep_min)/resolution; layer++){
        std::vector<double> means2;
        if(h_scint2.at(layer)->GetEntries() < 10){
            continue;
        }
        means2.push_back(edep_min+layer*resolution);
        means2.push_back(h_scint2.at(layer)->GetMean());
        means2.push_back(h_scint2.at(layer)->GetStdDev());
        means2.push_back(h_sipm2.at(layer)->GetMean());
        means2.push_back(h_sipm2.at(layer)->GetStdDev());

        scint_means2.push_back(means2);
    }

    TGraphErrors *g_scint = new TGraphErrors();
    TGraphErrors *g_sipm = new TGraphErrors();
    TGraphErrors *g_quotient = new TGraphErrors();
    double err = 0;
    for(size_t i = 0; i < scint_means.size(); i++){
        g_scint->SetPoint(i, scint_means.at(i).at(0), scint_means.at(i).at(1));
        g_scint->SetPointError(i, resolution/sqrt(12), scint_means.at(i).at(2));
        
        g_sipm->SetPoint(i, scint_means.at(i).at(0), scint_means.at(i).at(3));
        g_sipm->SetPointError(i, resolution/sqrt(12), scint_means.at(i).at(4));
        err = sqrt(pow(scint_means.at(i).at(4)/scint_means.at(i).at(1),2)+pow(scint_means.at(i).at(3)/(scint_means.at(i).at(1)*scint_means.at(i).at(1))*scint_means.at(i).at(2),2));
        g_quotient->SetPoint(i, scint_means.at(i).at(0), scint_means.at(i).at(3)/scint_means.at(i).at(1)*10000);
        g_quotient->SetPointError(i, 0, err);
    }

    TGraphErrors *g_scint2 = new TGraphErrors();
    TGraphErrors *g_sipm2 = new TGraphErrors();
    TGraphErrors *g_quotient2 = new TGraphErrors();

    for(size_t i = 0; i < scint_means2.size(); i++){
        g_scint2->SetPoint(i, scint_means2.at(i).at(0), scint_means2.at(i).at(1));
        g_scint2->SetPointError(i, resolution/sqrt(12), scint_means2.at(i).at(2));
        
        g_sipm2->SetPoint(i, scint_means2.at(i).at(0), scint_means2.at(i).at(3));
        g_sipm2->SetPointError(i, resolution/sqrt(12), scint_means2.at(i).at(4));
        err = sqrt(pow(scint_means2.at(i).at(4)/scint_means2.at(i).at(1),2)+pow(scint_means2.at(i).at(3)/(scint_means2.at(i).at(1)*scint_means2.at(i).at(1))*scint_means2.at(i).at(2),2));
        g_quotient2->SetPoint(i, scint_means2.at(i).at(0), scint_means2.at(i).at(3)/scint_means2.at(i).at(1)*10000);
        g_quotient2->SetPointError(i, 0, err);
    }

    double xmin_scint = g_scint->GetX()[0], xmax_scint = g_scint->GetX()[0];
    double ymin_scint = g_scint->GetY()[0], ymax_scint = g_scint->GetY()[0];
    
    double xmin_sipm = g_sipm->GetX()[0], xmax_sipm = g_sipm->GetX()[0];
    double ymin_sipm = g_sipm->GetY()[0], ymax_sipm = g_sipm->GetY()[0];

    double xmin_quotient = g_quotient->GetX()[0], xmax_quotient = g_quotient->GetX()[0];
    double ymin_quotient = g_quotient->GetY()[0], ymax_quotient = g_quotient->GetY()[0];

    TF1 *fit_scint = new TF1("fit_scint", "[0] + [1]*x", xmin_scint, xmax_scint);
    TF1 *fit_sipm = new TF1("fit_sipm", "[0] + [1]*x", xmin_sipm, xmax_sipm);
    TF1 *fit_quotient = new TF1("fit_quotient", "[0] + [1]*x", xmin_quotient, xmax_quotient);


    
    double xmin_scint2 = g_scint2->GetX()[0], xmax_scint2 = g_scint2->GetX()[0];
    double ymin_scint2 = g_scint2->GetY()[0], ymax_scint2 = g_scint2->GetY()[0];
    
    double xmin_sipm2 = g_sipm2->GetX()[0], xmax_sipm2 = g_sipm2->GetX()[0];
    double ymin_sipm2 = g_sipm2->GetY()[0], ymax_sipm2 = g_sipm2->GetY()[0];

    double xmin_quotient2 = g_quotient2->GetX()[0], xmax_quotient2 = g_quotient2->GetX()[0];
    double ymin_quotient2 = g_quotient2->GetY()[0], ymax_quotient2 = g_quotient2->GetY()[0];

    TF1 *fit_scint2 = new TF1("fit_scint", "[0] + [1]*x", xmin_scint2, xmax_scint2);
    TF1 *fit_sipm2 = new TF1("fit_sipm", "[0] + [1]*x", xmin_sipm2, xmax_sipm2);
    TF1 *fit_quotient2 = new TF1("fit_quotient", "[0] + [1]*x", xmin_quotient2, xmax_quotient2);



    fit_scint->SetLineColor(kBlue);
    fit_sipm->SetLineColor(kRed);
    fit_quotient->SetLineColor(kGreen+2);

    fit_scint->SetLineWidth(2);
    fit_sipm->SetLineWidth(2);
    fit_quotient->SetLineWidth(2);

    fit_scint->SetRange(0, 25);
    fit_sipm->SetRange(0, 25);
    
    g_scint->Fit(fit_scint, "R");
    g_sipm->Fit(fit_sipm, "R");
    g_quotient->Fit(fit_quotient, "R");

    fit_scint->SetRange(0, 60);
    fit_sipm->SetRange(0, 60);

    Int_t lightRed   = TColor::GetColor("#ff1994ff");
    Int_t lightGreen = TColor::GetColor("#91ff30ff");
    Int_t lightBlue  = TColor::GetColor("#3cbbffff");

    fit_scint2->SetLineColor(lightBlue);
    fit_sipm2->SetLineColor(lightRed);
    fit_quotient2->SetLineColor(lightGreen);

    fit_scint2->SetLineWidth(2);
    fit_sipm2->SetLineWidth(2);
    fit_quotient2->SetLineWidth(2);

    fit_scint2->SetRange(0, 25);
    fit_sipm2->SetRange(0, 25);

    g_scint2->Fit(fit_scint2, "R");
    g_sipm2->Fit(fit_sipm2, "R");
    g_quotient2->Fit(fit_quotient2, "R");

    fit_scint2->SetRange(0, 60);
    fit_sipm2->SetRange(0, 60);

    TCanvas *c2 = new TCanvas("c2", "Scintillator and SiPM Comparison", 1200, 800);
    c2->Divide(1, 1);  // Single pad
    c2->SetGrid();
    // Customize graphs
    g_scint->SetMarkerStyle(20);
    g_scint->SetMarkerColor(kBlue);
    g_scint->SetLineColor(kBlue);
    g_scint->SetTitle("Scintillator, SiPM, and Quotient; Energy Deposition / MeV; Counts");

    g_sipm->SetMarkerStyle(21);
    g_sipm->SetMarkerColor(kRed);
    g_sipm->SetLineColor(kRed);
    
    g_quotient->SetMarkerStyle(22);
    g_quotient->SetMarkerColor(kGreen+2);
    g_quotient->SetLineColor(kGreen+2);

    g_scint2->SetMarkerStyle(20);
    g_scint2->SetMarkerColor(lightBlue);
    g_scint2->SetLineColor(lightBlue);

    g_sipm2->SetMarkerStyle(21);
    g_sipm2->SetMarkerColor(lightRed);
    g_sipm2->SetLineColor(lightRed);
    
    g_quotient2->SetMarkerStyle(22);
    g_quotient2->SetMarkerColor(lightGreen);
    g_quotient2->SetLineColor(lightGreen);

    // Draw first graph with axes
    g_scint->GetXaxis()->SetRangeUser(0, 65);
    g_scint->GetYaxis()->SetRangeUser(0, 5000);
    g_scint->Draw("AP");

    // Overlay the others
    g_sipm->Draw("P SAME");
    g_quotient->Draw("P SAME");

    g_scint2->Draw("P SAME");
    g_sipm2->Draw("P SAME");
    g_quotient2->Draw("P SAME");


    fit_scint->Draw("SAME");
    fit_sipm->Draw("SAME");
    fit_quotient->Draw("SAME");
    
    fit_scint2->Draw("SAME");
    fit_sipm2->Draw("SAME");
    fit_quotient2->Draw("SAME");

    double R2_scint = computeR2(g_scint, fit_scint);
    double R2_sipm = computeR2(g_sipm, fit_sipm);
    double R2_quotient = computeR2(g_quotient, fit_quotient);

    // Add a legend
    TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
    leg->AddEntry(g_scint, "Scintillator signal", "lep");
    leg->AddEntry(g_sipm, "SiPM signal", "lep");
    leg->AddEntry(g_quotient, "SiPM/Scintillator ratio", "lep");
    leg->AddEntry(g_scint2, "Quenched: Scintillator signal", "lep");
    leg->AddEntry(g_sipm2, "Quenched: SiPM signal", "lep");
    leg->AddEntry(g_quotient2, "Quenched: SiPM/Scintillator ratio", "lep");
    leg->AddEntry(fit_scint, Form("y = %.2f + %.2f x", 
                                  fit_scint->GetParameter(0),
                                  fit_scint->GetParameter(1)), "l");
    leg->AddEntry(fit_sipm, Form("y = %.2f + %.2f x", 
                                  fit_sipm->GetParameter(0),
                                  fit_sipm->GetParameter(1)), "l");
    leg->AddEntry(fit_quotient, Form("y = %.2f + %.2f x", 
                                  fit_quotient->GetParameter(0),
                                  fit_quotient->GetParameter(1)), "l");
    
    leg->AddEntry(fit_scint2, Form("Quenched: y = %.2f + %.2f x", 
                                  fit_scint2->GetParameter(0),
                                  fit_scint2->GetParameter(1)), "l");
    leg->AddEntry(fit_sipm2, Form("Quenched: y = %.2f + %.2f x", 
                                  fit_sipm2->GetParameter(0),
                                  fit_sipm2->GetParameter(1)), "l");
    leg->AddEntry(fit_quotient2, Form("Quenched: y = %.2f + %.2f x", 
                                  fit_quotient2->GetParameter(0),
                                  fit_quotient2->GetParameter(1)), "l");
                                 
    leg->Draw();

    // Optional: grid and style
    c1->Update();

}