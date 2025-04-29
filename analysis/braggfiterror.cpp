#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"
#include "src/Plotter.cpp"

const char* datasets[2] = {"MIT_05_2024", "simulation"};
const char* in_data[3] = {"notarget", "homotarget", "heterotarget"};
const char* target_title[3] = {"without a target", "with the homogeneous target", "with the heterogeneous target"};

int dataset_select = 1;
int file_select = 2;

double CalcWETConv(double p_h2o, double p_m, double alpha_h2o, double alpha_m, double energy){
    return p_h2o/p_m*alpha_h2o/alpha_m*std::pow(energy, p_h2o-p_m);
}

double CalcWETConvSigma(double p_h2o, double p_m, double alpha_h2o, double alpha_m, double energy){
    return p_h2o/p_m*alpha_h2o/alpha_m*(p_h2o-p_m)*std::pow(energy, p_h2o-p_m-1);
}

void braggfiterror(){
    ROOT::EnableImplicitMT();
    Plotter plotter(20);
    Char_t file_curve[100];
    Char_t file_means[100];
    Char_t pdf_target[100];
    Char_t file[100];
    
    Char_t dataset[200];
    Char_t filename[200];
    Char_t title[100];

    sprintf(dataset, "%s",datasets[dataset_select]);
    sprintf(filename, "%s",in_data[file_select]);
    sprintf(title, "%s",target_title[file_select]);
    
    Float_t Phi0 = 0; //store tree values
    Float_t R0 = 0;
    Float_t sigma = 0;
    Float_t epsilon = 0;
    Float_t curve[4000];
    

    sprintf(file_curve, "../data/%s/%s/output/%s_fit.root", dataset, filename, filename);
    sprintf(file_means, "../data/%s/%s/output/%sMeans.root", dataset, filename, filename);

    TFile *input = new TFile(file_curve, "READ");
    TFile *input_error = new TFile(file_means, "READ");
    // Check if the file is successfully opened
    if (!input || input->IsZombie() || !input_error || input_error->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }

    //get data from tree
    TTree *datatree = (TTree*)input->Get("vtree");
    TTree *datatree_error = (TTree*)input_error->Get("meantree");

    // Check if the tree is successfully retrieved
    if (!datatree || !datatree_error) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        input_error->Close();
        return;
    }

    datatree->SetBranchAddress("Phi0", &Phi0);
    datatree->SetBranchAddress("R0", &R0);
    datatree->SetBranchAddress("sigma", &sigma);
    datatree->SetBranchAddress("epsilon", &epsilon);
    datatree->SetBranchAddress("curve", &curve);  
    
    Double_t mean;
    Double_t error;
    datatree_error->SetBranchAddress("mean", &mean);
    datatree_error->SetBranchAddress("error", &error);

    int max_heatmap = 70;
    Char_t h_title[200];
    sprintf(h_title, "Accumulated depth dose distribution %s", title);
    TH2D *h2_curve = new TH2D("h2_curve", h_title, 4000, 0, 40, 500, 0, max_heatmap); 
    sprintf(h_title, "Accumulated Phi0/N0 parameters %s", title);
    TH1D *h_Phi0 = new TH1D("h_Phi0", h_title, 100, 0, 3);
    sprintf(h_title, "Accumulated R0 parameters %s", title);
    TH1D *h_R0 = new TH1D("h_R0", h_title, 100, 28, 34);
    if (file_select == 1) {
        h_R0->SetBins(100, 24, 28); // modify the binning for "homo"
    }
    if (file_select == 2) {
        h_R0->SetBins(100, 22, 28); // modify the binning for "homo"
    }
    sprintf(h_title, "Accumulated sigma parameters %s", title);
    TH1D *h_sigma = new TH1D("h_sigma", h_title, 100, 0, 2);
    sprintf(h_title, "Accumulated epsilon parameters %s", title);
    TH1D *h_epsilon = new TH1D("h_epsilon", h_title, 100, -0.1, 1);

    int x[4001];
    for(int i = 0; i < 4001; i++){
        x[i] = i;
    }

    Double_t entries = datatree->GetEntries();
    datatree->GetEntry(0);
    Float_t Phi0_max = Phi0;
    Float_t Phi0_min = Phi0;
    Float_t R0_max = R0;
    Float_t R0_min = R0;
    Float_t sigma_max = sigma;
    Float_t sigma_min = sigma;
    Float_t epsilon_max = epsilon;
    Float_t epsilon_min = epsilon;

    for (double e = 0; e<entries; e++){
        datatree->GetEntry(e);
        h_Phi0->Fill(Phi0);
        h_R0->Fill(R0);
        h_sigma->Fill(sigma);
        h_epsilon->Fill(epsilon);

        if(Phi0 >= Phi0_max){
            Phi0_max = Phi0;
        }
        else if(Phi0 < Phi0_min){
            Phi0_min = Phi0;
        }

        if(R0 >= R0_max){
            R0_max = R0;
        }
        else if(R0 < R0_min){
            R0_min = R0;
        }

        if(sigma >= sigma_max){
            sigma_max = sigma;
        }
        else if(sigma < sigma_min){
            sigma_min = sigma;
        }

        if(epsilon >= epsilon_max){
            epsilon_max = epsilon;
        }
        else if(epsilon < epsilon_min){
            epsilon_min = epsilon;
        }


        for(int i = 0; i < 4001; i++){
            h2_curve->Fill(i*0.01,curve[i]);
        }
    }

    datatree->GetEntry(0);
    cout << "Phi0 Range: " << Phi0_min << " " << Phi0_max <<  " Best Fit: " << Phi0 << endl;
    cout << "R0 Range: " << R0_min << " " << R0_max <<  " Best Fit: " << R0 <<  endl;
    cout << "sigma Range: " << sigma_min << " " << sigma_max <<  " Best Fit: " << sigma <<  endl;
    cout << "epsilon Range: " << epsilon_min << " " << epsilon_max <<  " Best Fit: " << epsilon <<  endl;
    
    TGraphErrors *h_bestfit_error = new TGraphErrors();

    double alpha_h2o = 2.585e-2;
    double p_h2o = 1.738;
    double energy = 220;
    double alpha_pbwo4 = 7.275e-3;
    double p_pbwo4 = 1.690;
    double alpha_dsb = 1.030e-2;
    double p_dsb = 1.713;
    double alpha_ej256 = 2.393e-2;
    double p_ej256 = 1.742;
    double alpha_ej254 = 2.495e-2;
    double p_ej254 = 1.743;
    double alpha_ej212 = 2.483e-2;
    double p_ej212 = 1.743;
    double alpha_m;
    double p_m;

    p_m = p_pbwo4;
    alpha_m = alpha_pbwo4;
    double layerSize = 0.5;

    double wet_conv = CalcWETConv(p_h2o, p_m, alpha_h2o, alpha_m, energy);
    double wet_conv_sigma = CalcWETConvSigma(p_h2o, p_m, alpha_h2o, alpha_m, energy);

    for(int ch = 0; ch<15; ch++){
        datatree_error->GetEntry(ch);
        if(ch == 0){
            h_bestfit_error->SetPoint(ch, (layerSize/2+layerSize*ch)*wet_conv, mean);
        }
        else{
            h_bestfit_error->SetPoint(ch, (layerSize/2+layerSize*ch)*wet_conv, mean);
        }
        h_bestfit_error->SetPointError(ch, (layerSize/sqrt(12))*wet_conv,error);
    }

    // Define a TF1 function to represent the curve
    TGraph *h_bestfit = new TGraph();
    for(int i = 0; i<4000; i++){
        h_bestfit->SetPoint(i, i*0.01, curve[i]);
    }
    h_bestfit->SetLineColor(kRed);  // Set the curve color to red
    h_bestfit->SetLineWidth(1);
    h2_curve->SetStats(0);

    // Draw the curve over the histogram
    TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 4000, 2000);
    c1->Clear();
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    c1->Divide(3, 2);
    c1->cd(1);

    plotter.Histogram2D(h2_curve, "Depth / cm", "Norm. Energy Dose / MeV");
    h_bestfit->Draw("Same");

    h_bestfit_error->SetMarkerStyle(20);
    h_bestfit_error->SetMarkerSize(0.5);
    h_bestfit_error->SetMarkerColor(TColor::GetColor("#FF1493"));
    h_bestfit_error->SetLineColor(TColor::GetColor("#FF1493"));
    h_bestfit_error->Draw("P Same");

    TLine *line = new TLine(); // legend
    line->SetLineColor(kRed);

    c1->cd(2);
    plotter.Histogram1D(h_Phi0, "Phi0/N0 / cm^-2", "Counts");
    plotter.AddVerticalLine(Phi0, 0, h_Phi0->GetMaximum());
    plotter.Legend(h_Phi0, line, Phi0, "Phi0/N0", "cm^-2", 0.05);

    c1->cd(3);
    plotter.Histogram1D(h_R0, "R0 / cm", "Counts");
    plotter.AddVerticalLine(R0, 0, h_R0->GetMaximum());
    plotter.Legend(h_R0, line, R0, "R0", "cm", 0.05);

    c1->cd(4);
    plotter.Histogram1D(h_sigma, "sigma / cm", "Counts");
    plotter.AddVerticalLine(sigma, 0, h_sigma->GetMaximum());
    plotter.Legend(h_sigma, line, sigma, "sigma", "cm", 0.05);

    c1->cd(5);
    plotter.Histogram1D(h_epsilon, "epsilon", "Counts");
    gPad->SetLogy();
    plotter.AddVerticalLine(epsilon, 0, h_epsilon->GetMaximum());
    plotter.Legend(h_epsilon, line, epsilon, "epsilon", "", 0.05);

    c1->Update();
    gStyle->SetCanvasPreferGL(true);
    sprintf(pdf_target, "../data/%s/%s/output/pdf/output.pdf", dataset,filename);
    c1->SaveAs(pdf_target);

    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 2000);
    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plotter.Histogram2D(h2_curve, "Depth / cm", "Norm. Energy Dose / MeV");
    h_bestfit->Draw("Same");
    h_bestfit_error->Draw("P Same");
    sprintf(pdf_target, "../data/%s/%s/output/pdf/h2_curve.pdf", dataset, filename);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plotter.Histogram1D(h_Phi0, "Phi0/N0 / cm^-2", "Counts");
    plotter.AddVerticalLine(Phi0, 0, h_Phi0->GetMaximum());
    plotter.Legend(h_Phi0, line, Phi0, "Phi0/N0", "cm^-2", 0.05);
    sprintf(pdf_target, "../data/%s/%s/output/pdf/Phi0.pdf", dataset, filename);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plotter.Histogram1D(h_R0, "R0 / cm", "Counts");
    plotter.AddVerticalLine(R0, 0, h_R0->GetMaximum());
    plotter.Legend(h_R0, line, R0, "R0", "cm", 0.025);
    sprintf(pdf_target, "../data/%s/%s/output/pdf/R0.pdf", dataset, filename);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plotter.Histogram1D(h_sigma, "sigma / cm", "Counts");
    plotter.AddVerticalLine(sigma, 0, h_sigma->GetMaximum());
    plotter.Legend(h_sigma, line, sigma, "sigma", "cm", 0.025);
    sprintf(pdf_target, "../data/%s/%s/output/pdf/sigma.pdf", dataset, filename);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plotter.Histogram1D(h_epsilon, "epsilon", "Counts");
    gPad->SetLogy();
    plotter.AddVerticalLine(epsilon, 0, h_epsilon->GetMaximum());
    plotter.Legend(h_epsilon, line, epsilon, "epsilon", "", 0.05);
    sprintf(pdf_target, "../data/%s/%s/output/pdf/epsilon.pdf", dataset, filename);
    c2->SaveAs(pdf_target);
    c2->Close();
}