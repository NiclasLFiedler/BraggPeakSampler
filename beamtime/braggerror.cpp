#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"

const char* in_data[3] = {"beam2","homo-target", "hetero-target"};
const char* titles[3] = {"without a target","with the homogeneous target", "with the Heterogeneous target"};
const char* beamtime[1] = {"MIT_05_2024"};
int file_select = 1; 
int beamtime_select = 0; 

void draw_legend(TH1D* hist1){
    gPad->Update(); // Necessary to update the pad

    Char_t legendinfo[100];
    Double_t statsX1, statsX2, statsY1, statsY2;
    // Get the statistics box   
    TPaveStats *stats = (TPaveStats*)hist1->FindObject("stats");
    // Get the stats box dimensions
    statsX1 = stats->GetX1NDC();
    statsY1 = stats->GetY1NDC();
    statsX2 = stats->GetX2NDC();
    statsY2 = stats->GetY2NDC();
    
    // Create a legend
    TLegend *legend = new TLegend(statsX1, statsY2 - 0.2, statsX2, statsY2 - 0.2 - (statsY2 - statsY1));
    legend->AddEntry(hist1, "All Traces", "l");
    legend->Draw();
    return;
}

void draw_legend_fit(TH1D* hist, TLine* line, double bestfit, const char* param="param", const char* unit="unit", double offset=0){
    gPad->Update(); // Necessary to update the pad

    Double_t statsX1, statsX2, statsY1, statsY2;
    // Get the statistics box   
    TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
    // Get the stats box dimensions
    statsX1 = stats->GetX1NDC();
    statsY1 = stats->GetY1NDC();
    statsX2 = stats->GetX2NDC();
    statsY2 = stats->GetY2NDC();
    
    // Create a legend
    TLegend *legend = new TLegend(statsX1-offset, statsY2 - 0.2, statsX2, statsY2 - 0.15 - (statsY2 - statsY1));
    legend->SetTextSize(0.04);
    legend->AddEntry(line, Form("#splitline{Unperturbed fit:}{%s=%.3f %s}", param, bestfit, unit), "l");
    legend->Draw();
    return;
}

void addVerticalLine(double x, double ymin, double ymax) {
    TLine *line = new TLine(x, ymin, x, ymax);
    line->SetLineColor(kRed); // Set line color (optional)
    line->SetLineWidth(1); // Set line width (optional)
    line->Draw();
}

void plot_th2d(TH2D* hist, TString x_title, TString y_title){
    gPad->SetLogz();
    gPad->SetGrid();
    gStyle->SetTitleFontSize(0.07);
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.13);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->Update(); // Necessary to update the pad
    
    hist->GetYaxis()->SetTitleOffset(1.35);
    hist->GetXaxis()->SetTitle(x_title);
    hist->GetYaxis()->SetTitle(y_title);
    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelSize(0.07);
	hist->GetYaxis()->SetTitleSize(0.07);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelSize(0.07);
	hist->GetXaxis()->SetTitleSize(0.07);
    hist->GetZaxis()->SetLabelSize(0.07);
    hist->GetZaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->Draw("COLZ");
    gPad->Update(); // Necessary to update the pad
    
    //TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
    // Get the pointer to the stats box
    // Set the position of the stats box to top left
    // stats->SetX1NDC(0.15);
    // stats->SetX2NDC(0.45);
    // stats->SetY1NDC(0.65);
    // stats->SetY2NDC(0.85);
}

void plot_th1d(TH1D* hist, TString x_title, TString y_title){
    gPad->SetGrid();
    gStyle->SetTitleFontSize(0.07);
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);

    hist->GetYaxis()->SetNdivisions(506); // 5 primary ticks
    hist->GetYaxis()->SetTitle(y_title);
    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelSize(0.07);
	hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetTitleOffset(1);

    hist->GetXaxis()->SetTitle(x_title);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelSize(0.07);
	hist->GetXaxis()->SetTitleSize(0.07);
    hist->Draw("HIST");
}

void braggerror(){
    ROOT::EnableImplicitMT();
    Char_t histdesc[100];
    Char_t histname[100];
    Char_t title[100];
    Char_t file_curve[100];
    Char_t file_means[100];
    Char_t pdf_target[100];
    Char_t file[100];
    Char_t in_path[200];
    Char_t out_path[200];


    Float_t Phi0 = 0; //store tree values
    Float_t t = 0; //store tree values
    Float_t R0 = 0;
    Float_t sigma = 0;
    Float_t epsilon = 0;
    Float_t curve[4000];
    
    //init  histograms and file access
    //TFile *input = new TFile("bp-p-beam4.root", "READ");
    
    sprintf(file_curve, "%s/bp-p-%s/output/%s_fit.root", beamtime[beamtime_select], in_data[file_select], in_data[file_select]);
    sprintf(file_means, "%s/bp-p-%s/output/%sMeans.root", beamtime[beamtime_select], in_data[file_select], in_data[file_select]);

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

    if(file_select == 2){
        datatree->SetBranchAddress("t", &t);
        datatree->SetBranchAddress("sigma", &sigma);
        datatree->SetBranchAddress("curve", &curve);
    }
    else{
        datatree->SetBranchAddress("Phi0", &Phi0);
        datatree->SetBranchAddress("R0", &R0);
        datatree->SetBranchAddress("sigma", &sigma);
        datatree->SetBranchAddress("epsilon", &epsilon);
        datatree->SetBranchAddress("curve", &curve);  
    }
    Double_t mean;
    Double_t error;
    datatree_error->SetBranchAddress("mean", &mean);
    datatree_error->SetBranchAddress("error", &error);

    int max_heatmap = 120;
    Char_t h_title[200];
    sprintf(h_title, "Accumulated depth dose distribution %s", titles[file_select]);
    TH2D *h2_curve = new TH2D("h2_curve", h_title, 4000, 0, 40, 200, 0, max_heatmap); 
    sprintf(h_title, "Accumulated t parameters %s", titles[file_select]);
    TH1D *h_t = new TH1D("h_t", h_title, 100, 1, 2.5);
    sprintf(h_title, "Accumulated Phi0/N0 parameters %s", titles[file_select]);
    TH1D *h_Phi0 = new TH1D("h_Phi0", h_title, 100, 1, 2.5);
    sprintf(h_title, "Accumulated R0 parameters %s", titles[file_select]);
    TH1D *h_R0 = new TH1D("h_R0", h_title, 100, 28, 32);
    if (file_select == 1) {
        h_R0->SetBins(100, 24, 28); // modify the binning for "homo"
    }
    sprintf(h_title, "Accumulated sigma parameters %s", titles[file_select]);
    TH1D *h_sigma = new TH1D("h_sigma", h_title, 100, 0, 2);
    sprintf(h_title, "Accumulated epsilon parameters %s", titles[file_select]);
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
    double wet_conv = 4.603242332633419;
    for(int i = 0; i<15; i++){
        datatree_error->GetEntry(i);
        h_bestfit_error->SetPoint(i, (0.25+0.5*i)*wet_conv, mean);
        h_bestfit_error->SetPointError(i, (0.5/3.46410161514)*wet_conv, error); // 3.46410161514 = sqrt(12)
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

    plot_th2d(h2_curve, "Depth / cm", "Norm. Energy Dose / MeV");
    h_bestfit->Draw("Same");

    h_bestfit_error->SetMarkerStyle(20);
    h_bestfit_error->SetMarkerSize(0.5);
    h_bestfit_error->SetMarkerColor(TColor::GetColor("#FF1493"));
    h_bestfit_error->SetLineColor(TColor::GetColor("#FF1493"));
    h_bestfit_error->Draw("P Same");

    c1->cd(2);
    plot_th1d(h_Phi0, "Phi0/N0 / cm^-2", "Counts");
    TLine *line = new TLine(Phi0, 0, Phi0, h_Phi0->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_Phi0, line, Phi0, "Phi0/N0", "cm^-2", 0.05);

    c1->cd(3);
    plot_th1d(h_R0, "R0 / cm", "Counts");
    line = new TLine(R0, 0, R0, h_R0->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_R0, line, R0, "R0", "cm", 0.05);

    c1->cd(4);
    plot_th1d(h_sigma, "sigma / cm", "Counts");
    line = new TLine(sigma, 0, sigma, h_sigma->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_sigma, line, sigma, "sigma", "cm", 0.05);

    c1->cd(5);
    plot_th1d(h_epsilon, "epsilon", "Counts");
    gPad->SetLogy();
    line = new TLine(epsilon, 0, epsilon, h_epsilon->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_epsilon, line, epsilon, "epsilon", "", 0.05);

    c1->Update();
    gStyle->SetCanvasPreferGL(true);
    sprintf(pdf_target, "%s/bp-p-%s/output/output.pdf", beamtime[beamtime_select],in_data[file_select]);
    c1->SaveAs(pdf_target);

    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 2000);
    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plot_th2d(h2_curve, "Depth / cm", "Norm. Energy Dose / MeV");
    h_bestfit->Draw("Same");
    h_bestfit_error->Draw("P Same");
    sprintf(pdf_target, "%s/bp-p-%s/output/h2_curve.pdf", beamtime[beamtime_select], in_data[file_select]);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plot_th1d(h_Phi0, "Phi0/N0 / cm^-2", "Counts");
    line = new TLine(Phi0, 0, Phi0, h_Phi0->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_Phi0, line, Phi0, "Phi0/N0", "cm^-2", 0.05);
    sprintf(pdf_target, "%s/bp-p-%s/output/Phi0.pdf", beamtime[beamtime_select], in_data[file_select]);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plot_th1d(h_R0, "R0 / cm", "Counts");
    line = new TLine(R0, 0, R0, h_R0->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_R0, line, R0, "R0", "cm", 0.025);
    sprintf(pdf_target, "%s/bp-p-%s/output/R0.pdf", beamtime[beamtime_select], in_data[file_select]);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plot_th1d(h_sigma, "sigma / cm", "Counts");
    line = new TLine(sigma, 0, sigma, h_sigma->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_sigma, line, sigma, "sigma", "cm", 0.025);
    sprintf(pdf_target, "%s/bp-p-%s/output/sigma.pdf", beamtime[beamtime_select], in_data[file_select]);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plot_th1d(h_t, "t / cm", "Counts");
    line = new TLine(t, 0, t, h_t->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_t, line, t, "t", "cm", 0.025);
    sprintf(pdf_target, "%s/bp-p-%s/output/t.pdf", beamtime[beamtime_select], in_data[file_select]);
    c2->SaveAs(pdf_target);

    c2->Clear();
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);
    plot_th1d(h_epsilon, "epsilon", "Counts");
    gPad->SetLogy();
    line = new TLine(epsilon, 0, epsilon, h_epsilon->GetMaximum());
    line->SetLineColor(kRed);  // Match the color with your desired vertical line color
    line->Draw();
    draw_legend_fit(h_epsilon, line, epsilon, "epsilon", "", 0.05);
    sprintf(pdf_target, "%s/bp-p-%s/output/epsilon.pdf", beamtime[beamtime_select], in_data[file_select]);
    c2->SaveAs(pdf_target);
}