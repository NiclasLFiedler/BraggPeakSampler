#include "../include/Plotter.h"
#include <TPaveStats.h>
#include <TLegend.h>
#include <TCanvas.h>

Plotter::Plotter(int coincidence_time)
    : coincidence_time(coincidence_time) {}

void Plotter::Legend(TH1D* hist){
    gPad->Update(); // Necessary to update the pad

    Char_t legendinfo[100];
    Double_t statsX1, statsX2, statsY1, statsY2;
    TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");

    statsX1 = stats->GetX1NDC();
    statsY1 = stats->GetY1NDC();
    statsX2 = stats->GetX2NDC();
    statsY2 = stats->GetY2NDC();

    TLegend *legend = new TLegend(statsX1, statsY2 - 0.2, statsX2, statsY2 - 0.2 - (statsY2 - statsY1));
    legend->SetTextSize(0.05);
    
    sprintf(legendinfo, "#splitline{%i ns coinc.}{with ch. 0}", coincidence_time);
    legend->AddEntry(hist, legendinfo, "f1");
    legend->Draw();
    return;
}

void Plotter::Legend(TGraphErrors* graph){
    gPad->Update(); // Necessary to update the pad
    Char_t legendinfo[100];
    Double_t statsX1, statsX2, statsY1, statsY2;
    
    statsX1 = 0.1;
    statsY1 = 0.1;
    statsX2 = 0.2;
    statsY2 = 0.2;

    TLegend *legend = new TLegend(statsX1, statsY2 - 0.2, statsX2, statsY2 - 0.2 - (statsY2 - statsY1));
    legend->SetTextSize(0.05);
    
    sprintf(legendinfo, "#splitline{%i ns coinc.}{with ch. 0}", coincidence_time);
    legend->AddEntry(graph, legendinfo, "f1");
    legend->Draw();
    return;
}

void Plotter::Legend(TH1D* hist, TLine* line, double bestfit, const char* param="param", const char* unit="unit", double offset=0){
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

void Plotter::Histogram1D(TH1D* hist, TString x_title, TString y_title) {
    gPad->SetGrid();
    gPad->SetFrameLineColor(1);
    gPad->SetFrameLineWidth(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);

    hist->GetYaxis()->SetNdivisions(506);
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
    return;
}

void Plotter::Histogram2D(TH2D* hist, TString x_title, TString y_title){
    gPad->SetLogz();
    gPad->SetGrid();
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width    
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.13);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gStyle->SetTitleFontSize(0.07);
    
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
    gPad->Update();
    return;
}

void Plotter::GraphError(TGraphErrors* graph, TString x_title, TString y_title, TString title) {
    gStyle->SetEndErrorSize(5);
    gPad->SetGrid();
    
    graph->SetTitle(title);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->GetXaxis()->SetTitle(x_title);
    graph->GetYaxis()->SetTitle(y_title);
    graph->GetYaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetLabelSize(0.07);
	graph->GetYaxis()->SetTitleSize(0.07);
    graph->GetYaxis()->SetTitleOffset(1);
	graph->GetXaxis()->SetLabelFont(42);
	graph->GetXaxis()->SetLabelSize(0.07);
	graph->GetXaxis()->SetTitleSize(0.07);
    graph->SetLineColor(kOrange+1);
	graph->SetFillColor(kOrange+1);
    graph->Draw("AP");
    return;
}

void Plotter::AddVerticalLine(double x, double ymin, double ymax) {
    TLine *line = new TLine(x, ymin, x, ymax);
    line->SetLineColor(kRed);
    line->Draw();
    return;
}