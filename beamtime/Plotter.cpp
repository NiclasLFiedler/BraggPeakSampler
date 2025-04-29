#include "Plotter.h"
#include <TPaveStats.h>
#include <TLegend.h>
#include <TCanvas.h>

Plotter::Plotter(bool bhetero, int coincidence_time, double tot_edep_cut_off)
    : bhetero(bhetero), coincidence_time(coincidence_time), tot_edep_cut_off(tot_edep_cut_off) {}

void Plotter::draw_legend(TH1D* hist1, TH1D* hist2, TH1D* hist3, TH1D* hist4) {
    gPad->Update(); // Necessary to update the pad

    Char_t legendinfo[100];
    Double_t statsX1, statsX2, statsY1, statsY2;
    TPaveStats *stats = (TPaveStats*)hist1->FindObject("stats");

    statsX1 = stats->GetX1NDC();
    statsY1 = stats->GetY1NDC();
    statsX2 = stats->GetX2NDC();
    statsY2 = stats->GetY2NDC();

    TLegend *legend = new TLegend(statsX1, statsY2 - 0.2, statsX2, statsY2 - 0.2 - (statsY2 - statsY1));
    legend->SetTextSize(0.05);
    
    if (bhetero) {
        sprintf(legendinfo, "#splitline{%i ns coinc.}{with ch. 0}", coincidence_time);
        legend->AddEntry(hist1, legendinfo, "f1");
        legend->Draw();
        return;
    }

    sprintf(legendinfo, "#splitline{%i ns coinc. with ch. 0 &}{Tot_EDep >= %0.0f MeV}", coincidence_time, tot_edep_cut_off);
    legend->AddEntry(hist1, legendinfo, "f1");
    sprintf(legendinfo, "#splitline{%i ns coinc. with ch. 0 &}{Tot_EDep < %0.0f MeV}", coincidence_time, tot_edep_cut_off);
    legend->AddEntry(hist2, legendinfo, "f1");

    if (hist3) {
        legend->AddEntry(hist3, legendinfo, "f1");
    }

    if (hist4) {
        sprintf(legendinfo, "%i ns coinc. with ch. 1 & stopped particle", coincidence_time);
        legend->AddEntry(hist4, legendinfo, "f1");
    }

    legend->Draw();
}

void Plotter::draw_legend2(TH1D* hist1, TH1D* hist2) {
    gPad->Update();

    Char_t legendinfo[100];
    Double_t statsX1, statsX2, statsY1, statsY2;
    TPaveStats *stats = (TPaveStats*)hist1->FindObject("stats");

    statsX1 = stats->GetX1NDC();
    statsY1 = stats->GetY1NDC();
    statsX2 = stats->GetX2NDC();
    statsY2 = stats->GetY2NDC();

    TLegend *legend = new TLegend(statsX1, statsY2 - 0.2, statsX2, statsY2 - 0.2 - (statsY2 - statsY1));
    legend->SetTextSize(0.05);
    
    if (bhetero) {
        sprintf(legendinfo, "#splitline{%i ns coinc.}{with ch. 0}", coincidence_time);
        legend->AddEntry(hist1, legendinfo, "f1");
        legend->Draw();
        return;
    }

    sprintf(legendinfo, "#splitline{%i ns coinc. with ch. 0 &}{Tot_EDep >= %0.0f MeV}", coincidence_time, tot_edep_cut_off);
    legend->AddEntry(hist2, legendinfo, "f1");
    sprintf(legendinfo, "#splitline{%i ns coinc. with ch. 0 &}{Tot_EDep < %0.0f MeV}", coincidence_time, tot_edep_cut_off);
    legend->AddEntry(hist1, legendinfo, "f1");
    
    legend->Draw();
}

void Plotter::plot_th1d(TH1D* hist, TString x_title, TString y_title) {
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
}
