#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>

void plot_th1d(TH1D* hist, TString x_title, TString y_title){
    gPad->SetGrid();
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);


    hist->GetXaxis()->SetTitle(x_title);
    hist->GetYaxis()->SetTitle(y_title);
    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetTitleSize(0.04);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetXaxis()->SetTitleSize(0.04);
    hist->Draw("HIST");
}

void draw_legend(TH1D* hist1, TH1D* hist2, TH1D* hist3=nullptr){
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
    legend->AddEntry(hist1, "beam2", "f1");
    //legend->AddEntry(hist2, "beam1", "f1");
    legend->AddEntry(hist2, "homogeneous target", "f1");
    if(hist3){
        //sprintf(legendinfo, "%i ns coinc. w/ ch. 1 & max Tot_EDep of %0.0f MeV", coincidence_time, tot_edep_cut_off);
        legend->AddEntry(hist3, "heterogeneous target", "f1");
    }

    legend->Draw();
    return;
}

void drawHistogram() {
    // Open the ROOT file
    TFile *file1 = TFile::Open("MIT_05_2024/bp-p-beam2/output/MIT_05_2024_beam2.root");
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    TFile *file2 = TFile::Open("MIT_05_2024/bp-p-beam1/output/MIT_05_2024_beam1.root");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    TFile *file3 = TFile::Open("MIT_05_2024/bp-p-homo-target/output/MIT_05_2024_homo.root");
    if (!file3 || file3->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    TFile *file4 = TFile::Open("MIT_05_2024/bp-p-hetero-target/output/MIT_05_2024_hetero.root");
    if (!file4 || file4->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }


    // Retrieve the histogram
    TH1D *h_beam2[15]; //green
    TH1D *h_beam1[15]; //red
    TH1D *h_homo[15];  //orange
    TH1D *h_hetero[15]; //blue
    Char_t filename[200];

    
    for(int i = 0; i<15; i++){
        sprintf(filename, "h_total_edep_cut_%i", i);
        file1->GetObject(filename, h_beam2[i]);
        sprintf(filename, "Energy Deposition of CH%i", i);
        h_beam2[i]->SetTitle(filename);
        h_beam2[i]->SetLineColor(kGreen+1);
        
        sprintf(filename, "h_total_edep_cut_%i", i);
        file2->GetObject(filename, h_beam1[i]);
        sprintf(filename, "Energy Deposition of CH%i", i);
        h_beam1[i]->SetTitle(filename);
        h_beam1[i]->SetLineColor(kRed+1);

        sprintf(filename, "h_total_edep_cut_%i", i);
        file3->GetObject(filename, h_homo[i]);
        sprintf(filename, "Energy Deposition of CH%i", i);
        h_homo[i]->SetTitle(filename);
        h_homo[i]->SetLineColor(kOrange+1);

        sprintf(filename, "h_total_edep_cut_%i", i);
        file4->GetObject(filename, h_hetero[i]);
        sprintf(filename, "Energy Deposition of CH%i", i);
        h_hetero[i]->SetTitle(filename);
        h_hetero[i]->SetLineColor(kBlue+1);
    }

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Amplitude Spectra", 800, 600);
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    // Draw the histogram

    int start_ch = 8;
    c1->Divide(3, 2);
    for(int i = 0; i < 6; i++){
        c1->cd(i+1);
        plot_th1d(h_beam2[i+start_ch], "EDep [MeV]", "Counts");
        h_beam2[i+start_ch]->GetYaxis()->SetRangeUser(0, h_beam1[i+start_ch]->GetMaximum()*1.1);
        h_beam1[i+start_ch]->Draw("HIST SAME");
        draw_legend(h_beam2[i+start_ch],h_beam1[i+start_ch]);
        c1->Update();
    }
    
    // int start_ch = 7;
    // c1->Divide(3, 2);
    // for(int i = 0; i < 6; i++){
    //     c1->cd(i+1);
    //     plot_th1d(h_beam2[i+start_ch], "EDep [MeV]", "Counts");
    //     //h_beam2[i+start_ch]->GetYaxis()->SetRangeUser(0, h_beam1[i+start_ch]->GetMaximum()*1.1);
    //     h_homo[i+start_ch]->Draw("HIST SAME");
    //     h_hetero[i+start_ch]->Draw("HIST SAME");
    //     draw_legend(h_beam2[i+start_ch], h_homo[i+start_ch], h_hetero[i+start_ch]);
    //     c1->Update();
    // }

    Char_t file[200];
    sprintf(file, "/home/niclas_ubuntu/master/presentation/amp_beam1_beam2.root");
    //sprintf(file, "/home/niclas_ubuntu/master/presentation/amp_beam2_homo_hetero.root");
    TFile *hfile = new TFile(file, "RECREATE");

    c1->Write("All");
    hfile->Close();
}

int draw_hist() {
    drawHistogram();
    return 0;
}
