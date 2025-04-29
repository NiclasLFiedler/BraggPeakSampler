#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void substract_histograms() {
    // Open the ROOT files
    TFile* file1 = TFile::Open("co60_ch0_180s.root");
    TFile* file2 = TFile::Open("background_ch0_180s.root");

    // Check if files were opened successfully
    if (!file1 || !file2) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }
    const char* type[3] = {"h_ampl0", "h_integral0", "h_charge0"};
    // Get the histograms from the files
    TH1D* hist1 = dynamic_cast<TH1D*>(file1->Get(type[0]));
    TH1D* hist2 = dynamic_cast<TH1D*>(file2->Get(type[0]));

    // Check if histograms were retrieved successfully
    if (!hist1 || !hist2) {
        std::cerr << "Error retrieving histograms!" << std::endl;
        return;
    }

    // Create a histogram for the difference
    TH1D* histDiff = (TH1D*)hist1->Clone("histDiff");
    histDiff->SetTitle("Histogram Difference");
    histDiff->Add(hist2, -1); // Subtract hist2 from hist1

    // Create a canvas to draw the histograms
    TCanvas* canvas = new TCanvas("canvas", "Histogram Comparison", 800, 600);
    canvas->Divide(1, 3); // Divide canvas into 3 pads

    // Draw the first histogram
    canvas->cd(1);
    hist1->SetLineColor(kBlue);
    hist1->SetTitle("Histogram 1");
    hist1->Draw();

    // Draw the second histogram
    canvas->cd(2);
    hist2->SetLineColor(kRed);
    hist2->SetTitle("Histogram 2");
    hist2->Draw();

    // Draw the difference histogram
    canvas->cd(3);
    histDiff->SetLineColor(kGreen);
    histDiff->SetTitle("Difference (Hist1 - Hist2)");
    histDiff->Draw();

    // Update the canvas to show the histograms
    canvas->Update();
}

