#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

void plot_lightyield() {
    // Open the ROOT file
    const char* fileName = "BPS_CH10.root";
    
    // Open the ROOT file
    TFile *input = new TFile(fileName, "READ");
    if (!input || input->IsOpen() == false) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    //get data from tree
    TTree *datatree = (TTree*)input->Get("tree");

    // Check if the tree is successfully retrieved
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

    Int_t channel;
    datatree->Branch("channel", &channel);
    
    // Create a histogram
    TH1D* hist = new TH1D("hist", "Histogram from Branch Data", 2000, 0, 2000); // Adjust the bins and range as needed

    // Loop over entries in the tree
    Long64_t nEntries = datatree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        datatree->GetEntry(i);
        cout << "value " << channel << endl;
        hist->Fill(channel);
    }
    
    // Create a canvas to draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Histogram Canvas", 800, 600);

    // Draw the histogram
    hist->Draw("hist");

    // Save the canvas as an image file (optional)
    canvas->SaveAs("histogram_plot.png");
}
