#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>

struct trace_properties {    
    Double_t amp = 0;
    int amp_idx = 0;
    bool det_pileup = false;
    void clear(){
        amp = 0;
        amp_idx = 0;
        det_pileup = false;
    }
};

trace_properties get_max_trace(std::vector<double> trace){
    trace_properties trace_props;
    Double_t sum = 0;
    int pregate = 25;
    std::vector<Double_t>::iterator max_trace;
    for (int i = 0; i < pregate; ++i){
        sum += trace.at(i);
    }
    max_trace = std::min_element(trace.begin(), trace.end());
    
    trace_props.amp = static_cast<Double_t>(-*max_trace+sum/pregate);
    trace_props.amp_idx = static_cast<int>(std::distance(trace.begin(), max_trace));
    // trace_props.det_pileup = detectPileup(trace, 30);

    return trace_props;
}

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

void calib() {
    // Open the ROOT file
    Char_t filename[100];
    Char_t treeName[100];
    Char_t histdesc[100];
    Char_t histname[100];
    sprintf(filename, "na22_pos0.root");
    sprintf(treeName, "vtree");
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Get the TTree
    TTree* datatree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!datatree) {
        std::cerr << "Error: Could not find tree " << treeName << " in file " << filename << std::endl;
        file->Close();
        return;
    }
    Int_t nbofdetectors = 15;
    TH1D *h_ampl[nbofdetectors];


    for(int i = 0; i<nbofdetectors; i++){
        sprintf(histdesc, "layer %i", i);
        sprintf(histname, "h_ampl_%i", i);
        h_ampl[i] = new TH1D(histname, histdesc, 500, 0, 500);
    }

    Int_t event = 0; //store tree values
    Int_t bufferCounter = 0;
    Int_t channel = 0;
    uint32_t timestamp_ns = 0;
    Long64_t timestamp_ps = 0;
    std::vector<Double_t>* trace=nullptr;

    datatree->SetBranchAddress("EventCounter", &bufferCounter);
    datatree->SetBranchAddress("Channel", &channel);
    datatree->SetBranchAddress("TimeTag", &timestamp_ns);
    datatree->SetBranchAddress("TimeStampPico", &timestamp_ps);
    datatree->SetBranchAddress("Trace", &trace);

    Double_t entries = datatree->GetEntries();
    trace_properties trace_props;
    for (Long64_t i = 0; i < entries; ++i) {
        datatree->GetEntry(i);
        if(channel != 1) continue;
        trace_props = get_max_trace(*trace);
        h_ampl[channel]->Fill(trace_props.amp);
    }

    // Plot the histogram
    TCanvas *c1 = new TCanvas("c1", "Calibration", 10, 10, 1900, 1000);
    //set canvas settings
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    c1->Divide(4, 4);
    for(int i = 0; i<nbofdetectors; i++){
        c1->cd(i+1);
        plot_th1d(h_ampl[i], "Channel", "Counts");
    }
}

