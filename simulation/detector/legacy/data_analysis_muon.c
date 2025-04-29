#include "TH1D.h"


void set_standard_pad(TH1D* h_1, TH1D* h_2 = nullptr, TH1D* h_3 = nullptr, TH1D* h_4 = nullptr){
    gPad->SetGrid();
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);

    //gPad->SetTitle(title);

    h_1->GetXaxis()->SetTitle("Energy / MeV");
    h_1->GetYaxis()->SetTitle("Counts");
    h_1->GetYaxis()->SetLabelFont(42);
	h_1->GetYaxis()->SetLabelSize(0.04);
	h_1->GetYaxis()->SetTitleSize(0.04);
	h_1->GetXaxis()->SetLabelFont(42);
	h_1->GetXaxis()->SetLabelSize(0.04);
	h_1->GetXaxis()->SetTitleSize(0.04);
    h_1->Draw("HIST");
    if(h_2 != nullptr){
        h_2->SetLineColor(kRed);	
	    h_2->Draw("HIST SAME");
    }
    if(h_3 != nullptr){
        h_3->SetLineColor(kGreen+2);	
	    h_3->Draw("HIST SAME");
    }
    if(h_4 != nullptr){
        h_4->SetLineColor(kOrange+1);	
	    h_4->Draw("HIST SAME");
    }
}

void simple_plot(TH1D* hist, TString x_title, TString y_title){
    gPad->SetGrid();
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);

    //gPad->SetTitle(title);

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

void data_analysis_muon(){
    ROOT::EnableImplicitMT();
    //init variables
    Int_t event = 0, NDet = 0; //store tree values
    Double_t EDep = 0;
    Double_t Ekin = 0;
    Double_t Etot = 0;
    Int_t PType, TrackID;
    const int nbofdetectors = 15; //nbofdetectors
    Char_t histname[100]; // temp store of histnames
    Char_t histname2[100]; // temp store of histnames

    Double_t run_edep[nbofdetectors] = {0}, total_edep = 0; //energy storage for hist fill
    Double_t total_c12[nbofdetectors] = {0};
    Double_t means[nbofdetectors] = {0};
    Double_t conv[nbofdetectors] = {0};
    //init  histograms and file access
    TFile *input = new TFile("raw_data_muon.root", "READ");
    TH1D *h_edep[nbofdetectors];//, h_edep_full_c12[nbofdetectors];
    TH1D *h_total_edep = new TH1D("h_total_edep", "Total Energy Deposition", 200, 0, 50);

    //init of histograms of edep and edepcut
    for(int i = 0; i<nbofdetectors; i++){
            sprintf(histname, "Energy deposition in detector layer %i", i);
            sprintf(histname2, "h_edep_%i", i);
        	h_edep[i] = new TH1D(histname2, histname, 500, 0, 50);
    }

    // Check if the file is successfully opened
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }

    //get data from tree
    TTree *datatree = (TTree*)input->Get("braggsampler");

    // Check if the tree is successfully retrieved
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

	datatree->SetBranchAddress("event", &event);
    datatree->SetBranchAddress("track", &TrackID);
	datatree->SetBranchAddress("NDet", &NDet);
	datatree->SetBranchAddress("EDep", &EDep);
    datatree->SetBranchAddress("Ekin", &Ekin);
    datatree->SetBranchAddress("Etot", &Etot);	
    datatree->SetBranchAddress("PType", &PType);	
    
    datatree->GetEntry(0);
    int steps = 0;
    int prevEvent = event;

    double entries = datatree->GetEntries();
    for (double i = 0; i<entries; i++){
		datatree->GetEntry(i);

        if(event == prevEvent){
            run_edep[NDet]+=EDep;
            total_edep+=EDep;
        }
        else if(event != prevEvent){
            if(total_edep!=0){
                h_total_edep->Fill(total_edep);    
            }
            for(int ii = 0; ii < nbofdetectors; ii++){
                if(run_edep[ii]!=0){
                    h_edep[ii]->Fill(run_edep[ii]);        
                }
                run_edep[ii] = 0;
            }
            run_edep[NDet] = EDep; 
			total_edep = EDep;
        }
        prevEvent = event;
    }

    //init canvas
    TCanvas *c1 = new TCanvas("c1","Bragg Sampler Simulation Analysis", 10, 10, 1900, 1000);
    //set canvas settings
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    c1->Divide(4, 4); 

    //Draw EDep
    for(int i = 0; i < nbofdetectors; i++){
        c1->cd(i+1); //change pad for drawing
        set_standard_pad(h_edep[i]);
    }

    //Draw Total EDep
    c1->cd(nbofdetectors+1);
    set_standard_pad(h_total_edep); 

    TFile *	hfile = new TFile("analysed_data_muon.root","RECREATE");
    c1->Write("ALL");

    TCanvas *c2 = new TCanvas("c2","Bragg Sampler Simulation Analysis", 10, 10, 1900, 1000);
    //set canvas settings
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);

    for(int i = 0; i < nbofdetectors; i++){
	    set_standard_pad(h_edep[i]);
        c2->SetTitle("Energy deposition");
        sprintf(histname2, "edep_%i", i);
        c2->Write("edep");
    }
    c2->Clear();
    set_standard_pad(h_total_edep); 
    c1->Write("total edep");
    
    hfile->SetCompressionLevel(0);
    hfile->Close();
}
