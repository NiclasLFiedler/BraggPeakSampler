#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"
//I put everything together, take what you need. If you need anything else or if something is unclear I will provide additional information
//-----------------------------------------------------------------------------------//
//-----------------------------USE-OF-SCRIPT-----------------------------------------//
//-----------------------------------------------------------------------------------//
int coincidence_time = 20;              // set coincidence window in ns
double tot_edep_cut_off = 155;           // set total energy cut off  ----155 beam2 --75 beam1 -----130-hetero---140-homo
bool use_max = false;                    // enable (total_edep <= tot_edep_cut_off)
int discard_index = 55;                  // set index of prepeak step cut off
bool weighted = true;                   // enable weighted means
bool reversed = true;                    // enable reversed energy calibration
//-----------------------------------------------------------------------------------//
//-----------------------------DIR-STRUCTURE-----------------------------------------//
//-----------------------------------------------------------------------------------//
//----*some_dir*/edep_coincidence.cpp                                    // script
//----*some_dir*/*beamtimedir*/*beamrun_name*/input/*beamrun_name*.root  // input files
//----*some_dir*/*beamtimedir*/*beamrun_name*/output/                    // output dir
//-----------------------------------------------------------------------------------//
//-----------------------------INPUT-FILES-------------------------------------------//
const char* beamtime[3] = {"MIT_04_2024", "MIT_05_2024", "MIT_12_2024"};
//const char* in_data[9] = {"bp-p-beam1","bp-p-beam2", "bp-p-beam2.5","bp-p-beam3","bp-p-beam4","bp-p-target1","bp-p-target2", "bp-p-hetero-target", "bp-p-homo-target"}; // !!!!!!!!!!!!!! change to 2D array !!!!!!!!!!!!!!!!
const char* in_data[12] = {"notarget_low", "notarget_medium", "notarget_high", "homo_medium", "homo_high", "hetero_medium", "hetero_high"};
//-------0-------------1------------2------------3-------------4-------------5-------//
//   bp-p-beam1   bp-p-beam2.5  bp-p-beam3  bp-p-beam4   bp-p-target1   bp-p-target2
int file_select = 1;                    // select data file with index from above
int beamtime_select = 2;                // select beamtime
//-----------------------------------------------------------------------------------//
//-----------------------------REMARKS-----------------------------------------------//
//-----------------------------------------------------------------------------------//
// may need to change histogram range and binning
//
// histogram -> desc
// h_traceamp[i] -> all trace_amps of detector layer i
// h_traceamp_coinc[i] -> trace_amps with coincidence with ch0 of detector layer i
// h_total_edep_cut[i] -> trace_amps of h_traceamp_coinc with total energy requirements of < or > than 'tot_edep_cut_off' of detector layer i
// h_stopped[i] -> last trace_amps of h_traceamp_coinc - not necessarily stopped can also have escaped the detector of detector layer i
// h_means -> mean_values of h_traceamp[i]
// h_means_coinc -> mean_values of h_traceamp_coinc[i]
// h_means_edep_cut -> mean_values of h_total_edep_cut[i]
// h_total_edep -> total edep of h_traceamp_coinc[i]
// h_total_cut_edep -> total edep of h_total_edep_cut[i]
//-----------------------------------------------------------------------------------//

struct particle {    
    Double_t total_edep = 0; // total trace amp
    std::vector<Double_t> edep = std::vector<Double_t>(15, 0.0); //trace amp 
    std::vector<Double_t> time = std::vector<Double_t>(15, 0.0); // timestamp*2 -> ns
    bool det_prepeak_step = false; // flag for prepeak step detected
    bool det_pileup = false; // flag for pileup detection
    bool missing_buffer = false; // flag for missing buffer entry

    Double_t sum_edep(){
        total_edep = 0;
        for(int i=0; i<edep.size(); i++){
            total_edep += GetEDep(i);
        }
        return total_edep;
    }

    Double_t GetEDep(int channel){
        Double_t kB = 0.008694; //0.0333333
        return edep.at(channel)/(1-kB*edep.at(channel)/5);
    }

    void coincidence(Double_t event_time, Double_t event_edep, Double_t ampindex, bool detection, int channel){
        if (time.at(0) != 0) {
            if(std::abs(time.at(0) - event_time) <= coincidence_time){
                time.at(channel) = event_time;
                edep.at(channel) = event_edep;
                if(edep.at(channel-1)==0){
                    missing_buffer = true;
                }
                if(detection){
                    det_pileup = true;
                }
                if(ampindex>discard_index){
                    det_prepeak_step = true;
                }
            }
            return;
        }
    }
};

struct bufferdata {    
    std::vector<std::vector<Double_t>> edep = std::vector<std::vector<Double_t>>(15);
    std::vector<std::vector<int>> edep_ampindex = std::vector<std::vector<int>>(15);
    std::vector<std::vector<Double_t>> timestamp_ns = std::vector<std::vector<Double_t>>(15);
    std::vector<std::vector<Double_t>> timestamp_ps = std::vector<std::vector<Double_t>>(15);
    std::vector<std::vector<bool>> det_pileup = std::vector<std::vector<bool>>(15);
    void clear(){
        edep = std::vector<std::vector<Double_t>>(15); 
        edep_ampindex = std::vector<std::vector<int>>(15);
        timestamp_ns = std::vector<std::vector<Double_t>>(15);
        timestamp_ps = std::vector<std::vector<Double_t>>(15);
        det_pileup = std::vector<std::vector<bool>>(15);
    }
};

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

struct energy_ch{
    std::vector<Double_t> E = {};
    std::vector<Double_t> CH = {};

    Double_t GetEDep(int channel){
        Double_t kB = 0.008694; // 0.0333333
        return E.at(channel)/(1+kB*E.at(channel)/5);
    }
};

struct calibration{
    std::vector<Double_t> sl = std::vector<Double_t>(15, 0.0); // slope
    std::vector<Double_t> intrcpt = std::vector<Double_t>(15, 0.0); // intercept
};

void draw_legend(TH1D* hist1, TH1D* hist2, TH1D* hist3=nullptr, TH1D* hist4=nullptr){
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
    legend->AddEntry(hist1, "All Traces", "f1");
    sprintf(legendinfo, "%i ns coinc. w/ ch. 1 & w/o pileup, prepeak step", coincidence_time);
    legend->AddEntry(hist2, legendinfo, "f1");
    if(hist3){
        if(!use_max){
            sprintf(legendinfo, "%i ns coinc. w/ ch. 1 & min Tot_EDep of %0.0f MeV", coincidence_time, tot_edep_cut_off);
        }
        else{
            sprintf(legendinfo, "%i ns coinc. w/ ch. 1 & max Tot_EDep of %0.0f MeV", coincidence_time, tot_edep_cut_off);
        }
        legend->AddEntry(hist3, legendinfo, "f1");
    }
    if(hist4){
        sprintf(legendinfo, "%i ns coinc. w/ ch. 1 & stopped particle", coincidence_time);
        legend->AddEntry(hist4, legendinfo, "f1");
    }
    legend->Draw();
    return;
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

bool detectPileup(std::vector<double>& trace, double pulseDuration) {
    bool pileupDetected = false;
    double previousPeakTime = -1;

    for (int i = 1; i < trace.size() - 1; ++i) {
        if (trace[i] < trace[i-6] && trace[i] < trace[i+6] && (trace[i-6]-trace[i]) > 50 && (trace[i+6]-trace[i]) > 50){
            if (previousPeakTime != -1.0) {
                double timeSeparation = i - previousPeakTime;
                if (timeSeparation > pulseDuration) {
                    pileupDetected = true;
                    break;
                }
            }
            previousPeakTime = i;
        }
    }
    return pileupDetected;
}

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
    trace_props.det_pileup = detectPileup(trace, 30);

    return trace_props;
}

calibration energy_extrapolation(calibration calib){
    energy_ch muon;
    energy_ch co60;
    muon.CH = {970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8};
    muon.E = {4.886, 4.876, 4.873, 4.873, 4.87, 4.874, 4.876, 4.864, 4.886, 4.866, 4.869, 4.883, 4.893, 4.892, 4.876};

    co60.CH = {405, 410.2, 415.6, 447.3, 441.77, 438.24, 416.8, 48.9, 408.8, 422.34, 461.6, 51.59, 51.39, 49.26, 36.42};
    co60.E = {1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325, 1.3325};

    //Double_t kB = 0.0333333 //0.08694;
    //return edep.at(channel)/(1-kB*edep.at(channel)/5);

    for(int i = 0; i<muon.CH.size(); i++){
        if(reversed){
            calib.sl.at(i) = (muon.GetEDep(14-i) - co60.GetEDep(14-i)) / (muon.CH.at(14-i) - co60.CH.at(14-i));
            calib.intrcpt.at(i) = muon.GetEDep(14-i) - calib.sl.at(i) * muon.CH.at(14-i);
        }
        else{
            calib.sl.at(i) = (muon.GetEDep(i) - co60.GetEDep(i)) / (muon.CH.at(i) - co60.CH.at(i));
            calib.intrcpt.at(i) = co60.GetEDep(i) - calib.sl.at(i) * co60.CH.at(i);
        }
    }
    return calib;
}

double energy_calibration(double trace_amp, int channel, calibration calib){
    return trace_amp*calib.sl[channel];//+calib.intrcpt[channel];
}

void edep_old(){
    ROOT::EnableImplicitMT();
    const int nbofdetectors = 15;
    Char_t histdesc[100];
    Char_t histname[100];
    Char_t title[100];
    Char_t filename[100];
    Char_t file[100];
    Char_t in_path[200];
    Char_t out_path[200];

    Int_t event = 0; //store tree values
    Int_t bufferCounter = 0;
    Int_t channel = 0;
    uint32_t timestamp_ns = 0;
    Long64_t timestamp_ps = 0;
    std::vector<Double_t> *trace = 0;

    sprintf(filename, "%s", in_data[file_select]);
    sprintf(in_path, "%s/%s/input/", beamtime[beamtime_select], filename);
    sprintf(out_path, "%s/%s/output/", beamtime[beamtime_select], filename);
    sprintf(file, "%s%s.root", in_path, filename);
    
    //init  histograms and file access
    //TFile *input = new TFile("bp-p-beam4.root", "READ");
    TFile *input = new TFile(file, "READ");
    // Check if the file is successfully opened
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }

    //get data from tree
    TTree *datatree = (TTree*)input->Get("vtree");

    // Check if the tree is successfully retrieved
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

	//datatree->SetBranchAddress("EventCounter", &event);
    datatree->SetBranchAddress("EventCounter", &bufferCounter);
    datatree->SetBranchAddress("Channel", &channel);
    datatree->SetBranchAddress("TimeTag", &timestamp_ns);
    datatree->SetBranchAddress("TimeStampPico", &timestamp_ps);
    datatree->SetBranchAddress("Trace", &trace);	  

    TH1D *h_traceamp[nbofdetectors];
    TH1D *h_traceamp_coinc[nbofdetectors];
    TH1D *h_total_edep_cut[nbofdetectors];
    TH1D *h_stopped[nbofdetectors];
    
    if(weighted){
        sprintf(histdesc, "Weighted mean energy deposition per detector layer");
    }
    else{
        sprintf(histdesc, "Mean energy deposition per detector layer");
    }

    TH1D *h_means = new TH1D("h_means", histdesc, nbofdetectors+1, 0, nbofdetectors+1);
    TH1D *h_means_coinc = new TH1D("h_means_coinc", histdesc, nbofdetectors+1, 0, nbofdetectors+1);
    TH1D *h_means_edep_cut = new TH1D("h_means_edep_cut", histdesc, nbofdetectors+1, 0, nbofdetectors+1);
    sprintf(histdesc, "Total energy deposition");
    TH1D *h_total_edep = new TH1D("h_total_edep", histdesc, 270, 0, 270);
    TH1D *h_total_cut_edep = new TH1D("h_h_total_cut_edep", histdesc, 270, 0, 270);

    for(int i = 0; i<nbofdetectors; i++){
        sprintf(histdesc, "Energy deposition of CH%i and coincidence with first detector layer", i);
        sprintf(histname, "h_traceamp_%i", i);
        h_traceamp[i] = new TH1D(histname, histdesc, 240, 0, 60);
        // h_traceamp[i] = new TH1D(histname, histdesc, 1000, 0, 8000);
        sprintf(histname, "h_traceamp_coinc_%i", i);
        h_traceamp_coinc[i] = new TH1D(histname, histdesc, 240, 0, 60);
        // h_traceamp_coinc[i] = new TH1D(histname, histdesc, 1000, 0, 8000);
        sprintf(histname, "h_total_edep_cut_%i", i);
        sprintf(histdesc, "Energy deposition of CH%i with minimum total energy of %0.0f MeV", i, tot_edep_cut_off);
        h_total_edep_cut[i] = new TH1D(histname, histdesc, 240, 0, 60);
        // h_total_edep_cut[i] = new TH1D(histname, histdesc, 100, 0, 8000);
        sprintf(histname, "h_stopped_%i", i);
        sprintf(histdesc, "Energy deposition of stopped particle in CH%i", i);
        h_stopped[i] = new TH1D(histname, histdesc, 240, 0, 60);
        // h_stopped[i] = new TH1D(histname, histdesc, 1000, 0, 8000);
    }

    Double_t entries = datatree->GetEntries();

    std::vector<particle> particles;
    bufferdata buffer;
    trace_properties trace_props;

    calibration calib; // init calib
    calib = energy_extrapolation(calib); // get calib

    for(int i = 0; i < 15;i++){
        cout << i << " sl: " << calib.sl.at(i) << " int: " << calib.intrcpt.at(i) << endl;
    }

    Int_t prevbufferCounter = 0;
    Int_t prevbuffer = 0;
    bool epoch_skip = false;
    int detected_missing_buffer = 0;
    int detected_pileup = 0;
    int detected_prepeak_step = 0;
    Int_t start_channel = -1;


    cout << "Beamtime: " << beamtime[beamtime_select] << endl;
    cout << "Input file: " << filename << endl;
    for (double e = 0; e<entries; e++){
    	datatree->GetEntry(e);
        trace_props = get_max_trace(*trace);
        trace_props.amp = energy_calibration(trace_props.amp, channel-1, calib);
        h_traceamp[channel-1]->Fill(trace_props.amp/(1-0.008694*trace_props.amp/5));

        if(prevbufferCounter == 60 && prevbuffer != bufferCounter){
            prevbufferCounter = 0;
            for(int i = 0; i < buffer.edep.size(); i++){
                if(buffer.edep.at(i).size()>0){
                    particles.resize(buffer.edep.at(i).size());                    
                    start_channel = i;
                    break;
                }
            }
            for(int i=0; i<buffer.edep.at(start_channel).size();i++){//--------------- start data of epoch
                particles.at(i).edep.at(start_channel) = buffer.edep.at(start_channel).at(i);
                particles.at(i).time.at(start_channel) = buffer.timestamp_ns.at(start_channel).at(i);
            }        
            
            for(int i=start_channel+1; i<buffer.edep.size(); i++){ // -------------------- get coinc data
                for(int p=0; p<particles.size(); p++){
                    for(int j=0; j<buffer.edep.at(i).size(); j++){ //iterate over <double>max_amp & timestamp of one channel
                        particles.at(p).coincidence(buffer.timestamp_ns.at(i).at(j), buffer.edep.at(i).at(j), buffer.edep_ampindex.at(i).at(j), buffer.det_pileup.at(i).at(j), i);
                    }
                }
            }
            for(int p=0; p<particles.size(); p++){ //--------------------------------- fill h_total_edep & h_traceamp_coinc
                if(particles.at(p).missing_buffer == true){
                    detected_missing_buffer++;
                    continue;
                }
                if(particles.at(p).det_pileup == true){
                    detected_pileup++;
                    continue;
                }
                if(particles.at(p).det_prepeak_step == true){
                    detected_prepeak_step++;
                    continue;
                }

                if(particles.at(p).edep.at(start_channel)!=0){                      // fill total_edep
                    particles.at(p).sum_edep();
                    h_total_edep->Fill(particles.at(p).total_edep);
                    if(!use_max){
                        if(particles.at(p).total_edep>=tot_edep_cut_off){    // min_tot_edep                             
                            h_total_cut_edep->Fill(particles.at(p).total_edep);
                        }
                    }
                    else{
                        if(particles.at(p).total_edep<=tot_edep_cut_off){    // max_tot_edep                             
                            h_total_cut_edep->Fill(particles.at(p).total_edep);
                        }
                    }
                    for(int i=start_channel; i<buffer.edep.size(); i++){
                        if(particles.at(p).edep.at(i)!=0){
                            h_traceamp_coinc[i]->Fill(particles.at(p).GetEDep(i));
                            if(!use_max){
                                if(particles.at(p).total_edep>=tot_edep_cut_off){    // min_tot_edep 
                                    h_total_edep_cut[i]->Fill(particles.at(p).GetEDep(i));                                    
                                }
                            }
                            else{ 
                                if(particles.at(p).total_edep<=tot_edep_cut_off){    // max_tot_edep 
                                    h_total_edep_cut[i]->Fill(particles.at(p).GetEDep(i));                                    
                                }
                            }
                        }
                    }
                    for (int i = particles.at(p).edep.size() - 1; i >= 0; --i) { //--------- fill h_stopped
                        if(particles.at(p).edep.at(i) != 0){
                            h_stopped[i]->Fill(particles.at(p).GetEDep(i));
                            break;
                        }
                    }
                }
            }            
            particles.clear();
            buffer.clear();
        }
        buffer.edep.at(channel-1).push_back(trace_props.amp);
        buffer.edep_ampindex.at(channel-1).push_back(trace_props.amp_idx);
        buffer.timestamp_ns.at(channel-1).push_back(static_cast<double>(timestamp_ns)*2);
        buffer.timestamp_ps.at(channel-1).push_back(static_cast<double>(timestamp_ps));
        if(trace_props.det_pileup == true){
            buffer.det_pileup.at(channel-1).push_back(true);
        }
        else{
            buffer.det_pileup.at(channel-1).push_back(false);
        }
        if(prevbuffer != bufferCounter){
            prevbufferCounter++;
        }
        prevbuffer = bufferCounter;
        trace_props.clear();
    }
    
    //init canvas
    //TCanvas *c1 = new TCanvas("c1","Bragg Sampler MIT Analysis", 10, 10, 1000, 500);
    sprintf(title, "Bragg Sampler %s Analysis", filename);
    TCanvas *c1 = new TCanvas("c1", title, 10, 10, 1900, 1000);
    //set canvas settings
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    c1->Divide(5, 4);

    sprintf(histdesc, " w/ total EDep cutoff %0.0f MeV", tot_edep_cut_off);
    // Header
    cout << setfill(' ');
    cout << setw(10) << "Channel" << setw(2) << "|" << setw(15) << "Unfiltered" << setw(2) << "|"  << setw(25) << "Coinc. /w channel 1" << setw(2) << "|"  << setw(30) << histdesc << setw(2) << "|"  << setw(18) << " Stopped particle" << endl;
    cout << setfill('-') << setw(106) << "-" << endl;
    cout << setfill(' ');

    for(int i = 0; i<nbofdetectors; i++){
        c1->cd(i+1);
        plot_th1d(h_traceamp[i], "EDep [MeV]", "Counts");
        h_traceamp_coinc[i]->SetLineColor(kRed+1);
        h_traceamp_coinc[i]->Draw("HIST SAME");
        h_total_edep_cut[i]->SetLineColor(kGreen+1);
        h_total_edep_cut[i]->Draw("HIST SAME");
        h_stopped[i]->SetLineColor(kOrange+1);
        h_stopped[i]->Draw("HIST SAME");
        draw_legend(h_traceamp[i], h_traceamp_coinc[i], h_total_edep_cut[i], h_stopped[i]);

        cout << setw(10) << i << setw(2) << "|" << setw(15) << h_traceamp[i]->GetEntries() << setw(2) << "|" << setw(25) << h_traceamp_coinc[i]->GetEntries() << setw(2) << "|" << setw(30) << h_total_edep_cut[i]->GetEntries() << setw(2) << "|" << setw(18) <<  h_stopped[i]->GetEntries() << endl;
    }
    cout << setfill('-') << setw(106) << "-" << endl;
    cout << setfill(' ');
    cout << "Detected missing buffer entries: " << detected_missing_buffer << endl;
    cout << "Detected pileups: " << detected_pileup << endl;
    cout << "Detected prepeak steps: " << detected_prepeak_step << endl;

    c1->cd(nbofdetectors+1);

    Double_t mean_entries = 0;
    Double_t mean_coinc_entries = 0;
    Double_t mean_cut_entries = 0;
    for(int i = 1; i < nbofdetectors; i++){
        mean_entries+=h_traceamp[i]->GetMean();
        mean_coinc_entries+=h_traceamp_coinc[i]->GetMean();
        mean_cut_entries+=h_total_edep_cut[i]->GetMean();
    }
    for(int i = 1; i < nbofdetectors; i++){
        if(weighted){
            //weighted means
	        h_means->Fill(i, h_traceamp[i]->GetMean()*h_traceamp[i]->GetEntries());
            h_means_coinc->Fill(i, h_traceamp_coinc[i]->GetMean()*h_traceamp_coinc[i]->GetEntries());
            h_means_edep_cut->Fill(i, h_total_edep_cut[i]->GetMean()*h_total_edep_cut[i]->GetEntries());
        }
        else{
            h_means->Fill(i, h_traceamp[i]->GetMean());
            h_means_coinc->Fill(i, h_traceamp_coinc[i]->GetMean());
            h_means_edep_cut->Fill(i, h_total_edep_cut[i]->GetMean());
        }
    }

    gPad->SetGrid();
    h_means->GetXaxis()->SetTitle("Detector layer");
    h_means->GetYaxis()->SetTitle("Total Energy Dose [MeV]");
    if(h_means_coinc->GetMaximum() > h_means->GetMaximum()){
        h_means->GetYaxis()->SetRangeUser(0, h_means_coinc->GetMaximum()*1.1);
    }
    h_means->SetLineColor(kBlack);
	h_means->SetFillColor(kBlack);
    h_means->SetBarWidth(0.2); // Adjust width of bars
    h_means->SetBarOffset(-0.4); // Offset for first histogram
    h_means->Draw("bar hist");
    h_means_coinc->SetLineColor(kRed+1);
	h_means_coinc->SetFillColor(kRed+1);
    h_means_coinc->SetBarWidth(0.2); // Adjust width of bars
    h_means_coinc->SetBarOffset(-0.2); // Offset for second histogram
    h_means_coinc->Draw("bar hist same");
    h_means_edep_cut->SetLineColor(kGreen+1);
	h_means_edep_cut->SetFillColor(kGreen+1);
    h_means_edep_cut->SetBarWidth(0.2); // Adjust width of bars
    h_means_edep_cut->SetBarOffset(0); // Offset for second histogram
    h_means_edep_cut->Draw("bar hist same");
    draw_legend(h_means, h_means_coinc, h_means_edep_cut);

    c1->cd(nbofdetectors+2);
    plot_th1d(h_total_edep, "EDep [MeV]", "Counts");
    h_total_cut_edep->SetLineColor(kGreen+1);
    h_total_cut_edep->Draw("HIST SAME");
    
    //outputfile
    if(!use_max){
        sprintf(file, "%sedep_%s_min.root", out_path, filename);
    }
    else{
        sprintf(file, "%sedep_%s_max.root", out_path, filename);
    }
    TFile *	hfile = new TFile(file,"RECREATE");

    c1->Write("ALL");

    //save analysised data
    TCanvas *c2 = new TCanvas("c2", title, 10, 10, 1900, 1000);
    //set canvas settings
    c2->SetFillColor(0);
	c2->SetGrid();
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetFrameBorderMode(0);

    for(int i = 0; i < nbofdetectors; i++){
        c2->cd(i+1);
	    plot_th1d(h_traceamp[i], "EDep [MeV]", "Counts");
        if(i != 0){
            h_traceamp_coinc[i]->Draw("HIST SAME");
            h_total_edep_cut[i]->Draw("HIST SAME");
            h_stopped[i]->Draw("HIST SAME");
        }

        draw_legend(h_traceamp[i], h_traceamp_coinc[i], h_total_edep_cut[i], h_stopped[i]);

        sprintf(histname, "h_edep_%i", i);
        c2->Write(histname);
    }
    c2->cd(nbofdetectors+1);    
    h_means->Draw("bar hist");
    h_means_coinc->Draw("bar hist same");
    h_means_edep_cut->Draw("bar hist same");
    draw_legend(h_means, h_means_coinc, h_means_edep_cut);
    c2->Write("mean_edep");

    c2->cd(nbofdetectors+2);
    plot_th1d(h_total_edep, "EDep [MeV]", "Counts");
    h_total_cut_edep->Draw("HIST SAME");

    c2->Write("total_edep");

    hfile->Close();

    sprintf(file, "%sMIT_05_2024_%s.root", out_path, in_data[file_select]);
    hfile = new TFile(file,"RECREATE");

    //save analysised data
    TCanvas *c3 = new TCanvas("c3", title, 10, 10, 1900, 1000);
    //set canvas settings
    c3->SetFillColor(0);
    c3->SetGrid();
	c3->SetBorderMode(0);
	c3->SetBorderSize(2);
	c3->SetFrameBorderMode(0);

    h_means->Draw("bar hist");
    h_means_coinc->Draw("bar hist same");
    h_means_edep_cut->Draw("bar hist same");
    draw_legend(h_means, h_means_coinc, h_means_edep_cut);
    c3->Write("mean_edep");
        
    plot_th1d(h_total_edep, "EDep [MeV]", "Counts");
    h_total_cut_edep->Draw("HIST SAME");
    c3->Write("total_edep");

    c3->Clear();
    
    c3->SetFillColor(0);
	c3->SetGrid();
	c3->SetBorderMode(0);
	c3->SetBorderSize(2);
	c3->SetFrameBorderMode(0);
    c3->Divide(2, 2);
    int start = 10;
    int end = 1;
    for(int i = start; i < nbofdetectors-end; i++){
        c3->cd(i-start+1);
	    plot_th1d(h_traceamp[i], "EDep [MeV]", "Counts");
        if(i != 0){
            h_traceamp_coinc[i]->Draw("HIST SAME");
            h_total_edep_cut[i]->Draw("HIST SAME");
            h_stopped[i]->Draw("HIST SAME");
        }
        draw_legend(h_traceamp[i], h_traceamp_coinc[i], h_total_edep_cut[i], h_stopped[i]);

    }
    sprintf(histname, "edep_channels");
    c3->Write(histname);

    for(int i = 0; i < nbofdetectors; i++){
        h_total_edep_cut[i]->Write();
    }
    hfile->Close();

    //store hist means for bortfield fit
    const char* in_data2[9] = {"beam1","beam2", "beam2.5","beam3","beam4","target1","target2", "hetero-target", "homo-target"}; // !!!!!!!!!!!!!! change to 2D array !!!!!!!!!!!!!!!!
    sprintf(file, "%s%sMeans.root", out_path, in_data2[file_select]);
    hfile = new TFile(file, "RECREATE");
    
    // Create a TTree
    TTree *meantree = new TTree("meantree", "Tree storing histogram means");

    // Variable to store histogram means
    Int_t ch;
    Double_t mean;
    Double_t error;

    // Create a branch in the TTree
    meantree->Branch("ch", &ch);
    meantree->Branch("mean", &mean);
    meantree->Branch("error", &error);

    for (int i = 0; i < nbofdetectors; ++i) {
        ch = i;
        mean = h_total_edep_cut[i]->GetMean()*h_total_edep_cut[i]->GetEntries();
        error = h_total_edep_cut[i]->GetStdDev()*h_total_edep_cut[i]->GetEntries();
        
        meantree->Fill();        
    }
    meantree->Write();    
    hfile->Close();
}