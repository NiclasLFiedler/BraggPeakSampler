#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"

//-----------------------------------------------------------------------------------//
//-----------------------------USE-OF-SCRIPT-----------------------------------------//
//-----------------------------------------------------------------------------------//
int coincidence_time = 160;              // set coincidence window in ns
double tot_edep_cut_off = 155;           // set total energy cut off
bool use_max = false;                    // enable (total_edep <= tot_edep_cut_off)
int discard_index = 50;                  // set index of prepeak step cut off
bool weighted = false;                   // enable weighted means
bool reversed = false;                    // enable reversed energy calibration
//-----------------------------------------------------------------------------------//
//-----------------------------DIR-STRUCTURE-----------------------------------------//
//-----------------------------------------------------------------------------------//
//----*some_dir*/edep_coincidence.cpp                                    // script
//----*some_dir*/*beamtimedir*/*beamrun_name*/input/*beamrun_name*.root  // input files
//----*some_dir*/*beamtimedir*/*beamrun_name*/output/                    // output dir
//-----------------------------------------------------------------------------------//
//-----------------------------INPUT-FILES-------------------------------------------//
const char* beamtime[4] = {"MIT_04_2024", "MIT_05_2024", "MIT_11_2024", "MIT_12_2024"};
const char* in_data[12] = {"bp-p-beam1","bp-p-beam2", "bp-p-beam2.5","bp-p-beam3","bp-p-beam4","bp-p-target1","bp-p-target2", "bp-p-hetero-target", "bp-p-homo-target", "notarget_run1", "run3_notarget", "run6_notarget"}; // !!!!!!!!!!!!!! change to 2D array !!!!!!!!!!!!!!!!
// const char* in_data[12] = {"notarget_low", "notarget_medium", "notarget_high", "homo_medium", "homo_high", "hetero_medium", "hetero_high"};
//-------0-------------1------------2------------3-------------4-------------5-------//
//   bp-p-beam1   bp-p-beam2.5  bp-p-beam3  bp-p-beam4   bp-p-target1   bp-p-target2
int file_select = 1;                    // select data file with index from above
int beamtime_select = 1;                // select beamtime
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

struct trace_properties {   
    std::vector<Double_t> trace;
    Double_t amp = 0;
    Double_t energy = 0;
    Double_t energy_err = 0;
    Double_t time_ns = 0;
    Double_t time_ps = 0;
    int amp_idx = 0;
    bool det_pileup = false;
    
    void clear(){;
        amp = 0;
        energy = 0;
        energy_err = 0;
        amp_idx = 0;
        time_ns = 0;
        time_ps = 0;
        det_pileup = false;
    }
};

struct particle {
    Double_t total_edep = 0; 
    Double_t total_edep_err = 0; 
    std::vector<std::vector<Double_t>> trace = std::vector<std::vector<Double_t>>(15);
    std::vector<Double_t> amplitude = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> edep = std::vector<Double_t>(15, 0.0); 
    std::vector<Double_t> edep_err = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> time = std::vector<Double_t>(15, 0.0); //ns

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
        // Double_t kB = 0.008694; //0.0333333
        return edep.at(channel); // /(1-kB*edep.at(channel)/5);
    }

    void coincidence(trace_properties traceTest, int channel){
        if (time.at(0) != 0) {
            if(std::abs(time.at(0) - traceTest.time_ns) <= coincidence_time){
                time.at(channel) = traceTest.time_ns;
                edep.at(channel) = traceTest.energy;
                edep_err.at(channel) = traceTest.energy_err;
                amplitude.at(channel) = traceTest.amp;
                trace.at(channel) = traceTest.trace;
                det_pileup = traceTest.det_pileup;
                if(traceTest.amp_idx>discard_index){
                    det_prepeak_step = true;
                }
            }
            return;
        }
        return;
    }

    void TestBuffer(){
        bool emptyChannel = false;
        for(int channel = 0; channel<edep.size(); channel++){
            if(edep.at(channel) != 0){
                if(emptyChannel == false){
                    continue;
                }
                else{
                    missing_buffer = true;
                }
            }
            else{
                emptyChannel = true;
            }
        }
    }
};

struct bufferdata{
    std::vector<std::vector<trace_properties>> trace_props = std::vector<std::vector<trace_properties>>(15);
    void clear(){
        trace_props = std::vector<std::vector<trace_properties>>(15);
    }
};

struct energy_ch{
    Double_t E;
    Double_t o_E;
    std::vector<Double_t> CH = {};
    std::vector<Double_t> o_CH = {};
};

struct calibration{
    std::vector<Double_t> slope = std::vector<Double_t>(15, 0.0);
    Double_t epsilon_1 = 0;
    Double_t epsilon_2 = 0;
    Double_t o_epsilon_1 = 0;
    Double_t o_epsilon_2 = 0;
    std::vector<Double_t> chi_1 = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> chi_2 = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> o_chi_1 = std::vector<Double_t>(15, 0.0);
    std::vector<Double_t> o_chi_2 = std::vector<Double_t>(15, 0.0);
    double kB = 0.008694;
    double d = 5;

    double GetVar(int channel, double e_i){
        double ch = (e_i-epsilon_1)*(chi_1.at(channel)-chi_2.at(channel))/(epsilon_1-epsilon_2)+chi_1.at(channel);
        double sum1 = (epsilon_1-epsilon_2)*((chi_2.at(channel)-chi_1.at(channel))+(chi_1.at(channel)-ch))/(std::pow((chi_1.at(channel)-chi_2.at(channel)),2))*o_chi_1.at(channel);
        double sum2 = (epsilon_1-epsilon_2)*(ch-chi_1.at(channel))/(std::pow((chi_1.at(channel)-chi_2.at(channel)),2))*o_chi_2.at(channel);
        double sum3 = (1+(ch-chi_1.at(channel))/(chi_1.at(channel)-chi_2.at(channel)))*o_epsilon_1;
        double sum4 = (ch-chi_1.at(channel))/(chi_1.at(channel)-chi_2.at(channel))*o_epsilon_2;
        return std::pow(sum1,2)+std::pow(sum2,2)+std::pow(sum3,2)+std::pow(sum4,2);
    }
    double GetQuenchedVar(int channel, double energy){
        return  GetVar(channel, energy)/std::pow(1 - kB*energy/d, 4); // quenched_sigma^2 = sigma^2/(1-kb E /d)^4
    }

    double GetQuenchedEnergy(int channel, double trace_amp){
        return epsilon_1+slope.at(channel)*(trace_amp-chi_1.at(channel));
    }

    double GetEnergy(double QuenchedEnergy){
        return QuenchedEnergy/(1-kB*QuenchedEnergy/d);
    }

    double FromEnergyToQuenched(double Energy){
        return Energy/(1+kB*Energy/d);
    }

    double GetQuenchedStdDev(int channel, double Energy){
        return sqrt(GetVar(channel, FromEnergyToQuenched(Energy)))/(std::pow(1-kB*FromEnergyToQuenched(Energy)/d,2));
    }
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

void plot_th2d(TH2D* hist, TString x_title, TString y_title){    
    gPad->SetGrid();
    gPad->SetLogz();
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);

    //gPad->SetTitle(title);
    hist->GetYaxis()->SetTitleOffset(1.35);
    hist->GetXaxis()->SetTitle(x_title);
    hist->GetYaxis()->SetTitle(y_title);
    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetTitleSize(0.04);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetXaxis()->SetTitleSize(0.04);
    hist->Draw("colz");
    gPad->Update();
    TPaveStats *stats = (TPaveStats*) hist->FindObject("stats");
    // Set the position of the stats box: (x1, y1, x2, y2)
    stats->SetX1NDC(0.7); // Position X (0 to 1, relative to the canvas size)
    stats->SetY1NDC(0.15); // Position Y (0 to 1, relative to the canvas size)
    stats->SetX2NDC(0.9); // Width of the box
    stats->SetY2NDC(0.30); // Height of the box
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

    hist->GetYaxis()->SetTitleOffset(1.35);
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

calibration energy_extrapolation(calibration calib, energy_ch source1, energy_ch source2){
    int j = 0;
    if(reversed) j = 14;
    for(int i = 0; i<source1.CH.size(); i++){
            calib.epsilon_1 = source1.E;
            calib.o_epsilon_1 = source1.o_E;
            calib.epsilon_2 = source2.E;
            calib.o_epsilon_2 = source2.o_E;
            calib.chi_1.at(i) = source1.CH.at(std::abs(j-i));
            calib.o_chi_1.at(i) = source1.o_CH.at(std::abs(j-i));
            calib.chi_2.at(i) = source2.CH.at(std::abs(j-i));
            calib.o_chi_2.at(i) = source2.o_CH.at(std::abs(j-i));
            calib.slope.at(i) = (source2.E - source1.E) / (source2.CH.at(std::abs(j-i)) - source1.CH.at(std::abs(j-i)));
    }
    return calib;
}

void all_trace_analysis(){
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
    Int_t eventCounter = 0;
    Int_t channel = 0;
    uint32_t timestamp_ns = 0;
    Long64_t timestamp_ps = 0;
    std::vector<std::vector<Double_t>> *trace = 0;

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
    datatree->SetBranchAddress("BufferCounter", &eventCounter);
    datatree->SetBranchAddress("Channel", &channel);
    datatree->SetBranchAddress("TimeTag", &timestamp_ns);
    datatree->SetBranchAddress("TimeStampPico", &timestamp_ps);
    datatree->SetBranchAddress("Trace", &trace);	  

    TH2D *h2_traces[nbofdetectors];
    TH2D *h2_traces_clean1[nbofdetectors];
    TH2D *h2_traces_clean2[nbofdetectors];

    int x_samples[200];
    for(int i = 1; i < 201; i++){
        x_samples[i] = i;
    }
    for(int i = 0; i<nbofdetectors; i++){
        sprintf(histdesc,  "CH%i: Traces", i);
        sprintf(histname, "h2_traces_%i", i);
        h2_traces[i] = new TH2D(histname, histdesc, 200, 0, 200, 1000, -6000, 500);
        sprintf(histdesc, "CH%i: Traces w/o prepeak step", i);
        sprintf(histname, "h2_traces_clean1_%i", i);
        h2_traces_clean1[i] = new TH2D(histname, histdesc, 200, 0, 200, 1000, -6000, 500);
        sprintf(histdesc, "CH%i: Traces without pileups", i);
        sprintf(histname, "h2_traces_clean2_%i", i);
        h2_traces_clean2[i] = new TH2D(histname, histdesc, 200, 0, 200, 1000, -6000, 500);
    }

    Double_t entries = datatree->GetEntries();

    std::vector<particle> particles;
    bufferdata buffer;
    trace_properties trace_props;

    energy_ch muon;
    // muon.CH = {652.10, 633.39, 733.54, 667.98, 741.10, 712.55, 696.39, 715.85, 697.87, 725.30, 712.13, 387.25, 408.65, 404.82, 392.22};
    muon.CH = {970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8};
    muon.o_CH = {62.36, 57.54, 76.99, 67.60, 74.04, 77.66, 71.42, 70.80, 70.91, 75.69, 73.13, 35.84, 37.41, 37.17, 36.32};
    muon.E = 4.95;
    muon.o_E = 0.3433;
    
    energy_ch na22;
    na22.CH = {405, 410.2, 415.6, 447.3, 441.77, 438.24, 416.8, 48.9, 408.8, 422.34, 461.6, 51.59, 51.39, 49.26, 36.42};
    // na22.CH = {187.95, 182.31, 202.14, 185.35, 218.14, 175.91, 178.96, 187.66, 168.23, 173.47, 190.77, 116.02, 121.70, 114.04, 90.49};
    na22.o_CH = {33.10, 32.94, 42.41, 38.47, 46.53, 41.11, 38.38, 39.66, 39.22, 43.07, 41.87, 24.57, 25.99, 21.39, 18.50};
    // na22.E = 1.2745;
    na22.E = 1.3325;
    na22.o_E = 0;

    calibration calib; // init calib
    calib = energy_extrapolation(calib, na22, muon); // get calib

    Int_t prevbufferCounter = 0;
    Int_t prevbuffer = 0;
    bool epoch_skip = false;
    int missing_buffer_counter = 0;
    int pileup_counter = 0;
    int prepeak_step_counter = 0;

    cout << "Beamtime: " << beamtime[beamtime_select] << endl;
    cout << "Input file: " << filename << endl;
    for (double e = 0; e<entries; e++){
    	datatree->GetEntry(e);
        channel--;
        trace_props = get_max_trace(trace->at(0));
        trace_props.trace = trace->at(0);
        double quenched_energy = calib.GetQuenchedEnergy(channel, trace_props.amp);
        trace_props.energy = calib.GetEnergy(quenched_energy);
        trace_props.energy_err = calib.GetQuenchedStdDev(channel, quenched_energy);
        trace_props.time_ns = static_cast<double>(timestamp_ns)*2;
        trace_props.time_ps = static_cast<double>(timestamp_ps);

        if(prevbufferCounter == 20 && prevbuffer != eventCounter){
            prevbufferCounter = 0;
            for(int i = 0; i < buffer.trace_props.size(); i++){
                if(buffer.trace_props.at(i).size()>0){
                    particles.resize(buffer.trace_props.at(i).size());
                    break;
                }
            }
            for(int i=0; i<buffer.trace_props.at(0).size();i++){//--------------- start data of epoch
                particles.at(i).edep.at(0) = buffer.trace_props.at(0).at(i).energy;
                particles.at(i).time.at(0) = buffer.trace_props.at(0).at(i).time_ns;
                particles.at(i).amplitude.at(0) = buffer.trace_props.at(0).at(i).amp;
                particles.at(i).trace.at(0) = buffer.trace_props.at(0).at(i).trace;
            }        
            
            for(int ch=1; ch<buffer.trace_props.size(); ch++){ // -------------------- get coinc data
                for(int p=0; p<particles.size(); p++){
                    for(int tr=0; tr<buffer.trace_props.at(ch).size(); tr++){ //iterate over <double>max_amp & timestamp of one channel
                        particles.at(p).coincidence(buffer.trace_props.at(ch).at(tr), ch);
                    }
                }
            }  
            for(int p=0; p<particles.size(); p++){ //--------------------------------- fill h_total_edep & h_traceamp_coinc
                particles.at(p).TestBuffer();
                if(particles.at(p).edep.at(0)!=0){                      // fill total_edep                    
                    for(int i=0; i<buffer.trace_props.size(); i++){                        
                        if(particles.at(p).edep.at(i)!=0){                            
                            for(int ii = 0; ii < particles.at(p).trace.at(i).size(); ++ii){
                                h2_traces[i]->Fill(ii, particles.at(p).trace.at(i).at(ii));
                            }                  
                            if(particles.at(p).det_prepeak_step == false){    
                                for(int ii = 0; ii < particles.at(p).trace.at(i).size(); ++ii){                                                            
                                    h2_traces_clean1[i]->Fill(ii, particles.at(p).trace.at(i).at(ii));
                                }                                
                            }
                            if(particles.at(p).det_pileup == false && particles.at(p).det_prepeak_step == false){
                                for(int ii = 0; ii < particles.at(p).trace.at(i).size(); ++ii){
                                    h2_traces_clean2[i]->Fill(ii, particles.at(p).trace.at(i).at(ii));
                                }  
                            }
                        }
                    }
                }
            }            
            particles.clear();
            buffer.clear();
        }        
        buffer.trace_props.at(channel).push_back(trace_props);
        if(prevbuffer != eventCounter){
            prevbufferCounter++;
        }
        prevbuffer = eventCounter;
        trace_props.clear();
    }
    
    //init canvas
    //TCanvas *c1 = new TCanvas("c1","Bragg Sampler MIT Analysis", 10, 10, 1000, 500);
    sprintf(title, "Bragg Sampler %s Analysis", filename);
    TCanvas *c1 = new TCanvas("c1", title, 4000, 2000);
    //set canvas setting
    sprintf(file, "%sall_trace_%s.root", out_path, filename);
    TFile *	hfile = new TFile(file,"RECREATE");
    for(int ch=0; ch<nbofdetectors;ch++){
        c1->Clear();
        c1->SetFillColor(0);
	    c1->SetGrid();
	    c1->SetBorderMode(0);
	    c1->SetBorderSize(2);
	    c1->SetFrameBorderMode(0);
        c1->Divide(3, 1);
        c1->cd(1);
        plot_th2d(h2_traces[ch], "bins", "Trace Height[AU]");
        c1->Update();
        c1->cd(2);
        plot_th2d(h2_traces_clean1[ch], "bins", "Sample Height[AU]");
        c1->Update();
        c1->cd(3);
        plot_th2d(h2_traces_clean2[ch], "bins", "Sample Height[AU]");
        c1->Update();
        sprintf(file, "h_all_traces_ch%i", ch);
        c1->Write(file);

        c1->Clear();
        c1->SetFillColor(0);
	    c1->SetGrid();
	    c1->SetBorderMode(0);
	    c1->SetBorderSize(2);
	    c1->SetFrameBorderMode(0);
        plot_th2d(h2_traces[ch], "bins", "Sample Height[AU]");
        c1->Update();

        sprintf(file, "%sh_all_traces_ch%i.pdf", out_path, ch);
        c1->Update();
        c1->Print(file);
        sprintf(file, "h_all_traces_ch%i", ch);
        c1->Write(file);

        c1->Clear();
        c1->SetFillColor(0);
	    c1->SetGrid();
	    c1->SetBorderMode(0);
	    c1->SetBorderSize(2);
	    c1->SetFrameBorderMode(0);
        plot_th2d(h2_traces_clean1[ch], "bins", "Sample Height[AU]");
        c1->Update();
        sprintf(file, "h_all_traces_clean1_ch%i", ch);
        c1->Write(file);
        
        c1->Clear();
        c1->SetFillColor(0);
	    c1->SetGrid();
	    c1->SetBorderMode(0);
	    c1->SetBorderSize(2);
	    c1->SetFrameBorderMode(0);
        plot_th2d(h2_traces_clean2[ch], "bins   ", "Sample Height[AU]");
        c1->Update();
        sprintf(file, "%sh_all_traces_clean_ch%i.pdf", out_path, ch);
        c1->Update();
        c1->Print(file);
        sprintf(file, "h_all_traces_clean1_ch%i", ch);
        c1->Write(file);
    }

    hfile->Close();
}