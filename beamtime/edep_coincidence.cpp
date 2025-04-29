#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"
//#include "TraceProperties.h"
#include "TraceProperties.cpp"
#include "Particle.cpp"
#include "Plotter.cpp"
#include "EnergyCh.h"
#include "Calibration.h"

//I put everything together, take what you need. If you need anything else or if something is unclear I will provide additional information
//-----------------------------------------------------------------------------------//
//-----------------------------USE-OF-SCRIPT-----------------------------------------//
//-----------------------------------------------------------------------------------//
int coincidence_time = 20;              // set coincidence window in ns
int coincidence_layer = 10;
double tot_edep_cut_off = 0;           // set total energy cut off  ----155 beam2 --75 beam1 -----130-hetero---140-homo
int discard_index = 50;                  // set index of prepeak step cut off
bool reversed = true;                    // enable reversed energy calibration
bool bhetero = true;
//-----------------------------------------------------------------------------------//
//-----------------------------DIR-STRUCTURE-----------------------------------------//
//-----------------------------------------------------------------------------------//
//----*some_dir*/edep_coincidence.cpp                                    // script
//----*some_dir*/*beamtimedir*/*beamrun_name*/input/*beamrun_name*.root  // input files
//----*some_dir*/*beamtimedir*/*beamrun_name*/output/                    // output dir
//-----------------------------------------------------------------------------------//
//-----------------------------INPUT-FILES-------------------------------------------//
const char* beamtime[4] = {"MIT_04_2024", "MIT_05_2024", "MIT_11_2024", "MIT_12_2024"};
// const char* in_data[12] = {"bp-p-beam1","bp-p-beam2", "bp-p-beam2.5","bp-p-beam3","bp-p-beam4","bp-p-target1","bp-p-target2", "bp-p-hetero-target", "bp-p-homo-target", "notarget_run1", "run3_notarget", "run6_notargetcd no    "}; // 
const char* in_data[4] = {"bp-p-beam2", "bp-p-homo-target", "bp-p-hetero-target", "notarget_medium"};
const char* target_data[4] = {"without a target", "with the homogeneous target", "with the heterogeneous target", "without a target"};
//const char* in_data[12] = {"notarget_low", "notarget_medium", "notarget_high", "homo_medium", "homo_high", "hetero_medium", "hetero_high"};
//-------0-------------1------------2------------3-------------4-------------5-------//
//   bp-p-beam1   bp-p-beam2.5  bp-p-beam3  bp-p-beam4   bp-p-target1   bp-p-target2
int file_select = 2;                    // select data file with index from above
int beamtime_select = 1;                // select beamtime
// sprintf(target_title, "with the homogeneous Target");
// sprintf(target_title, "with the heterogeneous Target");
//-----------------------------------------------------------------------------------//
//-----------------------------REMARKS-----------------------------------------------//
//-----------------------------------------------------------------------------------//
// may need to change histogram range and binning
//
// histogram -> desc
// h_traceamp[i] -> all trace_amps of detector layer i
// h_energy[i] -> trace_amps with coincidence with ch0 of detector layer i
// h_energy_cut[i] -> trace_amps of h_energy with total energy requirements of < or > than 'tot_edep_cut_off' of detector layer i
// h_stopped[i] -> last trace_amps of h_energy - not necessarily stopped can also have escaped the detector of detector layer i
// h_means -> mean_values of h_traceamp[i]
// h_means_coinc -> mean_values of h_energy[i]
// h_means_edep_cut -> mean_values of h_energy_cut[i]
// h_total_edep -> total edep of h_energy[i]
// h_total_edep_cut -> total edep of h_energy_cut[i]
//-----------------------------------------------------------------------------------//

struct bufferdata{
    std::vector<std::vector<TraceProperties>> trace_props = std::vector<std::vector<TraceProperties>>(15);
    void clear(){
        trace_props = std::vector<std::vector<TraceProperties>>(15);
    }
};

void edep_coincidence(){
    ROOT::EnableImplicitMT();
    Plotter plotter(bhetero, coincidence_time, tot_edep_cut_off);
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
    std::vector<Double_t> *trace2 = 0;
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

    if(beamtime_select == 3){ 
        datatree->SetBranchAddress("EventCounter", &eventCounter);
    }
    else{
        datatree->SetBranchAddress("BufferCounter", &eventCounter);
    }
    datatree->SetBranchAddress("Channel", &channel);
    datatree->SetBranchAddress("TimeTag", &timestamp_ns);
    datatree->SetBranchAddress("TimeStampPico", &timestamp_ps);
    if(beamtime_select == 3){ 
        datatree->SetBranchAddress("Trace", &trace2);	  
    }
    else{
        datatree->SetBranchAddress("Trace", &trace);
    }
    TH1D *h_traceamp[nbofdetectors];
    TH1D *h_energy[nbofdetectors];
    TH1D *h_energy_cut[nbofdetectors];
    TH1D *h_energy_cut_low[nbofdetectors];
    TH1D *h_stopped[nbofdetectors];
    

    sprintf(histdesc, "Norm. energy depth dose distribution %s", target_data[file_select]);
    TGraphErrors *h_means = new TGraphErrors();
    h_means->SetTitle(histdesc);
    h_means->SetLineWidth(2);
    h_means->SetMarkerStyle(20);
    h_means->SetMarkerSize(1);

    sprintf(histdesc, "Total energy deposition %s", target_data[file_select]);
    TH1D *h_total_edep = new TH1D("h_total_edep", histdesc, 270, 0, 270);
    TH1D *h_total_edep_cut = new TH1D("h_total_edep_cut", histdesc, 270, 0, 270);

    for(int i = 0; i<nbofdetectors; i++){
        sprintf(histdesc, "Energy Deposition of CH%i %s", i, target_data[file_select]);
        sprintf(histname, "h_traceamp_%i", i);
        h_traceamp[i] = new TH1D(histname, histdesc, 240, 0, 60);        
        sprintf(histname, "h_energy_%i", i);
        h_energy[i] = new TH1D(histname, histdesc, 240, 0, 60);        
        sprintf(histname, "h_energy_cut_%i", i);
        sprintf(histdesc, "Energy Deposition of CH%i %s", i, target_data[file_select]);
        h_energy_cut[i] = new TH1D(histname, histdesc, 240, 0, 60);
        sprintf(histname, "h_energy_cut_low%i", i);
        sprintf(histdesc, "Energy deposition of CH%i with max total energy of %0.0f MeV", i, tot_edep_cut_off);
        h_energy_cut_low[i] = new TH1D(histname, histdesc, 240, 0, 60);     
        sprintf(histname, "h_stopped_%i", i);
        sprintf(histdesc, "Energy deposition of stopped particle in CH%i", i);
        h_stopped[i] = new TH1D(histname, histdesc, 240, 0, 60);        
    }

    std::vector<Particle> particles;
    bufferdata buffer;
    TraceProperties trace_props;

    energy_ch muon;
    muon.CH = {970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8};
    // muon.CH = {940.48, 928.89, 951.95, 1092.53, 1102.58, 1041.98, 946.03, 175.53, 946.45, 1083.89, 1180.65, 149.96, 164.90, 137.50, 90.82}; //old
    muon.o_CH = {100.52, 90.4, 100.65, 98.31, 107.62, 117.89, 98.19, 23.22, 99.15, 96.28, 116.67, 17.89, 19.18, 15.76, 9.33}; //old
    //muon.CH = {652.10, 633.39, 733.54, 667.98, 741.10, 712.55, 696.39, 715.85, 697.87, 725.30, 712.13, 387.25, 408.65, 404.82, 392.22}; new
    //muon.o_CH = {62.36, 57.54, 76.99, 67.60, 74.04, 77.66, 71.42, 70.80, 70.91, 75.69, 73.13, 35.84, 37.41, 37.17, 36.32}; new
    muon.E = 4.95;
    muon.o_E = 0.3433;
    
    energy_ch na22;
    na22.CH = {187.95, 182.31, 202.14, 185.35, 218.14, 175.91, 178.96, 187.66, 168.23, 173.47, 190.77, 116.02, 121.70, 114.04, 90.49};
    na22.o_CH = {33.10, 32.94, 42.41, 38.47, 46.53, 41.11, 38.38, 39.66, 39.22, 43.07, 41.87, 24.57, 25.99, 21.39, 18.50};
    na22.E = 1.1745;
    na22.o_E = 0;

    energy_ch co60;
    //co60.CH = {405, 410.2, 415.6, 447.3, 441.77, 438.24, 416.8, 48.9, 408.8, 422.34, 461.6, 51.59, 51.39, 49.26, 36.42}; old
    co60.CH = {402.56, 401.75, 421.52, 444.63, 437.76, 428.52, 415.91, 49.03, 403.4, 417.63, 457.43, 52.08, 52.11, 47.11, 35.55};
    co60.o_CH = {97.79, 97.82, 101.99, 99.31, 115.65, 102.12, 90.23, 11.05, 96.9, 97.69, 111.66, 20.22, 23.25, 19.32, 12.67};
    co60.E = 1.25275;
    co60.o_E = 0.068;

    Calibration calib; // init calib
    calib.energy_extrapolation(co60, muon); // get calib

    if(beamtime_select == 3){ 
        muon.CH = {652.10, 633.39, 733.54, 667.98, 741.10, 712.55, 696.39, 715.85, 697.87, 725.30, 712.13, 387.25, 408.65, 404.82, 392.22};
        muon.o_CH = {62.36, 57.54, 76.99, 67.60, 74.04, 77.66, 71.42, 70.80, 70.91, 75.69, 73.13, 35.84, 37.41, 37.17, 36.32};
        calib.energy_extrapolation(na22, muon); 
    }

    for(int i = 0; i < 15;i++){
        cout << i << " Epsilon: " << calib.epsilon_1 << " Slope: " << calib.slope.at(i) << " Chi: " << calib.chi_1.at(i) << endl;
    }

    Double_t entries = datatree->GetEntries();
    Int_t preveventCounter = 0;
    Int_t prevbuffer = 0;
    bool epoch_skip = false;
    int missing_buffer_counter = 0;
    int pileup_counter = 0;
    int prepeak_step_counter = 0;
    int coinc_layer_counter = 0;

    cout << "Beamtime: " << beamtime[beamtime_select] << endl;
    cout << "Input file: " << filename << endl;
    for (double e = 0; e<100000; e++){
    	datatree->GetEntry(e);
        if(beamtime_select == 3){ 
            trace_props.SetTrace(*trace2);
        }
        else{
            channel--;
            trace_props.SetTrace(trace->at(0));
        }
        trace_props.Process();
        double quenched_energy = calib.GetQuenchedEnergy(channel, trace_props.amp);
        trace_props.energy = calib.GetEnergy(quenched_energy);
        trace_props.energy_err = calib.GetQuenchedStdDev(channel, quenched_energy);
        trace_props.time_ns = static_cast<double>(timestamp_ns)*2;
        trace_props.time_ps = static_cast<double>(timestamp_ps);
        h_traceamp[channel]->Fill(quenched_energy);

        if(preveventCounter == 20 && prevbuffer != eventCounter){
            preveventCounter = 0;
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
            }
            
            for(int ch=1; ch<buffer.trace_props.size(); ch++){ // -------------------- get coinc data
                for(int p=0; p<particles.size(); p++){
                    for(int tr=0; tr<buffer.trace_props.at(ch).size(); tr++){ //iterate over <double>max_amp & timestamp of one channel
                        particles.at(p).coincidence(buffer.trace_props.at(ch).at(tr), ch);
                    }
                }
            }
            for(int p=0; p<particles.size(); p++){ //--------------------------------- fill h_total_edep & h_energy
                particles.at(p).TestBuffer();
                particles.at(p).TestCoincLayer();
                if(particles.at(p).missing_buffer == true){
                    missing_buffer_counter++;
                    continue;
                }
                if(particles.at(p).det_pileup == true){
                    pileup_counter++;
                    continue;
                }
                if(particles.at(p).det_prepeak_step == true){
                    prepeak_step_counter++;
                    continue;
                }
                if(particles.at(p).coinc_layer < coincidence_layer){
                    coinc_layer_counter++;
                    continue;
                }
                if(particles.at(p).edep.at(0)!=0){                      // fill total_edep
                    particles.at(p).sum_edep();
                    if(particles.at(p).total_edep<30) continue;
                    if(particles.at(p).total_edep>=tot_edep_cut_off){    // min_tot_edep                             
                        h_total_edep_cut->Fill(particles.at(p).total_edep);
                    }
                    else{
                        h_total_edep->Fill(particles.at(p).total_edep);
                    }
                    for(int i=0; i<buffer.trace_props.size(); i++){ // fill individual layers
                        if(particles.at(p).edep.at(i)!=0){
                            h_energy[i]->Fill(particles.at(p).GetEDep(i)); // all
                            if(particles.at(p).total_edep>=tot_edep_cut_off){    // min_tot_edep 
                                h_energy_cut[i]->Fill(particles.at(p).GetEDep(i));                                    
                            }
                            if(particles.at(p).total_edep<tot_edep_cut_off){    // max_tot_edep 
                                h_energy_cut_low[i]->Fill(particles.at(p).GetEDep(i));                                    
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
        buffer.trace_props.at(channel).push_back(trace_props);
        if(prevbuffer != eventCounter){
            preveventCounter++;
        }
        prevbuffer = eventCounter;
        trace_props.Clear();
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

    sprintf(histdesc, " with total EDep cutoff %0.0f MeV", tot_edep_cut_off);
    // Header
    cout << setfill(' ');
    cout << setw(10) << "Channel" << setw(2) << "|" << setw(15) << "Unfiltered" << setw(2) << "|"  << setw(25) << "Coinc. /w channel 1" << setw(2) << "|"  << setw(30) << histdesc << setw(2) << "|"  << setw(18) << " Stopped particle" << endl;
    cout << setfill('-') << setw(106) << "-" << endl;
    cout << setfill(' ');

    for(int i = 0; i<nbofdetectors; i++){
        c1->cd(i+1);
        plotter.plot_th1d(h_traceamp[i], "EDep [MeV]", "Counts");
        //h_energy[i]->SetLineColor(kRed+1);
        //h_energy[i]->Draw("HIST SAME");
        h_energy_cut[i]->SetLineColor(kOrange+1);
        h_energy_cut[i]->Draw("HIST SAME");
        h_energy_cut_low[i]->SetLineColor(kGreen+1);
        h_energy_cut_low[i]->Draw("HIST SAME");
        //h_stopped[i]->SetLineColor(kOrange+1);
        //h_stopped[i]->Draw("HIST SAME");
        plotter.draw_legend(h_traceamp[i], h_energy[i], h_energy_cut[i]);

        cout << setw(10) << i << setw(2) << "|" << setw(15) << h_traceamp[i]->GetEntries() << setw(2) << "|" << setw(25) << h_energy[i]->GetEntries() << setw(2) << "|" << setw(30) << h_energy_cut[i]->GetEntries() << setw(2) << "|" << setw(18) <<  h_stopped[i]->GetEntries() << endl;
    }
    cout << setfill('-') << setw(106) << "-" << endl;
    cout << setfill(' ');
    cout << "Detected missing buffer entries: " << missing_buffer_counter << endl;
    cout << "Detected pileups: " << pileup_counter << endl;
    cout << "Detected prepeak steps: " << prepeak_step_counter << endl;
    cout << "Number of lower coincidence particles: " << coinc_layer_counter << endl;
    c1->cd(nbofdetectors+1);

    double totalenergydose[nbofdetectors];
    double totalenergydose_quaderr[nbofdetectors];
    double x_sigma[nbofdetectors];
    double x[nbofdetectors];
    double alpha_pbwo4 = 7.275e-3;
    double p_pbwo4 = 1.690;
    double alpha_h2o = 2.585e-2;
    double p_h2o = 1.738;
    double energy = 220;
    double wet_coef = alpha_h2o*std::pow(energy, p_h2o)/(alpha_pbwo4*std::pow(energy,p_pbwo4));
    double wet_coef2 = alpha_h2o/alpha_pbwo4*(p_h2o-p_pbwo4)*std::pow(energy, p_h2o-p_pbwo4-1);
    for(int ch = 0; ch<nbofdetectors; ch++){
        x[ch] = (0.25 + 0.5*ch)*wet_coef;
        x_sigma[ch] = sqrt(std::pow((0.5/sqrt(12)*wet_coef),2) + std::pow(((0.25 + 0.5*ch)*wet_coef2*0.005*energy), 2));
    }
    for(int i = 0;  i < nbofdetectors; i++){
        totalenergydose[i] = 0;
        totalenergydose_quaderr[i] = 0;
	    for(int j = 1; j<=h_energy_cut[i]->GetNbinsX(); j++){
            double bin_content = h_energy_cut[i]->GetBinContent(j);
            if(bin_content == 0) continue;
            double energy = h_energy_cut[i]->GetBinCenter(j);
            double bin_width = h_energy_cut[i]->GetBinWidth(j);
            bin_width = 1;
            totalenergydose[i] += bin_content * energy * bin_width;
            totalenergydose_quaderr[i] += std::pow(energy * bin_width, 2) / bin_content + std::pow(bin_width * bin_content, 2) * calib.GetQuenchedVar(i, calib.FromEnergyToQuenched(energy));
        }
        std::cout << "Channel: " << i << " Total Energy Dose: " << totalenergydose[i] << " +- " << sqrt(totalenergydose_quaderr[i]) << " ~" << sqrt(totalenergydose_quaderr[i])/totalenergydose[i]*100 << "%" <<std::endl;
        h_means->SetPoint(i, x[i], totalenergydose[i]/h_energy_cut[0]->GetEntries());
        h_means->SetPointError(i, x_sigma[i], sqrt(totalenergydose_quaderr[i])/h_energy_cut[0]->GetEntries());
    }

    gPad->SetGrid();
    h_means->GetXaxis()->SetTitle("Depth [cm]");
    h_means->GetYaxis()->SetTitle("Norm. Energy Dose [MeV]");
    h_means->GetYaxis()->SetLabelFont(42);
	h_means->GetYaxis()->SetLabelSize(0.07);
	h_means->GetYaxis()->SetTitleSize(0.07);
    h_means->GetYaxis()->SetTitleOffset(1);
	h_means->GetXaxis()->SetLabelFont(42);
	h_means->GetXaxis()->SetLabelSize(0.07);
	h_means->GetXaxis()->SetTitleSize(0.07);
    h_means->SetLineColor(kOrange+1);
	h_means->SetFillColor(kOrange+1);
    h_means->Draw("AP");
    // draw_legend(h_means, h_means_coinc, h_means_edep_cut);

    c1->cd(nbofdetectors+2);
    h_total_edep->GetYaxis()->SetRangeUser(0, h_total_edep_cut->GetMaximum() * 1.1);
    h_total_edep->SetLineColor(kGreen+1);
    plotter.plot_th1d(h_total_edep, "EDep [MeV]", "Counts");
    h_total_edep_cut->SetLineColor(kOrange+1);
    h_total_edep_cut->Draw("HIST SAME");
    
    //outputfile
    sprintf(file, "%sedep.root", out_path);
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
	    plotter.plot_th1d( h_energy_cut[i], "EDep [MeV]", "Counts");
        if(i != 0){
            h_energy_cut_low[i]->Draw("HIST SAME");
        }
        plotter.draw_legend(h_energy_cut[i], h_energy_cut_low[i]);

        sprintf(histname, "h_edep_%i", i);
        c2->Write(histname);
    }
    c2->cd(nbofdetectors+1);    
    h_means->Draw("AP");
    // draw_legend(h_means, h_means_coinc, h_means);
    c2->Write("mean_edep");

    c2->cd(nbofdetectors+2);
    plotter.plot_th1d(h_total_edep, "EDep [MeV]", "Counts");
    h_total_edep_cut->Draw("HIST SAME");

    c2->Write("total_edep");

    hfile->Close();
    sprintf(file, "%sMIT_05_2024_%s.root", out_path, in_data[file_select]);
    hfile = new TFile(file,"RECREATE");

    const char* in_data2[12] = {"beam2","homo-target", "hetero-target","beam3","beam4","target1","target2", "hetero-target", "homo-target", "homo-target", "homo-target", "homo-target"}; // maybe change to 2D array
    //save analysised data
    TCanvas *c3 = new TCanvas("c3", title, 4000, 2000);
    for(int i = 0; i < nbofdetectors; i++){
        c3->SetFillColor(0);
        c3->SetGrid();
	    c3->SetBorderMode(0);
	    c3->SetBorderSize(2);
	    c3->SetFrameBorderMode(0);
        if(bhetero){
            plotter.plot_th1d(h_energy_cut[i], "EDep [MeV]", "Counts");
            gStyle->SetTitleFontSize(0.07); // Title font size for the canvas
            plotter.draw_legend(h_energy_cut[i], h_energy_cut[i]);
        }
        else{
	        plotter.plot_th1d(h_energy[i], "EDep [MeV]", "Counts");
            h_energy_cut[i]->Draw("hist same");
            h_energy_cut_low[i]->Draw("hist same");
            plotter.draw_legend(h_energy[i], h_energy_cut[i], h_energy_cut_low[i]);
        }
        sprintf(file, "%sh_edep_%i.pdf", out_path, i);
        c3->Update();
        c3->Print(file);
        c3->Clear();
    }

    c3->SetFillColor(0);
    c3->SetGrid();
	c3->SetBorderMode(0);
	c3->SetBorderSize(2);
	c3->SetFrameBorderMode(0);
    
    gStyle->SetTitleFontSize(0.07); // Title font size for the canvas
    h_means->Draw("AP");
    gStyle->SetEndErrorSize(5);
    c3->Update();
    sprintf(file, "%sh_means.pdf", out_path);
    c3->Print(file);

    if(bhetero){
        plotter.plot_th1d(h_total_edep_cut, "EDep [MeV]", "Counts");
        gStyle->SetTitleFontSize(0.07); // Title font size for the canvas
        plotter.draw_legend2(h_total_edep_cut, h_total_edep_cut);
    }
    else{
        plotter.plot_th1d(h_total_edep, "EDep [MeV]", "Counts");
        h_total_edep_cut->Draw("HIST SAME");
        plotter.draw_legend2(h_total_edep, h_total_edep_cut);
    }
    sprintf(file, "%sh_total_edep.pdf", out_path);
    c3->Print(file);

    for(int i = 0; i < nbofdetectors; i++){
        h_energy_cut[i]->Write();
    }
    hfile->Close();

    //store hist means for bortfield fit
    sprintf(file, "%s%sMeans.root", out_path, in_data2[file_select]);
    hfile = new TFile(file, "RECREATE");
    
    // Create a TTree
    TTree *meantree = new TTree("meantree", "Tree storing histogram means");

    // Variable to store histogram means
    Int_t layer;
    Double_t mean;
    Double_t error;

    // Create a branch in the TTree
    meantree->Branch("ch", &layer);
    meantree->Branch("mean", &mean);
    meantree->Branch("error", &error);

    for (int ch = 0; ch < nbofdetectors; ++ch) {
        layer = ch;
        mean = totalenergydose[ch]/h_energy_cut[0]->GetEntries();
        error = sqrt(totalenergydose_quaderr[ch])/h_energy_cut[0]->GetEntries();
        
        meantree->Fill();        
    }
    meantree->Write();    
    hfile->Close();
}