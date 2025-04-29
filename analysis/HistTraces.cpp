#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <TBranch.h>
#include "TROOT.h"
#include "src/TraceProperties.cpp"
#include "src/Particle.cpp"
#include "src/Plotter.cpp"
#include "include/EnergyCh.h"
#include "include/Calibration.h"


double CalcWETConv(double p_h2o, double p_m, double alpha_h2o, double alpha_m, double energy){
    return p_h2o/p_m*alpha_h2o/alpha_m*std::pow(energy, p_h2o-p_m);
}

double CalcWETConvSigma(double p_h2o, double p_m, double alpha_h2o, double alpha_m, double energy){
    return p_h2o/p_m*alpha_h2o/alpha_m*(p_h2o-p_m)*std::pow(energy, p_h2o-p_m-1);
}

void HistTraces(){
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
    TH2D *h2_traces_clean[nbofdetectors];

    int x_samples[200];
    for(int i = 1; i < 201; i++){
        x_samples[i] = i;
    }
    for(int i = 0; i<nbofdetectors; i++){
        sprintf(histdesc,  "CH%i: Traces", i);
        sprintf(histname, "h2_traces_%i", i);
        h2_traces[i] = new TH2D(histname, histdesc, 200, 0, 200, 1000, -6000, 500);
        sprintf(histdesc, "CH%i: Traces w/o pileup", i);
        sprintf(histname, "h2_traces_clean1_%i", i);
        h2_traces_clean[i] = new TH2D(histname, histdesc, 200, 0, 200, 1000, -6000, 500);
    }

    energy_ch muon;
    muon.CH = {970, 905, 968, 1072.4, 1094.3, 1078.8, 983.6, 177.9, 965.5, 1035.08, 1171.35, 162.08, 172.74, 145, 101.8};
    muon.o_CH = {100.52, 90.4, 100.65, 98.31, 107.62, 117.89, 98.19, 23.22, 99.15, 96.28, 116.67, 17.89, 19.18, 15.76, 9.33}; //old
    muon.E = 4.95;
    muon.o_E = 0.3433;
    
    energy_ch co60;
    co60.CH = {402.56, 401.75, 421.52, 444.63, 437.76, 428.52, 415.91, 49.03, 403.4, 417.63, 457.43, 52.08, 52.11, 47.11, 35.55};
    co60.o_CH = {97.79, 97.82, 101.99, 99.31, 115.65, 102.12, 90.23, 11.05, 96.9, 97.69, 111.66, 20.22, 23.25, 19.32, 12.67};
    co60.E = 1.25275;
    co60.o_E = 0.068;

    Calibration calib; // init calib
    calib.energy_extrapolation(co60, muon); // get calib

    calibration calib; // init calib
    calib = energy_extrapolation(calib, na22, muon); // get calib

    bool epoch_skip = false;
    int missing_buffer_counter = 0;
    int pileup_counter = 0;
    int prepeak_step_counter = 0;

    cout << "Beamtime: " << beamtime[beamtime_select] << endl;
    cout << "Input file: " << filename << endl;


    Double_t entries = datatree->GetEntries();
    Int_t prevEvent = -1;
    int missing_buffer_counter = 0;
    int pileup_counter = 0;
    int prepeak_step_counter = 0;
    int coinc_layer_counter = 0;
    bool discard = false;

    TraceProperties trace_props;
    std::vector<Particle> particles;
    std::vector<Particle> ScintParticles;
    Particle proton(nbofdetectors, coincidence_time, coincidence_layer, discard_index);
    std::unordered_map<double, TraceProperties> initialEvents;
    std::multimap<double, TraceProperties> postEvents;

    cout << "Dataset: " << dataset << endl;
    cout << "Input file: " << filename << endl;
    
    cout << "Measurement: Processing raw data from ROOT file." << endl;
    for (double e = 0; e<entries; e++){
    	datatree->GetEntry(e);
        channel--;
        trace_props.Clear();    
        trace_props.SetParameters(trace->at(0), channel, static_cast<double>(timestamp_ns)*2, static_cast<double>(timestamp_ps), calib);
        h_edep_all[channel]->Fill(trace_props.quenched_energy);
        if(trace_props.channel == 0){
            initialEvents.insert({trace_props.time_ps, trace_props});
        }
        else{
            postEvents.insert({trace_props.time_ps, trace_props});
        }        
    }
    for(const auto& [initialTime, initialTrace] : initialEvents) {
        proton.Clear();
        proton.InsertInitial(initialTrace);
        auto lowerTraces = postEvents.lower_bound(initialTime - coincidence_time);
        auto upperTraces = postEvents.upper_bound(initialTime + coincidence_time);            
        for (auto coincTrace = lowerTraces; coincTrace != upperTraces; ++coincTrace) {                                
            proton.Coincidence(coincTrace->second);
        }            
        particles.push_back(proton);
    }
    for(auto coinProton : particles){
        coinProton.Test();   
        coinProton.SumEDep();
        for(auto sampleTrace : coinProton.trace){
            
        }
        discard = false;
    }  


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
        plot_th2d(h2_traces_clean[ch], "bins", "Sample Height[AU]");
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