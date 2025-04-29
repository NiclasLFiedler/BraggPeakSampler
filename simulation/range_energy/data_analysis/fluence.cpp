#include "TH1D.h"

double resolution = 0.1; // resolution of ideal detector in mm
double rho_pwo = 8.28;
double rho_h2o = 0.997;
double rho_pen = 1.36;
double rho_air = 0.001225;

//-------------bethe bloch implementation------------


constexpr double pi = TMath::Pi();
constexpr double e = TMath::E();

constexpr double c = 299792458;
constexpr double electron_mass = 0.51099895069;  // Electron mass (MeV/c^2)
constexpr double e_charge = 1.602176634e-19;         // Elementary charge (C)
constexpr double avogadro = 6.02214076e23;      // Avogadro's number 1/mol
constexpr double epsilon0 = 8.854187817e-12; // Vacuum permittivity in C/(N*m^2) (farads per meter)

//----------------------------------------------------

double we_air_length = rho_air/rho_h2o * 3;
double we_pwo_length = rho_pwo/rho_h2o * 5;
double we_pen_length = rho_pen/rho_h2o * 5;

int coincidence_layer = 10;
int MeV = 1e6;

double I_pbwo4 = 542.3e-6;
double I_h2o = 75e-6;

int scaleFactor = 1/resolution;

double bethe1 = 4*pi/(electron_mass);// unit = 1/MeV
double bethe3 = std::pow(e_charge/(4*pi*epsilon0),2)*1e-6; // unit = MeV^2*mm^2

std::unordered_map<int, double> StoppingPowerMap;

double range_energy_relationship(double E){
    double a_pwo=7.275e-4;
    double p_pwo=1.690;
    return a_pwo* std::pow(E, p_pwo);
}

void ReadStoppingPower(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double E, dEdX_per_mass_unit;

        // Reading the key and value from the line
        ss >> E >> dEdX_per_mass_unit;
        int E_int = static_cast<int>(std::round(E * 1/resolution));
              
        // Store them in the map
        StoppingPowerMap[E_int] = rho_h2o*dEdX_per_mass_unit*1/10;
    }

    file.close();
}

double bethe_bloch(double beta=0, double I=0, double Z=0, double A=0, double rho=0, int z=0){
    
    // beta = v/c, I = mean excitation energy (MeV), Z = materials atomic number, A = material relative atomic mass, rho = density (g/cm^3), z = particle charge, M_u = molar mass (g/mol),
    //cout << "bethe1 " << bethe1 << endl;
    double bethe2 = avogadro*Z*rho*std::pow(z/beta, 2)/(A)*1/1000; // unit = 1/mm^3 
    //bethe2 = avogadro*Z*std::pow(z/beta, 2)/(A)*1/1000; // unit = 1/mm^3 
    //cout << "bethe2 " << bethe2 << endl;
    
    double bethe4 = (std::log(2*electron_mass*std::pow(beta,2)/(I*(1-std::pow(beta,2))))-std::pow(beta,2)); // unit = 1
    //cout << "bethe4 " << bethe4 << endl;
    //cout << "bethe bloch " << bethe1*bethe2*bethe3*bethe4 << endl;
    return bethe1*bethe2*bethe3*bethe4;
}

int ScaleToInt(double value){
    return static_cast<int>(std::round(value * scaleFactor));
}

double ScaleToDouble(double scaledValue) {
    return static_cast<double>(scaledValue) / scaleFactor;
}

double roundEnergy(double value) {
    return std::round(value * scaleFactor) / scaleFactor;
}

double round(double value, int roundFactor) {
    return std::round(value * roundFactor) / roundFactor;
}

struct dose{
    double fDose = 0;

    double GetDose(){
        return fDose;
    }

    void AddDose(double edep){
        fDose += edep;
        return;
    }

    void Clear(){
        fDose = 0;
        return;
    }
};

struct detector_theory
{
    int fnbofdetectors;
    int fnbofdetectorsPbWO4;
    std::vector<std::unordered_map<int, std::shared_ptr<dose>>> fTotalDose;
    std::unordered_map<int, std::shared_ptr<dose>> fWETotalDose;

    detector_theory(const int nbofdetectors, const int nbofdetectorsPbWO4){
        fnbofdetectors = nbofdetectors;
        fnbofdetectorsPbWO4 = nbofdetectorsPbWO4;
        fTotalDose = std::vector<std::unordered_map<int, std::shared_ptr<dose>>>(fnbofdetectors); // depth,dose
        fWETotalDose = std::unordered_map<int, std::shared_ptr<dose>>(fnbofdetectors); // water equivalent depth,dose
    }

    void SetDoses(std::vector<std::unordered_map<int, std::shared_ptr<dose>>> inDose, std::vector<std::unordered_map<int, std::shared_ptr<dose>>> inWEDose){
        int channel = 0;
        int offset = 0;
        for(const auto& doseMap : inDose){
            for(const auto& [depth, mdose]: doseMap){
                auto dose_pair = fTotalDose.at(channel).find(depth);
                if(dose_pair != fTotalDose.at(channel).end()){ //find existing pair
                    dose_pair->second->AddDose(mdose->GetDose());
                }
                else{
                    std::shared_ptr<dose> newDose = std::make_shared<dose>(); //create new pair
                    newDose->AddDose(mdose->GetDose());
                    fTotalDose.at(channel)[depth] = newDose;
                }
            }
            channel++;
        }    
        channel = 0;

        // for(int i = 0; i < fnbofdetectors; i++){
        //     offset = 5*i+3*i;
        //     for(const auto& [depth, mdose]: inDose.at(i)){
        //         double WEdepth = CalculateWaterEquivalentDepth(channel, depth-offset);
        //         auto WEdose_pair = fWETotalDose.at(channel).find(std::ceil(WEdepth*1/resolution));

        //         if(WEdose_pair != fWETotalDose.at(channel).end()){ //find existing pair
        //             WEdose_pair->second->AddDose(mdose->GetDose());///std::abs(WEdepth-WEdepth_prev));
        //         }
        //         else{                                                                
        //             std::shared_ptr<dose> newWEDose = std::make_shared<dose>(); //create new pair
        //             newWEDose->AddDose(mdose->GetDose());///std::abs(WEdepth-WEdepth_prev));
        //             fWETotalDose.at(channel)[std::ceil(WEdepth*1/resolution)] = newWEDose;
        //         }
        //     }
        // }
        for(const auto& doseMap : inWEDose){
            for(const auto& [depth, mdose] : doseMap){
                auto dose_pair = fWETotalDose.find(depth);
                if(dose_pair != fWETotalDose.end()){ //find existing pair
                    dose_pair->second->AddDose(mdose->GetDose());
                }
                else{
                    std::shared_ptr<dose> newDose = std::make_shared<dose>(); //create new pair
                    newDose->AddDose(mdose->GetDose());
                    fWETotalDose[depth] = newDose;
                }
            }
            channel++;
        }
        return;
    }

    double WaterEquivalentLength(double length, double density){
        return density/rho_h2o * length;
    }

    double CalculateWaterEquivalentDepth(int channel, double depth){
        double WEdepth = 0;
        if(channel < fnbofdetectorsPbWO4){
            WEdepth = channel*(we_air_length+we_pwo_length)+WaterEquivalentLength(depth, rho_pwo);
        }
        else{
            WEdepth = fnbofdetectorsPbWO4*(we_air_length+we_pwo_length)+(channel-fnbofdetectorsPbWO4)*(we_air_length+we_pen_length)+WaterEquivalentLength(depth, rho_pen);
        }
        return WEdepth;
    }

};

struct particle {
    Int_t fnbofdetectors;
    Int_t fnbofdetectorsPbWO4;
    Double_t TotalEDep;
    double offset;
    double WEdepth_prev = 0;
    double WEdepth = 0;
    double initialEnergy = 0;
    double range = 0;
    std::vector<double> edep;
    std::vector<int> edep_idx;
    std::vector<std::vector<double>> dE;
    std::vector<std::vector<double>> fbeta;
    std::vector<std::vector<double>> feKin;
    std::vector<std::vector<double>> fprojZ;
    std::vector<double> fmax_projZ;
    std::vector<std::vector<Double_t>> StepLength; 
    std::vector<std::vector<Double_t>> dEdX; 
    std::vector<std::unordered_map<int, std::shared_ptr<dose>>> fDose;
    std::vector<std::unordered_map<int, std::shared_ptr<dose>>> fWEDose;

    particle(const int& nbofdetectors, const int& nbofdetectorsPbWO4){
        fnbofdetectors = nbofdetectors;
        fnbofdetectorsPbWO4 = nbofdetectorsPbWO4;
        TotalEDep = 0;
        initialEnergy = 0;
        range = 0;
        edep = std::vector<double>(nbofdetectors);
        edep_idx = std::vector<int>(nbofdetectors);
        dE = std::vector<std::vector<double>>(nbofdetectors);
        feKin = std::vector<std::vector<double>>(nbofdetectors);
        fbeta = std::vector<std::vector<double>>(nbofdetectors);
        fprojZ = std::vector<std::vector<double>>(nbofdetectors);
        fmax_projZ = std::vector<double>(nbofdetectors);
        StepLength = std::vector<std::vector<Double_t>>(nbofdetectors); 
        dEdX = std::vector<std::vector<Double_t>>(nbofdetectors);
        fDose = std::vector<std::unordered_map<int, std::shared_ptr<dose>>>(nbofdetectors); // depth,dose
        fWEDose = std::vector<std::unordered_map<int, std::shared_ptr<dose>>>(nbofdetectors); // water equivalent depth,dose
    }

    Double_t GetTotalEDep(){
        TotalEDep = 0;
        for(int i=0; i<edep.size(); i++){
            TotalEDep += edep.at(i);
        }
        return TotalEDep;
    }

    void CalculateEDep(const int channel){
        for(const double& dE_i :dE.at(channel)){
            edep.at(channel) += dE_i;
        }
        return;
    }

    void CalculateEDep_idx(const int channel){
        double fEDep_idx = 0;
        for(const double& dE_i :dE.at(channel)){
            fEDep_idx += dE_i;
        }
        edep_idx.at(channel) = ScaleToInt(fEDep_idx);
        return;
    }

    void CalculateRange(){
        range = 0;
        for(const auto& StepLengthChannel : StepLength){
            for(const double& step : StepLengthChannel){
                range += step;
            }
        }
        return;
    }

    void CalculateMaxProjZ(){
        for(int channel = 0; channel<fprojZ.size(); channel++){
            auto maxElement = std::max_element(fprojZ.at(channel).begin(), fprojZ.at(channel).end());
            if (maxElement != fprojZ.at(channel).end()) {
                fmax_projZ.at(channel) = *maxElement;
            } else {
                std::cout << "The vector is empty." << std::endl;
            }
        }
        return;
    }

    void ProcessParticle(){
        CalculateRange();
        CalculateMaxProjZ();
        for(int i=0; i<edep.size(); i++){
            CalculateEDep(i);
            CalculateEDep_idx(i);
        }
        return;
    }

    Double_t GetStepLength(const int channel, const int index){
        return StepLength.at(channel).at(index);
    }

    std::vector<std::vector<double>> GetprojZ(){
        return fprojZ;
    }

    void SetInitialEnergy(double E){
        initialEnergy = E;
        return;
    }

    void SetdE(const int channel, double dEdep){
        dE.at(channel).push_back(dEdep);
        return;
    }
    
    void Setbeta(const int channel, double beta){
        fbeta.at(channel).push_back(beta);
        return;
    }
    
    void SeteKin(const int channel, double eKin){
        feKin.at(channel).push_back(eKin);
        return;
    }

    void SetprojZ(const int channel, double projZ){
        fprojZ.at(channel).push_back(projZ);
        return;
    }

    void SetStepLength(const int channel, double StepLength_i){
        StepLength.at(channel).push_back(StepLength_i);
        return;
    }

    void SetdEdX(const int channel, double dEdX_i){
        dEdX.at(channel).push_back(dEdX_i);
        return;
    }

    bool Coincidence(int layer){
        if(layer <= coincidence_layer){
            layer = coincidence_layer;
        }
        for(int i = 0; i < layer; i++){
            if(edep.at(i) < 0.1){
                return false;
            }
        }
        return true;
    }

    void Clear(){
        TotalEDep = 0;
        initialEnergy = 0;
        range = 0;
        edep = std::vector<double>(fnbofdetectors); 
        edep_idx = std::vector<int>(fnbofdetectors);
        dE = std::vector<std::vector<double>>(fnbofdetectors);
        feKin = std::vector<std::vector<double>>(fnbofdetectors);
        fbeta = std::vector<std::vector<double>>(fnbofdetectors);
        fprojZ = std::vector<std::vector<double>>(fnbofdetectors);
        StepLength = std::vector<std::vector<Double_t>>(fnbofdetectors); 
        dEdX = std::vector<std::vector<Double_t>>(fnbofdetectors); 
    }
};

void draw_legend(TH1D* hist1){
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
    gPad->Update(); // Necessary to update the pad

    hist->GetXaxis()->SetTitle(x_title);
    hist->GetYaxis()->SetTitle(y_title);
    hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetTitleSize(0.04);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetXaxis()->SetTitleSize(0.04);
    hist->Draw("COLZ");
    gPad->Update(); // Necessary to update the pad
    
    TPaveStats *stats = (TPaveStats*)hist->FindObject("stats");
    // Get the pointer to the stats box
    // Set the position of the stats box to top left
    stats->SetX1NDC(0.15);
    stats->SetX2NDC(0.45);
    stats->SetY1NDC(0.65);
    stats->SetY2NDC(0.85);
}

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

void plot_tgraph(TGraph* Graph, TString x_title, TString y_title, bool same = false){
    gPad->SetGrid();
    // Set border line color and width
    gPad->SetFrameLineColor(1); // White color
    gPad->SetFrameLineWidth(1); // Borderline width
    // Adjust margins
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.1);
    
    Graph->SetMarkerStyle(2);  // Crosses
    Graph->SetMarkerSize(1);
    Graph->GetXaxis()->SetTitle(x_title);
    Graph->GetYaxis()->SetTitle(y_title);
    Graph->GetYaxis()->SetLabelFont(42);
	Graph->GetYaxis()->SetLabelSize(0.04);
	Graph->GetYaxis()->SetTitleSize(0.04);
	Graph->GetXaxis()->SetLabelFont(42);
	Graph->GetXaxis()->SetLabelSize(0.04);
	Graph->GetXaxis()->SetTitleSize(0.04);
    if(same){
        Graph->Draw("P Same");
        return;
    }
    Graph->Draw("AP");
    return;
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

void fluence(int energy = 0){
    ROOT::EnableImplicitMT();
    Char_t filename[100]; // temp store of histnames
    cout << "energy " << energy << endl; 
    if(energy == 0){
        sprintf(filename, "raw_data_p.root");
    }
    else{
        sprintf(filename, "raw_data_%i.root", energy);
    }
    cout << "filename " << filename << endl;
    //init variables
    Int_t event = 0; //store tree values
    Double_t EDep = 0;
    Double_t Pos = 0;
    Double_t dEdX = 0;
    Double_t beta = 0;
    Double_t TrackID = 0;
    Double_t eKin = 0;
    Double_t eTot = 0;
    Double_t StepLength = 0;
    Char_t histname[100]; // temp store of histnames
    Char_t histName[100]; // temp store of histnames
    Char_t histTitle[100]; // temp store of histnames
    Char_t histname2[100]; // temp store of histnames

    Double_t edep = 0; //energy storage for hist fill

    //init  histograms and file access
    TFile *input = new TFile(filename, "READ");

    // Check if the file is successfully opened
    if (!input || input->IsZombie()) {
        cout << "Error: Failed to open file 'data.root'!" << endl;
        return;
    }

    sprintf(filename, "range_%i.root", energy);

    TFile *OutputFile = new TFile(filename, "RECREATE");
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Failed to create new ROOT file." << std::endl;
        return;
    }

    //get data from tree
    TTree *datatree = (TTree*)input->Get("braggsampler");
    
    TTree *outputtree = new TTree("range", "Tree of proton ranges");

    // Check if the tree is successfully retrieved
    if (!datatree) {
        cout << "Error: Failed to retrieve tree 'tree' from file!" << endl;
        input->Close();
        return;
    }

    Double_t tree_range;
    outputtree->Branch("range", &tree_range);

    TGraph *g_detectortheory = new TGraph();
    sprintf(histName, "g_detectortheory");
    sprintf(histTitle, "Simulation depth dose with 0.1 mm resolution");
    g_detectortheory->SetName(histName);
    g_detectortheory->SetTitle(histTitle);

	datatree->SetBranchAddress("event", &event);
	datatree->SetBranchAddress("eDep", &EDep);	
    datatree->SetBranchAddress("pos", &Pos);
    datatree->SetBranchAddress("dEdX", &dEdX);
    datatree->SetBranchAddress("beta", &beta);
    datatree->SetBranchAddress("trackid", &TrackID);
    datatree->SetBranchAddress("eKin", &eKin);
    datatree->SetBranchAddress("eTot", &eTot);
    datatree->SetBranchAddress("StepLength", &StepLength);

    double prevEDep = 0;
    datatree->GetEntry(0);
    int prevEvent = -1;
    double prev_WEpos = 0;
    double prev_pos = 0;
    double WEpos = 0;

    std::vector<particle> particles;
    particle proton(1, 1);
    
    //get Stopping Power Map
    ReadStoppingPower("pstar_stopping_power_water_protons.txt");
    
    int64_t entries = datatree->GetEntries();
	TH2D *h2_curve = new TH2D("h2_curve", "h2_edep", 500, 0, 500, 500, 0, 5);
    TH1D *h_range = new TH1D("h_range", "h_range", 310, 0, 310);
    TH1D *h_fluence = new TH1D("h_fluence", "h_fluence", 400, 0, 400);

    for (int64_t i = 0; i<entries; i++){
        if((i + 1) % (entries / 10) == 0) {
            int percentage = (i + 1) * 100 / entries;
            std::cout << "Progress: " << percentage << "% completed" << std::endl;
        }
		datatree->GetEntry(i);
        
        if(prevEvent != event){
            if(i != 0){
                particles.push_back(proton);
            }
            proton.Clear();
            //if(particles.size() == 1000) break;
            proton.SetInitialEnergy(eKin);
        }

        if(TrackID == 1){
            proton.SetdE(0, roundEnergy(EDep));            
            proton.SetStepLength(0, StepLength);
            proton.SetdEdX(0, roundEnergy(dEdX));
            proton.SetprojZ(0, Pos);
            proton.Setbeta(0, beta);
            proton.SeteKin(0, eKin);
        }
        prevEvent = event; 
    }

    cout << "Calculating energy depositions for coincidences." << endl;
    int nbOfProtons = particles.size();
    int particlesIdx = 0;
    double planeThickness = 5;
    detector_theory DetectorTheory(1, 1);
    double E = 250;
    cout << "Expected Range: " << range_energy_relationship(E) << endl;
    for(particle p : particles){
        if ((particlesIdx + 1) % (nbOfProtons / 10) == 0) {
            int percentage = (particlesIdx + 1) * 100 / nbOfProtons;
            std::cout << "Progress: " << percentage << "% completed" << std::endl;
        }
        p.ProcessParticle();
        for(int hitLayers = 0; hitLayers<p.range/planeThickness; hitLayers++){
            h_fluence->Fill(hitLayers*planeThickness);
        }
        p.GetTotalEDep();
        tree_range = p.range;
        h_range->Fill(tree_range);
        outputtree->Fill();
    }
    cout << "writing data " << endl;
    outputtree->Write();
    
    TFile *fluencefile = new TFile("fluence.root", "RECREATE");
    TTree *fluencetree = new TTree("fluence", "fluence");

    // Variables to hold bin center and bin content
    double binCenter = 0;
    double binContent = 0;

    // Set up the branches
    fluencetree->Branch("depth", &binCenter, "depth/D"); 
    fluencetree->Branch("counts", &binContent, "counts/D");

    for (int i = 1; i <= h_fluence->GetNbinsX(); i++) {
        if(h_fluence->GetBinContent(i) > 0){
            binCenter = h_fluence->GetBinCenter(i)-0.5;
            binContent = h_fluence->GetBinContent(i);
            fluencetree->Fill();
        }
    }
    fluencetree->Write();
    fluencefile->Close();

    // Draw the curve over the histogram
    TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 3000, 1000);
    c1->Clear();
    c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetFrameBorderMode(0);
    
    double j = 0;

    simple_plot(h_fluence, "Range / mm", "Counts");
}