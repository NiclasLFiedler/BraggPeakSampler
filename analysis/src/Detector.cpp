#include "../include/Detector.h"
#include "../include/DetectorProperties.h"


Detector::Detector(DetectorProperties* detectorProperties)
    : detectorProperties(detectorProperties) {

    Char_t histdesc[100];
    Char_t histname[100];
    

    int binningEdep[3] = {32768/8, 0, 32768};
    int binningPhotons[3] = {200, 0, 50};
    int binningEntry[3] = {200, -static_cast<int>(detectorProperties->GetLayerSizeX())/2, static_cast<int>(detectorProperties->GetLayerSizeX())/2};
    int binningAngle[3] = {400, 0, 4};

    sprintf(histdesc, "Total energy deposition %s", detectorProperties->GetTarget());
    h_total_edep = new TH1D("h_total_edep", histdesc, 230*100, 0, 230);
    
    g_means = new TGraphErrors();
    
    for(int layer = 0; layer<detectorProperties->GetNLayers(); layer++){
        sprintf(histdesc, "Energy Deposition of CH%i %s", layer, detectorProperties->GetTarget());
        sprintf(histname, "h_edep_all_%i", layer);
        h_edep_all.push_back(new TH1D(histname, histdesc, binningEdep[0], binningEdep[1], binningEdep[2]));
        sprintf(histname, "h_edep_coinc_%i", layer);
        h_edep_coinc.push_back(new TH1D(histname, histdesc, binningEdep[0], binningEdep[1], binningEdep[2]));
        
        sprintf(histname, "h_stopped_%i", layer);
        sprintf(histdesc, "Energy deposition of stopped particle in CH%i", layer);
        h_stopped.push_back(new TH1D(histname, histdesc, binningEdep[0], binningEdep[1], binningEdep[2]));

        sprintf(histname, "h_photons_%i", layer);
        sprintf(histdesc, "Number of photons in SiPM%i", layer);
        h_photons.push_back(new TH1D(histname, histdesc, binningPhotons[0], binningPhotons[1], binningPhotons[2]));
        sprintf(histname, "h_photons_coinc_%i", layer);
        h_photons_coinc.push_back(new TH1D(histname, histdesc, binningPhotons[0], binningPhotons[1], binningPhotons[2]));

        sprintf(histname, "h_entry_%i", layer);
        sprintf(histdesc, "Number of entry protons layer %i", layer);
        h_entry.push_back(new TH1D(histname, histdesc, binningEntry[0], binningEntry[1], binningEntry[2]));
        sprintf(histname, "h_exit_%i", layer);
        sprintf(histdesc, "Number of exit protons layer %i", layer);
        h_exit.push_back(new TH1D(histname, histdesc, binningEntry[0], binningEntry[1], binningEntry[2]));
        sprintf(histname, "h_angle_%i", layer);
        sprintf(histdesc, "Angle of protons in layer %i", layer);
        h_angle.push_back(new TH1D(histname, histdesc, binningAngle[0], binningAngle[1], binningAngle[2]));
    }

    crystals = std::vector<Crystal>(detectorProperties->GetNLayers());

    Rw = detectorProperties->materials["h2o"].R;
    alphaw = detectorProperties->materials["h2o"].alpha;
    pw = detectorProperties->materials["h2o"].p;

    Rm = detectorProperties->GetR();
    alpham = detectorProperties->GetAlpha();
    pm = detectorProperties->GetP();

    CalcWETConv();
    std::cout << "WET Conversion Factors: " << std::endl;
    std::cout << "WETConv: " << WETConv << std::endl;
    cout << "Detector Constructed" << endl;
}

TH1D* Detector::EnergyHist(const int channel){
    return h_edep_all.at(channel);
}

TH1D* Detector::CoincEnergyHist(const int channel){
    return h_edep_coinc.at(channel);
}

TH1D* Detector::StoppedEnergyHist(const int channel){
    return h_stopped.at(channel);
}

TH1D* Detector::PhotonHist(const int channel){
    return h_photons.at(channel);
}

TH1D* Detector::CoincPhotonHist(const int channel){
    return h_photons_coinc.at(channel);
}

TH1D* Detector::EntryHist(const int channel){
    return h_entry.at(channel);
}

TH1D* Detector::ExitHist(const int channel){
    return h_exit.at(channel);
}

TH1D* Detector::AngleHist(const int channel){
    return h_angle.at(channel);
}

TH1D* Detector::TotalEnergyHist(){
    return h_total_edep;
}

TGraphErrors* Detector::MeansGraph(){
    return g_means;
}

void Detector::CalcWETConv(){
    WETConv = detectorProperties->GetP("h2o")/detectorProperties->GetP()*detectorProperties->GetAlpha("h2o")/detectorProperties->GetAlpha()*std::pow(detectorProperties->GetBeamEnergy(), detectorProperties->GetP("h2o")-detectorProperties->GetP());
    WETConvTeflon = detectorProperties->GetP("h2o")/detectorProperties->GetP("teflon")*detectorProperties->GetAlpha("h2o")/detectorProperties->GetAlpha("teflon")*std::pow(detectorProperties->GetBeamEnergy(), detectorProperties->GetP("h2o")-detectorProperties->GetP("teflon"));
    WETConvAlu = detectorProperties->GetP("h2o")/detectorProperties->GetP("alu")*detectorProperties->GetAlpha("h2o")/detectorProperties->GetAlpha("alu")*std::pow(detectorProperties->GetBeamEnergy(), detectorProperties->GetP("h2o")-detectorProperties->GetP("alu"));
    return;
}

void Detector::Process(){
    for(int layer = 0; layer<detectorProperties->GetNLayers(); layer++){
        CalcDose(layer);
        CalcPosition(layer);
    }
    FillMeansGraph();
    return;
}

void Detector::CalcDose(int layer){
    double fbin_content = 0;
    double fenergy = 0;
    double fdeviationSq = 0;

    // for(int j = 1; j<=h_edep_coinc[layer]->GetNbinsX(); j++){
    //     fbin_content = h_edep_coinc[layer]->GetBinContent(j);
    //     if(fbin_content == 0) continue;
    //     fenergy = h_edep_coinc[layer]->GetBinCenter(j);
    //     crystals.at(layer).dose.dose += fbin_content * fenergy;
    //     if(detectorProperties->GetCalibrationStatus()){
    //         fdeviationSq += fenergy*fenergy / fbin_content + fbin_content* fbin_content * detectorProperties->GetCalibration()->GetQuenchedVar(layer, detectorProperties->GetCalibration()->FromEnergyToQuenched(layer, fenergy));
    //     }
    //     else{
    //         fdeviationSq += fenergy*fenergy / fbin_content;
    //     }
    // }
    // crystals.at(layer).dose.stddev = sqrt(fdeviationSq);
    
    
    // int layerEntries = h_edep_coinc[layer]->GetEntries();
    // double totalProtons = h_edep_coinc[0]->GetEntries();

    // // if (layerEntries < 0.01
    // //      * totalProtons) {
    // //     std::cout << "Layer " << layer << " has less than 0.01% of the total protons " << layerEntries << "  " << totalProtons << std::endl;
    // //     crystals.at(layer).dose.dose = 0; 
    // //     crystals.at(layer).dose.stddev = 0;
    // // }

    // // if(detectorProperties->GetNormStatus() && h_edep_coinc[layer]->GetEntries() > 0){
    // //     crystals.at(layer).dose.dose*=1/h_edep_coinc[layer]->GetEntries();
    // //     crystals.at(layer).dose.stddev*=1/h_edep_coinc[layer]->GetEntries();
    // // }

    // if(detectorProperties->GetNormStatus() && h_edep_coinc[0]->GetEntries() > 0){
    //     crystals.at(layer).dose.dose*=1/h_edep_coinc[0]->GetEntries();
    //     crystals.at(layer).dose.stddev*=1/h_edep_coinc[0]->GetEntries();
    // }

    crystals.at(layer).dose.dose = h_edep_all[layer]->GetMean();
    crystals.at(layer).dose.stddev = h_edep_all[layer]->GetStdDev();

    if(detectorProperties->GetNormStatus()){
        crystals.at(layer).dose.dose*=h_edep_all[layer]->GetEntries()/h_edep_all[0]->GetEntries();
        crystals.at(layer).dose.stddev*=h_edep_all[layer]->GetEntries()/h_edep_all[0]->GetEntries();
    }

    crystals.at(layer).dose.dose*=1/detectorProperties->GetLayerSizeZ(layer);
    crystals.at(layer).dose.stddev*=1/detectorProperties->GetLayerSizeZ(layer);

    return;
}

void Detector::CalcPosition(int layer){
    double depth = 0;
    if(layer != 0){
        depth = crystals.at(layer-1).pos.depth;
    }
    double aluthick = detectorProperties->GetAluThickness();
    double teflonthick = detectorProperties->GetTeflonThickness();
    double delta = 0;
    if(layer == 0){
        depth = CalcWET(depth, aluthick/2, detectorProperties->GetMaterial("alu"));
        depth = CalcWET(depth, teflonthick/2, detectorProperties->GetMaterial("teflon"));
        delta = CalcWET(depth, detectorProperties->GetLayerSizeZ(layer), detectorProperties->GetMaterial(detectorProperties->GetScintillator()));
        
        delta = delta-depth; 
    }
    else{
        depth += crystals.at(layer-1).size/2;
        depth = CalcWET(depth, teflonthick/2, detectorProperties->GetMaterial("teflon"));
        if(detectorProperties->GetAbsorberStatus() && layer == 1){
            depth = CalcWET(depth, detectorProperties->GetAbsorberSize(), detectorProperties->GetMaterial(detectorProperties->GetAbsorberType()));
        } 
        depth = CalcWET(depth, aluthick, detectorProperties->GetMaterial("alu"));
        depth = CalcWET(depth, teflonthick/2, detectorProperties->GetMaterial("teflon"));
        delta = CalcWET(depth, detectorProperties->GetLayerSizeZ(layer), detectorProperties->materials[detectorProperties->GetScintillator()]);
        delta = delta-depth;
    }
    crystals.at(layer).pos.depth = depth+delta/2;
    crystals.at(layer).size = delta;
    //std::cout << crystals.at(layer).pos.depth << std::endl;

    // crystals.at(layer).pos.depth = (detectorProperties->GetLayerSizeZ(layer-1)/2+detectorProperties->GetLayerSizeZ(layer)/2)*WETConv;
    // if(layer == 0){
    //     crystals.at(layer).pos.depth += (detectorProperties->GetTeflonThickness(layer)/2)*WETConvTeflon+(detectorProperties->GetAluThickness(layer)/2)*WETConvAlu;
    // }
    // if(layer > 0){
    //     crystals.at(layer).pos.depth += crystals.at(layer-1).pos.depth+(detectorProperties->GetTeflonThickness(layer)/2+detectorProperties->GetTeflonThickness(layer-1)/2)*WETConvTeflon+(detectorProperties->GetAluThickness(layer)/2+detectorProperties->GetAluThickness(layer-1)/2)*WETConvAlu;
    // }
    // if(detectorProperties->GetAbsorberStatus() && layer == 1){
    //     crystals.at(layer).pos.depth += absorberConv;
    // }

    crystals.at(layer).pos.stddev = sqrt(std::pow((crystals.at(layer).size/sqrt(12)),2) + std::pow(crystals.at(layer).pos.depth*0.005*(detectorProperties->GetP("h2o")-detectorProperties->GetP()), 2));
    return;
}

void Detector::FillMeansGraph(){
    for(int layer = 0; layer<detectorProperties->GetNLayers(); layer++){
        g_means->SetPoint(layer, crystals.at(layer).pos.depth, crystals.at(layer).dose.dose);
        g_means->SetPointError(layer, crystals.at(layer).pos.stddev, crystals.at(layer).dose.stddev);
    }
}

void Detector::Info(){
    std::cout << "Detector Information" << std::endl;
    std::cout << setfill('-') << setw(80) << "-" << std::endl;
    
    std::cout << "Scintillator Type: " << detectorProperties->GetScintillator() << std::endl;
    std::cout << "Beam energy: " << detectorProperties->GetBeamEnergy() << " MeV" << std::endl;
    std::cout << "Number of Layers: " << detectorProperties->GetNLayers() << std::endl;
    std::cout << "Crystal Size: ";
    PrintVector(detectorProperties->GetScintillatorDimensions());
    std::cout << "GapSize in z: " << detectorProperties->GetGapSizeZ() << " mm" <<  std::endl;
    std::cout << "Aluminum Thickness: " << detectorProperties->GetAluThickness() << " mm" <<  std::endl;
    std::cout << "Teflon Thickness: " << detectorProperties->GetTeflonThickness() << " mm" << std::endl;
    std::cout << "Absorber Type: " << detectorProperties->GetAbsorberType() <<  std::endl;
    std::cout << "Absorber Size: " << detectorProperties->GetAbsorberSize() << " mm" <<  std::endl;
    std::cout << "Number of Secondary Layers " << detectorProperties->GetNSecondaryLayers() << std::endl;
    std::cout << "Secondary Layer Size in z: " << detectorProperties->GetSecondaryLayerSizeZ() << " mm" <<  std::endl;
    std::cout << "Target: " << detectorProperties->GetTarget() << std::endl;
    std::cout << "Norm status: " << detectorProperties->GetNormStatus() << std::endl;
    std::cout << "Reversed status: " << detectorProperties->GetReversedStatus() << std::endl;
    std::cout << "Simulation status: " << detectorProperties->GetSimulationStatus() << std::endl;
    std::cout << "Calibration status: " << detectorProperties->GetCalibrationStatus() << std::endl;
    
    std::cout << setfill('-') << setw(80) << "-" << std::endl;
    std::cout << setfill(' ');
    return;
}

void Detector::PrintVector(const std::vector<double>& vec) {
    std::cout << "{";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "}" << std::endl;
    return;
}

double Detector::CalcWET(double z0, double delta, MaterialProperties mat){
    MaterialProperties h2o = detectorProperties->materials["h2o"];
    return h2o.alpha*std::pow((std::pow((z0/h2o.alpha), mat.p/h2o.p)+delta/mat.alpha), h2o.p/mat.p);
}
