#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TrackerSD.hh"
#include "SiPMSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorParameterisationColour.hh"
#include "HeteroParametrisation.cc"
#include "G4RegularNavigation.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"
#include "G4MultiFunctionalDetector.hh"
#include "Randomize.hh"

#include "G4EmCalculator.hh"
#include "G4Proton.hh"

//#include "G4LogicalBorderSurface.hh"
#include "nlohmann/json.hpp"
#include <filesystem>

namespace fs = std::filesystem;
using json = nlohmann::json;

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
{  
    std::ifstream configFile("../../../analysis/config.json");
    if (!configFile) {
        std::cerr << "Error opening config file!" << std::endl;
        return;
    }

    json allConfigs;
    configFile >> allConfigs;
    
    const auto& config                                = allConfigs["detectors"][int(allConfigs["detectorSelect"])];
    ftarget                                           = allConfigs["targetSelect"];

    detectorType                                      = config["detectorType"];
    //double beamEnergy                               = config["beamEnergy"]*MeV;
    fLayers                                           = config["nLayers"];
    std::vector<double> crystalSize                   = config["crystalSize"];
    gapSizeZ                                          = config["gapSizeZ"];
    secondaryLayerStatus                              = config["secondaryLayerStatus"];
    fLayersCut                                        = config["nSecondaryLayers"];
    absSizeZ                                          = config["secLayerSizeZ"];
    absorberStatus                                    = config["absorberStatus"];
    absorberSize                                      = config["absorberSize"];
    ThicknessTeflon                                   = config["teflonThickness"];
    ThicknessAlu                                      = config["aluThickness"];
    absorberType                                      = config["absorberType"];
    targetThickness                                   = config["targetThickness"];
    
    pmod                                              = allConfigs["pmod"];
    heteroThickness                                   = targetThickness;
    

    detSizeX = crystalSize.at(0)* mm;
    detSizeY = crystalSize.at(1)* mm;
    detSizeZ = crystalSize.at(2)* mm;

  switch(ftarget) {
    case 1:
      G4cout << "Simulation with homogeneous target" << G4endl;
      break;
    case 2:
      G4cout << "Simulation with heterogeneous target" << G4endl;
      break;
    default:
      G4cout << "Simulation without target" << G4endl;
  }
  fMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout << "Define materials" << G4endl;
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Material definition

  G4NistManager* nistManager = G4NistManager::Instance();

  G4Element *elH = nistManager->FindOrBuildElement("H");    // 1
  G4Element *elB = nistManager->FindOrBuildElement("B");    // 5
  G4Element *elC = nistManager->FindOrBuildElement("C");    // 6
  G4Element *elN = nistManager->FindOrBuildElement("N");    // 7
  G4Element *elO = nistManager->FindOrBuildElement("O");    // 8
  G4Element *elF = nistManager->FindOrBuildElement("F");    // 9
  G4Element *elNa = nistManager->FindOrBuildElement("Na");  // 11
  G4Element *elMg = nistManager->FindOrBuildElement("Mg");  // 12
  G4Element *elAl = nistManager->FindOrBuildElement("Al");  // 13
  G4Element *elSi = nistManager->FindOrBuildElement("Si");  // 14
  G4Element *elP = nistManager->FindOrBuildElement("P");    // 15
  G4Element *elS = nistManager->FindOrBuildElement("S");    // 16
  G4Element *elCl = nistManager->FindOrBuildElement("Cl");  // 17
  G4Element *elAr = nistManager->FindOrBuildElement("Ar");  // 18
  G4Element *elK = nistManager->FindOrBuildElement("K");    // 19
  G4Element *elCa = nistManager->FindOrBuildElement("Ca");  // 20
  G4Element *elFe = nistManager->FindOrBuildElement("Fe");  // 26
  G4Element *elZn = nistManager->FindOrBuildElement("Zn");  // 30
  G4Element *elBa = nistManager->FindOrBuildElement("Ba");  // 56
  G4Element *elCe = nistManager->FindOrBuildElement("Ce");  // 56
  G4Element *elGd = nistManager->FindOrBuildElement("Gd");  // 64
  G4Element *elW = nistManager->FindOrBuildElement("W");    // 74
  G4Element *elPb = nistManager->FindOrBuildElement("Pb");  // 82

  // Air defined using NIST Manager
  //worldMat = nistManager->FindOrBuildMaterial("G4_AIR"); 
  worldMat = nistManager->FindOrBuildMaterial("G4_Galactic"); 
  urethandimethacrylat = new G4Material("urethandimethacrylat", 1.11*g/cm3,4);
  urethandimethacrylat->AddElement(elC,23);
  urethandimethacrylat->AddElement(elH,38);
  urethandimethacrylat->AddElement(elN,2);
  urethandimethacrylat->AddElement(elO,8); //50 35 10 4 1

  methacrylatmonomere = new G4Material("methacrylatmonomere", 1.029*g/cm3,3);
  methacrylatmonomere->AddElement(elC,7);
  methacrylatmonomere->AddElement(elH,12);
  methacrylatmonomere->AddElement(elO,3);

  phosphinoxid = new G4Material("phosphinoxid", 1.19*g/cm3,4);
  phosphinoxid->AddElement(elC,26);
  phosphinoxid->AddElement(elH,27);
  phosphinoxid->AddElement(elO,3);
  phosphinoxid->AddElement(elP,1);

  resinMaterial = new G4Material("resinMaterial", 1.09*g/cm3, 3);  
  resinMaterial->AddMaterial(urethandimethacrylat,0.741);
  resinMaterial->AddMaterial(methacrylatmonomere,0.25);
  resinMaterial->AddMaterial(phosphinoxid,0.009);

  lungTissue = new G4Material("lungTissue", 1.05*g/cm3,13);
  lungTissue->AddElement(elH, 0.101278);
  lungTissue->AddElement(elC, 0.102310);
  lungTissue->AddElement(elN, 0.028650);
  lungTissue->AddElement(elO, 0.757072);
  lungTissue->AddElement(elNa, 0.001840);
  lungTissue->AddElement(elMg, 0.000730);
  lungTissue->AddElement(elP, 0.000800);
  lungTissue->AddElement(elS, 0.002250);
  lungTissue->AddElement(elCl, 0.002660);
  lungTissue->AddElement(elK, 0.001940);
  lungTissue->AddElement(elCa, 0.000090);
  lungTissue->AddElement(elFe, 0.000370);
  lungTissue->AddElement(elZn, 0.000010);

  water = nistManager->FindOrBuildMaterial("G4_WATER"); 
  water->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material *SiO2 = new G4Material("SiO2", 2.65*g/cm3, 2);
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);

  G4Material *BaO = new G4Material("BaO", 5.72*g/cm3, 2);
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);

  G4Material* Gd2O3 = new G4Material("Gd2O3", 7.41*g/cm3, 2);
  Gd2O3->AddElement(elGd, 2);
  Gd2O3->AddElement(elO, 3);

  G4Material* AlF3 = new G4Material("AlF3", 2.88*g/cm3, 2);
  AlF3->AddElement(elAl, 1);
  AlF3->AddElement(elF, 3);

  G4Material* Ce2O3 = new G4Material("Ce2O3", 6.2*g/cm3, 2);
  Ce2O3->AddElement(elCe, 2);
  Ce2O3->AddElement(elO, 3);

  G4Material* DSB_Gd = new G4Material("DSB_Gd", 4.3*g/cm3, 5);
  // DSB_Gd->AddMaterial(SiO2, 54.81*perCent);
  // DSB_Gd->AddMaterial(BaO, 18.70*perCent);
  // DSB_Gd->AddMaterial(AlF3, 2.45*perCent);
  // DSB_Gd->AddMaterial(Ce2O3, 1.32*perCent);
  // DSB_Gd->AddMaterial(Gd2O3, 22.72*perCent);
  //old comp
  DSB_Gd->AddMaterial(BaO, 26.8*perCent);
  DSB_Gd->AddMaterial(SiO2, 30.8*perCent);
  DSB_Gd->AddMaterial(AlF3, 1.9*perCent);
  DSB_Gd->AddMaterial(Gd2O3, 38.5*perCent);
  DSB_Gd->AddMaterial(Ce2O3, 2*perCent);

  G4Element* BoronEnriched = new G4Element("BoronEnriched", "B", 2);
  G4Isotope* B10 = new G4Isotope("B10", 5, 10, 10.012937*g/mole);
  G4Isotope* B11 = new G4Isotope("B11", 5, 11, 11.009305*g/mole);

  BoronEnriched->AddIsotope(B10, 20.0*perCent);
  BoronEnriched->AddIsotope(B11, 80.0*perCent);

  G4Material* EJ254 = new G4Material("EJ254", 1.026*g/cm3, 3); //5% boron loaded
  EJ254->AddElement(elC, 86.925*perCent);
  EJ254->AddElement(elH, 8.075*perCent);
  EJ254->AddElement(BoronEnriched, 5.0*perCent);

  EJ254->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);
  EJ254->GetIonisation()->SetBirksConstant(0.154*mm/MeV);

  G4int nEntries = 3;
  G4double PhotonEnergy[nEntries] = {2.48*eV, 2.88*eV, 3.1*eV}; // ~500-380 

  G4double RefractiveIndex[nEntries] = {1.58, 1.58, 1.58}; // Approximate for PVT
  G4double Scintillation[nEntries] = {0.01, 1.0, 0.1}; // Relative intensity
  G4double ScintillationYield = 7500 / MeV;  // Light yield (~10,000 photons/MeV)
  G4double ScintillationFastTime = 2.1 * ns;   // Fast decay component

  G4MaterialPropertiesTable* EJ254_MPT = new G4MaterialPropertiesTable();
  EJ254_MPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, nEntries);
  
  EJ254_MPT->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy, Scintillation, nEntries);
  EJ254_MPT->AddProperty("ABSLENGTH", {2.48*eV, 2.88*eV, 3.1*eV}, {100*cm, 100*cm, 100*cm});
  EJ254_MPT->AddConstProperty("SCINTILLATIONYIELD", ScintillationYield);
  EJ254_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Poisson statistics
  EJ254_MPT->AddConstProperty("SCINTILLATIONYIELD1", 1.0); // Single-component
  EJ254_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", ScintillationFastTime);

  EJ254->SetMaterialPropertiesTable(EJ254_MPT);
  
  G4Material* EJ256 = new G4Material("EJ256", 1.081*g/cm3, 3); //5% lead loaded 
  EJ256->AddElement(elC, 86.925*perCent);
  EJ256->AddElement(elH, 8.075*perCent);
  EJ256->AddElement(elPb, 5*perCent);

  EJ256->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);
  EJ256->GetIonisation()->SetBirksConstant(0.154*mm/MeV);
  
  G4double PhotonEnergy256[nEntries] = {2.48*eV, 2.88*eV, 3.1*eV}; // ~500-380 

  G4double RefractiveIndex256[nEntries] = {1.58, 1.58, 1.58}; // Approximate for PVT
  ScintillationYield = 5200 / MeV;  // Light yield (~10,000 photons/MeV)
  G4double Scintillation256[nEntries] = {0.01, 1.0, 0.1}; // Relative intensity
  ScintillationFastTime = 1.51 * ns;   // Fast decay component

  G4MaterialPropertiesTable* EJ256_MPT = new G4MaterialPropertiesTable();
  EJ256_MPT->AddProperty("RINDEX", PhotonEnergy256, RefractiveIndex256, nEntries);
  
  EJ256_MPT->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergy256, Scintillation256, nEntries);
  EJ256_MPT->AddProperty("ABSLENGTH", {2.48*eV, 2.88*eV, 3.1*eV}, {100*cm, 100*cm, 100*cm});
  EJ256_MPT->AddConstProperty("SCINTILLATIONYIELD", ScintillationYield);
  EJ256_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Poisson statistics
  EJ256_MPT->AddConstProperty("SCINTILLATIONYIELD1", 1.0); // Single-component
  EJ256_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", ScintillationFastTime);
  
  EJ256->SetMaterialPropertiesTable(EJ256_MPT);
  
  G4Material* EJ212 = new G4Material("EJ212", 1.023*g/cm3, 2);
  G4double numH = 5.17e+22;
  G4double numC = 4.69e+22;
  G4double massEJ212 = numH*elH->GetAtomicMassAmu()+numC*elC->GetAtomicMassAmu();
  
  EJ212->AddElement(elC, numC*elC->GetAtomicMassAmu()/massEJ212);
  EJ212->AddElement(elH, numH*elH->GetAtomicMassAmu()/massEJ212);
  
  EJ212->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);
  EJ212->GetIonisation()->SetBirksConstant(0.154*mm/MeV);

  nEntries = 3;
  G4double PhotonEnergyEJ212[nEntries] = {2.48*eV, 2.88*eV, 3.1*eV};
  G4double RefractiveIndexEJ212[nEntries] = {1.58, 1.58, 1.58};
  G4double ScintillationEJ212[nEntries] = {0.01, 1.0, 0.1};
  ScintillationYield = 10000 / MeV;
  ScintillationFastTime = 2.1 * ns;

  G4MaterialPropertiesTable* EJ212_MPT = new G4MaterialPropertiesTable();
  EJ212_MPT->AddProperty("RINDEX", PhotonEnergyEJ212, RefractiveIndexEJ212, nEntries);
  
  EJ212_MPT->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergyEJ212, ScintillationEJ212, nEntries);
  EJ212_MPT->AddProperty("ABSLENGTH", {2.48*eV, 2.88*eV, 3.1*eV}, {250*cm, 250*cm, 250*cm});
  EJ212_MPT->AddConstProperty("SCINTILLATIONYIELD", ScintillationYield);
  EJ212_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  EJ212_MPT->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  EJ212_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", ScintillationFastTime);

  EJ212->SetMaterialPropertiesTable(EJ212_MPT);

  G4Material* EJ200 = new G4Material("EJ200", 1.023*g/cm3, 2);
  numH = 5.17e+22;
  numC = 4.69e+22;
  G4double massEJ200 = numH*elH->GetAtomicMassAmu()+numC*elC->GetAtomicMassAmu();
  
  EJ200->AddElement(elC, numC*elC->GetAtomicMassAmu()/massEJ200);
  EJ200->AddElement(elH, numH*elH->GetAtomicMassAmu()/massEJ200);
  
  EJ200->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);
  EJ200->GetIonisation()->SetBirksConstant(0.154*mm/MeV);

  nEntries = 3;
  G4double PhotonEnergyEJ200[nEntries] = {2.48*eV, 2.88*eV, 3.1*eV};
  G4double RefractiveIndexEJ200[nEntries] = {1.58, 1.58, 1.58};
  G4double ScintillationEJ200[nEntries] = {0.01, 1.0, 0.1};
  ScintillationYield = 10000 / MeV;
  ScintillationFastTime = 2.1 * ns;

  G4MaterialPropertiesTable* EJ200_MPT = new G4MaterialPropertiesTable();
  EJ200_MPT->AddProperty("RINDEX", PhotonEnergyEJ200, RefractiveIndexEJ200, nEntries);
  
  EJ200_MPT->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergyEJ200, ScintillationEJ200, nEntries);
  EJ200_MPT->AddProperty("ABSLENGTH", {2.48*eV, 2.88*eV, 3.1*eV}, {380*cm, 380*cm, 380*cm});
  EJ200_MPT->AddConstProperty("SCINTILLATIONYIELD", ScintillationYield);
  EJ200_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  EJ200_MPT->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  EJ200_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", ScintillationFastTime);

  EJ200->SetMaterialPropertiesTable(EJ200_MPT);


  Teflon = new G4Material("Teflon", 2.2*g/cm3, 2);
  Teflon->AddElement(elC, 2);
  Teflon->AddElement(elF, 4);
  
  wrappingSurface = new G4OpticalSurface("WrappingSurface");
  wrappingSurface->SetType(dielectric_metal);
  wrappingSurface->SetModel(unified);
  wrappingSurface->SetFinish(polished);

  G4MaterialPropertiesTable* wrappingMPT = new G4MaterialPropertiesTable();
  std::vector<G4double> photonEnergies = {2.0*eV, 3.5*eV};
  std::vector<G4double> reflectivity = {0.98, 0.98};
  wrappingMPT->AddProperty("REFLECTIVITY", photonEnergies, reflectivity);

  wrappingSurface->SetMaterialPropertiesTable(wrappingMPT);

  G4MaterialPropertiesTable* pbwo4MPT = new G4MaterialPropertiesTable();

  nEntries = 2;
  photonEnergies = {2.0*eV, 3.5*eV}; // Energy range
  std::vector<G4double> refractiveIndex = {2.2, 2.2};  // Constant index
  std::vector<G4double> absorption = {100*cm, 100*cm}; //https://doi.org/10.1016/S0168-9002(98)00321-0
  G4double photonEnergy[nEntries] = { 2.0 * eV, 3.5 * eV };
  G4double scintillation[nEntries] = { 1.0, 1.0 };
  
  pbwo4MPT->AddProperty("RINDEX", photonEnergies, refractiveIndex);
  pbwo4MPT->AddProperty("ABSLENGTH", photonEnergies, absorption);
  pbwo4MPT->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintillation, nEntries);
  pbwo4MPT->AddConstProperty("SCINTILLATIONYIELD", 200./MeV); // Adjust yield
  pbwo4MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  pbwo4MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10.*ns);

  G4double ttt = 75.47+23.20+1.28;
  Air = new G4Material("Air", 0.001225*g/cm3,3);
  Air->AddElement(elN, 75.47/ttt);
  Air->AddElement(elO, 23.20/ttt);
  Air->AddElement(elAr, 1.28/ttt);

  // Air = nistManager->FindOrBuildMaterial("G4_Galactic"); 
  aluminumFoil = new G4Material("aluminumFoil", 2.71*g/cm3,1);
  aluminumFoil->AddElement(elAl,1);

  SiPMGlassMat = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4double refractiveIndexSiPMGlass[nEntries] = { 1.52, 1.52 };
  G4MaterialPropertiesTable* mptSiPMGlass = new G4MaterialPropertiesTable();
  mptSiPMGlass->AddProperty("RINDEX", photonEnergy, refractiveIndexSiPMGlass, nEntries);
  G4double quantumEnergies[5] = { 1.4*eV, 3.0*eV, 3.54*eV, 4.0*eV, 4.96*eV };
  G4double quantumEfficiencies[5] = { 0.06, 0.63, 0.5, 0.43, 0.06 };  // Example values
  mptSiPMGlass->AddProperty("EFFICIENCY", quantumEnergies, quantumEfficiencies, 5);

  dielectricSurface = new G4OpticalSurface("dielectricSurface");
  dielectricSurface->SetType(dielectric_dielectric);
  dielectricSurface->SetModel(unified);
  dielectricSurface->SetFinish(polished);

  SiPMGlassMat->SetMaterialPropertiesTable(mptSiPMGlass);

  SiPMSiliconMat = nistManager->FindOrBuildMaterial("G4_Si");
  G4double refractiveIndexSiPMSi[nEntries] = { 3.98, 3.98 };
  G4MaterialPropertiesTable* mptSiPMSilicon = new G4MaterialPropertiesTable();
  mptSiPMSilicon->AddProperty("RINDEX", photonEnergy, refractiveIndexSiPMSi, nEntries);
  SiPMSiliconMat->SetMaterialPropertiesTable(mptSiPMSilicon);

  G4double refractiveIndexWorld[nEntries] = { 1.0, 1.0 };
  G4MaterialPropertiesTable* mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", photonEnergy, refractiveIndexWorld, nEntries);
  worldMat->SetMaterialPropertiesTable(mptWorld);

  PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
  PMMA->AddElement(elC, 5);
  PMMA->AddElement(elH, 8);
  PMMA->AddElement(elO, 2);

  if(detectorType == "pbwo4"){
    G4cout << "Setting Materialproperties" << G4endl;
    detMaterial = nistManager->FindOrBuildMaterial("G4_PbWO4");
    detMaterial->GetIonisation()->SetMeanExcitationEnergy(600.7*eV);
    detMaterial->GetIonisation()->SetBirksConstant(0.008694);
    detMaterial->SetMaterialPropertiesTable(pbwo4MPT);
  }
  else if(detectorType == "dsb"){
    detMaterial = DSB_Gd;
  }
  else if(detectorType == "ej256"){
    detMaterial = EJ256;
  }
  else if(detectorType == "ej254"){
    detMaterial = EJ254;
  }
  else if(detectorType == "ej212"){
    detMaterial = EJ212;
  }
  else if(detectorType == "ej200"){
    detMaterial = EJ200;
  }
  else if(detectorType == "h2o"){
    detMaterial = water;
  }
  else{
    detMaterial = nistManager->FindOrBuildMaterial("G4_PbWO4");
  }

  
  if(absorberType == "pmma"){
    absorberMaterial = PMMA;
  }
  else if(absorberType == "h2o"){
    absorberMaterial = water;
  }
  else{
    absorberMaterial = water;
  }


  G4cout << "Mean excitation energy: " << detMaterial->GetIonisation()->GetMeanExcitationEnergy() << G4endl;
  G4cout << "Birks constant: " << detMaterial->GetIonisation()->GetBirksConstant() << G4endl;
}

void DetectorConstruction::fillDCMcontainer(){
  DetectorParameterisationColour* param = new DetectorParameterisationColour;

  param->SetVoxelDimensions(voxelXY/2, voxelXY/2, voxelZ/2);
  
  param->SetNoVoxels(nbofvoxelsX, nbofvoxelsY, NbOfSlices);

  std::vector<G4Material*> fMaterials = {Air, water};
  param->SetMaterials(fMaterials);

  param->SetMaterialIndices(fMateIDs);

  //logicalVoxel->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

  param->BuildContainerSolid(physContainer);

  param->CheckVoxelsFillContainer(solidContainer->GetXHalfLength(),
                                  solidContainer->GetYHalfLength(),
                                  solidContainer->GetZHalfLength());

  param->SetSkipEqualMaterials(true);
  phantom_phys = new G4PVParameterised("phantom", logicalVoxel, logicalContainer, kYAxis, nbofvoxelsX * nbofvoxelsY * NbOfSlices, param);

  phantom_phys->SetRegularStructureId(1);
}

std::vector<int> DetectorConstruction::readMatrixFromFile(const std::string &filename) {
    std::vector<std::vector<int>> matrix;
    std::ifstream infile(filename);
    std::vector<int> indices; 

    if (!infile.is_open()) {
        std::cerr << "Error opening file" << G4endl;
        return indices;
    }

    std::string line;
    while (getline(infile, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<int> row;

        while (getline(ss, item, ' ')) {
            row.push_back(stoi(item));
        }

        matrix.push_back(row);
    }
  
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            if (value == -999.0f) {
                indices.push_back(1);
            } else if (value == -1000.0f) {
                indices.push_back(0);
            } else {
                indices.push_back(static_cast<int>(value));
            }
        }
    }

    infile.close();
    return indices;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  G4double worldXY = 250*cm;
  G4double worldZ = 250*cm;
  
  G4double translation;
  G4ThreeVector physicalPosition;
  G4STL stl;
  stl.SetVerbosity(1);
  
  G4UserLimits* userLimits = new G4UserLimits();
  // userLimits->SetMaxAllowedStep(0.1*mm);
  
  solidworld = new G4Box("solidworld", worldXY / 2, worldXY / 2, worldZ / 2);
  logicalworld = new G4LogicalVolume(solidworld, worldMat, "logicalWorld");
  physworld = new G4PVPlacement(nullptr, G4ThreeVector(), logicalworld, "physworld", nullptr, false, 0, fCheckOverlaps);  
  

  G4Box* waterBox = new G4Box("Water", 50*cm, 50*cm, 50*cm);
  waterLog = new G4LogicalVolume(waterBox, water, "WaterLog");
  new G4PVPlacement(0, G4ThreeVector(0,0,-50*cm), waterLog, "WaterPhys", logicalworld, false, 0);
  
  G4double dBeamSpot = 0.1*mm;
  G4int nx, ny, nz;
  G4double cubeSizeX, cubeSizeY, cubeSizeZ;
  if (ftarget == 2){ //heterogenous
      auto waterVisAttr = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue
      waterVisAttr->SetVisibility(true);
      waterVisAttr->SetForceSolid(true);
    
      auto airVisAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // Gray
      airVisAttr->SetVisibility(true);
      airVisAttr->SetForceSolid(true);

      cubeSizeX = 0.3 * mm, cubeSizeY = 0.3 * mm;
      cubeSizeZ = static_cast<double>(pmod)/static_cast<double>(0.76923)*0.001;
      // cubeSizeZ = static_cast<double>(pmod)/static_cast<double>(0.5)*0.001;
      // cubeSizeZ = static_cast<double>(pmod)/static_cast<double>(0.7900)*0.001;
      // cubeSizeZ = static_cast<double>(pmod)/static_cast<double>(0.7355)*0.001;
      
      nx = 150;
      ny = nx;
      nz = static_cast<int>(std::round(heteroThickness/cubeSizeZ));
      // nz = 1;
      G4double tickness_real = nz * cubeSizeZ;
      G4double pmod_real = static_cast<double>(0.76923)*1000*tickness_real/nz;
      // G4double pmod_real = static_cast<double>(0.5)*1000*tickness_real/nz;
      std::cout << "pmod_real: " << pmod_real << std::endl;
      std::cout << "tickness_real: " << tickness_real << std::endl;
      auto voxelSolid = new G4Box("Voxel", cubeSizeX/2, cubeSizeY/2, cubeSizeZ/2);

      std::cout << "pmod: " << pmod << std::endl;
      std::cout << "heteroThickness: " << (nz*cubeSizeZ) << std::endl;
      std::cout << "heteroThickness/cubeSizeZ " << heteroThickness/cubeSizeZ << std::endl;
      std::cout << "std::round(heteroThickness/cubeSizeZ) " << std::round(heteroThickness/cubeSizeZ) << G4endl;

      auto voxelContainerSolid = new G4Box("VoxelContainer", (nx*cubeSizeX)/2, (ny*cubeSizeY)/2, (nz*cubeSizeZ)/2);
      auto voxelContainerLV = new G4LogicalVolume(voxelContainerSolid, Air, "VoxelContainerLV");
      new G4PVPlacement(nullptr, G4ThreeVector(0,0,heteroThickness/2+5*mm), voxelContainerLV, "VoxelContainer", logicalworld, false, 0);

      auto* parameterisation = new HeteroParameterisation(nx, ny, nz,
        cubeSizeX, cubeSizeY, cubeSizeZ,
        water, Air,
        waterVisAttr, airVisAttr);
    
      int nVoxels = nx * ny * nz;

      auto voxelLogic = new G4LogicalVolume(voxelSolid, Air, "Voxel");
      voxelLogic->SetVisAttributes(airVisAttr);

      new G4PVParameterised(
        "Voxel", voxelLogic, voxelContainerLV, kUndefined, nVoxels, parameterisation
      );

      G4Box* solidTracker = new G4Box("solidTracker", nx*cubeSizeX/2, ny*cubeSizeY/2, 1*um);
      logicalTracker = new G4LogicalVolume(solidTracker, worldMat, "logicalTracker");
      logicalTracker->SetVisAttributes(waterVisAttr);
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 4*mm), logicalTracker, "physTracker", logicalworld, false, 0);
  } 
  if(ftarget == 1){ //homogeneous
    solidHomo =  new G4Box("solidHomo", 50*cm, 50*cm, targetThickness/2);
    logicalHomo = new G4LogicalVolume(solidHomo, water, "logicalHomo");    
    physHomo = new G4PVPlacement(nullptr, G4ThreeVector(0.,0., targetThickness/2+5*mm), logicalHomo,"physHomo", logicalworld, false, 0, fCheckOverlaps);

    G4VisAttributes* visHomo = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.5)); // Red for the detector
    visHomo->SetVisibility(true);
    visHomo->SetForceSolid(true); // Ensure the detector is solid
    logicalHomo->SetVisAttributes(visHomo);
  }

  solidNozzle = new G4EllipticalTube("solidNozzle", FWHMNozzleX/2, FWHMNozzleY/2, dBeamSpot/2);
  logicalNozzle = new G4LogicalVolume(solidNozzle, worldMat, "logicalNozzle");
  logicalNozzle->SetVisAttributes(G4Color::Red());
  //logicalNozzle->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
  
  solidIsocentre = new G4EllipticalTube("solidIsocentre", FWHMIsocentreX/2, FWHMIsocentreY/2, dBeamSpot/2);
  logicalIsocentre = new G4LogicalVolume(solidIsocentre, worldMat, "logicalIsocentre");
  logicalIsocentre->SetVisAttributes(G4Color::Red());
  //logicalIsocentre->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

  physNozzle = new G4PVPlacement(nullptr, G4ThreeVector(x_off, y_off, d_NozzleIsocentre), logicalNozzle,"physNozzle", logicalworld, false, 0, fCheckOverlaps);
  physIsocentre = new G4PVPlacement(nullptr, G4ThreeVector(x_off, y_off, dBeamSpot/2), logicalIsocentre,"physIsocentre", logicalworld, false, 0, fCheckOverlaps);

  // logicalworld->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
  return physworld;
}

void DetectorConstruction::ConstructSDandField()
{  
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
  G4ScoringBox* mesh1 = new G4ScoringBox("waterMesh1");
  // G4ScoringBox* mesh2 = new G4ScoringBox("waterMesh2");
  // G4ScoringBox* mesh3 = new G4ScoringBox("waterMesh3");
  
  G4double meshSize = 25*cm;
  G4double sizeZ1 = meshSize;
  // G4double sizeZ1 = 8*cm;
  // G4double sizeZ2 = 9*cm;
  // G4double sizeZ3 = meshSize-sizeZ1-sizeZ2;

  G4double meshSize1[3] = {50*cm, 50*cm, sizeZ1};
  // G4double meshSize2[3] = {50*cm, 50*cm, sizeZ2};
  // G4double meshSize3[3] = {50*cm, 50*cm, sizeZ3};

  G4int meshSegNumber1[3] = {1, 1, fLayers};
  // G4int meshSegNumber2[3] = {1, 1, 800};
  // G4int meshSegNumber3[3] = {1, 1, 300};

  G4double meshCenter1[3] = {0,0, -sizeZ1};
  // G4double meshCenter2[3] = {0,0, -(sizeZ1*2+sizeZ2)};
  // G4double meshCenter3[3] = {0,0, -(sizeZ1*2+sizeZ2*2+sizeZ3)};
  
  mesh1->SetSize(meshSize1);
  mesh1->SetNumberOfSegments(meshSegNumber1);
  mesh1->SetCenterPosition(meshCenter1);
  mesh1->SetPrimitiveScorer(new G4PSEnergyDeposit("Edep"));
  
  // mesh2->SetSize(meshSize2);
  // mesh2->SetNumberOfSegments(meshSegNumber2);
  // mesh2->SetCenterPosition(meshCenter2);
  // mesh2->SetPrimitiveScorer(new G4PSEnergyDeposit("Edep"));

  // mesh3->SetSize(meshSize3);
  // mesh3->SetNumberOfSegments(meshSegNumber3);
  // mesh3->SetCenterPosition(meshCenter3);
  // mesh3->SetPrimitiveScorer(new G4PSEnergyDeposit("Edep"));



  scoringManager->SetVerboseLevel(1); // optional: for debug info
  scoringManager->RegisterScoringMesh(mesh1);
  // scoringManager->RegisterScoringMesh(mesh2);
  // scoringManager->RegisterScoringMesh(mesh3);
  
  return;
}

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
