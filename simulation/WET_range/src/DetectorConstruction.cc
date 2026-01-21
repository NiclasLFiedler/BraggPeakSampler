#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
//#include "G4LogicalBorderSurface.hh"
#include "nlohmann/json.hpp"
#include <filesystem>

namespace fs = std::filesystem;
using json = nlohmann::json;

using namespace B2;

namespace B2a
{

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
    heteroThickness                                   = allConfigs["heteroThickness"];
    

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
  // Define materials
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
  worldMat = nistManager->FindOrBuildMaterial("G4_AIR"); 
  // worldMat = nistManager->FindOrBuildMaterial("G4_Galactic"); 
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
  
  PTFEmembrane = new G4Material("PTFEmembrane", 0.35*g/cm3, 2);
  PTFEmembrane->AddElement(elC, 2);
  PTFEmembrane->AddElement(elF, 4);

  Teflon=PTFEmembrane;

  G4MaterialPropertiesTable* pbwo4MPT = new G4MaterialPropertiesTable();

  nEntries = 2;
  std::vector<G4double> photonEnergies = {2.0*eV, 3.5*eV}; // Energy range
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
  air = new G4Material("air", 0.001225*g/cm3,3);
  air->AddElement(elN, 75.47/ttt);
  air->AddElement(elO, 23.20/ttt);
  air->AddElement(elAr, 1.28/ttt);

  aluminumFoil = new G4Material("aluminumFoil", 2.71*g/cm3,1);
  aluminumFoil->AddElement(elAl,1);

  SiPMGlassMat = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4double refractiveIndexSiPMGlass[nEntries] = { 1.52, 1.52 };
  G4MaterialPropertiesTable* mptSiPMGlass = new G4MaterialPropertiesTable();
  mptSiPMGlass->AddProperty("RINDEX", photonEnergy, refractiveIndexSiPMGlass, nEntries);
  G4double quantumEnergies[5] = { 1.4*eV, 3.0*eV, 3.54*eV, 4.0*eV, 4.96*eV };
  G4double quantumEfficiencies[5] = { 0.06, 0.63, 0.5, 0.43, 0.06 };  // Example values
  mptSiPMGlass->AddProperty("EFFICIENCY", quantumEnergies, quantumEfficiencies, 5);

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
  auto printMat = [](G4Material* m)
  {
    G4cout << "Material: " << m->GetName() << G4endl;
    G4cout << "  Density = "
           << m->GetDensity()/(g/cm3) << " g/cm3" << G4endl;

    G4cout << "  I = "
           << m->GetIonisation()->GetMeanExcitationEnergy()/eV
           << " eV" << G4endl;

    G4cout << "  Radlen = "
           << m->GetRadlen()/cm << " cm" << G4endl;

    G4cout << "  Nuclear L = "
           << m->GetNuclearInterLength()/cm << " cm" << G4endl;

    const G4ElementVector* els = m->GetElementVector();
    const G4double* fr = m->GetFractionVector();

    for (size_t i = 0; i < m->GetNumberOfElements(); ++i) {
        G4cout << "    "
               << (*els)[i]->GetName()
               << " : " << fr[i] << G4endl;
    }
  };

  auto PbWO4_nist =
    nistManager->FindOrBuildMaterial("G4_PbWO4");
  G4cout << "=== NIST PbWO4 ===" << G4endl;
  printMat(PbWO4_nist);
  PbWO4_nist->GetIonisation()->SetMeanExcitationEnergy(600.7*eV); 

  auto PbWO4_custom = new G4Material(
    "PbWO4_custom",
    PbWO4_nist->GetDensity(),
    PbWO4_nist
  );
  
  G4cout << "=== NIST Custom PbWO4 ===" << G4endl;
  printMat(PbWO4_custom);
  

  G4cout << "Setting Materialproperties" << G4endl;

  detMaterial = PbWO4_custom;
    

  detMaterial->GetIonisation()->SetBirksConstant(0.008694);
  // detMaterial->SetMaterialPropertiesTable(pbwo4MPT);

  G4cout << "Mean excitation energy: " << detMaterial->GetIonisation()->GetMeanExcitationEnergy() << G4endl;
  G4cout << "Birks constant: " << detMaterial->GetIonisation()->GetBirksConstant() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  bool justEnergy = true;
  G4double worldX = 500*cm;
  G4double worldZ = 500*cm;
  G4double worldY = 500*cm;

  G4double WETSizeZ = 2000*mm; //det size in x
  G4double WETSizeX = 2000*mm; //det size in y
  G4double WETSizeY = 2000*mm; //det size in z

  G4double translation;
  G4ThreeVector physicalPosition;
  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World
  //G4UserLimits* userLimits = new G4UserLimits();
  //userLimits->SetMaxAllowedStep(1*mm);
  //logicaldetector->SetUserLimits(userLimits);
  
  solidworld = new G4Box("solidworld",              // its name
    worldX / 2, worldY / 2, worldZ / 2);               // its size

  logicalworld = new G4LogicalVolume(solidworld,    // its solid
    worldMat,                                                 // its material
    "logicalWorld");                                     // its name

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  physworld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    logicalworld,                                 // its logical volume
    "physworld",                         // its name
    nullptr,                                 // its mother volume
    false,                                   // no boolean operations
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps

  // Target

  solidWET = new G4Box("solidWET",     // its name
    WETSizeX/2, WETSizeY/2, WETSizeZ/2);             // its size
  
  logicalWET = new G4LogicalVolume(solidWET, water, "logicalWET");

  G4UserLimits* userLimits = new G4UserLimits();
  //userLimits->SetMaxAllowedStep(0.1*mm);
  logicalWET->SetUserLimits(userLimits);

  physWET = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0.,0., WETSizeZ/2+1*um),              // at (x,y,z)
    // G4ThreeVector(0.,0.,detSizeZ*i),  
    logicalWET,            // its logical volume
    "physWET",           // its name
    logicalworld,               // its mother volume
    false,                    // no boolean operations
    0,                        // copy number
    fCheckOverlaps);          // checking overlaps

  
  // if(!justEnergy){
  //   return physworld;
  // }
  G4double SiPMSizeYZ = 3.7*mm;
  if(SiPMSizeYZ>detSizeZ) SiPMSizeYZ = detSizeZ;
  G4double SiPMSizeX = ThicknessTeflon/2;
  G4Box* solidSiPM = new G4Box("solidSiPM", SiPMSizeX/2, SiPMSizeYZ/2, SiPMSizeYZ/2);
  logicalSiPM = new G4LogicalVolume(solidSiPM, SiPMGlassMat, "logicalSiPM");

  if(absSizeZ == 0) absSizeZ = 0.1*mm;
  solidAbsorber = new G4Box("solidAbsorber", detSizeX/2, detSizeY/2, absSizeZ/2);
  logicalAbsorber = new G4LogicalVolume(solidAbsorber, detMaterial, "logicalAbsorber");
  G4cout << "Secondary layer active? "<< secondaryLayerStatus << G4endl;
  G4Box* solidAluFoilAbs = new G4Box("solidAluFoilAbs", detSizeX/2+ThicknessTeflon/2+ThicknessAlu/2, detSizeY/2+ThicknessTeflon/2+ThicknessAlu/2, absSizeZ/2+ThicknessTeflon/2+ThicknessAlu/2);
  G4Box* solidTeflonFoilAbs = new G4Box("solidTeflonFoilAbs", detSizeX/2+ThicknessTeflon/2, detSizeY/2+ThicknessTeflon/2, absSizeZ/2+ThicknessTeflon/2);

  G4double offSet = 59*mm;
  if(fLayers <= fLayersCut){
    offSet = offSet + (solidAluFoilAbs->GetZHalfLength()*2+gapSizeZ)*fLayers;
  }
  else{
    offSet = offSet + (solidAluFoilAbs->GetZHalfLength()*2+gapSizeZ)*fLayersCut + ((detSizeZ/2+ThicknessAlu/2+ThicknessTeflon/2)*2+gapSizeZ)*(fLayers-fLayersCut);
  }


  if(secondaryLayerStatus){
    logicalAbsorber->SetUserLimits(userLimits);

    logicalTeflonFoilAbs = new G4LogicalVolume(solidTeflonFoilAbs, Teflon, "logicalTeflonFoilAbs");
    logicalAluFoilAbs = new G4LogicalVolume(solidAluFoilAbs, aluminumFoil, "logicalAluFoilAbs");

    physTeflonFoilAbs = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicalTeflonFoilAbs, "physTeflonFoilAbs", logicalAluFoilAbs, false, 0, fCheckOverlaps);
    
    physAbsorber = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicalAbsorber, "physAbsorber", logicalTeflonFoilAbs, false, 0, fCheckOverlaps);

    physSiPMAbs = new G4PVPlacement(nullptr, G4ThreeVector(solidAbsorber->GetXHalfLength()+solidSiPM->GetXHalfLength(), 0, 0), logicalSiPM, "physSiPMAbs", logicalTeflonFoilAbs, false, 0, fCheckOverlaps);

    G4VisAttributes* visAluAbs = new G4VisAttributes(G4Colour(1.0, 0.431, 0.0, 1));
    visAluAbs->SetVisibility(true);
    logicalAluFoilAbs->SetVisAttributes(visAluAbs);
    
    G4VisAttributes* visTeflonAbs = new G4VisAttributes(G4Colour(0.502, 0.0, 0.502, 1));
    visTeflonAbs->SetVisibility(true);
    logicalTeflonFoilAbs->SetVisAttributes(visTeflonAbs);

    for(G4int i=0; i<fLayersCut; i++){
      if(fLayers == i) break;
      translation = (solidAluFoilAbs->GetZHalfLength()*2+gapSizeZ)*(i)+d_IsocentreDetector;
      physicalPosition = G4ThreeVector(0.,0., translation-offSet);
      
      physAluFoilAbs = new G4PVPlacement(nullptr, physicalPosition, logicalAluFoilAbs, "physAluFoilAbs", logicalworld, false, i, fCheckOverlaps);
    }
  }

  solidDetector = new G4Box("solidDetector", detSizeX/2, detSizeY/2, detSizeZ/2);
  logicalDetector = new G4LogicalVolume(solidDetector, detMaterial, "logicalDetector");
  logicalDetector->SetUserLimits(userLimits);
  
  G4Box* solidAluFoil = new G4Box("solidAluFoil", detSizeX/2+ThicknessAlu/2 +ThicknessTeflon/2, detSizeY/2+ThicknessAlu/2+ThicknessTeflon/2, detSizeZ/2+ThicknessAlu/2+ThicknessTeflon/2);
  G4Box* solidTeflonFoil = new G4Box("solidTeflonFoil", detSizeX/2+ThicknessTeflon/2, detSizeY/2+ThicknessTeflon/2, detSizeZ/2+ThicknessTeflon/2);
  
  logicalAluFoil = new G4LogicalVolume(solidAluFoil , aluminumFoil, "logicalAluFoil");
  logicalTeflonFoil = new G4LogicalVolume(solidTeflonFoil, Teflon, "logicalTeflonFoil");
  
  physTeflonFoil = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicalTeflonFoil, "physTeflonFoil", logicalAluFoil, false, 0, fCheckOverlaps);
  
  physDetector = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), logicalDetector, "physDetector", logicalTeflonFoil, false, 0, fCheckOverlaps);
  
  physSiPM = new G4PVPlacement(nullptr, G4ThreeVector(solidDetector->GetXHalfLength()+solidSiPM->GetXHalfLength(), 0, 0), logicalSiPM, "physSiPM", logicalTeflonFoil, false, 0, fCheckOverlaps);

  G4double passiveFill = 0;
  for(G4int i=0; i<fLayers-fLayersCut; i++){
    if(absorberStatus) passiveFill = absorberSize+gapSizeZ;
    G4double offSetChannel = offSet + fLayersCut*(solidAluFoilAbs->GetZHalfLength()*2+gapSizeZ)+ (solidAluFoil->GetZHalfLength()*2+gapSizeZ)*(i);
    translation = (solidAluFoilAbs->GetZHalfLength()*2+gapSizeZ)*fLayersCut + (solidAluFoil->GetZHalfLength()*2+gapSizeZ)*(i)+d_IsocentreDetector+passiveFill;
    physicalPosition = G4ThreeVector(0.,0., translation-offSet);
    if(i == 0){
      physicalPosition = G4ThreeVector(0.,0., translation-passiveFill-offSet);
    }
    physAluFoil = new G4PVPlacement(nullptr, physicalPosition, logicalAluFoil, "physAluFoil", logicalworld, false, i+fLayersCut, fCheckOverlaps);
  }  

  G4VisAttributes* visDetector = new G4VisAttributes(G4Colour(1.0, 0.84, 0.0, 0.5));
  visDetector->SetVisibility(true);
  visDetector->SetForceWireframe(true);
  //visDetector->SetForceSolid(true);

  G4VisAttributes* visAbsorber = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.5));
  visAbsorber->SetVisibility(true);
  visAbsorber->SetForceWireframe(true);
  //visAbsorber->SetForceSolid(true);

  G4VisAttributes* visTeflonFoil = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.3));
  visTeflonFoil->SetVisibility(true);
  visTeflonFoil->SetForceWireframe(true);
  //visTeflonFoil->SetForceSolid(true);

  G4VisAttributes* visAluFoil = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  visAluFoil->SetVisibility(true);
  visAluFoil->SetForceWireframe(true);

  // Apply visualization attributes to logical volumes
  logicalDetector->SetVisAttributes(visDetector);
  logicalAbsorber->SetVisAttributes(visAbsorber);
  logicalAluFoil->SetVisAttributes(visAluFoil);
  logicalTeflonFoil->SetVisAttributes(visTeflonFoil);

  //Always return the physical world
  //logicalworld->SetVisAttributes(G4VisAttributes);
  return physworld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerDetectorSDname = "/TrackerDetectorSD";
  auto aTrackerSD = new TrackerSD(trackerDetectorSDname, "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name
  // of "Chamber_LV".
  SetSensitiveDetector("logicalWET", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}