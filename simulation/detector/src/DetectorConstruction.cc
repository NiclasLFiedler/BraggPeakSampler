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
#include "G4RegularNavigation.hh"

#include "Randomize.hh"

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
    
    const auto& config = allConfigs["detectors"][int(allConfigs["detectorSelect"])];
    
    ftarget                      = allConfigs["targetSelect"];

    detectorType            = config["detectorType"];
    //double beamEnergy                   = config["beamEnergy"]*MeV;
    fLayers                         = config["nLayers"];
    std::vector<double> crystalSize     = config["crystalSize"];
    
    detSizeX = crystalSize.at(0) * mm;
    detSizeY = crystalSize.at(1)* mm;
    detSizeZ = crystalSize.at(2)* mm;

    gapSizeZ                     = config["gapSizeZ"];
    secondaryLayerStatus           = config["secondaryLayerStatus"];
    fLayersCut                      = config["nSecondaryLayers"];
    absSizeZ                        = config["secLayerSizeZ"];
    
    absorberStatus                 = config["absorberStatus"];
    absorberSize                 = config["absorberSize"];
    std::vector<double> teflonThickness = config["teflonThickness"];
    std::vector<double> aluThickness    = config["aluThickness"];
    pmod                 = allConfigs["pmod"];
    heteroThickness                 = allConfigs["heteroThickness"];
    
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

  homoMaterial = nistManager->FindOrBuildMaterial("G4_WATER"); 
  homoMaterial->GetIonisation()->SetMeanExcitationEnergy(79.7*eV);

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
  EJ212->AddElement(elC, 91.5*perCent);
  EJ212->AddElement(elH, 8.5*perCent);
  
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
  heteroAir = new G4Material("heteroAir", 0.001225*g/cm3,3);
  heteroAir->AddElement(elN, 75.47/ttt);
  heteroAir->AddElement(elO, 23.20/ttt);
  heteroAir->AddElement(elAr, 1.28/ttt);

  heteroWater = new G4Material("heteroWater", 1.00*g/cm3,2);
  heteroWater->AddElement(elH, 2);
  heteroWater->AddElement(elO, 1);
  heteroWater->GetIonisation()->SetMeanExcitationEnergy(79.7*eV);

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
  else if(detectorType == "h2o"){
    detMaterial = heteroWater;
  }
  else{
    detMaterial = nistManager->FindOrBuildMaterial("G4_PbWO4");
  }
  // detMaterial=DSB_Gd;
  G4cout << "Mean excitation energy: " << detMaterial->GetIonisation()->GetMeanExcitationEnergy() << G4endl;
  G4cout << "Birks constant: " << detMaterial->GetIonisation()->GetBirksConstant() << G4endl;
}

void DetectorConstruction::fillDCMcontainer(){
  DetectorParameterisationColour* param = new DetectorParameterisationColour;

  param->SetVoxelDimensions(voxelXY/2, voxelXY/2, voxelZ/2);
  
  param->SetNoVoxels(nbofvoxelsX, nbofvoxelsY, NbOfSlices);

  std::vector<G4Material*> fMaterials = {heteroAir, heteroWater};
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
  
  G4double ThicknessTeflon = 1*nm;
  G4double ThicknessAlu = 1*nm;
  
  G4double SiPMSizeYZ = 3.7*mm;
  if(SiPMSizeYZ>detSizeZ) SiPMSizeYZ = detSizeZ;
  G4double SiPMSizeX = ThicknessTeflon+ThicknessAlu;
  G4Box* solidSiPM = new G4Box("solidSiPM", SiPMSizeX/2, SiPMSizeYZ/2, SiPMSizeYZ/2);
  
  G4double dz  = 100.0 * mm;
  G4double dx1 = detSizeZ;
  G4double dx2 = detSizeZ;
  G4double dy1 = detSizeY;
  G4double dx3 = SiPMSizeYZ;
  G4double dx4 = SiPMSizeYZ;
  G4double dy2 = SiPMSizeYZ;

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
  //rotationMatrix->rotateZ(90.0 * deg);

  G4LogicalVolume *logicalLGAlu;
  G4LogicalVolume *logicalLGTeflon;
  G4LogicalSkinSurface *skinLGTeflon;
  G4LogicalVolume *logicalLightGuide;
  
  // G4Trap* solidLGAlu = new G4Trap("solidLGAlu", dz/2, 0.0, 0.0, dy1/2+ThicknessAlu+ThicknessTeflon, dx1/2+ThicknessAlu+ThicknessTeflon, dx2/2+ThicknessAlu+ThicknessTeflon, 0.0, dy2/2+ThicknessAlu+ThicknessTeflon, dx3/2+ThicknessAlu+ThicknessTeflon, dx4/2+ThicknessAlu+ThicknessTeflon, 0.0);
  // G4Trap* solidLGTeflon = new G4Trap("solidLGTeflon", dz/2, 0.0, 0.0, dy1/2+ThicknessTeflon, dx1/2+ThicknessTeflon, dx2/2+ThicknessTeflon, 0.0, dy2/2+ThicknessTeflon, dx3/2+ThicknessTeflon, dx4/2+ThicknessTeflon, 0.0);
  // G4Trap* solidLightGuide = new G4Trap("solidLightGuide", dz/2, 0.0, 0.0, dy1/2, dx1/2, dx2/2, 0.0, dy2/2, dx3/2, dx4/2, 0.0);
  // logicalLGAlu = new G4LogicalVolume(solidLGAlu , aluminumFoil, "logicalLGAlu");
  // logicalLGTeflon = new G4LogicalVolume(solidLGTeflon , Teflon, "logicalLGTeflon");
  // skinLGTeflon = new G4LogicalSkinSurface("skinLGTeflon", logicalLGTeflon, wrappingSurface);
  // logicalLightGuide = new G4LogicalVolume(solidLightGuide , SiPMGlassMat, "logicalLightGuide");
  
  // G4TessellatedSolid* solidLGAlu = stl.Read("../constructs/lightguide_alu.stl");
  // G4TessellatedSolid* solidLGTeflon = stl.Read("../constructs/lightguide_teflon.stl");
  // G4TessellatedSolid* solidLightGuide = stl.Read("../constructs/lightguide.stl");
  
  G4Box* solidSiPMSub = new G4Box("solidSiPMSub", SiPMSizeX/2+5*mm, SiPMSizeYZ/2, SiPMSizeYZ/2);
  logicalSiPM = new G4LogicalVolume(solidSiPM, SiPMGlassMat, "logicalSiPM");
  G4Box* solidLGAluFull = new G4Box("solidLGAluFull", 0.1*mm/2+ThicknessAlu/2+ThicknessTeflon/2, detSizeY/2+ThicknessAlu+ThicknessTeflon, detSizeZ/2+ThicknessAlu+ThicknessTeflon);
  G4Box* solidLGTeflonFull = new G4Box("solidLGTeflonFull", 0.1*mm/2+ThicknessTeflon/2, detSizeY/2+ThicknessTeflon, detSizeZ/2+ThicknessTeflon);
  
  G4SubtractionSolid* solidLGAlu =  new G4SubtractionSolid("solidLGAlu", solidLGAluFull, solidSiPMSub, rotationMatrix, G4ThreeVector(0,0,0));
  G4SubtractionSolid* solidLGTeflon =  new G4SubtractionSolid("solidLGTeflon", solidLGTeflonFull, solidSiPMSub, rotationMatrix, G4ThreeVector(0,0,0));
  
  logicalLGAlu = new G4LogicalVolume(solidLGAlu , aluminumFoil, "logicalLGAlu");
  logicalLGTeflon = new G4LogicalVolume(solidLGTeflon , Teflon, "logicalLGTeflon");
  skinLGTeflon = new G4LogicalSkinSurface("skinLGTeflon", logicalLGTeflon, wrappingSurface);
  // logicalLightGuide = new G4LogicalVolume(solidLightGuide , SiPMGlassMat, "logicalLightGuide");
  // G4LogicalSkinSurface *skinLightGuide = new G4LogicalSkinSurface("skinLightGuide", logicalDetector, dielectricSurface);
  
  G4VisAttributes* visLGAlu = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 1));
  visLGAlu->SetVisibility(true);
  // visLGAlu->SetForceSolid(true);
  logicalLGAlu->SetVisAttributes(visLGAlu);
  G4VisAttributes* visLGTeflon = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1));
  visLGTeflon->SetVisibility(true);
  // visLGTeflon->SetForceSolid(true);
  logicalLGTeflon->SetVisAttributes(visLGTeflon);

  // G4VisAttributes* visLightGuide = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 1));
  // visLightGuide->SetVisibility(true);
  // visLightGuide->SetForceSolid(true);
  // logicalLightGuide->SetVisAttributes(visLightGuide);
  
  //G4Box
  G4PVPlacement* physLGTeflon = new G4PVPlacement(nullptr, G4ThreeVector(-(ThicknessAlu)/2, 0, 0), logicalLGTeflon, "physLGTeflon", logicalLGAlu, false, 0, fCheckOverlaps);
  //G4PVPlacement* physLightGuide = new G4PVPlacement(nullptr, G4ThreeVector(-(ThicknessTeflon)/2, 0, 0), logicalLightGuide, "physLightGuide", logicalLGTeflon, false, 0, fCheckOverlaps);
  
  //Lightguide
  // G4PVPlacement* physLGTeflon = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, ThicknessAlu), logicalLGTeflon, "physLGTeflon", logicalLGAlu, false, 0, fCheckOverlaps);
  // G4PVPlacement* physLightGuide = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, ThicknessTeflon), logicalLightGuide, "physLightGuide", logicalLGTeflon, false, 0, fCheckOverlaps);
  if(absSizeZ == 0) absSizeZ = 0.1*mm;
  solidAbsorber = new G4Box("solidAbsorber", detSizeX/2, detSizeY/2, absSizeZ/2);
  logicalAbsorber = new G4LogicalVolume(solidAbsorber, detMaterial, "logicalAbsorber");
  G4cout << "Secondary layer is active "<< secondaryLayerStatus << G4endl;
  if(secondaryLayerStatus){
    logicalAbsorber->SetUserLimits(userLimits);
    
    G4Box* solidLGAluFullAbs = new G4Box("solidLGAluFullAbs", 0.1*mm/2+ThicknessAlu/2+ThicknessTeflon/2, detSizeY/2+ThicknessAlu+ThicknessTeflon, absSizeZ/2+ThicknessAlu+ThicknessTeflon);
    G4Box* solidLGTeflonFullAbs = new G4Box("solidLGTeflonFullAbs", 0.1*mm/2+ThicknessTeflon/2, detSizeY/2+ThicknessTeflon, absSizeZ/2+ThicknessTeflon);
    G4SubtractionSolid* solidLGAluAbs =  new G4SubtractionSolid("solidLGAluAbs", solidLGAluFullAbs, solidSiPMSub, rotationMatrix, G4ThreeVector(0,0,0));
    G4SubtractionSolid* solidLGTeflonAbs =  new G4SubtractionSolid("solidLGTeflonAbs", solidLGTeflonFullAbs, solidSiPMSub, rotationMatrix, G4ThreeVector(0,0,0));

    G4LogicalVolume* logicalLGAluAbs = new G4LogicalVolume(solidLGAluAbs , aluminumFoil, "logicalLGAluAbs");
    G4LogicalVolume* logicalLGTeflonAbs = new G4LogicalVolume(solidLGTeflonAbs , Teflon, "logicalLGTeflonAbs");
    G4LogicalSkinSurface* skinLGTeflonAbs = new G4LogicalSkinSurface("skinLGTeflonAbs", logicalLGTeflonAbs, wrappingSurface);
    G4PVPlacement* physLGTeflonAbs = new G4PVPlacement(nullptr, G4ThreeVector(-(ThicknessAlu)/2, 0, 0), logicalLGTeflonAbs, "physLGTeflonAbs", logicalLGAluAbs, false, 0, fCheckOverlaps);

    G4Box* solidTeflonFoilAbs = new G4Box("solidTeflonFoilAbs", (detSizeX+ThicknessTeflon)/2, detSizeY/2+ThicknessTeflon, absSizeZ/2+ThicknessTeflon);
    G4Box* solidAluFoilAbs = new G4Box("solidAluFoilAbs", (detSizeX+ThicknessTeflon+ThicknessAlu)/2, detSizeY/2+ThicknessTeflon+ThicknessAlu, absSizeZ/2+ThicknessTeflon+ThicknessAlu);
    
    logicalTeflonFoilAbs = new G4LogicalVolume(solidTeflonFoilAbs, Teflon, "logicalTeflonFoilAbs");
    G4LogicalSkinSurface* skinTeflonAbs = new G4LogicalSkinSurface("skinTeflonAbs", logicalTeflonFoilAbs, wrappingSurface);
    logicalAluFoilAbs = new G4LogicalVolume(solidAluFoilAbs, aluminumFoil, "logicalAluFoilAbs");

    physTeflonFoilAbs = new G4PVPlacement(nullptr, G4ThreeVector((ThicknessAlu)/2, 0, 0), logicalTeflonFoilAbs, "physTeflonFoilAbs", logicalAluFoilAbs, false, 0, fCheckOverlaps);
    physAbsorber = new G4PVPlacement(nullptr, G4ThreeVector((ThicknessTeflon)/2, 0, 0), logicalAbsorber, "physAbsorber", logicalTeflonFoilAbs, false, 0, fCheckOverlaps);

    G4VisAttributes* visLGAluAbs = new G4VisAttributes(G4Colour(1.0, 0.431, 0.0, 1));
    visLGAluAbs->SetVisibility(true);
    // visLGAluAbs->SetForceSolid(true);
    logicalLGAluAbs->SetVisAttributes(visLGAluAbs);
    logicalAluFoilAbs->SetVisAttributes(visLGAluAbs);
    
    G4VisAttributes* visLGTeflonAbs = new G4VisAttributes(G4Colour(0.502, 0.0, 0.502, 1));
    visLGTeflonAbs->SetVisibility(true);
    // visLGTeflonAbs->SetForceSolid(true);
    logicalLGTeflonAbs->SetVisAttributes(visLGTeflonAbs);
    logicalTeflonFoilAbs->SetVisAttributes(visLGTeflonAbs);

    for(G4int i=0; i<fLayersCut; i++){
      if(fLayers == i) break;
      translation = absSizeZ*(i)+(i)*gapSizeZ+d_IsocentreDetector;
      physicalPosition = G4ThreeVector(0.,0., -translation);
      
      physAluFoilAbs = new G4PVPlacement(nullptr, physicalPosition, logicalAluFoilAbs, "physAluFoilAbs", logicalworld, false, i, fCheckOverlaps);
      G4PVPlacement* physLGAluAbs = new G4PVPlacement(rotationMatrix, physicalPosition+G4ThreeVector(solidAluFoilAbs->GetXHalfLength()+solidLGAluFullAbs->GetXHalfLength(), 0, 0), logicalLGAluAbs, "physLGAluAbs", logicalworld, false, i+fLayersCut, fCheckOverlaps);
      physSiPMAbs = new G4PVPlacement(nullptr, physicalPosition+G4ThreeVector(solidAluFoilAbs->GetXHalfLength()+solidSiPM->GetXHalfLength(), 0, 0), logicalSiPM, "physSiPMAbs", logicalworld, false, i, fCheckOverlaps);
    }
  }
  solidDetector = new G4Box("solidDetector", detSizeX/2, detSizeY/2, detSizeZ/2);
  logicalDetector = new G4LogicalVolume(solidDetector, detMaterial, "logicalDetector");
  logicalDetector->SetUserLimits(userLimits);

  G4Box* solidAluFoil = new G4Box("solidAluFoil", (detSizeX+ThicknessAlu+ThicknessTeflon)/2, detSizeY/2+ThicknessAlu+ThicknessTeflon, detSizeZ/2+ThicknessAlu+ThicknessTeflon);
  logicalAluFoil = new G4LogicalVolume(solidAluFoil , aluminumFoil, "logicalAluFoil");
  G4Box* solidTeflonFoil = new G4Box("solidTeflonFoil", (detSizeX+ThicknessTeflon)/2, detSizeY/2+ThicknessTeflon, detSizeZ/2+ThicknessTeflon);
  logicalTeflonFoil = new G4LogicalVolume(solidTeflonFoil, Teflon, "logicalTeflonFoil");
  G4LogicalSkinSurface *skinTeflon = new G4LogicalSkinSurface("skinTeflonFoil", logicalTeflonFoil, wrappingSurface);
  physTeflonFoil = new G4PVPlacement(nullptr, G4ThreeVector((ThicknessAlu)/2, 0, 0), logicalTeflonFoil, "physTeflonFoil", logicalAluFoil, false, 0, fCheckOverlaps);
  physDetector = new G4PVPlacement(nullptr, G4ThreeVector((ThicknessTeflon)/2, 0, 0), logicalDetector, "physDetector", logicalTeflonFoil, false, 0, fCheckOverlaps);

  //DetectorGuideBoarder = new G4LogicalBorderSurface("DetectorGuideBoarder", physDetector, physTeflonFoil, dielectricSurface);
  //GuideDetectorBoarder = new G4LogicalBorderSurface("GuideDetectorBoarder", physTeflonFoil, physDetector, dielectricSurface);
  
  G4double passiveFill = 0;
  std::cout << "AbsorberSize: " << absorberSize << std::endl;
  for(G4int i=0; i<fLayers-fLayersCut; i++){
    if(absorberStatus) passiveFill = absorberSize+gapSizeZ;
    translation = (absSizeZ+gapSizeZ)*fLayersCut + detSizeZ*(i)+(i)*gapSizeZ+d_IsocentreDetector+passiveFill;
    physicalPosition = G4ThreeVector(0.,0., -translation);
    if(i == 0){
      physicalPosition = G4ThreeVector(0.,0., -translation+passiveFill);
    }
    G4cout << "Translation: " << translation << G4endl;
    physAluFoil = new G4PVPlacement(nullptr, physicalPosition, logicalAluFoil,"physAluFoil", logicalworld, false, i+fLayersCut, fCheckOverlaps);
    physLGAlu = new G4PVPlacement(rotationMatrix, physicalPosition+G4ThreeVector(solidAluFoil->GetXHalfLength()+solidLGAluFull->GetXHalfLength(), 0, 0), logicalLGAlu, "physLGAlu", logicalworld, false, i+fLayersCut, fCheckOverlaps);
    physSiPM = new G4PVPlacement(nullptr, physicalPosition+G4ThreeVector(solidAluFoil->GetXHalfLength()+solidSiPM->GetXHalfLength(), 0, 0), logicalSiPM, "physSiPM", logicalworld, false, i+fLayersCut, fCheckOverlaps);
  
    //lightguide
    // physLGAlu = new G4PVPlacement(rotationMatrix, physicalPosition+G4ThreeVector(solidAluFoil->GetXHalfLength(), 0, -solidLGAlu->GetMaxZExtent()/2), logicalLGAlu, "physLGAlu", logicalworld, false, i+fLayersCut, fCheckOverlaps);
    // physSiPM = new G4PVPlacement(nullptr, physicalPosition+G4ThreeVector(solidAluFoil->GetXHalfLength()+solidLGAlu->GetMaxYExtent()+solidSiPM->GetXHalfLength()), logicalSiPM, "physSiPM", logicalworld, false, i+fLayersCut, fCheckOverlaps);
    //GuideSiPMBoarder = new G4LogicalBorderSurface("GuideSiPMBoarder", physLightGuide, physSiPM,  dielectricSurface);
    //SiPMGuideBoarder = new G4LogicalBorderSurface("SiPMGuideBoarder", physSiPM, physLightGuide, dielectricSurface);
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
  
  G4String holder="../constructs/halterung.stl";
  if(fLayers == 15){
    solidHolder = stl.Read(holder); // halterungsspacer 6 mm; spacer 2 mm; seitenanfang 15.25 mm & 6.75mm
    logicalHolder = new G4LogicalVolume(solidHolder, resinMaterial, "logicalHolder");
    G4cout << "Max Extent " << solidHolder->GetMaxZExtent() << G4endl;
    translation = d_IsocentreDetector+solidHolder->GetMaxZExtent()/2+16.725*mm+8.475*mm;
    physHolder = new G4PVPlacement(nullptr, G4ThreeVector(0.,0., -translation), logicalHolder,"physHolder", logicalworld, false, 0, fCheckOverlaps);    

    G4VisAttributes* visHolder = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.6));
    visHolder->SetVisibility(true);
    visHolder->SetForceSolid(true); // Ensure the detector is solid
    logicalHolder->SetVisAttributes(visHolder);
  }
  G4double dBeamSpot = 0.1*mm;
  

  if (ftarget == 2){ //heterogenous
    bool test = false;
    if(test){
      std::string path = "../constructs/dicom_conv/" + std::to_string(pmod) + "mu_" + std::to_string(heteroThickness) + "mm";
      int file_count = 0;

      try {
          for (const auto& entry : fs::directory_iterator(path)) {
              if (fs::is_regular_file(entry.status())) {
                  ++file_count;
              }
          }
          std::cout << "Number of files: " << file_count << std::endl;
      } catch (const fs::filesystem_error& e) {
          std::cerr << "Filesystem error: " << e.what() << std::endl;
      }

      voxelXY = 0.5*mm;
      NbOfSlices = file_count;
      voxelZ = heteroThickness/static_cast<double>(NbOfSlices);

      //voxelZ = static_cast<double>(pmod)/static_cast<double>(0.75)*0.001;
      //NbOfSlices = static_cast<int>(std::round(heteroThickness/voxelZ));

      std::cout << "pmod: " << pmod << std::endl;
      std::cout << "heteroThickness: " << heteroThickness << std::endl;
      std::cout << "voxelZ: " << voxelZ << std::endl;
      std::cout << "NbOfSlices: " << NbOfSlices << " count "<< file_count<< std::endl;


      nbofvoxelsX = 100;
      nbofvoxelsY = 100;
      solidContainer = new G4Box("solidContainer", voxelXY*nbofvoxelsX/2, voxelXY*nbofvoxelsY/2, voxelZ*NbOfSlices/2);

      logicalContainer = new G4LogicalVolume(solidContainer, worldMat, "logicalContainer",0,0,0);
      fMateIDs = new size_t[nbofvoxelsX*nbofvoxelsY*NbOfSlices];

      solidVoxel = new G4Box("solidVoxel", voxelXY/2, voxelXY/2, voxelZ/2);
      logicalVoxel = new G4LogicalVolume(solidVoxel, lungTissue, "logicalVoxel", 0, 0, 0);

      physContainer = new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer", logicalworld, false, 0, fCheckOverlaps);     
      // new G4PVPlacement(nullptr, G4ThreeVector(voxelXY*nbofvoxelsX,0.,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer1", logicalworld, false, 1, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(-voxelXY*nbofvoxelsX,0.,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer2", logicalworld, false, 2, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(0.,voxelXY*nbofvoxelsY,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer3", logicalworld, false, 3, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(0.,-voxelXY*nbofvoxelsY,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer4", logicalworld, false, 4, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(voxelXY*nbofvoxelsX,voxelXY*nbofvoxelsY,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer5", logicalworld, false, 5, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(-voxelXY*nbofvoxelsX,voxelXY*nbofvoxelsY,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer6", logicalworld, false, 6, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(voxelXY*nbofvoxelsX,-voxelXY*nbofvoxelsY,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer7", logicalworld, false, 7, fCheckOverlaps);
      // new G4PVPlacement(nullptr, G4ThreeVector(-voxelXY*nbofvoxelsX,-voxelXY*nbofvoxelsY,voxelZ*NbOfSlices/2+dBeamSpot/2), logicalContainer,"physContainer8", logicalworld, false, 8, fCheckOverlaps);

      for(int i = 0; i < NbOfSlices; i++){
        //std::vector<int> indices = readMatrixFromFile("../constructs/dicom_conv/DICOM_Waterphantom" + std::to_string(i+1) + ".txt");  
        std::vector<int> indices = readMatrixFromFile("../constructs/dicom_conv/" + std::to_string(100) + "mu_" + std::to_string(200) + "mm" + "/DICOM_Waterphantom" + std::to_string(i+1) + ".txt");
        for(int j=0; j<nbofvoxelsX*nbofvoxelsY; j++){                                  
          fMateIDs[i*nbofvoxelsX*nbofvoxelsY+j] = indices.at(j);
        }      
      }
      fillDCMcontainer();
    }
    else{
      auto waterVisAttr = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue
      waterVisAttr->SetVisibility(true);
      waterVisAttr->SetForceSolid(true);

      auto airVisAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // Gray
      airVisAttr->SetVisibility(true);
      airVisAttr->SetForceSolid(true);

      G4double cubeSizeX = 0.2 * mm, cubeSizeY = 0.2 * mm, cubeSizeZ = static_cast<double>(pmod)/static_cast<double>(0.75)*0.001;
      G4int nx = 100, ny = 100, nz = static_cast<int>(std::round(heteroThickness/cubeSizeZ));
      auto voxelSolid = new G4Box("Voxel", cubeSizeX/2, cubeSizeY/2, cubeSizeZ/2);

      std::cout << "pmod: " << pmod << std::endl;
      std::cout << "heteroThickness: " << heteroThickness << std::endl;
      std::cout << "voxelZ: " << cubeSizeZ << std::endl;
      std::cout << "NbOfSlices: " << nz << G4endl;

      // Place voxels in a grid
      for (G4int ix = 0; ix < nx; ix++) {
          for (G4int iy = 0; iy < ny; iy++) {
              for (G4int iz = 0; iz < nz; iz++) {
                  // Randomly assign material
                  G4Material* mat = (G4UniformRand() < 0.25) ? homoMaterial : worldMat;

                  // Define voxel logical volume with assigned material
                  auto voxelLogic = new G4LogicalVolume(voxelSolid, mat, "Voxel");
                  if (mat == homoMaterial) {
                    voxelLogic->SetVisAttributes(waterVisAttr);
                  } else {
                    voxelLogic->SetVisAttributes(airVisAttr);
                  }
                  G4ThreeVector position(
                      (ix - nx/2) * cubeSizeX,
                      (iy - ny/2) * cubeSizeY,
                      (iz - nz/2) * cubeSizeZ + cubeSizeZ*nz/2 + 2*cm
                  );

                  // Unique copy number if needed
                  G4int copyNo = ix * ny * nz + iy * nz + iz;

                  new G4PVPlacement(nullptr, position, voxelLogic, "Voxel", logicalworld, false, copyNo);
              }
          }
      }
    }
  } 
  if(ftarget == 1 || absorberStatus){ //homogeneous
    G4double targetX = solidAluFoil->GetXHalfLength()*2;
    G4double targetY = solidAluFoil->GetYHalfLength()*2;

    solidHomo =  new G4Box("solidHomo", targetX/2, targetY/2, heteroThickness/2);
    logicalHomo = new G4LogicalVolume(solidHomo, homoMaterial, "logicalHomo");
    
    if(absorberSize == 0) absorberSize = 0.1;
    G4Box *solidAbsorberPlate =  new G4Box("solidAbsorberPlate", targetX/2, targetY/2, absorberSize/2);
    G4LogicalVolume* logicalAbsorberPlate = new G4LogicalVolume(solidAbsorberPlate, homoMaterial, "logicalAbsorberPlate");
    
    if(absorberStatus){
      G4PVPlacement* physAbsorberPlate = new G4PVPlacement(nullptr, G4ThreeVector(0.,0., -d_IsocentreDetector-absorberSize/2-detSizeZ/2-gapSizeZ), logicalAbsorberPlate,"physAbsorberPlate", logicalworld, false, 0, fCheckOverlaps);
      //physHomo = new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,absorberSize/2+dBeamSpot/2), logicalHomo,"physHomo", logicalworld, false, 0, fCheckOverlaps);
    }
    else{
      physHomo = new G4PVPlacement(nullptr, G4ThreeVector(0.,0., heteroThickness/2+dBeamSpot/2), logicalHomo,"physHomo", logicalworld, false, 0, fCheckOverlaps);
    }
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
  physIsocentre = new G4PVPlacement(nullptr, G4ThreeVector(x_off, y_off, 0.), logicalIsocentre,"physIsocentre", logicalworld, false, 0, fCheckOverlaps);

  // logicalworld->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
  return physworld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  std::string trackerDetectorSDname = "/TrackerDetectorSD";
  auto aTrackerSD = new TrackerSD(trackerDetectorSDname, "TrackerHitsCollection", fLayers);
  auto aSiPMSD = new SiPMSD("/SiPMSD", "TrackerHitsCollectionSiPM", fLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  G4SDManager::GetSDMpointer()->AddNewDetector(aSiPMSD);

  SetSensitiveDetector("logicalDetector", aTrackerSD, true);
  SetSensitiveDetector("logicalAbsorber", aTrackerSD, true);
  SetSensitiveDetector("logicalSiPM", aSiPMSD, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
