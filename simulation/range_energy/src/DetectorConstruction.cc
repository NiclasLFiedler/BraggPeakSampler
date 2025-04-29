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

using namespace B2;

namespace B2a
{

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
{
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

  // detMaterial = nistManager->FindOrBuildMaterial("G4_PbWO4");
  // detMaterial->GetIonisation()->SetBirksConstant(0.008694); // mm/MeV
  
  detMaterial = nistManager->FindOrBuildMaterial("G4_WATER");
  detMaterial->GetIonisation()->SetBirksConstant(0.052*mm/MeV);

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

  G4Material* DSB = new G4Material("DSB", 3.8*g/cm3, 3);
  DSB->AddElement(elBa, 1);
  DSB->AddElement(elSi, 2);
  DSB->AddElement(elO, 5);

  G4Material* Al2O3 = new G4Material("Al2O3", 3.94*g/cm3, 2);
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  
  G4Material* EJ256 = new G4Material("EJ256", 1.081*g/cm3, 3); //5% lead loaded 
  EJ256->AddElement(elC, 86.925*perCent);
  EJ256->AddElement(elH, 8.075*perCent);
  EJ256->AddElement(elPb, 5*perCent);

  EJ256->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);

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

  G4Material* EJ212 = new G4Material("EJ212", 1.023*g/cm3, 2);
  EJ212->AddElement(elC, 91.5*perCent);
  EJ212->AddElement(elH, 8.5*perCent);
  
  EJ212->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);
  EJ212->GetIonisation()->SetBirksConstant(0.154*mm/MeV);

  G4int nEntries = 3;
  G4double PhotonEnergyEJ212[nEntries] = {2.48*eV, 2.88*eV, 3.1*eV};
  G4double RefractiveIndexEJ212[nEntries] = {1.58, 1.58, 1.58};
  G4double ScintillationEJ212[nEntries] = {0.01, 1.0, 0.1};
  G4double ScintillationYield = 10000 / MeV;
  G4double ScintillationFastTime = 2.1 * ns;

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

  Alu = new G4Material("Alu", 2.71*g/cm3,1);
  Alu->AddElement(elAl,1);

  G4double density = 1.19 * g/cm3;
  G4int ncomponents = 3;
  
  G4Material* PMMA = new G4Material("PMMA", density, ncomponents);
  PMMA->AddElement(elC, 5);
  PMMA->AddElement(elH, 8);
  PMMA->AddElement(elO, 2);

  //detMaterial = PMMA;

  G4cout <<"Mean excitation Energy: " << detMaterial->GetIonisation()->GetMeanExcitationEnergy() << G4endl;
  G4cout <<"Density: " << detMaterial->GetDensity()/(g/cm3) << G4endl;
  G4cout <<"Radiation length: " << detMaterial->GetRadlen()/cm << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  G4double worldX = 500*cm;
  G4double worldZ = 500*cm;
  G4double worldY = 500*cm;

  G4double detSizeZ = 2000*mm; //det size in x
  G4double detSizeX = 2000*mm; //det size in y
  G4double detSizeY = 2000*mm; //det size in z


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

  soliddetector = new G4Box("soliddetector",     // its name
    detSizeX/2, detSizeY/2, detSizeZ/2);             // its size
  
  logicaldetector = new G4LogicalVolume(soliddetector, detMaterial, "logicaldetector");

  G4UserLimits* userLimits = new G4UserLimits();
  //userLimits->SetMaxAllowedStep(0.1*mm);
  logicaldetector->SetUserLimits(userLimits);

  physdetector = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0.,0., detSizeZ/2+1*um),              // at (x,y,z)
    // G4ThreeVector(0.,0.,detSizeZ*i),  
    logicaldetector,            // its logical volume
    "physdetector",           // its name
    logicalworld,               // its mother volume
    false,                    // no boolean operations
    0,                        // copy number
    fCheckOverlaps);          // checking overlaps

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
  SetSensitiveDetector("logicaldetector", aTrackerSD, true);

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