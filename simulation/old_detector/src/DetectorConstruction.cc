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
#include "DetectorParameterisationColour.hh"
#include "G4RegularNavigation.hh"

using namespace B2;

namespace B2a
{

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
{
  voxelXY = 0.5*mm;
  voxelZ = 0.2932*mm; //0.278 for pmod = 209um
  NbOfSlices = 781; //
  nbofvoxelsX = 100;
  nbofvoxelsY = 100;

  fNbOfDetectors = 15;
  fNbOfDetectorsPbWO4 = 15;
  fNbOfDetectorsPEN = fNbOfDetectors-fNbOfDetectorsPbWO4;

  #ifdef _WIN32
      std::cout << "Running on Windows, getting environmental variables" << std::endl;
      const char* cNbOfDetectors = std::getenv("ENV_NBDETECTORS");
      const char* cNbOfDetectorsPbWO4 = std::getenv("ENV_NBDETECTORSPWO");
  #elif __linux__
      std::cout << "Running on Linux, getting environmental variables" << std::endl;
      const char* cNbOfDetectors = std::getenv("ENV_NBDETECTORS");
      const char* cNbOfDetectorsPbWO4 = std::getenv("ENV_NBDETECTORSPWO");
  #elif __APPLE__
      std::cout << "Running on macOS, getting environmental variables" << std::endl;
      const char* cNbOfDetectors = std::getenv("ENV_NBDETECTORS");
      const char* cNbOfDetectorsPbWO4 = std::getenv("ENV_NBDETECTORSPWO");
  #else
      std::cout << "Unknown operating system, not getting environmental variables using defaults" << std::endl;
  #endif

  if (cNbOfDetectors && cNbOfDetectorsPbWO4) {
    try {
        fNbOfDetectors = std::stoi(cNbOfDetectors);
        fNbOfDetectorsPbWO4 = std::stoi(cNbOfDetectorsPbWO4);
        std::cout << "ENV_NBDETECTORS: " << fNbOfDetectors << std::endl;
        std::cout << "ENV_NBDETECTORSPWO: " << fNbOfDetectorsPbWO4 << std::endl;
        fNbOfDetectorsPEN = fNbOfDetectors-fNbOfDetectorsPbWO4;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
        return;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << std::endl;
        return;
    }
  } 
  else {
    fNbOfDetectors = 15;
    fNbOfDetectorsPbWO4 = 15;
    fNbOfDetectorsPEN = fNbOfDetectors-fNbOfDetectorsPbWO4;
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
  G4Element *elC = nistManager->FindOrBuildElement("C");    // 6
  G4Element *elN = nistManager->FindOrBuildElement("N");    // 7
  G4Element *elO = nistManager->FindOrBuildElement("O");    // 8
  G4Element *elNa = nistManager->FindOrBuildElement("Na");  // 11
  G4Element *elMg = nistManager->FindOrBuildElement("Mg");  // 12
  G4Element *elAl = nistManager->FindOrBuildElement("Al");  // 13
  G4Element *elP = nistManager->FindOrBuildElement("P");    // 15
  G4Element *elS = nistManager->FindOrBuildElement("S");    // 16
  G4Element *elCl = nistManager->FindOrBuildElement("Cl");  // 17
  G4Element *elAr = nistManager->FindOrBuildElement("Ar");  // 18
  G4Element *elK = nistManager->FindOrBuildElement("K");    // 19
  G4Element *elCa = nistManager->FindOrBuildElement("Ca");  // 20
  G4Element *elFe = nistManager->FindOrBuildElement("Fe");  // 26
  G4Element *elZn = nistManager->FindOrBuildElement("Zn");  // 30
  G4Element *elW = nistManager->FindOrBuildElement("W");    // 74
  G4Element *elPb = nistManager->FindOrBuildElement("Pb");  // 82

  // Air defined using NIST Manager
  worldMat = nistManager->FindOrBuildMaterial("G4_AIR"); 
  // https://www.3dmensionals.de/formlabs-standard-resin-v4-2574?category=121&attrib=50-2560#attr=10677,2862,10678,2850,2852,2856,2861,21812,23715
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

  mirrorSurface = new G4OpticalSurface("mirrorSurface");
  mirrorSurface->SetType(dielectric_metal);
  mirrorSurface->SetFinish(ground);
  mirrorSurface->SetModel(unified);

  G4MaterialPropertiesTable *mptMirror = new G4MaterialPropertiesTable();
  G4double energy_reflective[2] = {1.239841929*eV/0.9,1.239841929*eV/0.2};
  G4double reflectivity[2] = {1,1};
  mptMirror->AddProperty("REFLECTIVITY", energy_reflective, reflectivity, 2);
  mirrorSurface->SetMaterialPropertiesTable(mptMirror);

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

  detMaterial = nistManager->FindOrBuildMaterial("G4_PbWO4");
  detMaterial->GetIonisation()->SetBirksConstant(0.008694); // mm/MeV
  G4cout << "Mean excitation energy of PbWO4: " << detMaterial->GetIonisation()->GetMeanExcitationEnergy() << G4endl;
  G4cout << "Birks constant PbWO4: " << detMaterial->GetIonisation()->GetBirksConstant() << G4endl;

  heteroAir = new G4Material("heteroAir", 0.001225*g/cm3,3);
  G4double ttt = 75.47+23.20+1.28;
  heteroAir->AddElement(elN, 75.47/ttt);
  heteroAir->AddElement(elO, 23.20/ttt);
  heteroAir->AddElement(elAr, 1.28/ttt);

  heteroWater = new G4Material("heteroWater", 1.00*g/cm3,2);
  heteroWater->AddElement(elH, 2);
  heteroWater->AddElement(elO, 1);
  heteroWater->GetIonisation()->SetMeanExcitationEnergy(79.7*eV);

  aluminiumFoil = new G4Material("aluminiumFoil", 2.71*g/cm3,1);
  aluminiumFoil->AddElement(elAl,1);

  penMaterial = new G4Material("penMaterial", 1.36*g/cm3, 3);
  penMaterial->AddElement(elC, 14);
  penMaterial->AddElement(elH, 10);
  penMaterial->AddElement(elO, 4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
        std::cerr << "Error opening file" << std::endl;
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
  
    // Iterate through each row in the original matrix
    for (const auto& row : matrix) {
        // Iterate through each float in the row
        for (const auto& value : row) {
            // Convert float to int based on the conditions
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
  G4double worldXY = 50*cm;
  G4double worldZ = 250*cm;

  G4double detSizeZ = 5.0*mm; //det size in x
  G4double detSizeX = 50.0*mm; //det size in y
  G4double detSizeY = 5.0*mm; //det size in z

  G4String holder;
  G4double translation;
  G4double offset;
  G4ThreeVector physicalPosition;

  G4bool use_pen = false;
  if(fNbOfDetectors>15) use_pen = true;

  if(!use_pen){
    holder="../constructs/halterung.stl";
    translation = 7*cm-15.75*mm;
    offset = 1.825*cm-15.75*mm;
  }
  else {
    holder="../constructs/halterung_extended_origin.stl";
    translation = 10.175*cm-7.5*mm-0.0999*mm;
    offset = 1*cm-7.5*mm-0.0999*mm;
  }

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World
  G4UserLimits* userLimits = new G4UserLimits();
  // userLimits->SetMaxAllowedStep(0.1*mm);
  
  solidworld = new G4Box("solidworld",              // its name
    worldXY / 2, worldXY / 2, worldZ / 2);               // its size

  logicalworld = new G4LogicalVolume(solidworld,    // its solid
    worldMat,                                                 // its material
    "logicalWorld");                                     // its name
  // logicalworld->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
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
  
  logicaldetectorPbWO4 = new G4LogicalVolume(soliddetector, detMaterial, "logicaldetectorPbWO4");
  logicaldetectorPEN = new G4LogicalVolume(soliddetector, penMaterial, "logicaldetectorPEN");
  
  //G4LogicalSkinSurface *skin = new G4LogicalSkinSurface("skin", logicaldetectorPbWO4, mirrorSurface);
  logicaldetectorPbWO4->SetUserLimits(userLimits);
  logicaldetectorPEN->SetUserLimits(userLimits);

  G4double wrappingThickness = 20*um;
  G4Box* wrappingSubtraction = new G4Box("wrappingSubtraction", detSizeX/2+wrappingThickness, detSizeY/2+wrappingThickness, detSizeZ/2+wrappingThickness);
  solidAluFoil = new G4SubtractionSolid("solidAluFoil", wrappingSubtraction, soliddetector);
  
  logicalAluFoil = new G4LogicalVolume(solidAluFoil, aluminiumFoil, "logicalAluFoil");

  for(G4int i=0; i<fNbOfDetectorsPbWO4; i++){
    if(i>11 && use_pen) physicalPosition = G4ThreeVector(0.,0.,detSizeZ*(i)+(i)*3*mm+offset+2*mm);
    else physicalPosition = G4ThreeVector(0.,0.,detSizeZ*(i)+offset+(i)*3*mm); // 

    physdetectorPbWO4 = new G4PVPlacement(nullptr,  // no rotation
      physicalPosition,      
      logicaldetectorPbWO4,            // its logical volume
      "physdetectorPbWO4",           // its name
      logicalworld,               // its mother volume
      false,                    // no boolean operations
      i,                        // copy number
      fCheckOverlaps);          // checking overlaps
  }

  G4int idx_pos = fNbOfDetectorsPbWO4;

  for(G4int i=0; i<fNbOfDetectorsPEN; i++){
    if(i+idx_pos==12 && use_pen) physicalPosition = physicalPosition+G4ThreeVector(0.,0., detSizeZ+3*mm+2*mm);
    else physicalPosition = physicalPosition+G4ThreeVector(0.,0., detSizeZ+3*mm);
    
    physdetectorPEN = new G4PVPlacement(nullptr,  // no rotation
      physicalPosition,      
      logicaldetectorPEN,            // its logical volume
      "physdetectorPEN",           // its name
      logicalworld,               // its mother volume
      false,                    // no boolean operations
      i+idx_pos,                        // copy number
      fCheckOverlaps);          // checking overlaps
  }
  
  idx_pos = fNbOfDetectorsPbWO4+fNbOfDetectorsPEN;
  
  for(G4int i=0; i<idx_pos; i++){
    if(i>11 && use_pen) physicalPosition = G4ThreeVector(0.,0.,detSizeZ*(i)+(i)*3*mm+offset+2*mm);
    else physicalPosition = G4ThreeVector(0.,0.,detSizeZ*(i)+(i)*3*mm+offset);

    physAluFoil = new G4PVPlacement(nullptr,  // no rotation
      physicalPosition,      
      logicalAluFoil,            // its logical volume
      "physAluFoil",           // its name
      logicalworld,               // its mother volume
      false,                    // no boolean operations
      i+idx_pos,                        // copy number
      fCheckOverlaps);          // checking overlaps
  }

  idx_pos += idx_pos;

  G4VisAttributes* visDetectorPbWO4 = new G4VisAttributes(G4Colour(1.0, 0.84, 0.0, 1));
  visDetectorPbWO4->SetVisibility(true);
  visDetectorPbWO4->SetForceSolid(true); // Ensure the detector is solid

  G4VisAttributes* visDetectorPEN = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 1));
  visDetectorPEN->SetVisibility(true);
  visDetectorPEN->SetForceSolid(true); // Ensure the detector is solid

  G4VisAttributes* visAluFoil = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAluFoil->SetVisibility(true);
  visAluFoil->SetForceWireframe(true); // Wireframe for the aluminum foil to see the detector inside
  
  // Apply visualization attributes to logical volumes
  logicaldetectorPbWO4->SetVisAttributes(visDetectorPbWO4);
  logicaldetectorPEN->SetVisAttributes(visDetectorPEN);
  logicalAluFoil->SetVisAttributes(visAluFoil);
  
  G4STL stl;
  stl.SetVerbosity(1);
  solidHolder = stl.Read(holder); // halterungsspacer 6 mm; spacer 2 mm; seitenanfang 15.25 mm & 6.75mm
  logicalHolder = new G4LogicalVolume(solidHolder, resinMaterial, "logicalHolder");
  physHolder = new G4PVPlacement(nullptr,  // no rotation
       G4ThreeVector(0.,0., translation),              // at (x,y,z)
       logicalHolder,            // its logical volume
       "physHolder",           // its name
       logicalworld,               // its mother volume
       false,                    // no boolean operations
       idx_pos,                        // copy number
       fCheckOverlaps);          // checking overlaps
  idx_pos+=1;
  if(use_pen){
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateX(180.*deg);
    physHolder = new G4PVPlacement(rotation, 
         G4ThreeVector(0.,0., translation+1*cm-7.5*mm),//+103.25*mm),              // at (x,y,z)
         logicalHolder,            // its logical volume
         "physHolder",           // its name
         logicalworld,               // its mother volume
         false,                    // no boolean operations
         idx_pos,                        // copy number
         fCheckOverlaps);          // checking overlaps  
  }
  G4VisAttributes* visHolder = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.6));
  visHolder->SetVisibility(true);
  visHolder->SetForceSolid(true); // Ensure the detector is solid
  logicalHolder->SetVisAttributes(visHolder);
  // logicalHolder->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
  
  if(ftarget=="heterogeneous_lung"){
    solidLung = stl.Read("../constructs/small_lung.stl"); // "very_small_lung_cut.stl, "small_lung.stl", "Lung_phantom_hires_clean_scaled"
    logicalLung = new G4LogicalVolume(solidLung, lungTissue, "logicalLung");

    // Create a rotation matrix
    G4RotationMatrix* rotate_matrix = new G4RotationMatrix();
    rotate_matrix->rotateY(90 * degree); // Rotate by 90 degrees around y-axis

    physLung = new G4PVPlacement(rotate_matrix,  // no rotation
         G4ThreeVector(16.1632*mm,-74.955*mm,-20*cm),              // at (x,y,z) (-51.6565
         logicalLung,            // its logical volume
         "physLung",           // its name
         logicalworld,               // its mother volume
         false,                    // no boolean operations
         idx_pos,                        // copy number
         fCheckOverlaps);          // checking overlaps
    idx_pos += 1;
    G4VisAttributes* visLung = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 1)); // Red for the detector
    visLung->SetVisibility(true);
    visLung->SetForceSolid(true); // Ensure the detector is solid
    logicalLung->SetVisAttributes(visLung);
  }
  else if (ftarget == "heterogeneous"){
    solidContainer = new G4Box("solidContainer", voxelXY*nbofvoxelsX/2, voxelXY*nbofvoxelsY/2, voxelZ*NbOfSlices/2);

    logicalContainer = new G4LogicalVolume(solidContainer, worldMat, "logicalContainer",0,0,0);
    fMateIDs = new size_t[nbofvoxelsX*nbofvoxelsY*NbOfSlices];

    solidVoxel = new G4Box("solidVoxel", voxelXY/2, voxelXY/2, voxelZ/2);
    logicalVoxel = new G4LogicalVolume(solidVoxel, lungTissue, "logicalVoxel", 0, 0, 0);

    // G4UserLimits* userLimits = new G4UserLimits();
    // userLimits->SetMaxAllowedStep(1*mm);
    // logicalVoxel->SetUserLimits(userLimits);

    physContainer = new G4PVPlacement(nullptr,  // no rotation
      G4ThreeVector(0.,0.,-25*cm),              // at (x,y,z) 
      logicalContainer,            // its logical volume
      "physContainer",           // its nameLungInhale
      logicalworld,               // its mother volume
      false,                    // no boolean operations
      idx_pos,                        // copy number
      fCheckOverlaps);          // checking overlaps

    idx_pos += 1;

    for(int i = 0; i < NbOfSlices; i++){
      std::vector<int> indices = readMatrixFromFile("../constructs/dicom_conv/DICOM_Waterphantom" + std::to_string(i+1) + ".txt");  
      for(int j=0; j<nbofvoxelsX*nbofvoxelsY; j++){                                  
        fMateIDs[i*nbofvoxelsX*nbofvoxelsY+j] = indices.at(j);
      }      
    }
    fillDCMcontainer();
  }
  else if(ftarget == "homogeneous"){
    G4double targetDis = 15*cm;
    G4double targetX = 25*cm;
    G4double targetY = 25*cm;
    G4double targetZ = 5*cm;

    solidHomo =  new G4Box("solidHomo", targetX/2, targetY/2, targetZ/2);
    logicalHomo = new G4LogicalVolume(solidHomo, homoMaterial, "logicalHomo");
    physHomo = new G4PVPlacement(nullptr,  // no rotation
      G4ThreeVector(0.,0.,-targetDis),              // at (x,y,z)
      logicalHomo,            // its logical volume
      "physHomo",           // its name
      logicalworld,               // its mother volume
      false,                    // no boolean operations
      idx_pos,                        // copy number
      fCheckOverlaps);          // checking overlaps

      idx_pos += 1;

      G4VisAttributes* visHomo = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.5)); // Red for the detector
      visHomo->SetVisibility(true);
      visHomo->SetForceSolid(true); // Ensure the detector is solid
      logicalHomo->SetVisAttributes(visHomo);
  }

  if(true){
    if(fbeam_particle == "proton"){
      solidNozzleFWHM = new G4EllipticalTube("solidNozzleFWHM", beamFWHM_p_nozzle_x/2, beamFWHM_p_nozzle_y/2, 0.5*mm);
      solidFilmFWHM = new G4EllipticalTube("solidFilmFWHM", beamFWHM_p_film_x/2, beamFWHM_p_film_y/2, 0.5*mm);
    }
    else if(fbeam_particle == "c12"){
      solidNozzleFWHM = new G4EllipticalTube("solidNozzleFWHM", beamFWHM_c12_nozzle_x/2, beamFWHM_c12_nozzle_y/2, 0.5*mm);
      solidFilmFWHM = new G4EllipticalTube("solidFilmFWHM", beamFWHM_c12_film_x/2, beamFWHM_c12_film_y/2, 0.5*mm);
    }
    logicalNozzleFWHM = new G4LogicalVolume(solidNozzleFWHM, worldMat, "logicalNozzleFWHM");
    logicalFilmFWHM = new G4LogicalVolume(solidFilmFWHM, worldMat, "logicalFilmFWHM");

    logicalNozzleFWHM->SetVisAttributes(G4Color::Red());
    logicalFilmFWHM->SetVisAttributes(G4Color::Red());
    logicalNozzleFWHM->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));
    logicalFilmFWHM->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

    physNozzleFWHM = new G4PVPlacement(nullptr,  // no rotation
       G4ThreeVector(x_off,0., -5.975*cm-dis_film_nozzle),// at (x,y,z) //-5.875*cm end
       logicalNozzleFWHM,            // its logical volume
       "physNozzleFWHM",           // its name
       logicalworld,               // its mother volume
       false,                    // no boolean operations
       idx_pos,                        // copy number
       fCheckOverlaps);          // checking overlaps

    idx_pos += 1;

    physFilmFWHM = new G4PVPlacement(nullptr,  // no rotation
       G4ThreeVector(x_off,0., -5.975*cm-dis_det_film),              // at (x,y,z) //-5.875*cm end
       logicalFilmFWHM,            // its logical volume
       "physFilmFWHM",           // its name
       logicalworld,               // its mother volume
       false,                    // no boolean operations
       idx_pos,                        // copy number
       fCheckOverlaps);          // checking overlaps
    idx_pos += 1;
  }
  //Always return the physical world
  logicalworld->SetVisAttributes(new G4VisAttributes(G4VisAttributes::GetInvisible()));

  return physworld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  std::string trackerDetectorSDname = "/TrackerDetectorSD";
  auto aTrackerSD = new TrackerSD(trackerDetectorSDname, "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name
  // of "Chamber_LV".
  SetSensitiveDetector("logicaldetectorPbWO4", aTrackerSD, true);
  SetSensitiveDetector("logicaldetectorPEN", aTrackerSD, true);

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

}