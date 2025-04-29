#ifndef B2aDetectorConstruction_h
#define B2aDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "tls.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4STL.hh"

#include "G4PhantomParameterisation.hh"
#include "G4PVParameterised.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;
class DicomPhantomZSliceMerged;

namespace B2a
{

class DetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Set methods
    void SetCheckOverlaps(G4bool );
    //G4double 
    G4double voxelXY;
    G4double voxelZ;

    int NbOfSlices; //781 max
    int nbofvoxelsX;
    int nbofvoxelsY;

    G4double beamFWHM_p_nozzle_x = 4.2*mm; //4.2
    G4double beamFWHM_p_nozzle_y = 4.2*mm; //4.2
    G4double beamFWHM_p_film_x = 8.1*mm;
    G4double beamFWHM_p_film_y = 8.1*mm;

    G4double beamFWHM_c12_nozzle_x = 6.6*mm;
    G4double beamFWHM_c12_nozzle_y = 6.6*mm;
    G4double beamFWHM_c12_film_x = 3.6*mm;
    G4double beamFWHM_c12_film_y = 3.6*mm;
    
    G4double dis_film_nozzle = 103*cm;
    G4double dis_det_film = 1*cm;
    G4double x_off = 0*cm;
    
    G4String ftarget = "heterogeneous"; // "homogeneous", "heterogeneous", "heterogeneous_lung", "n/a";
    
    size_t* fMateIDs;
    G4bool fuse_cone = true;
    G4String fbeam_particle = "proton"; //proton or c12
    void SetBeamShape(G4bool use_cone) { G4cout << "Use Cone: " << use_cone << G4endl; fuse_cone = use_cone; };
    void SetBeamParticle(G4String beam_particle) { G4cout << "Particle Name: " << beam_particle << G4endl; fbeam_particle = beam_particle; };

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    void ReadPhantomData();
    void ReadPhantomDataFile(const G4String& fname);
    void ConstructPhantomContainer();
    void fillDCMcontainer();
    std::vector<int> readMatrixFromFile(const std::string &filename);
    // static data members
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                         // magnetic field messenger
    // data members
    G4int fNbOfDetectors = 0;
    G4int fNbOfDetectorsPbWO4 = 0;
    G4int fNbOfDetectorsPEN = 0;
    G4Material *detMaterial = nullptr;
    G4Material *worldMat = nullptr;
    G4Material *urethandimethacrylat = nullptr;
    G4Material *methacrylatmonomere = nullptr;
    G4Material *phosphinoxid = nullptr;
    G4Material *resinMaterial = nullptr;
    G4Material *lungTissue = nullptr;
    G4Material *homoMaterial = nullptr;
    G4Material *heteroAir = nullptr;
    G4Material *heteroWater = nullptr;
    G4Material *aluminiumFoil = nullptr;
    G4Material *penMaterial =nullptr;

    G4TessellatedSolid *solidHolder, *solidLung;  
    G4Box *soliddetector, *solidworld, *solidHomo, *solidVoxel, *solidContainer;
    G4EllipticalTube *solidNozzleFWHM, *solidFilmFWHM;
    G4SubtractionSolid *solidAluFoil;
    G4LogicalVolume *logicaldetectorPbWO4, *logicaldetectorPEN, *logicalworld, *logicalHolder, *logicalNozzleFWHM, *logicalFilmFWHM, *logicalHomo, *logicalLung, *logicalVoxel, *logicalContainer,*logicalAluFoil;
    G4VPhysicalVolume *physdetectorPbWO4, *physdetectorPEN, *physworld, *physHolder, *physNozzleFWHM, *physFilmFWHM, *physHomo, *physLung, *physContainer, *physAluFoil;
    G4PVParameterised* phantom_phys; 

    G4OpticalSurface *mirrorSurface;

    DetectorMessenger* fMessenger = nullptr; // messenger

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

}

#endif
