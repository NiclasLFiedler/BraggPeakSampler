#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
//#include "tls.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4STL.hh"

#include "G4PhantomParameterisation.hh"
#include "G4PVParameterised.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;
class DicomPhantomZSliceMerged;

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
    std::string detectorType;
    G4int NbOfSlices; //781 max
    G4int nbofvoxelsX;
    G4int nbofvoxelsY;
    G4int fLayers = 0;
    G4int ftarget = 0;
    G4bool secondaryLayerStatus = false;
    G4bool absorberStatus = false;
    G4int fLayersCut = 0;
    G4bool useAbsorber = false;
    G4double detSizeZ = 0;
    G4double absSizeZ = 0;
    G4double detSizeX = 0;
    G4double detSizeY = 0;
    G4double gapSizeZ = 0;
    G4double absorberSize = 0;
    G4int heteroThickness = 0;
    G4int pmod = 0;

    G4double FWHMNozzleX = 4.2*mm; //4.2
    G4double FWHMNozzleY = 4.2*mm; //4.2
    G4double FWHMIsocentreX = 8.1*mm;
    G4double FWHMIsocentreY = 8.1*mm;
    
    G4double d_NozzleIsocentre = 104*cm;
    G4double d_IsocentreDetector = 6*cm;
    G4double x_off = 0;//-12.5*mm;
    G4double y_off = 0;//-12.5*mm;
    
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
    G4Material *aluminumFoil = nullptr;
    G4Material *penMaterial = nullptr;
    G4Material *Teflon = nullptr;
    G4Material *SiPMGlassMat = nullptr;
    G4Material *SiPMSiliconMat = nullptr;
    G4Material* PMMA = nullptr;
    
    G4TessellatedSolid *solidHolder, *solidLung;  
    G4Box *solidDetector, *solidAbsorber, *solidworld, *solidHomo, *solidVoxel, *solidContainer;
    G4SubtractionSolid *solidAluFoil, *solidAluFoilAbs, *solidTeflonFoil, *solidTeflonFoilAbs;
    G4EllipticalTube *solidNozzle, *solidIsocentre;
    G4LogicalVolume *logicalDetector, *logicalAbsorber, *logicalworld, *logicalHolder, *logicalNozzle, *logicalIsocentre, *logicalHomo, *logicalLung, *logicalVoxel, *logicalContainer,*logicalAluFoil,*logicalAluFoilAbs, *logicalTeflonFoil,*logicalTeflonFoilAbs, *logicalSiPM;
    G4VPhysicalVolume *physDetector, *physAbsorber, *physworld, *physHolder, *physNozzle, *physIsocentre, *physHomo, *physLung, *physContainer, *physAluFoil, *physAluFoilAbs, *physTeflonFoil, *physTeflonFoilAbs, *physSiPM, *physSiPMAbs, *physLGAlu;
    G4PVParameterised* phantom_phys; 

    G4OpticalSurface *wrappingSurface, *dielectricSurface;
    G4LogicalBorderSurface *DetectorGuideBoarder, *GuideDetectorBoarder, *GuideSiPMBoarder, *SiPMGuideBoarder;
    DetectorMessenger* fMessenger = nullptr; // messenger

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

#endif
