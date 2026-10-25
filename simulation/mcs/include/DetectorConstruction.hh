#ifndef B2aDetectorConstruction_h
#define B2aDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4TessellatedSolid.hh"
#include "tls.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4STL.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVParameterised.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

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

    void SetBeamShape(G4bool use_cone) { G4cout << "Use Cone: " << use_cone << G4endl; fuse_cone = use_cone; };
    void SetBeamParticle(G4String beam_particle) { G4cout << "Particle Name: " << beam_particle << G4endl; fbeam_particle = beam_particle; };

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // static data members
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                         // magnetic field messenger
    // data members
    G4int fNbOfDetectors = 0;
    G4Material *detMaterial = nullptr;
    G4Material *worldMat = nullptr;
    G4Material *urethandimethacrylat = nullptr;
    G4Material *methacrylatmonomere = nullptr;
    G4Material *phosphinoxid = nullptr;
    G4Material* Teflon = nullptr;
    G4Material* Alu = nullptr;
    G4Material *resinMaterial = nullptr;
    G4Material *lungTissue = nullptr;
    G4Material *homoMaterial = nullptr;
    G4Material *PTFEmembrane = nullptr;

    G4TessellatedSolid *solidHolder, *solidLung;  
    G4Box *soliddetector, *solidworld, *solidHomo;
    G4EllipticalTube *solidNozzleFWHM, *solidFilmFWHM;
    G4LogicalVolume *logicaldetector, *logicalworld, *logicalHolder, *logicalNozzleFWHM, *logicalFilmFWHM, *logicalHomo, *logicalLung;
    G4VPhysicalVolume *physdetector, *physworld, *physHolder, *physNozzleFWHM, *physFilmFWHM, *physHomo, *physLung;

    G4OpticalSurface *mirrorSurface;

    DetectorMessenger* fMessenger = nullptr; // messenger

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

}

#endif
