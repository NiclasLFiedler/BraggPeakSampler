#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* ) override;

    G4ParticleGun* GetParticleGun() {return fParticleGun;}

  private:
    G4ParticleGun* fParticleGun = nullptr;
    G4ThreeVector beamPos = G4ThreeVector(0.,0.,0.);
    G4ThreeVector beamDis = G4ThreeVector(0.,0.,0.);
    G4double beamFWHM = 0;
    G4double sigma = 0;
    G4double x0 = 0;
    G4double y0 = 0;
    G4double beamZ = 0;
    G4double beamY = 0;
    G4double beamX = 0;
};

#endif
