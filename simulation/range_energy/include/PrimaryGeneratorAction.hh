#ifndef B2PrimaryGeneratorAction_h
#define B2PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UImessenger.hh"

class G4ParticleGun;
class G4Event;

namespace B2
{

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* ) override;

    void SetBeamEnergy(G4double);

    G4ParticleGun* GetParticleGun() {return fParticleGun;}


  private:
    G4ParticleGun* fParticleGun = nullptr; // G4 particle gun
    G4ThreeVector beamPos = G4ThreeVector(0.,0.,0.);
    G4ThreeVector beamDis = G4ThreeVector(0.,0.,0.);
    G4double beamFWHM = 0;
    G4double sigma = 0;
    G4double x0 = 0;
    G4double y0 = 0;
    G4double beamZ = 0;
    G4double beamY = 0;
    G4double beamX = 0;

    G4double fBeamEnergy;
    G4UIcmdWithADoubleAndUnit* fSetBeamEnergyCmd;

    PrimaryGeneratorMessenger* fPrimaryGeneratorMessenger = nullptr; // messenger
};

}

#endif
