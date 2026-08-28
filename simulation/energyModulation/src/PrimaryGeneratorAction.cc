#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ChargedGeantino.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserLimits.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr), fBeamEnergy(200.0*MeV)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleEnergy(fBeamEnergy);
  
  fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PrimaryGeneratorAction::SetBeamEnergy(G4double energy) {
    return;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double energy = 90*MeV;
  fBeamEnergy = G4RandGauss::shoot(energy, energy*0.005);
  fParticleGun->SetParticleEnergy(fBeamEnergy);

  G4double x0 = G4RandGauss::shoot(0.0, 8*mm);
  G4double y0 = G4RandGauss::shoot(0.0, 8*mm);

  beamPos = G4ThreeVector(x0, y0, -40*cm);
  
  fParticleGun->SetParticlePosition(beamPos);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1.));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

}

