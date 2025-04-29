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
 : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr), fBeamEnergy(10.0*MeV)
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
    fBeamEnergy = energy;
    // G4cout << "Changed beamenergy to " << energy << " MeV" << G4endl;
    G4cout << "Energy change " << energy << G4endl;

    G4UserLimits* userLimits = new G4UserLimits();
    if(energy < 5){
      userLimits->SetMaxAllowedStep(0.1*mm);
      G4LogicalVolume* logicaldetector = G4LogicalVolumeStore::GetInstance()->GetVolume("logicaldetector");
      logicaldetector->SetUserLimits(userLimits);
    }
    fParticleGun->SetParticleEnergy(fBeamEnergy);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  beamPos = G4ThreeVector(0,0,0);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1.));
  
  // Apply beam position and direction
  fParticleGun->SetParticlePosition(beamPos);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

}

