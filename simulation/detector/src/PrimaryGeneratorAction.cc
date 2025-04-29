#include "PrimaryGeneratorAction.hh"

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
#include "G4Gamma.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  G4LogicalVolume* logicalworld
    = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalWorld");
  G4Box* solidworld = nullptr;
  if ( logicalworld ) solidworld = dynamic_cast<G4Box*>(logicalworld->GetSolid());
  if ( solidworld ){
    beamZ = solidworld->GetZHalfLength();
    beamY = solidworld->GetYHalfLength();
    beamX = solidworld->GetXHalfLength();
  }
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be placed in the center." << G4endl;
  }

  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    //c12
    G4int Z = 6, A = 12;
    G4double ionCharge   = 6.*eplus;
    G4double excitEnergy = 0.*keV;

    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A, excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    fParticleGun->SetParticleEnergy(A*430*MeV);
    beamFWHM = 6.6*mm;
    sigma = beamFWHM / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    x0 = G4RandGauss::shoot(0.0, sigma);
    y0 = G4RandGauss::shoot(0.0, sigma);
    beamPos = G4ThreeVector(x0, y0, -beamZ+1*um);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  }
  else if(particle->GetParticleName() == "proton"){
    //G4ParticleDefinition* particle = G4Gamma::GammaDefinition();
    fParticleGun->SetParticleDefinition(particle);
    G4double energy = 220*MeV;
    //G4double energy = 0.662*MeV;
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(energy, energy*0.005));
    fParticleGun->SetParticleEnergy(energy);
    G4LogicalVolume* logicalWorld = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalWorld");
    G4LogicalVolume* logicalNozzle = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalNozzle");
    G4LogicalVolume* logicalIsocentre = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalIsocentre");
    
    G4EllipticalTube* solidNozzle = nullptr;
    G4EllipticalTube* solidIsocentre = nullptr;
    G4VPhysicalVolume *physNozzle = nullptr;
    G4VPhysicalVolume *physIsocentre = nullptr;
    solidNozzle = dynamic_cast<G4EllipticalTube*>(logicalNozzle->GetSolid());
    solidIsocentre = dynamic_cast<G4EllipticalTube*>(logicalIsocentre->GetSolid());
    
    physNozzle = logicalWorld->GetDaughter(logicalWorld->GetNoDaughters()-2);
    physIsocentre = logicalWorld->GetDaughter(logicalWorld->GetNoDaughters()-1);
    G4double FWHMNozzleX = solidNozzle->GetDx()*2;
    G4double FWHMNozzleY = solidNozzle->GetDy()*2;
    G4double FWHMIsocentreX = solidIsocentre->GetDx()*2;
    G4double FWHMIsocentreY = solidIsocentre->GetDy()*2;
    
    G4ThreeVector PosNozzle = physNozzle->GetObjectTranslation();
    G4ThreeVector Pos_Isocentre = physIsocentre->GetObjectTranslation();
    G4double sigma_nozzle_x = FWHMNozzleX / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    G4double sigma_nozzle_y = FWHMNozzleY / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    G4double sigma_film_x = FWHMIsocentreX / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    G4double sigma_film_y = FWHMIsocentreY / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    
    // beam path
    G4double x0 = G4RandGauss::shoot(0.0, sigma_nozzle_x);
    G4double y0 = G4RandGauss::shoot(0.0, sigma_nozzle_y);
    G4double x1 = G4RandGauss::shoot(0.0, sigma_film_x);
    G4double y1 = G4RandGauss::shoot(0.0, sigma_film_y);      
    G4ThreeVector beamPosNozzle = G4ThreeVector(x0+PosNozzle.x(), y0+PosNozzle.y(), PosNozzle.z());
    G4ThreeVector beamPosIsocentre = G4ThreeVector(x1+Pos_Isocentre.x(), y1+Pos_Isocentre.y(), Pos_Isocentre.z());
    //beamPosNozzle = G4ThreeVector(0, 0, PosNozzle.z());
    //beamPosIsocentre = G4ThreeVector(0, 0, 0);
    beamPosIsocentre = G4ThreeVector(x0+PosNozzle.x(), y0+PosNozzle.y(), Pos_Isocentre.z());
    fParticleGun->SetParticleMomentumDirection(beamPosIsocentre-beamPosNozzle);
    
    beamPos = beamPosNozzle;
  }
  else if(particle->GetParticleName() == "mu-"){
    fParticleGun->SetParticleEnergy(3*GeV);
    beamX = beamX/2*(1 - 2*G4UniformRand());
    beamZ  = beamZ/2*(1 - 2*G4UniformRand());
    beamPos = G4ThreeVector(beamX, beamY-1*um, beamZ);

    G4double theata = pow(cos(M_PI/2*G4UniformRand()), 2);
    G4double phi = 2*M_PI*(1-2*G4UniformRand());
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theata)*sin(phi),-cos(theata), sin(theata)*cos(phi)));
    // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
  }

  // Apply beam position and direction
  fParticleGun->SetParticlePosition(beamPos);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  return;
}

