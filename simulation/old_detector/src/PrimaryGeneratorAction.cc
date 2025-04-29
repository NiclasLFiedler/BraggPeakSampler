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

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
    G4double energy = 220*MeV;
    // G4cout << "E " << G4RandGauss::shoot(energy, energy*0.005) << G4endl;
    fParticleGun->SetParticleEnergy(G4RandGauss::shoot(energy, energy*0.005));
    G4LogicalVolume* logicalWorld = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalWorld");
    if(true){//logicalWorld->GetNoDaughters()>=16)
      G4LogicalVolume* logicalNozzleFWHM = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalNozzleFWHM");
      G4LogicalVolume* logicalFilmFWHM = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalFilmFWHM");
      
      G4EllipticalTube* solidNozzleFWHM = nullptr;
      G4EllipticalTube* solidFilmFWHM = nullptr;

      G4VPhysicalVolume *physNozzleFWHM = nullptr;
      G4VPhysicalVolume *physFilmFWHM = nullptr;

      solidNozzleFWHM = dynamic_cast<G4EllipticalTube*>(logicalNozzleFWHM->GetSolid());
      solidFilmFWHM = dynamic_cast<G4EllipticalTube*>(logicalFilmFWHM->GetSolid());
      
      physNozzleFWHM = logicalWorld->GetDaughter(logicalWorld->GetNoDaughters()-2);
      physFilmFWHM = logicalWorld->GetDaughter(logicalWorld->GetNoDaughters()-1);

      G4double beamFWHM_nozzle_x = solidNozzleFWHM->GetDx()*2;
      G4double beamFWHM_nozzle_y = solidNozzleFWHM->GetDy()*2;
      G4double beamFWHM_film_x = solidFilmFWHM->GetDx()*2;
      G4double beamFWHM_film_y = solidFilmFWHM->GetDy()*2;
      
      G4ThreeVector trans_FWHM_nozzle = physNozzleFWHM->GetObjectTranslation();
      G4ThreeVector trans_FWHM_film = physFilmFWHM->GetObjectTranslation();

      G4double sigma_nozzle_x = beamFWHM_nozzle_x / (2.0 * std::sqrt(2.0 * std::log(2.0)));
      G4double sigma_nozzle_y = beamFWHM_nozzle_y / (2.0 * std::sqrt(2.0 * std::log(2.0)));
      G4double sigma_film_x = beamFWHM_film_x / (2.0 * std::sqrt(2.0 * std::log(2.0)));
      G4double sigma_film_y = beamFWHM_film_y / (2.0 * std::sqrt(2.0 * std::log(2.0)));
      
      G4double x0_nzl = G4RandGauss::shoot(0.0, sigma_nozzle_x);
      G4double y0_nzl = G4RandGauss::shoot(0.0, sigma_nozzle_y);
      G4double x0_flm = G4RandGauss::shoot(0.0, sigma_film_x);
      G4double y0_flm = G4RandGauss::shoot(0.0, sigma_film_y);      
      G4ThreeVector beamPos_nozzle = G4ThreeVector(x0_nzl+trans_FWHM_nozzle.x(), y0_nzl, trans_FWHM_nozzle.z());
      G4ThreeVector beamPos_film = G4ThreeVector(x0_flm+trans_FWHM_film.x(), y0_flm, trans_FWHM_film.z());
      G4ThreeVector beamdir = beamPos_film - beamPos_nozzle;
      //beamdir = G4ThreeVector(0., 0., 1.);
      beamPos = beamPos_nozzle;
      fParticleGun->SetParticleMomentumDirection(beamdir);
    }
    else{
      beamFWHM = 4.2*mm;
      sigma = beamFWHM / (2.0 * std::sqrt(2.0 * std::log(2.0)));
      x0 = G4RandGauss::shoot(0.0, sigma);
      y0 = G4RandGauss::shoot(0.0, sigma);
      beamPos = G4ThreeVector(x0, y0, -beamZ+1*um);
      //beamPos = G4ThreeVector(20*mm,0., -beamZ+1*um);
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1.));
    }
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
}
}

