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

#include "G4EmCalculator.hh"
#include "G4Proton.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
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

    G4NistManager* nist = G4NistManager::Instance();

    // Define material (you can use G4_Al, G4_WATER, G4_Si, etc.)
    G4Material* material = nist->FindOrBuildMaterial("G4_WATER");
  
    // Particle: proton
    G4ParticleDefinition* particle = G4Proton::ProtonDefinition();
  
    // EM calculator
    G4EmCalculator emCal;
  
    std::cout << "Energy [MeV]\tStopping Power [MeV mm^2 / g]" << std::endl;
    std::vector<G4double> energies = {0.001*MeV, 0.0015*MeV, 0.002*MeV, 0.0025*MeV, 0.003*MeV, 0.004*MeV, 0.005*MeV, 0.006*MeV, 0.007*MeV, 0.008*MeV, 0.009*MeV, 0.01*MeV, 0.0125*MeV, 0.015*MeV, 0.0175*MeV, 0.02*MeV, 0.0225*MeV, 0.025*MeV, 0.0275*MeV, 0.03*MeV, 0.035*MeV, 0.04*MeV, 0.045*MeV, 0.05*MeV, 0.055*MeV, 0.06*MeV, 0.065*MeV, 0.07*MeV, 0.075*MeV, 0.08*MeV, 0.085*MeV, 0.09*MeV, 0.095*MeV, 0.1*MeV, 0.125*MeV, 0.15*MeV, 0.175*MeV, 0.2*MeV, 0.225*MeV, 0.25*MeV, 0.275*MeV, 0.3*MeV, 0.35*MeV, 0.4*MeV, 0.45*MeV, 0.5*MeV, 0.55*MeV, 0.6*MeV, 0.65*MeV, 0.7*MeV, 0.75*MeV, 0.8*MeV, 0.85*MeV, 0.9*MeV, 0.95*MeV, 1.0*MeV, 1.25*MeV, 1.5*MeV, 1.75*MeV, 2.0*MeV, 2.25*MeV, 2.5*MeV, 2.75*MeV, 3.0*MeV, 3.5*MeV, 4.0*MeV, 4.5*MeV, 5.0*MeV, 5.5*MeV, 6.0*MeV, 6.5*MeV, 7.0*MeV, 7.5*MeV, 8.0*MeV, 8.5*MeV, 9.0*MeV, 9.5*MeV, 10.0*MeV, 12.5*MeV, 15.0*MeV, 17.5*MeV, 20.0*MeV, 25.0*MeV, 27.5*MeV, 30.0*MeV, 35.0*MeV, 40.0*MeV, 45.0*MeV, 50.0*MeV, 55.0*MeV, 60.0*MeV, 65.0*MeV, 70.0*MeV, 75.0*MeV, 80.0*MeV, 85.0*MeV, 90.0*MeV, 95.0*MeV, 100.0*MeV, 125.0*MeV, 150.0*MeV, 175.0*MeV, 200.0*MeV, 225.0*MeV, 250.0*MeV, 275.0*MeV, 300.0*MeV, 350.0*MeV, 400.0*MeV, 450.0*MeV, 500.0*MeV, 550.0*MeV, 600.0*MeV, 650.0*MeV, 700.0*MeV, 750.0*MeV, 800.0*MeV, 850.0*MeV, 900.0*MeV, 950.0*MeV, 1000.0*MeV, 1500.0*MeV, 2000.0*MeV, 2500.0*MeV, 3000.0*MeV, 4000.0*MeV, 5000.0*MeV, 6000.0*MeV, 7000.0*MeV, 8000.0*MeV, 9000.0*MeV, 10000.0*MeV};

    for (G4double energy : energies) {
        G4double dedx = emCal.GetDEDX(energy, particle, material)*10; // in MeV/mm
        G4double rho = material->GetDensity(); // in g/cm3
        rho = 1;
  
        G4double massStoppingPower = dedx / rho; // in MeV cm^2 / g
  
        std::cout << massStoppingPower << std::endl;
    }

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

