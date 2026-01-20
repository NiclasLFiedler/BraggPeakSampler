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
#include "nlohmann/json.hpp"
#include <filesystem>

namespace fs = std::filesystem;
using json = nlohmann::json;

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

  std::ifstream configFile("../../../analysis/config.json");
    if (!configFile) {
        std::cerr << "Error opening config file!" << std::endl;
        return;
    }

    json allConfigs;
    configFile >> allConfigs;
    
    const auto& config                                = allConfigs["detectors"][int(allConfigs["detectorSelect"])];
    ftarget                                           = allConfigs["targetSelect"];

    detectorType                                      = config["detectorType"];
    fBeamEnergy                                       = config["beamEnergy"];
    fBeamEnergy = fBeamEnergy*MeV;
    fLayers                                           = config["nLayers"];
    std::vector<double> crystalSize                   = config["crystalSize"];
    gapSizeZ                                          = config["gapSizeZ"];
    secondaryLayerStatus                              = config["secondaryLayerStatus"];
    fLayersCut                                        = config["nSecondaryLayers"];
    absSizeZ                                          = config["secLayerSizeZ"];
    absorberStatus                                    = config["absorberStatus"];
    absorberSize                                      = config["absorberSize"];
    ThicknessTeflon                                   = config["teflonThickness"];
    ThicknessAlu                                      = config["aluThickness"];
    absorberType                                      = config["absorberType"];
    targetThickness                                   = config["targetThickness"];
    
    pmod                                              = allConfigs["pmod"];
    heteroThickness                                   = allConfigs["heteroThickness"];
    

    detSizeX = crystalSize.at(0)* mm;
    detSizeY = crystalSize.at(1)* mm;
    detSizeZ = crystalSize.at(2)* mm;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PrimaryGeneratorAction::SetBeamEnergy(G4double energy) {
    // G4cout << "Changed beamenergy to " << energy << " MeV" << G4endl;
    G4cout << "Energy change " << fBeamEnergy << G4endl;

    // G4UserLimits* userLimits = new G4UserLimits();
    // if(energy < 5){
    //   userLimits->SetMaxAllowedStep(0.1*mm);
    //   G4LogicalVolume* logicaldetector = G4LogicalVolumeStore::GetInstance()->GetVolume("logicaldetector");
    //   logicaldetector->SetUserLimits(userLimits);
    // }
    fParticleGun->SetParticleEnergy(fBeamEnergy);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  if(fLayers > 29){
    fParticleGun->SetParticleEnergy(250*MeV); 
  }
  else if (fLayers > 27){
    fParticleGun->SetParticleEnergy(240*MeV); 
  }
  else if (fLayers > 25){
    fParticleGun->SetParticleEnergy(230*MeV); 
  }
  else{
    fParticleGun->SetParticleEnergy(fBeamEnergy);
  }
  G4double offSet = 63*mm;
  offSet = 0*mm;
  if(fLayers <= fLayersCut){
    offSet = offSet + ((absSizeZ/2+ThicknessTeflon/2+ThicknessAlu/2)*2+gapSizeZ)*fLayers;
  }
  else{
    offSet = offSet + ((absSizeZ/2+ThicknessTeflon/2+ThicknessAlu/2)*2+gapSizeZ)*fLayersCut + ((detSizeZ/2+ThicknessAlu/2+ThicknessTeflon/2)*2+gapSizeZ)*(fLayers-fLayersCut);
  }

  beamPos = G4ThreeVector(0,0,-offSet-6*mm);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1.));
  
  // Apply beam position and direction
  fParticleGun->SetParticlePosition(beamPos);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

}

