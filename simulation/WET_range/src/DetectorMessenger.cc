#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"

namespace B2a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
 : fDetectorConstruction(det)
{
  fBeamshapDirectory = new G4UIdirectory("/beamshape/");
  fBeamshapDirectory->SetGuidance("Enable cone beanshape");

  fBeamShape = new G4UIcmdWithABool("/beamshape/setBeamCone", this);
  fBeamShape->SetGuidance("Set beamshape to cone");
  fBeamShape->SetParameterName("setBeamCone", false);

  fBeamParticle = new G4UIcmdWithAString("/beamshape/setBeamParticle", this);
  fBeamParticle->SetGuidance("Set beamParticle");
  fBeamParticle->SetParameterName("setBeamParticle", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fBeamshapDirectory;
  delete fBeamShape;
  delete fBeamParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fBeamParticle ) {
    fDetectorConstruction->SetBeamParticle(newValue);
  }
  if( command == fBeamShape ) {
    fDetectorConstruction->SetBeamShape(fBeamShape->GetNewBoolValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
