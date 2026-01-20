#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prime)
 : fPrimaryGeneratorAction(prime)
{
  fPrimaryGeneratorDirectory = new G4UIdirectory("/set/");
  fPrimaryGeneratorDirectory->SetGuidance("Set beam energy");

  fSetBeamEnergy = new G4UIcmdWithADoubleAndUnit("/set/beamenergy", this);
  fSetBeamEnergy->SetGuidance("Set setBeamEnergy");
  fSetBeamEnergy->SetParameterName("setBeamEnergy", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimaryGeneratorDirectory;
  delete fSetBeamEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fSetBeamEnergy ) {
    fPrimaryGeneratorAction->SetBeamEnergy(fSetBeamEnergy->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
