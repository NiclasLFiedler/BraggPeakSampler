#ifndef B2aPrimaryGeneratorMessenger_h
#define B2aPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

namespace B2
{

class PrimaryGeneratorAction;


class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction* );
    ~PrimaryGeneratorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    PrimaryGeneratorAction*  fPrimaryGeneratorAction = nullptr;

    G4UIdirectory*         fPrimaryGeneratorDirectory = nullptr;
    G4UIcmdWithADoubleAndUnit* fSetBeamEnergy = nullptr;


};

}

#endif
