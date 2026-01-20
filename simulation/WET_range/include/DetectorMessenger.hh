#ifndef B2aDetectorMessenger_h
#define B2aDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

namespace B2a
{

class DetectorConstruction;

/// Messenger class that defines commands for DetectorConstruction.
///
/// It implements commands:
/// - /B2/det/stepMax value unit

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    DetectorConstruction*  fDetectorConstruction = nullptr;

    G4UIdirectory*         fBeamshapDirectory = nullptr;

    G4UIcmdWithABool* fBeamShape = nullptr;
    G4UIcmdWithAString* fBeamParticle = nullptr;
};

}

#endif
