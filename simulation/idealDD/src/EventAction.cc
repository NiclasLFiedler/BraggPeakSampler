#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "SiPMSD.hh"
#include "TrackerSD.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include "HeteroParametrisation.cc"
#include "G4RegularNavigation.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"
#include "G4MultiFunctionalDetector.hh"
#include "Randomize.hh"

void EventAction::BeginOfEventAction(const G4Event*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  return;
}


