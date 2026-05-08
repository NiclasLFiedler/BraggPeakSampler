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
#include <map>

void EventAction::BeginOfEventAction(const G4Event*)
{
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();	
  if (!analysisManager->OpenFile("../data/temp/data.root")) {
    G4cout << "Error: Unable to open file!" << G4endl;
  }
  G4int eventID = event->GetEventID();
  if ( eventID % 1000 == 0){
    G4cout << ">>> Event: " << eventID  << G4endl;
  }
  G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
  
  TrackerSD* trackerSD = (TrackerSD*)SDmanager->FindSensitiveDetector("TrackerDetectorSD");

  std::vector<G4double> fEdep = trackerSD->GetEdep();
  std::map<G4int, G4double> fWetAccum = trackerSD->fWetAccum;
  G4double wetAccumValue = 0;
  for(int i = 0; i < fEdep.size(); i++) {
    if(fWetAccum[i] != 0) {
      wetAccumValue = fWetAccum[i];
    }
    // std::cout << "wetAccumValue: " << wetAccumValue << " mm" << G4endl;
    if(fEdep.at(i) >= 0) {
      analysisManager->FillNtupleIColumn(0, eventID);
      analysisManager->FillNtupleIColumn(1, i);
      analysisManager->FillNtupleDColumn(2, fEdep.at(i)/MeV);
      analysisManager->FillNtupleDColumn(3, wetAccumValue/mm);
      analysisManager->AddNtupleRow(0);
    }
  }
  analysisManager->Write();  

  trackerSD->ClearHits(); 

  return;
}

