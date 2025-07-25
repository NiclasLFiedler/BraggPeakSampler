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

void EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();	
  if (!analysisManager->OpenFile("../data/data.root")) {
    G4cout << "Error: Unable to open file!" << G4endl;
  }
  G4int eventID = event->GetEventID();
  if ( eventID % 1000 == 0){
    G4cout << ">>> Event: " << eventID  << G4endl;
  }
  G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
  
  SiPMSD* siPMSD = (SiPMSD*)SDmanager->FindSensitiveDetector("SiPMSD");
  TrackerSD* trackerSD = (TrackerSD*)SDmanager->FindSensitiveDetector("TrackerDetectorSD");

  if(siPMSD->fHit == 1 || siPMSD->target != 2) {
    std::vector<G4double> fEdep = trackerSD->GetEdep();
    
    for(int i = 0; i < fEdep.size(); i++) {
      if(fEdep.at(i) > 0) {

        analysisManager->FillNtupleIColumn(0, eventID);
        analysisManager->FillNtupleIColumn(1, i);
        analysisManager->FillNtupleDColumn(2, fEdep.at(i)/MeV);
        analysisManager->AddNtupleRow(0);
      }
    }
    analysisManager->Write();  
  }
  trackerSD->ClearHits(); 
  siPMSD->ClearHits();
  return;
}


