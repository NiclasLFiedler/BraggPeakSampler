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
  G4int hitsCollectionID = SDmanager->GetCollectionID("TrackerHitsCollection");


  TrackerHitsCollection* hitsCollection(0);
  // -- Collection ID is correct, we get the pointer of the collection:
  if (hitsCollectionID>=0) hitsCollection = (TrackerHitsCollection*)(event->GetHCofThisEvent()->GetHC(hitsCollectionID));
  else G4cout << "Collection `hitsCollection' not found!" << G4endl;
  
  if(hitsCollection){
    G4int n_hit = hitsCollection->entries();
    for (G4int i = 0; i < n_hit; i++){
      G4int TrackID = (*hitsCollection)[i]->GetTrackID();
      G4int NDet = (*hitsCollection)[i]->GetNDet();
	    G4double partEdep = (*hitsCollection)[i]->GetEdep();

      analysisManager->FillNtupleIColumn(0, 0, eventID);
      analysisManager->FillNtupleIColumn(0, 1, TrackID);
	    analysisManager->FillNtupleIColumn(0, 2, NDet);
	    analysisManager->FillNtupleDColumn(0, 3, partEdep/MeV);
 	    analysisManager->AddNtupleRow(0);
    }
  }

  SiPMSD* sipmSD = (SiPMSD*)SDmanager->FindSensitiveDetector("SiPMSD");
  TrackerSD* trackerSD = (TrackerSD*)SDmanager->FindSensitiveDetector("TrackerDetectorSD");
  if (sipmSD) {      
    std::vector<G4int> hitMap = sipmSD->hitMap;
    // std::vector<G4int> hitMap = trackerSD->hitMap;
    std::vector<std::vector<double>> entryPosMap = trackerSD->entryPosMap;
    if (!hitMap.empty()){
      for(int channel = 0; channel<hitMap.size(); channel++) {
        if(hitMap.at(channel) == 0) continue;
        analysisManager->FillNtupleIColumn(1, 0, eventID);
        analysisManager->FillNtupleIColumn(1, 1, channel);
	      analysisManager->FillNtupleIColumn(1, 2, hitMap.at(channel));
        analysisManager->FillNtupleDColumn(1, 3, entryPosMap.at(channel).at(0));
        analysisManager->FillNtupleDColumn(1, 4, entryPosMap.at(channel).at(1));
        analysisManager->FillNtupleDColumn(1, 5, entryPosMap.at(channel).at(2));
        analysisManager->AddNtupleRow(1);
      }
    }
    sipmSD->ClearHits();
    trackerSD->ClearHits();
  }
  else {
    G4cout << "SiPMSensitiveDetector not found!" << G4endl;
  }
  analysisManager->Write();   
}


