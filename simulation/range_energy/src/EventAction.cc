#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4SystemOfUnits.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();	
 
  G4int eventID = event->GetEventID();
  //G4cout << "EventID " << eventID << G4endl;
  if ( eventID % 10000 == 0){
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
    // for (G4int i = 0; i < n_hit; i++){
	  //   G4double partEdep = (*hitsCollection)[i]->GetEdep();
    //   G4double dEdX = (*hitsCollection)[i]->GetdEdX();
    //   G4double eTot = (*hitsCollection)[i]->GetEtot();
    //   G4double TrackID = (*hitsCollection)[i]->GetTrackID();
    //   G4double eKin = (*hitsCollection)[i]->GetEkin();
    //   G4double StepLength = (*hitsCollection)[i]->GetStepLength();
    //   G4ThreeVector Pos = (*hitsCollection)[i]->GetPos();
      // G4double Theta = (*hitsCollection)[i]->GetTheta();

      G4double partEdep = (*hitsCollection)[n_hit-1]->GetEdep();
      G4double dEdX = (*hitsCollection)[n_hit-1]->GetdEdX();
      G4double eTot = (*hitsCollection)[n_hit-1]->GetEtot();
      G4double TrackID = (*hitsCollection)[n_hit-1]->GetTrackID();
      G4double eKin = (*hitsCollection)[n_hit-1]->GetEkin();
      G4double StepLength = (*hitsCollection)[n_hit-1]->GetStepLength();
      G4ThreeVector Pos = (*hitsCollection)[n_hit-1]->GetPos();
      
      analysisManager->FillNtupleIColumn(0, eventID);	
      analysisManager->FillNtupleDColumn(1, Pos.z());
	    analysisManager->FillNtupleDColumn(2, partEdep/MeV);
      analysisManager->FillNtupleDColumn(3, dEdX);
      analysisManager->FillNtupleDColumn(4, TrackID);
      analysisManager->FillNtupleDColumn(5, eKin);
      analysisManager->FillNtupleDColumn(6, eTot);
      analysisManager->FillNtupleDColumn(7, StepLength);
      // analysisManager->FillNtupleDColumn(8, Theta);
      
 	    analysisManager->AddNtupleRow();
    // }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}


