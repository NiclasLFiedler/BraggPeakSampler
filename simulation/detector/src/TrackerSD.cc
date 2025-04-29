#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4OpticalPhoton.hh"
#include "G4Proton.hh"
#include "G4VProcess.hh"


TrackerSD::TrackerSD(const G4String& name, const G4String& hitsCollectionName, G4double layers)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
  fLayers = layers;
  entryPosMap = std::vector<std::vector<double>>(fLayers, std::vector<double>(3, 0.0));
  hitMap = std::vector<G4int>(fLayers, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection
    = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4Track* track = aStep->GetTrack();
  G4int NDet = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(2);
  if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
    //hitMap.at(NDet)++;
    //aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return true; 
  }
  auto newHit = new TrackerHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  G4double eDep = aStep->GetTotalEnergyDeposit();
  if (track->GetDefinition() == G4Proton::ProtonDefinition()) {
    // G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
    if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
      G4ThreeVector entryPosition = aStep->GetPreStepPoint()->GetPosition();
      // G4cout << "Layer: " << NDet << " EDep: " << eDep << " Process: " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
      if(entryPosMap.at(NDet).at(0) != 0){
        // G4cout << "Already filled" << G4endl;
      }
      else{
        entryPosMap.at(NDet) = {entryPosition.x(), entryPosition.y(), entryPosition.z()};
      }
    }
    // newHit->SetNDet(NDet);
    // newHit->SetEdep(energy);
    // fHitsCollection->insert(newHit);
  }
  // aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  newHit->SetNDet(NDet);
  newHit->SetEdep(eDep);
  fHitsCollection->insert(newHit);
  return true;
}

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>2 ) {
     std::size_t nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

void TrackerSD::ClearHits()
{
  hitMap = std::vector<G4int>(fLayers, 0);
  entryPosMap = std::vector<std::vector<double>>(fLayers, std::vector<double>(3, 0.0));
  return;
}