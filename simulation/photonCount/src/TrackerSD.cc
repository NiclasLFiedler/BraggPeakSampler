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
  fEdep.resize(fLayers, 0.0);
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
  G4int layerID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);

  if (aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::Definition())
    return true; // skip optical photons

  const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
  if (!secondaries) return true;

  for (auto sec : *secondaries) {
      if (sec->GetDefinition() == G4OpticalPhoton::Definition()) {
          if (sec->GetCreatorProcess() && 
              sec->GetCreatorProcess()->GetProcessName() == "Scintillation") {
              hitMap.at(layerID)++;
          }
      }
  }

  G4double eDep = aStep->GetTotalEnergyDeposit();

  if (eDep > 0 && (track->GetDefinition() == G4Proton::ProtonDefinition())) {
    AddEdep(layerID, eDep);
  }
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
  fEdep = std::vector<G4double>(fLayers, 0);
  return;
}

void TrackerSD::AddEdep(G4int layerID, G4double edep) {
  if (layerID >= 0 && layerID < fEdep.size()){
      fEdep.at(layerID) += edep;
  }
  return;
}