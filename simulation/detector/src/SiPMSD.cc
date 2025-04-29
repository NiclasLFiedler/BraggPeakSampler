#include "SiPMSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4OpticalPhoton.hh"
#include "G4Proton.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"
#include "G4RandomTools.hh"
#include <vector>

SiPMSD::SiPMSD(const G4String& name, const G4String& hitsCollectionName, G4double layers)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
  fLayers = layers;
  hitMap = std::vector<G4int>(fLayers, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiPMSD::Initialize(G4HCofThisEvent* hce)
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

G4bool SiPMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{	
  G4Track* track = aStep->GetTrack();
  if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    G4MaterialPropertiesTable* mpt = aStep->GetPreStepPoint()->GetMaterial()->GetMaterialPropertiesTable();
    if (!mpt){
      G4cout << "No MPT found!" << G4endl;
      return false;
    }
    G4MaterialPropertyVector* efficiency = mpt->GetProperty("EFFICIENCY");
    if (!efficiency){
      G4cout << "No MPT efficiency found!" << G4endl;
      return false;
    }
    if (G4UniformRand() < efficiency->Value(track->GetTotalEnergy())) {
      G4int NDet = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
      hitMap.at(NDet)++;
    }
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return true;
  }

  return true;
}

void SiPMSD::ClearHits()
{
    hitMap = std::vector<G4int>(fLayers, 0);
    return;
}

void SiPMSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>2 ) {
     std::size_t nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}