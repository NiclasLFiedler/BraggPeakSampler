#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
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

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory*)
{
  auto newHit = new TrackerHit();
	
  G4double eDep = aStep->GetTotalEnergyDeposit();
  G4ThreeVector PPos = aStep->GetPostStepPoint()->GetPosition();
  G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();
  G4double eKin = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ThreeVector momentumDirection = aStep->GetTrack()->GetMomentumDirection();
  
  G4double depth = aStep->GetPostStepPoint()->GetPosition().z();
  G4ThreeVector momentum = aStep->GetPostStepPoint()->GetMomentum();
  G4double angle = std::atan2(momentumDirection.x(), momentumDirection.z());
  
  if (fabs(depth - fNextDepthBoundary) < 0.001) {
      fNextDepthBoundary += 0.1;      
      newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
      newHit->SetEkin(eKin);
      newHit->SetEdep(eDep);
      newHit->SetPos(PPos);
      newHit->SetEtot(energy);
      newHit->SetTheta(angle);
      fHitsCollection->insert( newHit );
  }
      
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

// hadd textfile.root textfile*.root

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

