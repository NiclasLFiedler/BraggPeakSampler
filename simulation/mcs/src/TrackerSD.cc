#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <algorithm>

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
  
  G4double trackid = aStep->GetTrack()->GetTrackID();
  if (trackid == 1){
    G4StepPoint* pre  = aStep->GetPreStepPoint();
    G4StepPoint* post = aStep->GetPostStepPoint();

    if(pre->GetStepStatus() == fGeomBoundary) {
      auto newHit = new TrackerHit();
      G4int layerID = pre->GetTouchable()->GetCopyNumber();
      G4ThreeVector PPos = aStep->GetPostStepPoint()->GetPosition();
      // G4double depth = PPos.z();
      G4double layerThickness = 10;
      G4double depth = std::round(pre->GetPosition().z() / layerThickness) * layerThickness/10;
      G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();
      G4double eKin = aStep->GetPreStepPoint()->GetKineticEnergy();     
      
      G4ThreeVector dirOut = post->GetMomentumDirection();
      G4double thetaOut = std::atan2(dirOut.x(), dirOut.z());

      G4ThreeVector preMom = aStep->GetPreStepPoint()->GetMomentumDirection();
      G4double cosDeflectionAngle = preMom.dot(dirOut);
      G4double thetaXPre  = std::atan2(preMom.x(), preMom.z());
      G4double thetaXPost = std::atan2(dirOut.x(), dirOut.z());

      G4double deltaThetaX = thetaXPost - thetaXPre;
      // cosDeflectionAngle = std::clamp(cosDeflectionAngle, -1.0, 1.0);
      // G4double deflectionAngle = std::acos(cosDeflectionAngle);

      newHit->SetTrackID(trackid);
      newHit->SetEkin(eKin);
      newHit->SetEtot(energy);
      newHit->SetDepth(depth);
      newHit->SetCumVariance(thetaOut*thetaOut);
      newHit->SetScatteringAngle(deltaThetaX);
      newHit->SetLayerID(layerID);

      fHitsCollection->insert(newHit);
    }
    else{
      return true;
    }
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

