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

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                              G4TouchableHistory*)
{
    G4double trackid = aStep->GetTrack()->GetTrackID();

    if (trackid != 1)
        return true;

    G4StepPoint* pre  = aStep->GetPreStepPoint();
    G4StepPoint* post = aStep->GetPostStepPoint();
    
    // Entering a layer
    if (pre->GetStepStatus() == fGeomBoundary)
    {
        layerID = pre->GetTouchable()->GetCopyNumber();
        prepos = pre->GetPosition().z();
        const G4ThreeVector& dirIn =
            pre->GetMomentumDirection();

        thetaXIn = std::atan2(dirIn.x(), -dirIn.z());
    }

    // Leaving a layer
    if (post->GetStepStatus() == fGeomBoundary)
    {
        // Check that this is the layer we stored at entrance
        if (layerID != pre->GetTouchable()->GetCopyNumber())
        { 
            layerID = -1;
            thetaXIn = 0.;
            return true;
        }

        auto newHit = new TrackerHit();

        const G4ThreeVector& dirOut =
            post->GetMomentumDirection();

        thetaXOut =
            std::atan2(dirOut.x(), -dirOut.z());


        G4double deltaThetaX =
            thetaXOut - thetaXIn;
    
        G4double layerThickness = 1.; // mm
        
        // G4double depth = std::abs(std::round(post->GetPosition().z() / layerThickness)* layerThickness / 10.);
        G4double depth = (layerID+1)* layerThickness / 10.;
        postpos = post->GetPosition().z();
        G4double energy =
            pre->GetTotalEnergy();

        G4double eKin =
            pre->GetKineticEnergy();
        if (prepos == postpos){
          layerID = -1;
          thetaXIn = 0.;
          thetaXOut = 0;
          return true;
        }

        newHit->SetTrackID(trackid);
        newHit->SetEkin(eKin);
        newHit->SetEtot(energy);
        newHit->SetDepth(depth);

        // Cumulative angle relative to original z direction
        newHit->SetCumVariance(thetaXOut);

        // Angular change accumulated over this entire layer
        newHit->SetScatteringAngle(deltaThetaX);

        newHit->SetLayerID(layerID);

        fHitsCollection->insert(newHit);

        // Reset
        prepos = 0;
        postpos = 0;
        layerID = -1;
        thetaXIn = 0.;
        thetaXOut = 0;
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

// hadd textfile.root textfile*.root

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

