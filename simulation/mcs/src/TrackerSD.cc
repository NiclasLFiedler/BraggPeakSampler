#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <algorithm>
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

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
    layerID = -1;
    thetaXIn = 0.0;
    thetaXOut = 0.0;
    prepos = 0.0;
    postpos = 0.0;
    fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection( hcID, fHitsCollection );
}

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                              G4TouchableHistory*)
{
    const G4int trackid = aStep->GetTrack()->GetTrackID();

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

        /////////
        // Use the PRE-step touchable: it describes the layer being left.
        const auto touchable = pre->GetTouchableHandle();

        const auto* box = dynamic_cast<const G4Box*>(touchable->GetVolume()->GetLogicalVolume()->GetSolid());

        if (!box){
            G4Exception( "TrackerSD::ProcessHits", "UnexpectedSolid", FatalException, "Expected a G4Box layer.");
            return true;
        }

        const auto globalToLocal = touchable->GetHistory()->GetTopTransform();
        const G4ThreeVector localExit = globalToLocal.TransformPoint(post->GetPosition());
        const G4double halfZ = box->GetZHalfLength();

        const G4double tolerance = 10.0 * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

        // Beam travels along -z, so the downstream face is local z = -halfZ.
        const G4bool downstreamExit = std::abs(localExit.z() + halfZ) <= tolerance;

        if (!downstreamExit){
            layerID = -1; thetaXIn = 0.0; thetaXOut = 0.0; return true;
        }

        // Nominal downstream-face position, identical for every hit in this layer.
        const G4ThreeVector downstreamWorld = globalToLocal.Inverse().TransformPoint(G4ThreeVector(0, 0, -halfZ));
        
        // Phantom entrance is at global z = 0.
        const G4double depth = -downstreamWorld.z() / cm;

        /////////

        const G4ThreeVector& dirOut = post->GetMomentumDirection();
        thetaXOut = std::atan2(dirOut.x(), -dirOut.z());


        G4double deltaThetaX = thetaXOut - thetaXIn;
        
        postpos = post->GetPosition().z();
        G4double deltaX = post->GetPosition().x();
        const G4double energy = post->GetTotalEnergy();
        const G4double eKin = post->GetKineticEnergy();

        auto newHit = new TrackerHit();
        newHit->SetTrackID(trackid);
        newHit->SetEkin(eKin);
        newHit->SetEtot(energy);
        newHit->SetDepth(depth);

        newHit->SetCumAngle(thetaXOut);

        newHit->SetScatteringAngle(deltaThetaX);
        newHit->SetLayerID(layerID);
        newHit->SetDeltaX(deltaX/cm);
        
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

