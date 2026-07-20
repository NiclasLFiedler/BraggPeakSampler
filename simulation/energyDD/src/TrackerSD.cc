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
  fEtot = 0.0;
  if (!fEmCalc) fEmCalc = new G4EmCalculator();
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
  fWetAccum.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4Track* track = aStep->GetTrack();

  G4double eDep = aStep->GetTotalEnergyDeposit();
  G4int layerID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
  if (track->GetTrackID() == 1 && track->GetDefinition() == G4Proton::Proton()){

    G4double stepLength = aStep->GetStepLength();
    G4double kineticE   = aStep->GetPreStepPoint()->GetKineticEnergy();

    G4Material* material = track->GetMaterial();
    G4Material* water    = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    G4double dedx_mat   = fEmCalc->ComputeTotalDEDX(kineticE, G4Proton::Proton(), material);
    G4double dedx_water = fEmCalc->ComputeTotalDEDX(kineticE, G4Proton::Proton(), water);
    G4double wet_step = (dedx_mat / dedx_water) * stepLength;

    fWetAccum[layerID] += wet_step; 
    fEtot = kineticE;
  }

  if (eDep >= 0) {
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
  fEdep = std::vector<G4double>(fLayers, 0);
  fEtot = 0;
  fWetAccum.clear();
  return;
}

void TrackerSD::AddEdep(G4int layerID, G4double edep) {
  if (layerID >= 0 && layerID < fEdep.size()){
      fEdep.at(layerID) += edep;
  }
  return;
}