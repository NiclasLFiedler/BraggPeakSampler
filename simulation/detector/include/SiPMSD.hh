#ifndef SiPMSD_h
#define SiPMSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;

class SiPMSD : public G4VSensitiveDetector
{
  public:
    SiPMSD(const G4String& name, const G4String& hitsCollectionName, G4double layers);
    ~SiPMSD() override = default;

    // methods from base class
    void   Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void   EndOfEvent(G4HCofThisEvent* hitCollection) override;

    void ClearHits(); // Reset the map at the end of each event

    std::vector<G4int> hitMap;
    G4int fLayers = 0;
  private:
    TrackerHitsCollection* fHitsCollection = nullptr;
};

#endif
