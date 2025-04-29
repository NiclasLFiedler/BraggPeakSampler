#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class TrackerSD : public G4VSensitiveDetector
{
  public:
    TrackerSD(const G4String& name, const G4String& hitsCollectionName, G4double layers);
    ~TrackerSD() override = default;

    // methods from base class
    void   Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void   EndOfEvent(G4HCofThisEvent* hitCollection) override;
    G4int GetHitsCount() const { return hitsCount; }
    void ClearHits();
    void Reset() { hitsCount = 0; }

    std::vector<std::vector<double>> entryPosMap;
    std::vector<G4int> hitMap;
    G4int fLayers = 0;
  private:
    G4int hitsCount = 0;
    TrackerHitsCollection* fHitsCollection = nullptr;
};

#endif