#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include <vector>
#include <map>


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
    void AddEdep(G4int layerID, G4double edep);
    const std::vector<G4double>& GetEdep() const { return fEdep; }
    G4double GetEtot() const { return fEtot; }

    G4int fLayers = 0;
    std::map<G4int, G4double> fWetAccum;
  private:
    G4int hitsCount = 0;
    TrackerHitsCollection* fHitsCollection = nullptr;
    std::vector<G4double> fEdep;
    G4double fEtot;
    G4EmCalculator* fEmCalc = nullptr;  // lazy init
    
};

#endif