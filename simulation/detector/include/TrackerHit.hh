#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class TrackerHit : public G4VHit
{
  public:
    TrackerHit() = default;
    TrackerHit(const TrackerHit&) = default;
    ~TrackerHit() override = default;

    // operators
    TrackerHit& operator=(const TrackerHit&) = default;
    G4bool operator==(const TrackerHit&) const;

    // Set methods
    void SetTrackID  (G4int track)      { fTrackID = track; };
    void SetNDet  (G4int ndet)      { fNDet = ndet; };
    void SetNPhotons  (G4int nPhotons)      { fNPhotons = nPhotons; };
    void SetEdep     (G4double de)      { fEdep = de; };
    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetNDet() const   { return fNDet; };
    G4int GetNPhotons() const   { return fNPhotons; };
    G4double GetEdep() const     { return fEdep; };
  
    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    void Draw() override;
    void Print() override;

  private:
    G4int         fTrackID = -1;
    G4int         fNDet = -1;
    G4int         fNPhotons = 0;
    G4double      fEdep = 0.;
  
};

using TrackerHitsCollection = G4THitsCollection<TrackerHit>;

extern G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t)
{
  if(!TrackerHitAllocator)
      TrackerHitAllocator = new G4Allocator<TrackerHit>;
  return (void *) TrackerHitAllocator->MallocSingle();
}

inline void TrackerHit::operator delete(void *hit)
{
  TrackerHitAllocator->FreeSingle((TrackerHit*) hit);
}

#endif
