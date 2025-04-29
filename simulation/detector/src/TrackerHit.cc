#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator = nullptr;

G4bool TrackerHit::operator==(const TrackerHit& right) const
{
  return ( this == &right ) ? true : false;
}

void TrackerHit::Draw()
{

}

void TrackerHit::Print()
{
    G4cout
     << "  trackID: " << fTrackID << G4endl;
}