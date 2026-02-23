#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"
#include "G4RootAnalysisManager.hh"
#include "G4SDManager.hh"
#include "TrackerHit.hh"

class EventAction : public G4UserEventAction
{
  public:
    EventAction() = default;
    ~EventAction() override = default;

    void  BeginOfEventAction(const G4Event* ) override;
    void  EndOfEventAction(const G4Event* ) override;
};


#endif