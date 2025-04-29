//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2/B2a/include/TrackerHit.hh
/// \brief Definition of the B2::TrackerHit class

#ifndef B2TrackerHit_h
#define B2TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

namespace B2
{

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

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
    void SetEdep     (G4double de)      { fEdep += de; };
    void SetEkin     (G4double eKin)      { feKin += eKin; };
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; };
    void SetdEdX      (G4double dEdX){ fdEdX = dEdX; };
    void Setbeta      (G4double beta){ fbeta = beta; };
    void SetEtot      (G4double eTot){ feTot = eTot; };
    void SetStepLength      (G4double StepLength){ fStepLength = StepLength; };
    void SetTrackID   (G4double TrackID){ fTrackID = TrackID; };
    // Get methods
    G4double GetEdep() const     { return fEdep; };
    G4double GetEkin() const     { return feKin; };
    G4double GetdEdX() const     { return fdEdX; };
    G4double Getbeta() const     { return fbeta; };
    G4double GetEtot() const     { return feTot; };
    G4double GetStepLength() const     { return fStepLength; };
    G4double GetTrackID() const     { return fTrackID; };
    G4ThreeVector GetPos() const { return fPos; };
    
    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    void Draw() override;
    void Print() override;

  private:
    G4double      fEdep = 0.;
    G4double      feKin = 0.;
    G4double      fdEdX = 0.;
    G4double      fbeta = 0.;
    G4double      feTot = 0.;
    G4double      fStepLength = 0.;
    G4ThreeVector fPos;
    G4double      fTrackID = 0.;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using TrackerHitsCollection = G4THitsCollection<TrackerHit>;

extern G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TrackerHit::operator new(size_t)
{
  if(!TrackerHitAllocator)
      TrackerHitAllocator = new G4Allocator<TrackerHit>;
  return (void *) TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerHit::operator delete(void *hit)
{
  TrackerHitAllocator->FreeSingle((TrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

#endif
