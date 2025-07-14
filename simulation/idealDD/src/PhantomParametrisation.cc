#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" // For G4UniformRand()
#include <unordered_map>

class PhantomParameterisation : public G4VPVParameterisation {
    public:
        PhantomParameterisation(G4int nz, G4double dz,
                               G4Material* mat1, G4VisAttributes* vis1)
            : fNz(nz), fDz(dz), fMat1(mat1), fVis1(vis1) {}
    
        void ComputeTransformation(G4int copyNo, G4VPhysicalVolume* physVol) const override {
            int iz = copyNo;
    
            G4ThreeVector position(0, 0, (-iz + fNz * 0.5 - 0.5) * fDz);

            // G4cout << "CopyNo = " << copyNo
            // << ", iz = " << iz
            // << ", Z-pos (local) = " << position.z()/mm << " mm"
            // << G4endl;
            physVol->SetTranslation(position);
            return;
        }
    
        G4Material* ComputeMaterial(G4int copyNo, G4VPhysicalVolume*, const G4VTouchable*) {
            return fMat1;
            }

    private:
        G4int fNz;
        G4double fDz;
        G4Material* fMat1;
        G4VisAttributes* fVis1;
        mutable std::unordered_map<G4int, G4Material*> fMaterials;
    };
    