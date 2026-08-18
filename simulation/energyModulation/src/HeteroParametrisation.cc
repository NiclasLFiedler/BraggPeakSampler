#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" // For G4UniformRand()
#include <unordered_map>

class HeteroParameterisation : public G4VPVParameterisation {
    public:
        HeteroParameterisation(G4int nx, G4int ny, G4int nz, G4double dx, G4double dy, G4double dz,
                               G4Material* mat1, G4Material* mat2,
                               G4VisAttributes* vis1, G4VisAttributes* vis2)
            : fNx(nx), fNy(ny), fNz(nz), fDx(dx), fDy(dy), fDz(dz),
              fMat1(mat1), fMat2(mat2), fVis1(vis1), fVis2(vis2) {
                G4double rhoLung = 1.05;
                G4double rhoLunginf = 0.26;
                G4double rhoH2O = 1;
                G4double rhoair = 0.0012;
                probability = (rhoLunginf-rhoair)/(rhoLung-rhoair);
                probability = 0.2251;
                probability = 0.2467;
              }
    
        void ComputeTransformation(G4int copyNo, G4VPhysicalVolume* physVol) const override {
            int ix = copyNo / (fNy * fNz);
            int iy = (copyNo / fNz) % fNy;
            int iz = copyNo % fNz;
    
            G4ThreeVector position(
                -(0.5+ix)*fDx+(fNx * fDx)/2,
                -(0.5+iy)*fDy+(fNy * fDy)/2,
                -(0.5+iz)*fDz+(fNz * fDz)/2
            );
            physVol->SetTranslation(position);
        }
    
        G4Material* ComputeMaterial(G4int copyNo, G4VPhysicalVolume*, const G4VTouchable*) {
            if (fMaterials.count(copyNo) == 0) {
                G4Material* mat = (G4UniformRand() < probability) ? fMat1 : fMat2;
                fMaterials[copyNo] = mat;
            }
            return fMaterials[copyNo];
            }

    private:
        G4double probability;
        G4int fNx, fNy, fNz;
        G4double fDx, fDy, fDz;
        G4Material* fMat1;
        G4Material* fMat2;
        G4VisAttributes* fVis1;
        G4VisAttributes* fVis2;
        mutable std::unordered_map<G4int, G4Material*> fMaterials;
    };
    