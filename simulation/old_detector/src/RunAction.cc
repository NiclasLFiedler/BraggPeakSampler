#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();

  analysisManager->SetFileName("../data/data.root");
  analysisManager->CreateNtuple("vtree", "track, NDet and EDep");
  analysisManager->CreateNtupleIColumn("event");	
  analysisManager->CreateNtupleIColumn("track");	 
  analysisManager->CreateNtupleIColumn("NDet");
  analysisManager->CreateNtupleDColumn("EDep");

  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  const G4String fileName = "../data/data";
	analysisManager->OpenFile(fileName);
  G4cout << "Using " << analysisManager->GetType() << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
	analysisManager->SetFileName("../data/data");	
	analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}

