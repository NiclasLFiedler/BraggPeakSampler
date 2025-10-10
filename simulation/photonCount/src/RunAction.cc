#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

RunAction::RunAction()
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
  //analysisManager->SetNtupleMerging(true);
  analysisManager->SetNtupleFileName("../data/data.root");
  analysisManager->CreateNtuple("RawData", "track, NDet, EDep and PhotonCount");
  analysisManager->CreateNtupleIColumn(0, "event");	
  analysisManager->CreateNtupleIColumn(0, "NDet");
  analysisManager->CreateNtupleDColumn(0, "EDep");	 
  analysisManager->CreateNtupleIColumn(0, "PhotonScint");	 
  analysisManager->CreateNtupleIColumn(0, "PhotonSiPM");	 
  analysisManager->FinishNtuple();  
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->OpenFile("../data/data.root");
}

void RunAction::EndOfRunAction(const G4Run* )
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
	analysisManager->SetFileName("../data/data.root");
  analysisManager->CloseFile();
}