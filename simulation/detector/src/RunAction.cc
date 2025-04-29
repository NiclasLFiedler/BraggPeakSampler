#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

RunAction::RunAction()
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
  //analysisManager->SetNtupleMerging(true);
  analysisManager->SetNtupleFileName("../data/data.root");
  analysisManager->CreateNtuple("vtree", "track, NDet and EDep");
  analysisManager->CreateNtupleIColumn(0, "event");	
  analysisManager->CreateNtupleIColumn(0, "track");	 
  analysisManager->CreateNtupleIColumn(0, "NDet");
  analysisManager->CreateNtupleDColumn(0, "EDep");
  analysisManager->FinishNtuple(0);  

  analysisManager->CreateNtuple("ptree", "eventPhotons, NDetPhotons and NPhotons");
  analysisManager->CreateNtupleIColumn(1, "eventPhotons");
  analysisManager->CreateNtupleIColumn(1, "NDetPhotons");
  analysisManager->CreateNtupleIColumn(1, "NPhotons");
  analysisManager->CreateNtupleDColumn(1, "EntryPosX");
  analysisManager->CreateNtupleDColumn(1, "EntryPosY");
  analysisManager->CreateNtupleDColumn(1, "EntryPosZ");
  analysisManager->FinishNtuple(1);
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

