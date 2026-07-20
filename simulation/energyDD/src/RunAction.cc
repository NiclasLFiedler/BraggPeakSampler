#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ScoringBox.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"

RunAction::RunAction()
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
  //analysisManager->SetNtupleMerging(true);
  analysisManager->SetNtupleFileName("../data/temp/data.root");
  analysisManager->CreateNtuple("vtree", "track, NDet, EDep and WetAccum");
  analysisManager->CreateNtupleIColumn(0, "event");	
  analysisManager->CreateNtupleIColumn(0, "NDet");
  analysisManager->CreateNtupleDColumn(0, "EDep");
  analysisManager->CreateNtupleDColumn(0, "WetAccum");
  analysisManager->FinishNtuple();  
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->OpenFile("../data/temp/data.root");
}

void RunAction::EndOfRunAction(const G4Run* )
{
  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->SetFileName("../data/temp/data.root");
  analysisManager->CloseFile();
}