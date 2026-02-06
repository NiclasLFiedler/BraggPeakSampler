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
  analysisManager->SetNtupleFileName("../data/data.root");
  analysisManager->CreateNtuple("vtree", "track, NDet and EDep");
  analysisManager->CreateNtupleIColumn(0, "event");	
  analysisManager->CreateNtupleIColumn(0, "NDet");
  analysisManager->CreateNtupleDColumn(0, "EDep");	 
  analysisManager->FinishNtuple();  
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->OpenFile("../data/data.root");
}


void RunAction::EndOfRunAction(const G4Run*)
{
    auto sm = G4ScoringManager::GetScoringManager();
    auto mesh = sm->GetMesh(0);
    auto scoreMap = mesh->GetScoreMap();
    auto hitsMap = (scoreMap)["Edep"];

  hitsMap->PrintAllHits();
// for(const auto& it : *(hitsMap->GetMap()))
// {
    // G4int bin = it.first;
    // it->PrintAllHits();
// }

    sm->DumpAllQuantitiesToFile("waterMesh", "BraggPeak.txt");
}
