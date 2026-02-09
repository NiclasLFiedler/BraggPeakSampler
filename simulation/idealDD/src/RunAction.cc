#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ScoringBox.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"

RunAction::RunAction()
{

}

void RunAction::BeginOfRunAction(const G4Run*)
{

}


void RunAction::EndOfRunAction(const G4Run*)
{
    auto sm = G4ScoringManager::GetScoringManager();


    sm->DumpAllQuantitiesToFile("waterMesh1", "meshEdep1.txt");
    // sm->DumpAllQuantitiesToFile("waterMesh2", "meshEdep2.txt");
    // sm->DumpAllQuantitiesToFile("waterMesh3", "meshEdep3.txt");
}


    // auto mesh1 = sm->GetMesh(0);
    // auto mesh2 = sm->GetMesh(1);
    // auto mesh3 = sm->GetMesh(2);
    // auto scoreMap = mesh1->GetScoreMap();
    // auto hitsMap = (scoreMap)["Edep"];
    // hitsMap->PrintAllHits();
    // for(const auto& it : *(hitsMap->GetMap()))
    // {
        // G4int bin = it.first;
        // it->PrintAllHits();
    // }