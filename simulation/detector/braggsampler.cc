#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "QBBC.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4PhysListFactory.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  //
  auto runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction());
  //auto physicsList = new QBBC();
  G4PhysListFactory factory;
  auto physicsList = factory.GetReferencePhysList("QGSP_BIC_EMY");
  

  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());


  G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
  stepLimitPhys->SetApplyToAll(true); // activates step limit for ALL particles
  physicsList->RegisterPhysics(stepLimitPhys);
  
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  
  opticalPhysics->SetVerboseLevel(1);
  G4OpticalParameters* opticalParams = G4OpticalParameters::Instance();
  
  //opticalParams->SetScintillationYieldFactor(1.0);
  //opticalParams->SetScintByParticleType(true);
	opticalParams->SetScintTrackSecondariesFirst(true);	

  //physicsList->RegisterPhysics(opticalPhysics);

  runManager->SetUserInitialization(physicsList);

  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization());

  // Initialize visualization with the default graphics system
  auto visManager = new G4VisExecutive();

  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
    if(ui){
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }

  // Job termination
  delete visManager;
  delete runManager;
  return 0;
}