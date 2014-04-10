/*
 *  $Id: mirror.cc, 2014-02-20 14:32:44 liudongzy $
 *  Author(s):
 *    Liu Dong (dliu13@mail.ustc.edu.cn) 20/02/2014
*/ 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UImanager.hh"
#include "G4UItcsh.hh"
//#include "QGSP_BERT.hh"
//replace physical list
//#include "FTFP_BERT.hh" 

#include "MirrorDetectorConstruction.hh"
#include "EpamPhysicsList.hh"
#include "EpamPrimaryGeneratorAction.hh"

#include "EpamRunAction.hh"
#include "EpamEventAction.hh"
#include "EpamTrackingAction.hh"
#include "EpamSteppingAction.hh"
#include "EpamStackingAction.hh"
#include "EpamSteppingVerbose.hh"

#include "G4RunManager.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{   
  // User Verbose output class
  //
//  G4VSteppingVerbose* verbosity = new MirrorSteppingVerbose;
//  G4VSteppingVerbose::SetInstance(verbosity);
  G4String physName = "QGSP_BERT_EMV";

  G4int seed = time(0)+getpid();
  if (argc  > 2) seed = atoi(argv[argc-1]);  //"atoi" change string to integer, the last argument will be the seed.

  // Choose the Random engine

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

  // My Verbose output class

  G4VSteppingVerbose::SetInstance(new EpamSteppingVerbose);

  // Run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  MirrorDetectorConstruction* detector = new MirrorDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new EpamPhysicsList(physName);
  runManager->SetUserInitialization(physics);

  // User Action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = new EpamPrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);

  EpamRunAction* runAction = new EpamRunAction();
  EpamEventAction* eventAction = new EpamEventAction(runAction);

  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction( new EpamTrackingAction() );
  runManager->SetUserAction( new EpamSteppingAction(detector) );
  runManager->SetUserAction( new EpamStackingAction() );

  // Initialize G4 kernel
  //
  runManager->Initialize();
      
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();     
#endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
	    cout<<"start batch mode"<<endl;
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
     } 
  else           // interactive mode : define UI session
    { 
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis2.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
     }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete runManager;
//  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

