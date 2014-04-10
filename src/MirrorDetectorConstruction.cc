// $Id$
//

#include "MirrorDetectorConstruction.hh"
#include "EpamPhotonDetSD.hh"
#include "EpamMaterials.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4SDManager.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"

EpamPhotonDetSD* MirrorDetectorConstruction::mppcSD = NULL;

MirrorDetectorConstruction::MirrorDetectorConstruction()
 :  solidWorld(0), logicWorld(0),
    physiWorld(0), solidMirror(0),
    logicMirror(0), physiMirror(0)
{
  fWorldLength = 1. *m;
  MirrorLength = 320. *mm;
  MirrorThickness = 25. *mm;

  wlsfiberR = 0.5 *mm;
  wlsfiberZ = 320. *mm;
  clad1R = wlsfiberR * 1.01;
  clad1Z = wlsfiberZ;
  clad2R = wlsfiberR * 1.02;
  clad2Z = wlsfiberZ;

  holeRadius = wlsfiberR * 1.03;
  holeLength = wlsfiberZ;
  coatingThickness = wlsfiberR * 0.01;
  coatingRadius = wlsfiberR * 0.01;

  numOfCladLayers = 2;
  surfaceRoughness = 1.0;

  mppcPolish = 1.;
  mppcReflectivity = 0.;
  mppcShape = "Circle";
  mppcHalfL = wlsfiberR;
  mppcDist  = 0.00*mm;
  mppcTheta = 0.0*deg;
  mppcZ     = 0.05*mm;
  clrfiberZ  = mppcZ + 10.*nm ;
  clrfiberHalfL = mppcHalfL;

  G4double  worldSizeX = clad2R   + mppcDist + mppcHalfL + 1.*cm;
  G4double  worldSizeY = clad2R   + mppcDist + mppcHalfL + 1.*cm;
  G4double  worldSizeZ = wlsfiberZ + mppcDist + mppcHalfL + 1.*cm;

  coupleRX = worldSizeX;
  coupleRY = worldSizeY;
  coupleZ  = (worldSizeZ - wlsfiberZ) / 2;
  coupleOrigin = wlsfiberZ + coupleZ; 
}

MirrorDetectorConstruction::~MirrorDetectorConstruction()
{
}

G4VPhysicalVolume* MirrorDetectorConstruction::Construct()
{

  //------------------------------------------------------ materials
  //G4NistManager * nistManager = G4NistManager::Instance();
  //nistManager->SetVerbose(1);
  materials = EpamMaterials::GetInstance();

  //G4Material* Air = nistManager->FindOrBuildMaterial("G4_AIR");
  //G4Material* Vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
  //G4Material* Glass = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)
  solidWorld = new G4Box( "world", fWorldLength/2, fWorldLength/2, fWorldLength/2);
  logicWorld = new G4LogicalVolume( solidWorld, FindMaterial("G4_AIR"), "World");
  physiWorld = new G4PVPlacement( 0, G4ThreeVector(), logicWorld, "world", 0, false, 0);

  //----------------scintillator-----------
  solidMirror = new G4Box( "Mirror", MirrorLength/2, MirrorLength/2, MirrorThickness/2);
  logicMirror = new G4LogicalVolume( solidMirror, FindMaterial("Polystyrene"), "Mirror");
  physiMirror = new G4PVPlacement( 0, G4ThreeVector(), logicMirror, "mirror", logicWorld, false, 0);

  //--------set surface of scintillator------
  
  G4OpticalSurface* MirrorSurface = new G4OpticalSurface("MirrorSurface",  // Surface Name
		                                 glisur,                   // SetModel
						 				 ground,                   // SetFinish
						                 dielectric_dielectric,    // SetType
						                 1.);                      // SetPolish
  G4MaterialPropertiesTable* MirrorSurfaceProperty = new G4MaterialPropertiesTable();
  G4double p_Mirror[2] = { 2.00*eV, 3.47*eV};
  G4double refl_Mirror[2] = { 1., 1.};
  G4double effi_Mirror[2] = { 0, 0};
  MirrorSurfaceProperty->AddProperty("REFLECTIVITY", p_Mirror, refl_Mirror, 2);
  MirrorSurfaceProperty->AddProperty("EFFICIENCY", p_Mirror, effi_Mirror, 2);
  MirrorSurface->SetMaterialPropertiesTable( MirrorSurfaceProperty);
  new G4LogicalSkinSurface("MirrorSurface", logicMirror, MirrorSurface);

  //construct fibers
  G4double firstFiberPositionX = -148.9 *mm; //-154.5 * mm;// change to -148.9*mm??
  G4double fiberGap = 19.82 *mm; //20 *mm;// change to 19.82*mm??
  G4double firstFiberCenterPosition = 12.0 *mm; //11.5 *mm;// change to 12*mm??
  G4double fiberPitch = 1.07*mm;
  
  G4ThreeVector pTranslate;
  G4LogicalVolume* pCurrentLogical;
  G4LogicalVolume* pMotherLogical;
  G4String pName;
  G4int pMany, pCopyNumber;
  G4RotationMatrix *pRotateMatrix;

  G4LogicalVolume* logicFiber = MirrorDetectorConstruction::ConstructFiber();
  G4PVPlacement* physiFiber;
  for (G4int j =0; j<2; j++)
   { 
     if (j==0)
      {
	    pRotateMatrix = new G4RotationMatrix(0,90*deg,0);
        for (G4int i=0; i<16; i++)
        {
	       for(G4int k=0; k!=9; ++k){// 10 to 9
	       //G4RotationMatrix rm;
	       //rm.rotateX(90*deg);
	       physiFiber = new G4PVPlacement(
					    pRotateMatrix,
					    //pTransform = G4Transform3D(rm, G4ThreeVector((firstFiberPositionX + i*fiberGap) , 0, firstFiberCenterPosition- k*fiberPitch)),
					    pTranslate = G4ThreeVector((firstFiberPositionX + i*fiberGap + k*fiberPitch) , 0, firstFiberCenterPosition ),
					    pCurrentLogical = logicFiber,
					    pName = "fiber",
					    pMotherLogical = logicMirror,
					    pMany = false,
					    pCopyNumber = (i*9 + k) );
	       } 
	     }
      }
      else
      { 
	    pRotateMatrix = new G4RotationMatrix(90*deg,90*deg,0);
        for (G4int i=0; i<16; i++)
	    { 
	       for(G4int k=0; k!=9; ++k)
           { 
	         physiFiber = new G4PVPlacement( pRotateMatrix,
			              pTranslate = G4ThreeVector(0, (firstFiberPositionX + i*fiberGap + k*fiberPitch), -firstFiberCenterPosition ),
			              pCurrentLogical = logicFiber,     //should be checked
			              pName = "fiber",
		                  pMotherLogical = logicMirror,
		                  pMany = false,
		                  pCopyNumber = (j*144 + i*9 + k) );
	      }
        }
      } 
   }

  //----------set visibility-----
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  return physiWorld;


}

G4LogicalVolume* MirrorDetectorConstruction::ConstructFiber()
{
  //construct hole
  G4VSolid* solidHole = new G4Tubs("Hole",
                               0.0*cm,
                               holeRadius,
                               holeLength/2,
                               0.*deg,
                               360.*deg);
  logicHole = new G4LogicalVolume(solidHole,
                               FindMaterial("G4_AIR"),
                               "Hole");
/*
  physiHole = new G4PVPlacement(0,
                               G4ThreeVector(),
                               logicHole,
                               "Hole",
                               logicScintillator,
                               false,
                               0);
*/
/*
  if (!(logicHole) || !(physiHole) ) {
     std::ostringstream o;
     o << "The Fiber Hole has not been constructed";
     G4Exception("WLSDetectorConstruction::ConstructFiber","",
                  FatalException,o.str().c_str());
  }
*/
  // Pointers to the most recently constructed volume
  G4LogicalVolume* logicPlacement = new G4LogicalVolume(solidHole,
                               FindMaterial("G4_AIR"),
                               "Holein");
  G4VPhysicalVolume* physiPlacement = new G4PVPlacement(0,
                               G4ThreeVector(),
                               logicPlacement,
                               "Hole",
                               logicHole,
                               false,
                               0);

  //--------------------------------------------------
  // Fiber Construction
  //-------------------------------------------------- 

  // Boundary Surface Properties
  G4OpticalSurface* OpSurface = NULL;
 
  if (surfaceRoughness < 1.)
     OpSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                      glisur,                  // SetModel
                                      ground,                  // SetFinish
                                      dielectric_dielectric,   // SetType
                                      surfaceRoughness);       // SetPolish

  G4LogicalVolume   *logicClad1, *logicClad2;
  G4VPhysicalVolume *physiClad1, *physiClad2;

  // Determine the number of cladding layers to be built
  switch ( numOfCladLayers ) {
 
    case 2:

     //--------------------------------------------------
     // Cladding 2
     //--------------------------------------------------

     G4VSolid* solidClad2;
 

     solidClad2 = new G4Tubs("Clad2",0.,clad2R,clad2Z/2,0.0*rad,twopi*rad);

     logicClad2  = new G4LogicalVolume(solidClad2,
                                       FindMaterial("FPethylene"),
                                       "Clad2");

     physiClad2 = new G4PVPlacement(0,
                                    G4ThreeVector(0.0,0.0,0.0),
                                    logicClad2,
                                    "Clad2",
                                    logicPlacement,
                                    false,
                                    0);

     // Place the rough surface only if needed
     if (OpSurface) {
       new G4LogicalBorderSurface("surfaceClad2Out",
                                  physiClad2,
                                  physiPlacement,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceClad2In",
                                  physiPlacement,
                                  physiClad2,
                                  OpSurface);
     }

     logicPlacement = logicClad2;
     physiPlacement = physiClad2;

    case 1:

     //--------------------------------------------------
     // Cladding 1
     //--------------------------------------------------

     G4VSolid* solidClad1;


     solidClad1 = new G4Tubs("Clad1",0.,clad1R,clad1Z/2,0.0*rad,twopi*rad);

     logicClad1 = new G4LogicalVolume(solidClad1,
                                      FindMaterial("Pethylene"),
                                      "Clad1");

     physiClad1 = new G4PVPlacement(0,
                                    G4ThreeVector(0.0,0.0,0.0),
                                    logicClad1,
                                    "Clad1",
                                    logicPlacement,
                                    false,
                                    0);

     // Place the rough surface only if needed
     if (OpSurface) {
       new G4LogicalBorderSurface("surfaceClad1Out",
                                  physiClad1,
                                  physiPlacement,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceClad1In",
                                  physiPlacement,
                                  physiClad1,
                                  OpSurface);
     }

     logicPlacement = logicClad1;
     physiPlacement = physiClad1;

    default:

     //--------------------------------------------------
     // WLS Fiber
     //--------------------------------------------------

     G4VSolid* solidWLSfiber;

     solidWLSfiber = new G4Tubs("WLSFiber",0.,wlsfiberR,wlsfiberZ/2,0.0*rad,twopi*rad);

     G4LogicalVolume*   logicWLSfiber = new G4LogicalVolume(solidWLSfiber,
                                                         FindMaterial("PMMA"),
                                                         "WLSFiber");

     logicWLSfiber->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

     G4VPhysicalVolume* physiWLSfiber = new G4PVPlacement(0,
                                       G4ThreeVector(0.0,0.0,0.0),
                                       logicWLSfiber,
                                       "WLSFiber",
                                       logicPlacement,
                                       false,
                                       0);

     // Place the rough surface only if needed
     if (OpSurface) {
       new G4LogicalBorderSurface("surfaceWLSOut",
                                  physiWLSfiber,
                                  physiPlacement,
                                  OpSurface);
       new G4LogicalBorderSurface("surfaceWLSIn",
                                  physiPlacement,
                                  physiWLSfiber,
                                  OpSurface);
     }
  }

  //--------------------------------------------------
  // Mirror for reflection at one of the end
  //--------------------------------------------------

  // Place the mirror only if the user wants the mirror
  /*
  if (mirrorToggle) {  

     G4VSolid* solidMirror = new G4Box("Mirror",
                                       mirrorRmax,
                                       mirrorRmax,
                                       mirrorZ);
 
     G4LogicalVolume* logicMirror = new G4LogicalVolume(solidMirror,
                                                        FindMaterial("G4_Al"),
                                                        "Mirror");

     G4OpticalSurface* MirrorSurface = new G4OpticalSurface("MirrorSurface",
                                                             glisur,
                                                             ground,
                                                             dielectric_metal,
                                                             mirrorPolish);

     G4MaterialPropertiesTable* MirrorSurfaceProperty =
                                              new G4MaterialPropertiesTable();

     G4double p_mirror[2] = {2.00*eV, 3.47*eV};
     G4double refl_mirror[2] = {mirrorReflectivity,mirrorReflectivity};
     G4double effi_mirror[2] = {0, 0};

     MirrorSurfaceProperty->AddProperty("REFLECTIVITY",p_mirror,refl_mirror,2);
     MirrorSurfaceProperty->AddProperty("EFFICIENCY",p_mirror,effi_mirror,2);

     MirrorSurface -> SetMaterialPropertiesTable(MirrorSurfaceProperty);

     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,mirrorOrigin),
                       logicMirror,
                       "Mirror",
                       logicWorld,
                       false,
                       0);

     new G4LogicalSkinSurface("MirrorSurface",logicMirror,MirrorSurface);
  }
  */
  //--------------------------------------------------
  // Coupling at the read-out end
  //--------------------------------------------------  

  // Clear Fiber (Coupling Layer)
  
  G4VSolid* solidCouple = new G4Box("Couple",coupleRX,coupleRY,coupleZ);

  G4LogicalVolume*   logicCouple = new G4LogicalVolume(solidCouple,
                                                       FindMaterial("G4_AIR"),
                                                       "Couple");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0,0.0,coupleOrigin),
                    logicCouple,
                    "Couple",
                    logicHole,
                    false,
                    0);
  
  //--------------------------------------------------
  // A logical layer in front of PhotonDet
  //--------------------------------------------------  

  // Purpose: Preventing direct dielectric to metal contact  

  // Check for valid placement of PhotonDet
  
  if (mppcTheta > std::atan(mppcDist / mppcHalfL)) {

     mppcTheta = 0;
     mppcOriginX  = std::sin(mppcTheta) * (mppcDist + clrfiberZ);
     mppcOriginZ  = -coupleZ + std::cos(mppcTheta) * (mppcDist + clrfiberZ);
     G4cerr << "Invalid alignment.  Alignment Reset to 0" << G4endl;     
  }
 
  // Clear Fiber (Coupling Layer)
  G4VSolid* solidClrfiber;
 
  //if ( mppcShape == "Square" )
  //  solidClrfiber = 
  //              new G4Box("ClearFiber",clrfiberHalfL,clrfiberHalfL,clrfiberZ);
  //else
  solidClrfiber =
       new G4Tubs("ClearFiber",0.,clrfiberHalfL,clrfiberZ,0.0*rad,twopi*rad);

  G4LogicalVolume*   logicClrfiber =
                                          new G4LogicalVolume(solidClrfiber,
                                                       FindMaterial("G4_AIR"),
                                                       "ClearFiber");

  new G4PVPlacement(new G4RotationMatrix(CLHEP::HepRotationY(-mppcTheta)),
                    G4ThreeVector(mppcOriginX,0.0,mppcOriginZ),
                    logicClrfiber,
                    "ClearFiber",
                    logicCouple,
                    false,
                    0); 
  
  //--------------------------------------------------
  // PhotonDet (Sensitive Detector)
  //--------------------------------------------------  
  
  // Physical Construction
  G4VSolid* solidPhotonDet;

  //if ( mppcShape == "Square" )
    //solidPhotonDet = new G4Box("PhotonDet",mppcHalfL,mppcHalfL,mppcZ);
  //else
  solidPhotonDet =
                 new G4Tubs("PhotonDet",0.,mppcHalfL,mppcZ,0.0*rad,twopi*rad);

  G4LogicalVolume*   logicPhotonDet =
                                    new G4LogicalVolume(solidPhotonDet,
                                                        FindMaterial("G4_Al"),
                                                        "PhotonDet");

  new G4PVPlacement(0,
                    G4ThreeVector(0.0,0.0,0.0),
                    logicPhotonDet,
                    "PhotonDet",
                    logicClrfiber,
                    false,
                    0);

  // PhotonDet Surface Properties
  G4OpticalSurface* PhotonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       mppcPolish);

  G4MaterialPropertiesTable* PhotonDetSurfaceProperty =
                                               new G4MaterialPropertiesTable();

  G4double p_mppc[2] = {2.00*eV, 3.47*eV};
  G4double refl_mppc[2] = {mppcReflectivity,mppcReflectivity};
  G4double effi_mppc[2] = {1, 1};
 
  PhotonDetSurfaceProperty -> AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
  PhotonDetSurfaceProperty -> AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);

  PhotonDetSurface -> SetMaterialPropertiesTable(PhotonDetSurfaceProperty);
 
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,PhotonDetSurface); 

  if (!mppcSD) {
     G4String mppcSDName = "WLS/PhotonDet";
     mppcSD = new EpamPhotonDetSD(mppcSDName);

     G4SDManager* SDman = G4SDManager::GetSDMpointer();
     SDman->AddNewDetector(mppcSD);
  }

  // Setting the detector to be sensitive
  logicPhotonDet->SetSensitiveDetector(mppcSD);
  
  return logicHole;

}

G4double MirrorDetectorConstruction::GetWLSFiberEnd()
{
  return wlsfiberZ;
}

G4bool MirrorDetectorConstruction::IsPerfectFiber()
{ 
  return     surfaceRoughness == 1. ;
}

G4Material* MirrorDetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}

