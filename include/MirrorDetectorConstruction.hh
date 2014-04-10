//
// $Id$
//

#ifndef MirrorDetectorConstruction_H
#define MirrorDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

//-------------for mirror---------
class G4Box;
class G4Tubs;
class G4EllipticalTube;
//class G4Material;
//-------------for mirror---------
class EpamMaterials;
class G4Material;

class EpamPhotonDetSD;

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class MirrorDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    MirrorDetectorConstruction();
    ~MirrorDetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4LogicalVolume* ConstructFiber();

  private:

    EpamMaterials* materials;

    //---world parameters---
    G4Box*		solidWorld;
    G4LogicalVolume*	logicWorld;
    G4VPhysicalVolume*	physiWorld;
    //---world---

    //---Mirror parameters------
    G4Box*		solidMirror;
    G4LogicalVolume*	logicMirror;
    G4VPhysicalVolume*	physiMirror;

    G4LogicalVolume* logicHole;
    G4VPhysicalVolume* physiHole;

    G4double fWorldLength;
    G4double MirrorLength;
    G4double MirrorThickness;
    //----Mirror------

    G4double           wlsfiberR;
    G4double           wlsfiberZ;

    G4double           clad1R;
    G4double           clad1Z;
    G4double           clad2R;
    G4double           clad2Z;

    G4double           clrfiberHalfL;
    G4double           clrfiberZ;
    //G4double           barLength;
    //G4double           barBase;
    G4double           holeRadius;
    G4double           holeLength;
    G4double           coatingThickness;
    G4double           coatingRadius;

    G4int              numOfCladLayers;
    G4double           surfaceRoughness;

    G4double           coupleRX;
    G4double           coupleRY;
    G4double           coupleZ;
    G4double           coupleOrigin;
 
    G4String           mppcShape;
    G4double           mppcHalfL;
    G4double           mppcZ;
    G4double           mppcDist;
    G4double           mppcTheta;
    G4double           mppcOriginX;
    G4double           mppcOriginZ;
    G4double           mppcPolish;
    G4double           mppcReflectivity;

    static EpamPhotonDetSD* mppcSD;

    G4Material* FindMaterial(G4String);

public:  
    G4double GetWLSFiberEnd();

    G4bool   IsPerfectFiber();

};

#endif

