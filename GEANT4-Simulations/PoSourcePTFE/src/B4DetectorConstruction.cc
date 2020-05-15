//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class 

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


#include "G4Region.hh"



// CADMESH //
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4Tet.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  G4NistManager* nistManager = G4NistManager::Instance();

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 

  // Get materials
  //G4Material* G4_water = nistManager->FindOrBuildMaterial("G4_WATER");
 
  G4Material* G4_air = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* G4_NaI =   nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4Material* G4_W =   nistManager->FindOrBuildMaterial("G4_W");

  G4Material* G4_Steel = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material* G4_Copper = nistManager->FindOrBuildMaterial("G4_Cu");
  G4Material* G4_vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material* G4_Gold = nistManager->FindOrBuildMaterial("G4_Au");
  G4Material* G4_PTFE = nistManager->FindOrBuildMaterial("G4_TEFLON");
  
  density = 19.3*g/cm3;
  G4Material* sourceMaterial = new G4Material("sourceMaterial", density, 2);
  a = 196.96657*g/mole;
  G4Element* elAu = new G4Element ("Gold","Au",z = 79.0,a);
  a = 210.0*g/mole;
  G4Element* elPo = new G4Element("Polonium","Po",z = 84.,a);  
  
  G4double fraction;
  sourceMaterial->AddElement(elPo, fraction=0.000006);
  sourceMaterial->AddElement(elAu, fraction=0.999994);
  

  param_num = 0;

  // cm
  G4double sourceThickness = 0.00005;
  G4double AuThickness = 0.00015;
  G4double PTFEThickness = 1.0;
  G4double worldThickness = 10.0;
  G4double binWidth_cm = 0.0001;

   

  G4Box* world_solid
    = new G4Box("World", (worldThickness/2.0)*cm, (worldThickness/2.0)*cm, (worldThickness/2.0)*cm); // its size
                         
  G4LogicalVolume* world_logical
    = new G4LogicalVolume(
                 world_solid,           // its solid
                 G4_air,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* world_physical
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 world_logical,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  G4Box* SourceLayer_solid
    = new G4Box("SourceLayer", (2.0/2.0)*cm, (1.0/2.0)*cm, (sourceThickness/2.0)*cm); // its size

  G4LogicalVolume* SourceLayer_logical
    = new G4LogicalVolume(
                 SourceLayer_solid,           // its solid
                 sourceMaterial,  // its material
                 "SourceLayer");         // its name
                                   
  G4VPhysicalVolume* SourceLayer_physical
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0.0*m, 0.0*m, 0.0*cm),  
                 SourceLayer_logical,   // its logical volume                         
                 "SourceLayer",          // its name
                 world_logical,   // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps */

  
  G4Box* AuLayer_solid
    = new G4Box("AuLayer", (2.0/2.0)*cm, (1.0/2.0)*cm, (AuThickness/2.0)*cm); // its size

  G4LogicalVolume* AuLayer_logical
    = new G4LogicalVolume(
                 AuLayer_solid,           // its solid
                 G4_Gold,  // its material
                 "AuLayer");         // its name
                                   
  G4VPhysicalVolume* AuLayer_physical
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0.0*m, 0.0*m, (sourceThickness/2.0 + AuThickness/2.0)*cm),  
                 AuLayer_logical,   // its logical volume                         
                 "AuLayer",          // its name
                 world_logical,   // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps */

  // 1 micron thick
  G4Box* PTFE_solid
    = new G4Box("PTFE", (0.3/2.0)*cm, (0.3/2.0)*cm, (binWidth_cm/2.0)*cm); // its size

  G4LogicalVolume* PTFE_logical
    = new G4LogicalVolume(
                 PTFE_solid,           // its solid
                 G4_PTFE,  // its material
                 "PTFE");         // its name
                                   
G4VPhysicalVolume* PTFE_physical_0 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 1*binWidth_cm)*cm), PTFE_logical, "0", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_1 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 2*binWidth_cm)*cm), PTFE_logical, "1", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_2 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 3*binWidth_cm)*cm), PTFE_logical, "2", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_3 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 4*binWidth_cm)*cm), PTFE_logical, "3", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_4 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 5*binWidth_cm)*cm), PTFE_logical, "4", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_5 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 6*binWidth_cm)*cm), PTFE_logical, "5", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_6 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 7*binWidth_cm)*cm), PTFE_logical, "6", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_7 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 8*binWidth_cm)*cm), PTFE_logical, "7", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_8 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 9*binWidth_cm)*cm), PTFE_logical, "8", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_9 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 10*binWidth_cm)*cm), PTFE_logical, "9", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_10 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 11*binWidth_cm)*cm), PTFE_logical, "10", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_11 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 12*binWidth_cm)*cm), PTFE_logical, "11", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_12 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 13*binWidth_cm)*cm), PTFE_logical, "12", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_13 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 14*binWidth_cm)*cm), PTFE_logical, "13", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_14 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 15*binWidth_cm)*cm), PTFE_logical, "14", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_15 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 16*binWidth_cm)*cm), PTFE_logical, "15", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_16 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 17*binWidth_cm)*cm), PTFE_logical, "16", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_17 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 18*binWidth_cm)*cm), PTFE_logical, "17", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_18 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 19*binWidth_cm)*cm), PTFE_logical, "18", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_19 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 20*binWidth_cm)*cm), PTFE_logical, "19", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_20 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 21*binWidth_cm)*cm), PTFE_logical, "20", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_21 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 22*binWidth_cm)*cm), PTFE_logical, "21", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_22 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 23*binWidth_cm)*cm), PTFE_logical, "22", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_23 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 24*binWidth_cm)*cm), PTFE_logical, "23", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_24 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 25*binWidth_cm)*cm), PTFE_logical, "24", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_25 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 26*binWidth_cm)*cm), PTFE_logical, "25", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_26 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 27*binWidth_cm)*cm), PTFE_logical, "26", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_27 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 28*binWidth_cm)*cm), PTFE_logical, "27", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_28 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 29*binWidth_cm)*cm), PTFE_logical, "28", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_29 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 30*binWidth_cm)*cm), PTFE_logical, "29", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_30 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 31*binWidth_cm)*cm), PTFE_logical, "30", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_31 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 32*binWidth_cm)*cm), PTFE_logical, "31", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_32 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 33*binWidth_cm)*cm), PTFE_logical, "32", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_33 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 34*binWidth_cm)*cm), PTFE_logical, "33", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_34 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 35*binWidth_cm)*cm), PTFE_logical, "34", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_35 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 36*binWidth_cm)*cm), PTFE_logical, "35", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_36 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 37*binWidth_cm)*cm), PTFE_logical, "36", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_37 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 38*binWidth_cm)*cm), PTFE_logical, "37", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_38 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 39*binWidth_cm)*cm), PTFE_logical, "38", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_39 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 40*binWidth_cm)*cm), PTFE_logical, "39", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_40 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 41*binWidth_cm)*cm), PTFE_logical, "40", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_41 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 42*binWidth_cm)*cm), PTFE_logical, "41", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_42 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 43*binWidth_cm)*cm), PTFE_logical, "42", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_43 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 44*binWidth_cm)*cm), PTFE_logical, "43", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_44 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 45*binWidth_cm)*cm), PTFE_logical, "44", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_45 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 46*binWidth_cm)*cm), PTFE_logical, "45", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_46 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 47*binWidth_cm)*cm), PTFE_logical, "46", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_47 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 48*binWidth_cm)*cm), PTFE_logical, "47", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_48 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 49*binWidth_cm)*cm), PTFE_logical, "48", world_logical, false, 0, fCheckOverlaps);
G4VPhysicalVolume* PTFE_physical_49 = new G4PVPlacement(0, G4ThreeVector(0.0*m, 0.0*m, (0.5*sourceThickness + AuThickness + 0.5*binWidth_cm + 50*binWidth_cm)*cm), PTFE_logical, "49", world_logical, false, 0, fCheckOverlaps);
 

  G4Colour colorRed(G4Colour::Red());
  G4Colour colorGreen(G4Colour::Green());
  G4Colour colorBlue(G4Colour::Blue());
  G4Colour colorYellow(G4Colour::Yellow());
  G4Colour colorBrown(G4Colour::Brown());


  G4VisAttributes* detVisAtt1= new G4VisAttributes(colorRed);
  SourceLayer_logical->SetVisAttributes(detVisAtt1);

  G4VisAttributes* detVisAtt2= new G4VisAttributes(colorYellow);
  AuLayer_logical->SetVisAttributes(detVisAtt2);

  G4VisAttributes* detVisAtt3= new G4VisAttributes(colorGreen);
  PTFE_logical->SetVisAttributes(detVisAtt3);


  return world_physical;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}


 void B4DetectorConstruction::SetParameters(G4int p) 
    {
        param_num = p;
    }

















//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
