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
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "../include/DetectorConstruction.hh"

#include "G4RunManager.hh" 
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

DetectorConstruction::~DetectorConstruction(){}

G4VPhysicalVolume* DetectorConstruction::Construct(){

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Dimensions of BGO and plastic scintillators
  G4double shape_bgoXZ = 2.2*cm, shape_bgoY = 18*cm;
  G4double shape_plasticX = 5*cm, shape_plasticY = 1*cm, shape_plasticZ = 10*cm;
  G4double distance_scintillators = 0.5*cm;
  G4double distance_BGOscintillators = 10*cm;
  
  // minimal radius such that the experiment is circumscribed in a sphere
  // we willl use minimal radius + 20%
  G4double H = 2 * shape_plasticZ + distance_BGOscintillators + distance_scintillators;
  G4double minimal_radius = (4 * H*H + shape_bgoY*shape_bgoY) / (8*H);
  G4double radius_sphere = minimal_radius * 1.2;
  
  // Envelope material
  G4Material* envelope_material = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4Material* world_material = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Sphere* solidWorld = new G4Sphere("World", 0, 3*radius_sphere, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_material,      //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                   
  //     
  // Envelope
  //  
  G4Sphere* solidEnv = new G4Sphere("World", 0, radius_sphere, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        envelope_material,   //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //     
  // BGO Crystal bgo_material
  //  
  // the position in choosen in order to minimize the envelope sphere
  G4double pos_BGO = minimal_radius - 2 * shape_plasticZ - distance_BGOscintillators - distance_scintillators - shape_bgoXZ / 2;
  G4Material* bgo_material = nist->FindOrBuildMaterial("G4_BGO");
  G4ThreeVector bgo_position = G4ThreeVector(0, 0, pos_BGO);
  
        
  // BGO shape
  G4Box* BGOShape = new G4Box("BGO Box", 0.5*shape_bgoXZ, 0.5*shape_bgoY, 0.5*shape_bgoXZ);
                      
  G4LogicalVolume* BGO_LogicalVolume =                         
    new G4LogicalVolume(BGOShape,            //its solid
                        bgo_material,        //its material
                        "BGO Crystal");      //its name
               
  new G4PVPlacement(0,                       //no rotation
                    bgo_position,            //at position
                    BGO_LogicalVolume,       //its logical volume
                    "BGO Crystal",           //its name
                    logicEnv,                //its mother volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  //     
  // Plastic scintillator
  //
  // the position in choosen in order to minimize the envelope sphere
  G4double pos_scint1 = minimal_radius - distance_scintillators - 1.5 * shape_plasticZ;
  G4double pos_scint2 = minimal_radius - shape_plasticZ / 2;
  G4Material* plastic_material = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        
  // Plastic scintillators shape
  G4Box* PlasticShape = new G4Box("BGO Box", 0.5*shape_plasticX, 0.5*shape_plasticY, 0.5*shape_plasticZ);

  // Plastic scintillators positions
  G4ThreeVector trigger_1_position = G4ThreeVector(0, 0, pos_scint1);
  G4ThreeVector trigger_2_position = G4ThreeVector(0, 0, pos_scint2);

  G4LogicalVolume* Plastic_LogicalVolume =                         
    new G4LogicalVolume(PlasticShape,             //its solid
                        plastic_material,         //its material
                        "Plastic Scintillator");  //its name
               
  new G4PVPlacement(0,                       //no rotation
                    trigger_1_position,      //at position
                    Plastic_LogicalVolume,   //its logical volume
                    "Plastic Scintillator",  //its name
                    logicEnv,                //its mother volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    trigger_2_position,      //at position
                    Plastic_LogicalVolume,   //its logical volume
                    "Plastic Scintillator",  //its name
                    logicEnv,                //its mother volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  fScoringVolume = BGO_LogicalVolume;

  //
  //always return the physical World
  //
  return physWorld;
}
