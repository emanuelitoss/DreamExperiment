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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "../include/PrimaryGeneratorAction.hh"

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(ofstream* outfile)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeSphere(0)
{
  this->setOutput(outfile);
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  G4double envSizeR = 0;

  if (!fEnvelopeSphere)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeSphere = dynamic_cast<G4Sphere*>(envLV->GetSolid());
  }

  if ( fEnvelopeSphere ) {
    envSizeR = fEnvelopeSphere->GetOuterRadius();
  }
  else {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  ParticleKinematicsGenerator();

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::ParticleKinematicsGenerator(){

  	// this function sets kinematic and specifics of the incoming ray

  	// default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="mu-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(3.*GeV);

  // generation of radnomic angles
  G4double theta = acos(pow(G4UniformRand(),1./3));
  G4double phi = G4UniformRand() * 2 * M_PI;
  
  // direction of the beam
  G4ThreeVector* Direction_Beam = new G4ThreeVector(0, 0, fEnvelopeSphere->GetOuterRadius());
  Direction_Beam->rotateY(theta);
  Direction_Beam->rotateZ(phi);
  fParticleGun->SetParticleMomentumDirection(*Direction_Beam);

  // tangent plane position generation
  G4double position_x = (G4UniformRand() - 0.5) * 2 * fEnvelopeSphere->GetOuterRadius();
  G4double position_y = (G4UniformRand() - 0.5) * 2 * fEnvelopeSphere->GetOuterRadius();
  
  // METODO EMA/FEDE
  G4ThreeVector* Position_Beam = new G4ThreeVector(position_x, position_y, -fEnvelopeSphere->GetOuterRadius());
  
  std::cout << "Initial Vector = (" << Position_Beam->getX() << ", " << Position_Beam->getY() << ", " << Position_Beam->getZ() << ")"
  << "\n\tTheta = " << Position_Beam->theta() << "\tPhi = " << Position_Beam->phi() << std::endl;

  Position_Beam->rotateY(theta);

  std::cout << "Vector After \'rotateY()\'= (" << Position_Beam->getX() << ", " << Position_Beam->getY() << ", " << Position_Beam->getZ() << ")"
  << "\n\tTheta = " << Position_Beam->theta() << "\tPhi = " << Position_Beam->phi() << std::endl;

  Position_Beam->rotateZ(phi);

  std::cout << "Vector After \'rotateZ()\'= (" << Position_Beam->getX() << ", " << Position_Beam->getY() << ", " << Position_Beam->getZ() << ")"
  << "\n\tTheta = " << Position_Beam->theta() << "\tPhi = " << Position_Beam->phi() << std::endl;
  

  ostream* reference = (this->getOutput());
	// print informations
	*reference << "ciao" << endl;
  /*"\nx \t\t y \t\t theta \t\t phi \t\t dir x \t\t dir y \t\t dir z \t\t pos x \t\t pos y \t\t pos z\n" <<
	position_x << "\t" << position_y << "\t" << theta << "\t" << phi << "\t" <<
  Direction_Beam->getX() << "\t" << Direction_Beam->getY() << "\t" << Direction_Beam->getZ() << "\t" <<
  Position_Beam->getX() << "\t" << Position_Beam->getY() << "\t" << Position_Beam->getZ() <<  endl;*/

  // Set position of the particle
  fParticleGun->SetParticlePosition(*Position_Beam);
  delete Position_Beam;
  delete Direction_Beam;

}