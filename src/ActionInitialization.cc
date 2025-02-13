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
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "../include/ActionInitialization.hh"
#include "../include/PrimaryGeneratorAction.hh"
#include "../include/RunAction.hh"
#include "../include/EventAction.hh"
#include "../include/SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "G4MTRunManager.hh"

#include <fstream>
#include <iostream>
using namespace std;

ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
 fDetConstruction(detConstruction)
{}

ActionInitialization::~ActionInitialization(){}

void ActionInitialization::BuildForMaster() const {
  RunAction* runAction = new RunAction;
  SetUserAction(runAction);
}

void ActionInitialization::Build() const {

  //Writing first line in an output file
  G4String outputFile = "../analysisDreamSimulation/positions.txt";
  //this->InitializeOutputFile(outputFile);

  SetUserAction(new PrimaryGeneratorAction(outputFile));

  RunAction* runAction = new RunAction;
  SetUserAction(runAction);

  EventAction* eventAction = new EventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new SteppingAction(eventAction, fDetConstruction));

}

void ActionInitialization::InitializeOutputFile(G4String outputFile) const {

  ofstream* outfile = new ofstream();

  outfile->open(outputFile);
  if(!outfile) cout << OBOLDRED << "ERROR: Could not open the file" << ORESET << endl;

  *outfile << "#This file contains points and directions generated (tangent plane method)\n"
    << "#x\ty\ttheta\tphi\tdir_x\tdir_y\tdir_z\tpos_x\tpos_y\tpos_z\n" << endl;

  outfile->close();

  delete outfile;
}