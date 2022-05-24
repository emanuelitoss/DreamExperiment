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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "../include/SteppingAction.hh"
#include "../include/EventAction.hh"
#include "../include/DetectorConstruction.hh"
#include "../include/RunData.hh"
#include "../include/OutputColors.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

SteppingAction::SteppingAction(EventAction* eventAction, const DetectorConstruction* detConstruction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  fDetConstruction(detConstruction)
{}

SteppingAction::~SteppingAction(){}

void SteppingAction::UserSteppingAction(const G4Step* step){

  if (!fScoringVolume) { 
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get global volume
  // auto global_volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4VPhysicalVolume* physicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();

  auto runData = static_cast<RunData*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  //std::cout << ORED << physicalVolume->GetName() << ORESET << std::endl;
  //std::cout << OCYAN << fDetConstruction->GetBGOcrystal()->GetName() << ORESET << std::endl;

  if ( physicalVolume == fDetConstruction->GetBGOcrystal() ) {
    runData->Add(kBGO, edepStep);
    fEventAction->AddEdepBGO(edepStep);
  }

  if ( physicalVolume == fDetConstruction->GetPlasticScintillator_1() ) {
    fEventAction->PassedThroughScint1();
    runData->Add(kScint1, edepStep);
    fEventAction->AddEdepScint1(edepStep);
  }

  if ( physicalVolume == fDetConstruction->GetPlasticScintillator_2() ) {
    fEventAction->PassedThroughScint2();
    runData->Add(kScint2, edepStep);
    fEventAction->AddEdepScint2(edepStep);
  }
  
  // useful print to check the right behaviour:
  //fEventAction->PrintStatus();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  fEventAction->AddEdep(edepStep);

}
