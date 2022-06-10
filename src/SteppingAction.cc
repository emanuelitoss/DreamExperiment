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

#include "G4Cerenkov.hh"
#include "G4EventManager.hh"
#include "G4Scintillation.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <vector>

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

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4VPhysicalVolume* physicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// SCORE ENERGY ////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  auto runData = static_cast<RunData*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // for each PV:
  // modify a boolean value to check the passage of particle through PV and add energy
  if ( physicalVolume == fDetConstruction->GetBGOcrystal() ) {
    fEventAction->PassedThroughBGO();
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

  // Cherenkov and scintillation photons in the two PMTs
  
  if ( physicalVolume == fDetConstruction->GetBGOcrystal()){
    
    // loop over secondaries of this particle
    const std::vector <const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
    for(auto sec : *secondaries)
    {
      G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
      if(creator_process.compare("Cerenkov") == 0)
      {
        G4double cher_photon_energy = sec->GetKineticEnergy();
        runData->Add(kBGO_Cherenkov, cher_photon_energy);
        runData->Add(kNum_Cerenkov, 1);
        fEventAction->AddEdepBGOCerenkov(cher_photon_energy);
      }
      
      else if(creator_process.compare("Scintillation") == 0)
      {
        std::cout << OYELLOW << "check Scintillation" << ORESET << std::endl;
        G4double scint_photon_energy = sec->GetKineticEnergy();
        runData->Add(kBGO_Scintillation, scint_photon_energy);
        runData->Add(kNum_Scint, 1);
        fEventAction->AddEdepBGOScint(scint_photon_energy);
      }
    }
  }
  /*
  if ( physicalVolume == fDetConstruction->GetCerenkovVolume() || physicalVolume == fDetConstruction->GetScintillatorVolume())
  {
    // this particle
    G4Track* primary = step->GetTrack();
    G4String creator_process_thisparticle = primary->GetCreatorProcess()->GetProcessName();

    if(creator_process_thisparticle.compare("Cerenkov") == 0 && physicalVolume == fDetConstruction->GetCerenkovVolume()){
      G4double cher_photon_energy = primary->GetKineticEnergy();
      runData->Add(kBGO_Cherenkov, cher_photon_energy);
      runData->Add(kNum_Cerenkov, 1);
      fEventAction->AddEdepBGOCerenkov(cher_photon_energy);
      primary->SetKineticEnergy(0.);
    }

    if(creator_process_thisparticle.compare("Scintillation") == 0 && physicalVolume == fDetConstruction->GetScintillatorVolume()){
      std::cout << OYELLOW << "check Scintillation" << ORESET << std::endl;
      G4double scint_photon_energy = primary->GetKineticEnergy();
      runData->Add(kBGO_Scintillation, scint_photon_energy);
      runData->Add(kNum_Scint, 1);
      fEventAction->AddEdepBGOScint(scint_photon_energy);
      primary->SetKineticEnergy(0.);
    }
  
    // loop over secondaries of this particle
    const std::vector <const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
    for(auto sec : *secondaries)
    {
      G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
      if(creator_process.compare("Cerenkov") == 0 && physicalVolume == fDetConstruction->GetCerenkovVolume())
      {
        G4double cher_photon_energy = sec->GetKineticEnergy();
        runData->Add(kBGO_Cherenkov, cher_photon_energy);
        runData->Add(kNum_Cerenkov, 1);
        fEventAction->AddEdepBGOCerenkov(cher_photon_energy);
      }
      
      else if(creator_process.compare("Scintillation") == 0 && physicalVolume == fDetConstruction->GetScintillatorVolume())
      {
        std::cout << OYELLOW << "check Scintillation" << ORESET << std::endl;
        G4double scint_photon_energy = sec->GetKineticEnergy();
        runData->Add(kBGO_Scintillation, scint_photon_energy);
        runData->Add(kNum_Scint, 1);
        fEventAction->AddEdepBGOScint(scint_photon_energy);
      }
    }
    

    // once you register secondary particles, you delete them
    // step->DeleteSecondaryVector();

  }
  */

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  fEventAction->AddEdep(edepStep);

}
