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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "../include/EventAction.hh"
#include "../include/RunAction.hh"
#include "../include/RunData.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{}

EventAction::~EventAction(){}

void EventAction::PrintEventStatistics(G4double EdepBGO, G4double EdepTRG1, G4double EdepTRG2) const {
 // print event statistics
  G4cout
    << "   BGO: total energy: "
    << std::setw(7) << G4BestUnit(EdepBGO, "Energy")
    << G4endl
    << "   PlasticTrigger1: total energy: "
    << std::setw(7) << G4BestUnit(EdepTRG1, "Energy")
    << G4endl
    << "   PlasticTrigger2: total energy: "
    << std::setw(7) << G4BestUnit(EdepTRG2, "Energy")
    << G4endl;
}

void EventAction::BeginOfEventAction(const G4Event*){

  fEdep = 0.;
  fEdep_BGO = 0.;
  fEdep_PMT1 = 0.;
  fEdep_PMT2 = 0.;

  auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset();

}

void EventAction::EndOfEventAction(const G4Event* event){
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->FillPerEvent();

  //print per event (modulo n)
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;
  
    PrintEventStatistics(
      runData->GetEdep(kBGO),
      runData->GetEdep(kScint1),
      runData->GetEdep(kScint2));
  }
}

void EventAction::AddEdep(G4double edep){
  fEdep += edep;
  std::cout << "ENERGIA TOTALE: " << fEdep << std::endl;
}

void EventAction::AddEdepBGO(G4double edep){
  fEdep_BGO += edep;
  std::cout << "ENERGIA NEL BGO: " << fEdep_BGO << std::endl;
}

void EventAction::AddEdepPMT1(G4double edep){
  fEdep_PMT1 += edep;
}

void EventAction::AddEdepPMT2(G4double edep){
  fEdep_PMT2 += edep;
}
