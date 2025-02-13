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

void EventAction::BeginOfEventAction(const G4Event*){

  fEdep = 0.;
  fEdep_BGO = 0.;
  fEdep_Scint1 = 0.;
  fEdep_Scint2 = 0.;
  fEdep_BGO_Cherenkov = 0.;
  fEdep_BGO_Scintillation = 0.;
  Nphotons_Cerenkov = 0;
  Nphotons_Scint = 0;

  IsInBGO = false;
  IsInTrg1 = false;
  IsInTrg2 = false;
  PMT1detection = false;
  PMT2detection = false;

  Nproduced_Cerenkov = 0;
  Nproduced_Scintillation = 0;

  auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset();

}

void EventAction::EndOfEventAction(const G4Event* event){
  
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  //print per event (modulo n)
  auto eventID = event->GetEventID();
  if (( eventID % 10000 == 0 )) {
    G4cout << "-------> End of event: " << eventID << G4endl;
  }

  // output printing <=> particle pass thorugh both the plastic scintillators
  if ( IsInTrg1 && IsInTrg2 && IsInBGO ){

    std::cout << "... detecting a particle ..." << std::endl;
    
    auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    runData->FillPerEvent();

    // adding a particle
    fRunAction->addDetectedParticle();

    // output
    ofstream output;
    output.open("../datasets/July2/p30degrees.txt", ios::app);
    if(!output) cout << OBOLDRED << "ERROR: Could not open the file" << ORESET << endl;
    output << setw(7)
      << fEdep << "\t"
      << fEdep_BGO << "\t"
      << fEdep_Scint1 << "\t"
      << fEdep_Scint2 << "\t"
      << fEdep_BGO_Cherenkov << "\t"
      << fEdep_BGO_Scintillation << "\t"
      << Nphotons_Cerenkov << "\t"
      << Nphotons_Scint << endl;
    output.close();

    output.open("../datasets/July2/numbersP30.txt", ios::app);
    if(!output) cout << OBOLDRED << "ERROR: Could not open the file" << ORESET << endl;
    output << setw(7)
      << Nphotons_Cerenkov << "\t"
      << Nproduced_Cerenkov << "\t"
      << Nphotons_Scint << "\t"
      << Nproduced_Scintillation << endl;
    output.close();
    
  }
  /*
  ofstream output;
  output.open("../analysisDreamSimulation/good_angles.txt", ios::app);
  if(!output) cout << OBOLDRED << "ERROR: Could not open the file" << ORESET << endl;
  if(IsInTrg1 && IsInTrg2 && IsInBGO) output << 1 << endl;
  else output << 0 << endl;
  output.close();
  */
}

void EventAction::PrintStatus(){
  std::cout << OCYAN << "Status of the Event:\n"
  << " ------------- Status of passage in detectors -------------\n"
  << "\tBGO crystal:\t" << IsInBGO << "\n"
  << "\tTrigger1:\t" << IsInTrg1 << "\n"
  << "\tTrigger2:\t" << IsInTrg2 << "\n"
  << " ------------- Energy deposited  -------------\n"
  << "\tTotal:\t" << fEdep << "\n"
  << "\tBGO:\t" << fEdep_BGO << ", in particular:\n"
  << "\t\tCherenkov radiation:\t" << fEdep_BGO_Cherenkov << "\t with " << Nphotons_Cerenkov << "photons\n"
  << "\t\tScintillation:\t" << fEdep_BGO_Scintillation << "\t with " << Nphotons_Scint << "photons\n"
  << "\tTrigger1:\t" << fEdep_Scint1 << "\n"
  << "\tTrigger1:\t" << fEdep_Scint2 << "\n"
  << ORESET << std::endl;
}

void EventAction::AddEdep(G4double edep){
  fEdep += edep;
}

void EventAction::AddEdepBGO(G4double edep){
  fEdep_BGO += edep;
}

void EventAction::AddEdepScint1(G4double edep){
  fEdep_Scint1 += edep;
}

void EventAction::AddEdepScint2(G4double edep){
  fEdep_Scint2 += edep;
}

void EventAction::AddEdepBGOCerenkov(G4double edep){ 
  fEdep_BGO_Cherenkov += edep;
  Nphotons_Cerenkov ++;
}

void EventAction::AddEdepBGOScint(G4double edep){
  fEdep_BGO_Scintillation += edep;
  Nphotons_Scint ++;
}
