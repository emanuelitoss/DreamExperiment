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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


class RunAction;

/// Event action class
class EventAction : public G4UserEventAction
{

  public:

  EventAction(RunAction* runAction);
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  void PassedThroughScint1();
  void PassedThroughScint2();

  void AddEdep(G4double edep);
  void AddEdepBGO(G4double edep);
  void AddEdepScint1(G4double edep);
  void AddEdepScint2(G4double edep);

  void PrintStatus();

  private:

  RunAction* fRunAction;
  
  // deposited energies in: ScoringVolume, BGO crystal and two plastic Scintillators
  G4double fEdep;
  G4double fEdep_BGO;
  G4double fEdep_Scint1;
  G4double fEdep_Scint2;
  
  // boolean variables to check if the particle pass thorugh physical volumes
  G4bool IsInTrg1 = false;
  G4bool IsInTrg2 = false;
  
  // private method
  void PrintEventStatistics(G4double EdepBGO, G4double EdepTRG1, G4double EdepTRG2) const;

};

inline void EventAction::PassedThroughScint1(){
  IsInTrg1 = true;
}

inline void EventAction::PassedThroughScint2(){
  IsInTrg2 = true;
}

#endif
