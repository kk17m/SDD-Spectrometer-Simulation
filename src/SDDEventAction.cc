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
/// \file SDDEventAction.cc
/// \brief Implementation of the SDDEventAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SDDEventAction.hh"
#include "SDDRunAction.hh"
#include "SDDdigitize.hh"
#include "SDDDetectorHit.hh"
#include "SDDHistoManager.hh"

#include "G4DigiManager.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDEventAction::SDDEventAction(SDDRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  collID(-1)
{
    SDDdigitize * sddDigiModule = new SDDdigitize("SDDdigitize");
    G4DigiManager::GetDMpointer()->AddNewModule(sddDigiModule);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDEventAction::~SDDEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDEventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(collID == -1)
  {
    collID = SDman->GetCollectionID("SDsdd/DC");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDEventAction::EndOfEventAction(const G4Event* evt)
{   
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  //
  // SDD
  //
  G4HCofThisEvent* eventHC = evt->GetHCofThisEvent();
  DetectorHitsCollection* detHitColl = 0;
  G4DigiManager * pointerDM = G4DigiManager::GetDMpointer();

  if (eventHC)
  {
      detHitColl = (DetectorHitsCollection*)(eventHC->GetHC(collID));
      if (detHitColl)
      {
          G4int nHits = detHitColl->entries();
          if (nHits > 0)
          {
              SDDdigitize * poiterDigiModule = (SDDdigitize*)pointerDM->FindDigitizerModule("SDDdigitize");
              poiterDigiModule->Digitize();
          }
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
