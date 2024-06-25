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
// GEANT4 tag $Name: SDD-01-00
//
// Author: Kunal Kumar (kunal.kumar@ovgu.de)
//
// Citation: Kumar K, Fachet M, Hoeschen C. High-Spatial-Resolution Benchtop X-ray Fluorescence Imaging through Bragg-Diffraction-Based Focusing with Bent Mosaic Graphite Crystals: A Simulation Study. Int J Mol Sci. 2024 Apr 26;25(9):4733. doi: 10.3390/ijms25094733. PMID: 38731956; PMCID: PMC11083219.
//
// DOI: 10.3390/ijms25094733; PMID: 38731956; PMCID: PMC11083219
//
// History:
// -----------
// 15 Nov 2023 - Kunal Kumar:     Created. Built using B1 base example
// -------------------------------------------------------------------
//
/// \file SDDDetectorSD.cc
/// \brief Implementation of the SDDDetectorSD class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SDDDetectorSD.hh"
#include "SDDDetectorConstruction.hh"
#include "SDDDetectorHit.hh"

#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "SDDHistoManager.hh"
#include "G4UnitsTable.hh"

SDDDetectorSD::SDDDetectorSD(G4String name, SDDDetectorConstruction* detCon)
    : G4VSensitiveDetector(name),
      HitsColl(0),
      detConstruction(0)
{
    detConstruction = detCon;
    collectionName.insert("DC");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SDDDetectorSD::~SDDDetectorSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SDDDetectorSD::Initialize(G4HCofThisEvent* HCE)
{
    HitsColl = new DetectorHitsCollection(SensitiveDetectorName, collectionName[0]);
    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    HCE->AddHitsCollection(HCID, HitsColl);

    G4RunManager *detConst = G4RunManager::GetRunManager();
    SDDDetectorConstruction *userDetCon = (SDDDetectorConstruction *)(detConst->GetUserDetectorConstruction());
    setUnitTrackWeight = userDetCon->GetSetUnitTrackWeight();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool SDDDetectorSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4double eDep = aStep->GetTotalEnergyDeposit();
    if ((eDep/eV < 1.0))
        return false;

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(2, eDep); // Track weight =1, change if biasing

    // Track weight
    G4double trackWeight;
    if (setUnitTrackWeight == true)
        trackWeight = 1;
    else
        trackWeight = aStep->GetTrack()->GetWeight();

    G4TouchableHistory* theTouchable  = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4int detID = theTouchable->GetVolume(1)->GetCopyNo();

    // Global to local transformation
    G4ThreeVector pos = aStep->GetTrack()->GetPosition();
    const G4AffineTransform transformation =  aStep->GetTrack()->GetTouchable()->GetHistory()->GetTransform(1);
    G4ThreeVector vecLocal = transformation.TransformPoint(pos);

    // Get hit accounting data for this cell
    auto hitSDD = new  SDDDetectorHit();

    if ( ! hitSDD ) {
        G4ExceptionDescription msg;
        msg << "ERROR: in Hit processing";
        G4Exception("DetectorSD::ProcessHits()", "MyCode0004", FatalException, msg);
    }

    // Hit info
    hitSDD->Add(eDep, vecLocal, detID, trackWeight);

    HitsColl->insert(hitSDD);

    if (verboseLevel>0)
    {
        G4RunManager *detConst = G4RunManager::GetRunManager();
        const G4Event *currentEvent = detConst->GetCurrentEvent();

        G4String particle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();

        G4double gtime = aStep->GetPostStepPoint()->GetGlobalTime()/ns;
        G4cout << "### DEBUG(SD): NEW HIT ---> ["
               << "Event: " << currentEvent->GetEventID()
               << ", Particle: " << particle
               << ", Slot: " << detID
               << ", G-Time: " << gtime
               << ", Energy: " << eDep/keV
               << "]" << G4endl;
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SDDDetectorSD::EndOfEvent(G4HCofThisEvent*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
