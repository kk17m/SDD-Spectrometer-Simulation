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
/// \file SDDSteppingAction.cc
/// \brief Implementation of the SDDSteppingAction class

#include "SDDSteppingAction.hh"
#include "SDDEventAction.hh"
#include "SDDDetectorConstruction.hh"
#include "SDDHistoManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDSteppingAction::SDDSteppingAction(SDDEventAction* eventAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fScoringVolume(0)
{
    const SDDDetectorConstruction* userDetCon
            = static_cast<const SDDDetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    setUnitTrackWeight = userDetCon->GetSetUnitTrackWeight();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDSteppingAction::~SDDSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDSteppingAction::UserSteppingAction(const G4Step* step)
{
    //
    // Check next volume
    const G4VPhysicalVolume* nextVol =  step->GetTrack()->GetNextVolume();
    if(nextVol == nullptr)
        return;

    if (!fScoringVolume)
    {
        const SDDDetectorConstruction* detectorConstruction
                = static_cast<const SDDDetectorConstruction*>
                (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume = detectorConstruction->GetScoringVolume();
    }

    //
    // Essentials
    //
    const G4StepPoint* preStepPoint = step->GetPreStepPoint();
    const G4StepPoint* postStepPoint = step->GetPostStepPoint();
    const G4String particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    const G4LogicalVolume* volume = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    const G4String preStepVol = volume->GetName();
    const G4String postStepVol = postStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
    const G4int stepNo =  step->GetTrack()->GetCurrentStepNumber();

    //
    // Set track for special case
    //
    if(stepNo == 1 && particle == "gamma" && step->GetTrack()->GetParentID() >= 1)
    {
        if (volume->GetMaterial()->GetName() != "CONTRAST" && preStepVol == "phantom")
            step->GetTrack()->SetTrackStatus(fStopAndKill);
    }

    //
    // Analysis
    //
    G4double tW;
    if (setUnitTrackWeight == true)
        tW = 1; // Disregard track weight
    else
        tW = step->GetTrack()->GetWeight(); // track weight considering biasing (e.g., from source)

    if ((preStepVol == "World" || preStepVol == "MLC" || preStepVol == "ChipSubstrate") && postStepVol == "Pixel" && particle == "gamma")
    {
        // Incident total sum spectrum in all detectors
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        const G4double ekin = step->GetTrack()->GetKineticEnergy();
        analysisManager->FillH1(3, ekin, tW);
    }

    //
    // Dose calculation. Kept at last!
    //
    if (volume != fScoringVolume) return;

    const G4double edepStep = step->GetTotalEnergyDeposit();
    if (edepStep > 0.)
        fEventAction->AddEdep(edepStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

