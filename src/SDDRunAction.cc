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
/// \file SDDRunAction.cc
/// \brief Implementation of the SDDRunAction class

#include "SDDRunAction.hh"
#include "SDDPrimaryGeneratorAction.hh"
#include "SDDDetectorConstruction.hh"
#include "SDDRun.hh"
#include "SDDHistoManager.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDRunAction::SDDRunAction()
    : G4UserRunAction(),
      fEdep(0.),
      fEdep2(0.),
      fTimer(nullptr),
      fHistoManager(0),
      fRun(0)
{ 
    // add new units for dose
    //
    const G4double milligray = 1.e-3*gray;
    const G4double microgray = 1.e-6*gray;
    const G4double nanogray  = 1.e-9*gray;
    const G4double picogray  = 1.e-12*gray;

    new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
    new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
    new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
    new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

    // Register accumulable to the accumulable manager
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();

    accumulableManager->RegisterAccumulable(fEdep);
    accumulableManager->RegisterAccumulable(fEdep2);

    fHistoManager = new SDDHistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDRunAction::~SDDRunAction()
{
    delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* SDDRunAction::GenerateRun()
{
    fRun = new SDDRun();
    return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDRunAction::BeginOfRunAction(const G4Run* run)
{ 
    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    // reset accumulables to their initial values
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Reset();

    if (IsMaster())
    {
        // Timer start
        fTimer = new G4Timer();
        fTimer->Start();
        G4cout << "### Run " << run->GetRunID() << " starts (master)." << G4endl;
    }

    //histograms
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->IsActive() )
    {
        analysisManager->OpenFile();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDRunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;

    // Merge accumulables
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Merge();

    // Compute dose = total energy deposit in a run and its variance
    //
    const SDDDetectorConstruction* detectorConstruction
            = static_cast<const SDDDetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4double edep  = fEdep.GetValue();
    G4double edep2 = fEdep2.GetValue();

    G4double rms = edep2 - edep*edep/nofEvents;
    if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

    G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
    G4double dose = edep/mass;
    G4double rmsDose = rms/mass;

    G4cout
            << G4endl
            << " Cumulated dose per run, in scoring volume : "
            << G4endl
            << " Sample dose : "
            << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
            << " mass = " << G4BestUnit(mass,"Mass")
            << G4endl
            << "------------------------------------------------------------"
            << G4endl
            << G4endl;

    //save histograms
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->IsActive() )
    {
        analysisManager->Write();
        analysisManager->CloseFile();
    }

    // Print
    //
    if (IsMaster())
    {
        G4cout
                << G4endl
                << "--------------------End of Global Run-----------------------";

        // Timer stop
        fTimer->Stop();
        if(!((G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::sequentialRM)))
        {
            G4cout << "\n" << "Number of events processed in run " << run->GetRunID() << ":  " << run->GetNumberOfEventToBeProcessed() << G4endl;
            G4cout << "Master thread time:  "  << *fTimer << "\n" << G4endl;
        }
        delete fTimer;
    }
    else
    {
        G4cout
                << G4endl
                << "--------------------End of Local Run------------------------";
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDRunAction::PrintTableRow(const std::initializer_list<std::string>& columns)
{
    const int columnWidth = 35;

    G4cout << std::fixed << std::setprecision(2);

    for (const auto& col : columns) {
        G4cout << std::setw(columnWidth) << std::left << col;
    }

    G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDRunAction::AddEdep(G4double edep)
{
    fEdep  += edep;
    fEdep2 += edep*edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

