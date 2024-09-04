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
/// \file SDDdigitize.cc
/// \brief Implementation of the SDDdigitize class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SDDDetectorConstruction.hh"
#include "SDDdigitize.hh"
#include "SDDDetectorHit.hh"
#include "G4DigiManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "SDDHistoManager.hh"
#include "stdio.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SDDdigitize::SDDdigitize(G4String name)
    : G4VDigitizerModule(name)
{
    G4RunManager *pointerRunMan = G4RunManager::GetRunManager();
    SDDDetectorConstruction *detConstruction = (SDDDetectorConstruction *)(pointerRunMan->GetUserDetectorConstruction());
    sensorPz = detConstruction->GetSensorThickness();
    nEH = 3.61;   // avg. energy to produce an electron-hole pair in Silicon (eV)
    fano = 0.115; // fano factor
    DelE = 52.; // ref. Amptek 160mm2 FastSDD, 1us T_peak, 250K temp (https://www.amptek.com/-/media/ametekamptek/documents/resources/application-notes/high-sensitivity-detectors-for-xrf.pdf?la=en&revision=9d04dd37-c2ea-4f89-ad58-55579a8574b1)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SDDdigitize::~SDDdigitize()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SDDdigitize::Digitize()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    G4DigiManager* pointerDM = G4DigiManager::GetDMpointer();
    G4int HCID;
    HCID = pointerDM->GetHitsCollectionID("DC");
    DetectorHitsCollection* detHitColl = 0;
    detHitColl = (DetectorHitsCollection*)(pointerDM->GetHitsCollection(HCID));

    if (detHitColl)
    {
        G4double effEnergy_SDD = 0.;
        G4double depEnergy = 0.;

        G4int fNhits = detHitColl->entries();
        for (G4int i=0; i<fNhits; i++)
        {
            G4double energy = (*detHitColl)[i]->GetEnergy();
            G4double effCharge_Cloud_SDD = NoiseFano((energy/eV)/nEH);
            effEnergy_SDD += effCharge_Cloud_SDD * nEH * eV * (*detHitColl)[i]->GetTrackWeight();
            depEnergy += energy;
        }
        G4double energy_SDD = effEnergy_SDD;
        G4double trackWeight = (*detHitColl)[0]->GetTrackWeight();
        G4double firstZPos = (*detHitColl)[0]->GetPosition().z();

        analysisManager->FillH2(1, depEnergy, sensorPz/2. - firstZPos, trackWeight);

        if (energy_SDD > 0.)
            analysisManager->FillH1(1, energy_SDD, trackWeight);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SDDdigitize::NoiseFano(G4double eCharge)
{
    G4double sigma_eV = std::sqrt(fano * nEH * (eCharge * nEH) + pow(DelE, 2));
    G4double sigma = sigma_eV/nEH; // units: charge
    G4double resEcharge = CLHEP::RandGauss::shoot(eCharge, sigma);
    return resEcharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SDDdigitize::NoiseGaussian(G4double ene)
{
    G4double esigma = DelE*eV/(2.355); // units: eV
    G4double resEne = CLHEP::RandGauss::shoot(ene, esigma);
    return resEne;
}
