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
/// \file SDDSteppingAction.hh
/// \brief Definition of the SDDSteppingAction class

#ifndef SDDSteppingAction_h
#define SDDSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

class SDDEventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class SDDSteppingAction : public G4UserSteppingAction
{
  public:
    SDDSteppingAction(SDDEventAction* eventAction);
    virtual ~SDDSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    SDDEventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
    G4bool setUnitTrackWeight;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
