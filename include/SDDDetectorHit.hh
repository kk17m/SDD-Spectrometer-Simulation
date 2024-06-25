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
/// \file SDDDetectorHit.hh
/// \brief Definition of the SDDDetectorHit class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SDDDetectorHit_HH
#define SDDDetectorHit_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class SDDDetectorHit : public G4VHit
{
public:
    SDDDetectorHit();
    ~SDDDetectorHit();

    SDDDetectorHit(const SDDDetectorHit&);
    const SDDDetectorHit& operator=(const SDDDetectorHit&);
    int operator==(const SDDDetectorHit&) const;
    inline void* operator new(size_t);
    inline void  operator delete(void*);

public:
    inline void Add(G4double eDep, G4ThreeVector pos, G4int detID, G4double trackWeight)
    {
        fEdep += eDep;
        fPosition = pos;
        fTrackWeight = trackWeight;
        fdetID = detID;
    }

    inline G4double GetEnergy() {return fEdep;}
    inline G4ThreeVector GetPosition() {return fPosition;}
    inline G4int GetdetID() {return fdetID;}
    inline G4double GetTrackWeight() {return fTrackWeight;}

private:
    G4double fEdep;
    G4ThreeVector fPosition;
    G4double fTrackWeight;
    G4int fdetID;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<SDDDetectorHit> DetectorHitsCollection;

extern G4ThreadLocal G4Allocator<SDDDetectorHit> *AllocHit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* SDDDetectorHit::operator new(size_t)
{
  if (!AllocHit)
    AllocHit = new G4Allocator<SDDDetectorHit>;
  return (void*) AllocHit->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void SDDDetectorHit::operator delete(void* aHit)
{
  AllocHit->FreeSingle((SDDDetectorHit*) aHit);
}

#endif
