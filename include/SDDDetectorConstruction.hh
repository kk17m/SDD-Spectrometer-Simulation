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
/// \file SDDDetectorConstruction.hh
/// \brief Definition of the SDDDetectorConstruction class

#ifndef SDDDetectorConstruction_h
#define SDDDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <set>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Box;
class G4Material;

/// Detector construction class to define materials and geometry.

class SDDDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    SDDDetectorConstruction(G4bool visualize);
    virtual ~SDDDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    void ConstructSDandField();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4bool GetVisBool() const { return visBool; }
    G4bool GetVisBoolSpecific() const { return visBoolSpecific; }
    G4double GetfOffsetZ() const { return fOffsetZ; }

    std::vector<G4double> GetTargetPos() const { return targetPos; }

    // Detector
    G4bool GetSetUnitTrackWeight() const {return setUnitTrackweight;}
    G4double GetPixelSize() const {return pixelDia;}
    G4double GetSensorThickness() const {return sensorThickness;}
    G4int GetNbSensors() const {return nSensorsPerRow;}

protected:

    void PhantomUniform(); // Uniform phantom containing sample material

    //
    // Detector
    //
    void ConstructDetector();

protected:
    G4bool visBool;
    G4bool visBoolSpecific;
    G4Material* Fluo_Solution;

    G4LogicalVolume* fWorld_logic;
    G4VPhysicalVolume* fWorld_phys;

    G4LogicalVolume*  fScoringVolume;

    G4double fOffsetZ;

    //
    // Target sphere
    //
    G4int tcx;
    G4int tcy;
    G4int tcz;
    std::vector<G4double> targetPos;

    //
    // Detector
    //

    // Track weight
    G4bool setUnitTrackweight;

    // Sensor configuration (Amptek XR100 FastSDD)
    G4double sensorThickness;
    G4double pixelDia;
    G4double activeDia;
    G4double substrateThickness;
    G4double botElectrodeHeight;
    G4double tecThickness;
    G4double headerThickness;
    G4double coverHeight;
    G4double windowThickness;
    G4double studThickness;
    G4double extenderThickness;
    G4double elecBoxLength;

    // Sensor
    G4LogicalVolume* logicPixel;
    G4double sensorY;
    G4double sensorX;

    std::vector<G4double> polarAngle;
    std::vector<G4double> detAngle;
    std::vector<G4double> detSensorZ;
    G4int nSensorsPerRow;
    G4double setupBase_offset;

    G4double sampleDia;
    G4double samplePz;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

