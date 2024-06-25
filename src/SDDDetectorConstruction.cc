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
/// \file SDDDetectorConstruction.cc
/// \brief Implementation of the SDDDetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SDDDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDDetectorConstruction::SDDDetectorConstruction(G4bool visualize)
    : G4VUserDetectorConstruction(),
      visBool(false),
      visBoolSpecific(false),
      Fluo_Solution(0),
      fWorld_logic(0),
      fWorld_phys(0),
      fScoringVolume(0),
      fOffsetZ(0.),
      tcx(0),
      tcy(0),
      tcz(0),
      setUnitTrackweight(false),
      sensorThickness(0.),
      pixelDia(0),
      activeDia(0),
      substrateThickness(0),
      botElectrodeHeight(0),
      tecThickness(0),
      headerThickness(0),
      coverHeight(0),
      windowThickness(0),
      studThickness(0),
      extenderThickness(0),
      elecBoxLength(0),
      logicPixel(0),
      sensorY(0),
      sensorX(0),
      nSensorsPerRow(0),
      setupBase_offset(0),
      sampleDia(0),
      samplePz(0)
{
    visBool = visualize;
    visBoolSpecific = false;

    targetPos.resize(3);
    std::fill(targetPos.begin(), targetPos.end(), 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDDDetectorConstruction::~SDDDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SDDDetectorConstruction::Construct()
{
    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();

    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4double world_sizeXY = 600.*mm;
    G4double world_sizeZ  = 0.5*m;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box* solidWorld = new G4Box("World",
                                  0.5*world_sizeXY,
                                  0.5*world_sizeXY,
                                  0.5*world_sizeZ);

    fWorld_logic = new G4LogicalVolume(solidWorld,
                                       world_mat,
                                       "World");

    fWorld_phys = new G4PVPlacement(0,
                                    G4ThreeVector(),
                                    fWorld_logic,
                                    "World",
                                    0,
                                    false,
                                    0,
                                    checkOverlaps);


    // Build phantom and fluorescence material
    PhantomUniform();

    // SDD Detector
    ConstructDetector();

    return fWorld_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDDetectorConstruction::PhantomUniform()
{
    // Tissue material
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* materialTissue = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");

    // Fluorescence material
    G4double density;
    G4Material* Mo = nist->FindOrBuildMaterial("G4_Mo");
    G4Material* I  = nist->FindOrBuildMaterial("G4_I");
    G4double nistDensity_Mo = Mo->GetDensity();
    G4double nistDensity_I = I->GetDensity();
    G4double nistDensity_phantom = materialTissue->GetDensity();

    G4double fContrastFmassMo = 0.01;
    G4double fContrastFmassI = 0.01;
    G4double fmass_Tissue = 1. - fContrastFmassMo - fContrastFmassI;
    density = fContrastFmassMo * nistDensity_Mo
            + fContrastFmassI * nistDensity_I
            + fmass_Tissue * nistDensity_phantom;

    Fluo_Solution = new G4Material("CONTRAST", density, 3);
    Fluo_Solution->AddMaterial(materialTissue, fmass_Tissue);
    Fluo_Solution->AddMaterial(Mo, fContrastFmassMo);
    Fluo_Solution->AddMaterial(I, fContrastFmassI);

    //
    // Sample
    //
    sampleDia = 20.*mm;
    samplePz = 70.*mm;
    G4Tubs* solidSample = new G4Tubs("phantom", 0, sampleDia/2, samplePz/2,  0, twopi);
    G4LogicalVolume* logicSample = new G4LogicalVolume(solidSample, materialTissue, "phantom");

    G4Region* sample_Region = new G4Region("phantom");
    sample_Region->AddRootLogicalVolume(logicSample);

    G4double shiftSampleY = 0.*mm;
    fOffsetZ = 0.*mm;
    G4ThreeVector sampleTrans = G4ThreeVector(0.,0.,0.);
    sampleTrans.setZ(fOffsetZ);
    sampleTrans.setY(shiftSampleY);
    sampleTrans.setX(0.);

    G4RotationMatrix rot2 = G4RotationMatrix();
    rot2.rotateX(90.*deg);

    G4Transform3D transform = G4Translate3D(sampleTrans)*G4Rotate3D(rot2);
    new G4PVPlacement(transform,
                      logicSample,
                      "phantom",
                      fWorld_logic,
                      true,
                      0,
                      true
                      );

    G4VisAttributes* vis  = new G4VisAttributes();
    vis->SetForceWireframe(false);
    vis->SetColour(0.9, 0., 0., 0.9);
    vis->SetForceSolid(true);
    logicSample->SetVisAttributes(vis);

    //
    // Sphere target
    //
    G4double Dmax = 5.*mm;
    G4Sphere* solidTarget = new G4Sphere("Target", 0., Dmax/2., 0, twopi, 0, twopi);
    G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, Fluo_Solution, "Target");

    G4Region* target_Region = new G4Region("Target");
    target_Region->AddRootLogicalVolume(logicTarget);

    G4RotationMatrix rotTarget = G4RotationMatrix();
    rotTarget.rotateX(-90.*deg);

    G4ThreeVector targetTrans = G4ThreeVector(0.,0.,0.);

    G4Transform3D transformTarget = G4Translate3D(targetTrans)*G4Rotate3D(rotTarget);

    new G4PVPlacement(transformTarget,
                      logicTarget,
                      "Target",
                      logicSample,
                      true,
                      0,
                      true
                      );

    G4VisAttributes* targetColor = new G4VisAttributes(G4Colour(0., 1., 0., 0.8));
    targetColor->SetForceSolid(true);
    logicTarget->SetVisAttributes(targetColor);

    fScoringVolume = logicSample;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDDDetectorConstruction::ConstructDetector()
{
    // Sensor configuration (Amptek XR100 FastSDD)
    sensorThickness = 1000*um;                  // Sensor thickness XR100 FastSDD
    pixelDia = 9.44*mm;                         // dia for 70mm2 XR100 FastSDD
    activeDia = 7.11*mm;                        // dia for 40mm2 active area collimated XR100 FastSDD
    substrateThickness = 750.*um;               // Sensor substrate
    botElectrodeHeight = 50.*nm;                // Anode
    tecThickness = 4.5*mm;                      // Thermoelectric cooler
    headerThickness = 1.5*mm;                   // TO-8 header (Kovar)
    coverHeight = 8.636*mm - headerThickness;   // Sensor external cover
    windowThickness = 12.5*um;                  // Beryllium window
    studThickness = 6.096*mm;                   // Lower stud
    extenderThickness = 38.1*mm;                // Sensor extender
    elecBoxLength = 30.94*mm;                   // Electronics box

    // Distance params
    sensorY = 0.*mm;   // Sensor y-position
    sensorX = 0.*mm;   // Sensor x-position

    // Placement of detectors (here two detectors are used)
    G4double azimuthDet = 90.*deg;
    detAngle = {azimuthDet, -azimuthDet};
    polarAngle = {30.*deg, 30.*deg};
    detSensorZ = {-32.*mm, -32.*mm};
    nSensorsPerRow = static_cast<G4int>(detAngle.size());

    G4bool checkOverlaps = true;
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* worldMaterial = nist->FindOrBuildMaterial("G4_AIR");

    //
    // Sensor
    //
    G4Material* sensorMaterial = nist->FindOrBuildMaterial("G4_Si");

    // Create pixel. Should only be nulled in the constructor otherwise no hits !!
    G4Tubs* solidPixel = new G4Tubs("Pixel",
                                    0.,
                                    pixelDia/2.,
                                    sensorThickness/2.,
                                    0.,
                                    twopi
                                    );

    logicPixel = new G4LogicalVolume(solidPixel,
                                     sensorMaterial,
                                     "Pixel"
                                     );

    G4Region* Pixel_Region = new G4Region("Pixel");
    Pixel_Region->AddRootLogicalVolume(logicPixel);

    //
    // ANODE Amptek
    //
    G4Tubs* solidAnode = new G4Tubs("Anode",
                                    0.,
                                    activeDia/2.,
                                    botElectrodeHeight/2.,
                                    0.,
                                    twopi
                                    );

    G4Material* anodeMaterial = nist->FindOrBuildMaterial("G4_Al");

    G4LogicalVolume* logicAnode = new G4LogicalVolume(solidAnode,
                                                      anodeMaterial,
                                                      "Anode"
                                                      );

    //
    // Sensor mother volume. Mother needs to be nulled here !!
    //
    G4Tubs* solidSensorMother = new G4Tubs("SensorMother",
                                           0.,
                                           pixelDia/2.,
                                           sensorThickness/2.,
                                           0.,
                                           twopi
                                           );

    G4LogicalVolume* logicSensorMother = new G4LogicalVolume(solidSensorMother,
                                                             worldMaterial,
                                                             "SensorMother"
                                                             );

    //
    // MULTILAYER COLLIMATOR SDD
    //
    G4double pz_MLCmother = 100.*um + 35.*um + 15.*um + 75.*um;
    G4Tubs* solidMLCmother = new G4Tubs("MLC",
                                        0.,
                                        (pixelDia + 0.5*mm)/2.,
                                        pz_MLCmother/2.,
                                        0.,
                                        twopi
                                        );

    G4LogicalVolume* logicMLCmother = new G4LogicalVolume(solidMLCmother,
                                                          worldMaterial,
                                                          "MLC"
                                                          );

    // W
    G4Tubs* solidMLCW = new G4Tubs("MLC_W",
                                   7.11*mm/2.,
                                   (pixelDia + 0.5*mm)/2.,
                                   100.*um/2.,
                                   0.,
                                   twopi
                                   );

    G4Material* materialMLCW = nist->FindOrBuildMaterial("G4_W");
    G4LogicalVolume* logicMLCW = new G4LogicalVolume(solidMLCW,
                                                     materialMLCW,
                                                     "MLC_W"
                                                     );

    // Cr
    G4Tubs* solidMLCCr = new G4Tubs("MLC_Cr",
                                    7.11*mm/2.,
                                    (pixelDia + 0.5*mm)/2.,
                                    35.*um/2.,
                                    0.,
                                    twopi
                                    );

    G4Material* materialMLCCr = nist->FindOrBuildMaterial("G4_Cr");
    G4LogicalVolume* logicMLCCr = new G4LogicalVolume(solidMLCCr,
                                                      materialMLCCr,
                                                      "MLC_W"
                                                      );

    // Ti
    G4Tubs* solidMLCTi = new G4Tubs("MLC_Ti",
                                    7.11*mm/2.,
                                    (pixelDia + 0.5*mm)/2.,
                                    15.*um/2.,
                                    0.,
                                    twopi
                                    );

    G4Material* materialMLCTi = nist->FindOrBuildMaterial("G4_Ti");
    G4LogicalVolume* logicMLCTi = new G4LogicalVolume(solidMLCTi,
                                                      materialMLCTi,
                                                      "MLC_Ti"
                                                      );

    // Al
    G4Tubs* solidMLCAl = new G4Tubs("MLC_Al",
                                    7.11*mm/2.,
                                    (pixelDia + 0.5*mm)/2.,
                                    75.*um/2.,
                                    0.,
                                    twopi
                                    );

    G4Material* materialMLCAl = nist->FindOrBuildMaterial("G4_Al");
    G4LogicalVolume* logicMLCAl = new G4LogicalVolume(solidMLCAl,
                                                      materialMLCAl,
                                                      "MLC_Al"
                                                      );


    //
    // Anode mother volume. Mother needs to be nulled here !!
    //
    G4Tubs* solidElectrodeMother = new G4Tubs("AnodeMother",
                                              0.,
                                              pixelDia/2.,
                                              botElectrodeHeight/2.,
                                              0.,
                                              twopi
                                              );

    G4LogicalVolume* logicAnodeMother = new G4LogicalVolume(solidElectrodeMother,
                                                            worldMaterial,
                                                            "AnodeMother"
                                                            );

    //
    // SUBSTRATE Amptek XR100
    //
    G4Tubs* solidSubstrate = new G4Tubs("ChipSubstrate",
                                        0.,
                                        (pixelDia + 2.*mm)/2.,
                                        substrateThickness/2.,
                                        0.,
                                        twopi
                                        );

    G4Material* substrateMaterial = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    G4LogicalVolume* substrateChip = new G4LogicalVolume(solidSubstrate,
                                                         substrateMaterial,
                                                         "ChipSubstrate");

    //
    // THERMOELECTRIC COOLER Amptek XR100
    //

    // Si columns as a box
    G4Box* solidTECcolumn = new G4Box("columnTEC",
                                      (pixelDia - 2.*mm)/2.,
                                      (pixelDia - 2.*mm)/2.,
                                      tecThickness/2.
                                      );

    G4Material* tecColumnMaterial = nist->FindOrBuildMaterial("G4_Si");
    G4LogicalVolume* logicTECcolumn = new G4LogicalVolume(solidTECcolumn,
                                                          tecColumnMaterial,
                                                          "columnTEC");

    // Tin-Antimony solder as plate inside mother (columnTEC)
    G4Box* solidTECsolder = new G4Box("solderTEC",
                                      6.*mm/2.,
                                      6.*mm/2.,
                                      0.001*mm/2.
                                      );

    G4LogicalVolume* logicTECsolder = new G4LogicalVolume(solidTECsolder,
                                                          worldMaterial,
                                                          "solderTEC");


    //
    // HEADER KOVAR Amptek XR100
    //
    G4Tubs* solidHeader = new G4Tubs("header", 0., 15.24*mm/2., headerThickness/2., 0, twopi);

    G4Element* elC = new G4Element("Carbon", "C", 6, 12.01*g/mole);
    G4Element* elSi = new G4Element("Silicon", "Si", 14, 28.0855*g/mole);
    G4Element* elMn = new G4Element("Manganese", "Mn", 25, 54.938044*g/mole);
    G4Element* elFe = new G4Element("Iron", "Fe", 26, 55.845*g/mole);
    G4Element* elCo = new G4Element("Cobalt", "Co", 27, 58.933195*g/mole);
    G4Element* elNi = new G4Element("Nickel", "Ni", 28, 58.6934*g/mole);
    G4double densityKovar = 8.36*g/cm3;
    G4Material* headerMaterial = new G4Material("Kovar", densityKovar, 6);
    headerMaterial->AddElement(elC, 0.01*perCent);
    headerMaterial->AddElement(elSi, 0.2*perCent);
    headerMaterial->AddElement(elMn, 0.3*perCent);
    headerMaterial->AddElement(elFe, 53.49*perCent);
    headerMaterial->AddElement(elCo, 17.*perCent);
    headerMaterial->AddElement(elNi, 29.*perCent);

    G4LogicalVolume* logicHeader = new G4LogicalVolume(solidHeader,
                                                       headerMaterial,
                                                       "header");

    //
    // COVER Amptek XR100
    //
    G4Tubs* solidCoverBase = new G4Tubs("coverBase", 0., 13.97*mm/2., coverHeight/2., 0, twopi);
    G4Tubs* solidCoverS1 = new G4Tubs("coverS1", 0., 13.47*mm/2., coverHeight/2., 0, twopi);
    G4Tubs* solidCoverS2 = new G4Tubs("coverS2", 0., 9.0*mm/2., coverHeight/2., 0, twopi);

    G4SubtractionSolid* solidCoverSub1= new G4SubtractionSolid("cover", solidCoverBase, solidCoverS1, 0 , G4ThreeVector(0,0,-0.25*mm));
    G4SubtractionSolid* solidCover= new G4SubtractionSolid("cover", solidCoverSub1, solidCoverS2, 0 , G4ThreeVector(0,0,coverHeight/2.));

    G4Material* coverMaterial = nist->FindOrBuildMaterial("G4_Ni");

    G4LogicalVolume* logicCover = new G4LogicalVolume(solidCover,
                                                      coverMaterial,
                                                      "cover");

    //
    // Be WINDOW Amptek XR100
    //
    G4Tubs* solidWindow = new G4Tubs("window", 0., 9.5*mm/2., windowThickness/2., 0, twopi);

    G4Material* windowMaterial = nist->FindOrBuildMaterial("G4_Be");

    G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,
                                                       windowMaterial,
                                                       "window");

    //
    // SOLDER Be WINDOW Amptek XR100
    //
    G4Tubs* solidSolder = new G4Tubs("solder", 9.5*mm/2., 10.2*mm/2., windowThickness/2., 0, twopi);

    G4Material* solderMaterial = nist->FindOrBuildMaterial("G4_Sn");

    G4LogicalVolume* logicSolder = new G4LogicalVolume(solidSolder,
                                                       solderMaterial,
                                                       "solder");

    //
    // MOUNTING STUD Amptek XR100
    //
    G4Tubs* solidStud = new G4Tubs("stud", 0., 2.845*mm/2., studThickness/2., 0, twopi);

    G4Element* elP = new G4Element("Phosphorus", "P", 15, 30.973762*g/mole);
    G4Element* elS = new G4Element("Sulfer", "S", 16, 32.065*g/mole);
    G4Element* elCr = new G4Element("Chromium", "Cr", 24, 51.9961*g/mole);
    G4double densityStud = 8.0*g/cm3;
    G4Material* studMaterial = new G4Material("studMaterial", densityStud, 8);
    studMaterial->AddElement(elC, 0.04*perCent);
    studMaterial->AddElement(elSi, 0.5*perCent);
    studMaterial->AddElement(elP, 0.02*perCent);
    studMaterial->AddElement(elS, 0.02*perCent);
    studMaterial->AddElement(elCr, 19.*perCent);
    studMaterial->AddElement(elMn, 1.*perCent);
    studMaterial->AddElement(elFe, 70.16*perCent);
    studMaterial->AddElement(elNi, 9.26*perCent);

    G4LogicalVolume* logicStud= new G4LogicalVolume(solidStud,
                                                    studMaterial,
                                                    "stud");

    //
    // EXTENDER TUBE Amptek XR100
    //
    G4Tubs* solidExtendM = new G4Tubs("ExtenderMain", 0., 17.8*mm/2., extenderThickness/2., 0, twopi);
    G4Tubs* solidExtendS = new G4Tubs("ExtenderMain", 0., 2.845*mm/2., extenderThickness/1.5, 0, twopi);
    G4SubtractionSolid* solidExtend= new G4SubtractionSolid("Extender", solidExtendM, solidExtendS, 0 , G4ThreeVector(0,0,0));

    G4Material* extender_mat = nist->FindOrBuildMaterial("G4_Al");
    G4LogicalVolume* logicExtend= new G4LogicalVolume(solidExtend,
                                                      extender_mat,
                                                      "Extender");

    G4double pX_ElecBox = 8.79*mm;
    G4Box* solidElecBox = new G4Box("ElecBox",
                                    pX_ElecBox/2.,
                                    25.4*mm/2.,
                                    elecBoxLength/2.
                                    );

    G4Material* ElecBox_mat = nist->FindOrBuildMaterial("G4_Al");
    G4LogicalVolume* logicElecBox = new G4LogicalVolume(solidElecBox,
                                                        ElecBox_mat,
                                                        "ElecBox");


    //
    // Physical volumes
    //

    // Pixels
    G4ThreeVector shiftPixels = G4ThreeVector(0.,0.,0.);
    G4Transform3D transformPixel = G4Translate3D(shiftPixels);

    new G4PVPlacement(transformPixel,
                      logicPixel,
                      "Pixel",
                      logicSensorMother,
                      false,
                      0,
                      checkOverlaps
                      );

    // MLC
    G4ThreeVector shiftMLCAl = G4ThreeVector(0.,0.,0.);
    shiftMLCAl.setZ(-pz_MLCmother/2 + 75.*um/2);
    G4Transform3D transformMLCAl = G4Translate3D(shiftMLCAl);

    new G4PVPlacement(transformMLCAl,
                      logicMLCAl,
                      "MLC_Al",
                      logicMLCmother,
                      false,
                      0,
                      checkOverlaps
                      );

    G4ThreeVector shiftMLCTi = G4ThreeVector(0.,0.,0.);
    shiftMLCTi.setZ(75.*um/2 + 15.*um/2);
    G4Transform3D transformMLCTi = transformMLCAl*G4Translate3D(shiftMLCTi);

    new G4PVPlacement(transformMLCTi,
                      logicMLCTi,
                      "MLC_Ti",
                      logicMLCmother,
                      false,
                      0,
                      checkOverlaps
                      );

    G4ThreeVector shiftMLCCr = G4ThreeVector(0.,0.,0.);
    shiftMLCCr.setZ(15.*um/2 + 35.*um/2);
    G4Transform3D transformMLCCr = transformMLCTi*G4Translate3D(shiftMLCCr);

    new G4PVPlacement(transformMLCCr,
                      logicMLCCr,
                      "MLC_Cr",
                      logicMLCmother,
                      false,
                      0,
                      checkOverlaps
                      );

    G4ThreeVector shiftMLCW = G4ThreeVector(0.,0.,0.);
    shiftMLCW.setZ(35.*um/2 + 100.*um/2);
    G4Transform3D transformMLCW = transformMLCCr*G4Translate3D(shiftMLCW);

    new G4PVPlacement(transformMLCW,
                      logicMLCW,
                      "MLC_W",
                      logicMLCmother,
                      false,
                      0,
                      checkOverlaps
                      );


    // Anode
    G4ThreeVector shiftAnode = G4ThreeVector(0.,0.,0.);
    G4Transform3D transformAnode = G4Translate3D(shiftAnode);

    new G4PVPlacement(transformAnode,
                      logicAnode,
                      "Anode",
                      logicAnodeMother,
                      true,
                      0,
                      checkOverlaps
                      );

    // solder in TEC
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      logicTECsolder,
                      "solderTEC",
                      logicTECcolumn,
                      true,
                      0,
                      checkOverlaps);

    G4ThreeVector phantomTrans = G4ThreeVector(0.,0.,0.);
    phantomTrans.setZ(fOffsetZ);
    G4int sensorCpy = 0;

    for(G4int sensorNumber=0; sensorNumber<nSensorsPerRow; sensorNumber++)
    {
        G4RotationMatrix rot2 = G4RotationMatrix();
        rot2.rotateY(detAngle.at(sensorNumber));

        G4RotationMatrix rot3 = G4RotationMatrix();
        rot3.rotateX(polarAngle.at(sensorNumber));

        G4ThreeVector sensorTrans = G4ThreeVector(0.,0.,0.);
        sensorTrans.setZ(detSensorZ.at(sensorNumber));

        G4Transform3D transformSensor = G4Translate3D(phantomTrans)*G4Rotate3D(rot2)*G4Rotate3D(rot3)*G4Translate3D(sensorTrans);

        new G4PVPlacement(transformSensor,
                          logicSensorMother,
                          "SensorMother",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps
                          );

        // Bottom electrode
        G4ThreeVector translateAnodeMother = G4ThreeVector(0.,0.,0.);
        translateAnodeMother.setZ(-(sensorThickness/2.0 + botElectrodeHeight/2.0));

        G4Transform3D transformAnodeMother = transformSensor*G4Translate3D(translateAnodeMother);

        new G4PVPlacement(transformAnodeMother,
                          logicAnodeMother,
                          "AnodeMother",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps
                          );

        // MLC
        G4ThreeVector translateMLCmother = G4ThreeVector(0.,0.,0.);
        translateMLCmother.setZ(sensorThickness/2.0 + pz_MLCmother/2.);

        G4Transform3D transformMLCmother = transformSensor*G4Translate3D(translateMLCmother);

        new G4PVPlacement(transformMLCmother,
                          logicMLCmother,
                          "MLCmother",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps
                          );

        // Substrate
        G4ThreeVector substrateTranslation = G4ThreeVector(0.,0.,0.);
        substrateTranslation.setZ(- sensorThickness /  2.0 - botElectrodeHeight - substrateThickness / 2.);

        G4Transform3D transformSubstrate = transformSensor*G4Translate3D(substrateTranslation);

        new G4PVPlacement(transformSubstrate,
                          substrateChip,
                          "ChipSubstrate",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // TEC column
        G4ThreeVector translationTECcolumn = G4ThreeVector(0.,0.,0.);
        translationTECcolumn.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness/2.);

        G4Transform3D transformTECcolumn = transformSensor*G4Translate3D(translationTECcolumn);

        new G4PVPlacement(transformTECcolumn,
                          logicTECcolumn,
                          "columnTEC",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // Header
        G4ThreeVector translationHeader = G4ThreeVector(0.,0.,0.);
        translationHeader.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness - headerThickness/2.);

        G4Transform3D transformHeader = transformSensor*G4Translate3D(translationHeader);

        new G4PVPlacement(transformHeader,
                          logicHeader,
                          "header",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // Cover
        G4ThreeVector translationCover = G4ThreeVector(0.,0.,0.);
        translationCover.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness + coverHeight/2.);

        G4Transform3D transformCover = transformSensor*G4Translate3D(translationCover);

        new G4PVPlacement(transformCover,
                          logicCover,
                          "cover",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // Window
        G4ThreeVector translationWin = G4ThreeVector(0.,0.,0.);
        translationWin.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness + coverHeight - 0.25*mm - windowThickness/2);

        G4Transform3D transformWin = transformSensor*G4Translate3D(translationWin);

        new G4PVPlacement(transformWin,
                          logicWindow,
                          "window",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // Solder
        new G4PVPlacement(transformWin,
                          logicSolder,
                          "solder",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // Stud
        G4ThreeVector translationStud = G4ThreeVector(0.,0.,0.);
        translationStud.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness - headerThickness - studThickness/2.);

        G4Transform3D transformStud = transformSensor*G4Translate3D(translationStud);

        new G4PVPlacement(transformStud,
                          logicStud,
                          "stud",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // EXTENDER
        G4ThreeVector translationExtender = G4ThreeVector(0.,0.,0.);
        translationExtender.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness - headerThickness - extenderThickness/2.);

        G4Transform3D transformExtender= transformSensor*G4Translate3D(translationExtender);

        new G4PVPlacement(transformExtender,
                          logicExtend,
                          "Extender",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        // ELECTRONIC BOX
        G4ThreeVector translationElecBox = G4ThreeVector(0.,0.,0.);
        translationElecBox.setZ(- sensorThickness/2. - botElectrodeHeight - substrateThickness - tecThickness - headerThickness - extenderThickness - elecBoxLength/2.);

        G4Transform3D transformElecBox = transformSensor*G4Translate3D(translationElecBox);

        new G4PVPlacement(transformElecBox,
                          logicElecBox,
                          "ElecBox",
                          fWorld_logic,
                          false,
                          sensorCpy,
                          checkOverlaps);

        sensorCpy++;
    }

    if (visBool)
    {
        G4VisAttributes* coverColor = new G4VisAttributes(G4Colour(1., 0., 1., 0.1));
        coverColor->SetForceSolid(true);
        logicCover->SetVisAttributes(coverColor);

        G4VisAttributes* extenderColor = new G4VisAttributes(G4Colour(1., 1., 1., 1));
        extenderColor->SetForceSolid(true);
        logicExtend->SetVisAttributes(extenderColor);
        logicElecBox->SetVisAttributes(extenderColor);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4SDManager.hh"
#include "SDDDetectorSD.hh"

void SDDDetectorConstruction::ConstructSDandField()
{
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    SDDDetectorSD* SD = new SDDDetectorSD("/SDsdd", this);
    SDman->AddNewDetector(SD);
    logicPixel->SetSensitiveDetector(SD);
}
