#ifndef DetectorSD_HH
#define DetectorSD_HH

#include "G4VSensitiveDetector.hh"
#include "SDDDetectorConstruction.hh"
#include "SDDDetectorHit.hh"

class SDDDetectorConstruction;

class SDDDetectorSD : public G4VSensitiveDetector
{

public:
    SDDDetectorSD(G4String, SDDDetectorConstruction*detCon);
    ~SDDDetectorSD();

    void Initialize(G4HCofThisEvent*HCE);
    void EndOfEvent(G4HCofThisEvent*);

protected:
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);

private:
    DetectorHitsCollection *HitsColl;
    SDDDetectorConstruction* detConstruction;
    G4bool setUnitTrackWeight;
};

#endif
