#ifndef DetectorHit_HH
#define DetectorHit_HH

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
