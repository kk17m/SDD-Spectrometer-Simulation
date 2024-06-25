#ifndef SDDdigitize_HH
#define SDDdigitize_HH

#include "G4VDigitizerModule.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "boost/tuple/tuple.hpp"

class SDDdigitize : public G4VDigitizerModule
{
public:
    SDDdigitize(G4String name);
    ~SDDdigitize();

    void Digitize();

private:
    G4double sensorPz;
    G4double nEH;
    G4double fano;
    G4double DelE;

private:

    // Analytical model
    G4double NoiseFano(G4double eCharge);
    G4double NoiseGaussian(G4double ene);
};

#endif
