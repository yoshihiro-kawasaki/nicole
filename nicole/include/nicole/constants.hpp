#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cmath>
#define _USE_MATH_DEFINES

#include "types.hpp"

namespace nicole {
    namespace constants {
        // Fundamental constants
        constexpr Real kSpeedOfLight            = 2.99792458e10;        // Speed of light in cm/s
        constexpr Real kGravitationalConstant   = 6.67408e-08;          // Gravitational constant in cm^3 g^-1 s^-2
        constexpr Real kPlanckConstant          = 6.6260755e-27;        // Planck constant in erg s
        constexpr Real kDiracConstant           = kPlanckConstant / (2.0 * M_PI); // Dirac constant (h-bar) in erg s
        constexpr Real kBoltzmannConstant       = 1.38064852e-16;       // Boltzmann constant in erg/K
        constexpr Real kStefanBoltzmannConstant = 5.6705e-5;            // Stefan-Boltzmann constant in erg cm^-2 s^-1 K^-4
        constexpr Real kAtomicMassUnit          = 1.660539067e-24;      // Atomic mass unit in grams
        constexpr Real kElectronMass            = 9.109383702e-28;      // Electron mass in grams
        constexpr Real kProtonMass              = 1.672621924e-24;      // Proton mass in grams
        constexpr Real kNeutronMass             = 1.674927498e-24;      // Neutron mass in grams
        constexpr Real kChargeUnit              = 4.803204673e-10;      // Elementary charge in statcoulombs
        constexpr Real kElectronVolt            = 1.60218e-12;          // Electron volt in erg
        constexpr Real kBohrRadius              = 5.2917720859e-9;      // Bohr radius in cm
        constexpr Real kAvogadroConstant        = 6.0221e23;            // Avogadro constant in mol^-1
        constexpr Real kGasConstantMol          = 8.3145e7;             // Gas constant in erg/mol K
        constexpr Real kRadiationConstant       = 7.5646e-15;           // Radiation constant in erg cm^-3 K^-4
        constexpr Real kAstronomicalUnit        = 1.495979e+13;         // Astronomical unit in cm
        constexpr Real kParsec                  = 3.085677e18;          // Parsec in cm
        constexpr Real kLightYear               = 9.460730473e17;       // Light year in cm
        constexpr Real kSolarYear               = 3.1556925e7;          // Solar year in seconds
        constexpr Real kSolarMass               = 1.9884e33;            // Solar mass in grams
        constexpr Real kSolarRadius             = 6.955080e10;          // Solar radius in cm
        constexpr Real kSolarLuminosity         = 3.828e33;             // Solar luminosity in erg/s
        constexpr Real kMeanMolecularWeight     = 2.34;                 // Mean molecular weight (for a typical gas)
        constexpr Real kGasMolecularMass        = kMeanMolecularWeight * kProtonMass; // Molecular mass of the gas in grams
        constexpr Real kH2CrossSection          = 2.0e-15;              // H2 cross-section in cm^2
        constexpr Real kInvMg                   = 1.0 / kGasMolecularMass; // Inverse of the gas molecular mass

        // Mathematical constants
        constexpr Real kSqrt2Pi                 = 2.5066282746310002;    // sqrt(2 * pi)
        constexpr Real kInv3Pi                  = 1.0 / (3.0 * M_PI);    // 1 / (3 * pi)
        constexpr Real kInv4Pi                  = 1.0 / (4.0 * M_PI);    // 1 / (4 * pi)
        constexpr Real kSqrt8Pi3                = 2.8944050182330705;    // sqrt(8 * pi / 3)
        constexpr Real kPi8                     = M_PI / 8.0;            // pi / 8
        constexpr Real kSqrt8Pi                 = 1.5957691216057308;    // sqrt(8 / pi)
        constexpr Real kSqrt2overPi             = 0.7978845608028654;    // sqrt(2 / pi)
        constexpr Real kInvSqrt2Pi              = 0.3989422804014327;    // 1 / sqrt(2 * pi)
    }
}

#endif /* CONSTANTS_HPP */
