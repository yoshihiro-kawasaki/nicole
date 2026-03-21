/**
 * @file constants.hpp
 * @brief A header file containing physical constants for use in astrophysical calculations.
 * 
 * This file defines a collection of important physical constants in various domains, such as:
 * - Fundamental constants (e.g., speed of light, gravitational constant)
 * - Atomic and molecular constants (e.g., electron mass, atomic mass unit)
 * - Astronomical constants (e.g., astronomical unit, parsec, light year)
 * - Mathematical constants (e.g., Pi, square roots of Pi)
 * 
 * These constants are defined as `constexpr` to allow for compile-time evaluation.
 * The constants are grouped in the `constants` namespace for better organization and to avoid naming conflicts.
 * 
 * @date 2025-02-11
 * @author Y. Kawasaki
 */
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cmath>
#define _USE_MATH_DEFINES

namespace constants 
{
    // Fundamental constants
    constexpr double kSpeedOfLight            = 2.99792458e10;        // Speed of light in cm/s
    constexpr double kGravitationalConstant   = 6.67408e-08;          // Gravitational constant in cm^3 g^-1 s^-2
    constexpr double kPlanckConstant          = 6.6260755e-27;        // Planck constant in erg s
    constexpr double kDiracConstant           = kPlanckConstant / (2.0 * M_PI); // Dirac constant (h-bar) in erg s
    constexpr double kBoltzmannConstant       = 1.38064852e-16;       // Boltzmann constant in erg/K
    constexpr double kStefanBoltzmannConstant = 5.6705e-5;            // Stefan-Boltzmann constant in erg cm^-2 s^-1 K^-4
    constexpr double kAtomicMassUnit          = 1.660539067e-24;      // Atomic mass unit in grams
    constexpr double kElectronMass            = 9.109383702e-28;      // Electron mass in grams
    constexpr double kProtonMass              = 1.672621924e-24;      // Proton mass in grams
    constexpr double kNeutronMass             = 1.674927498e-24;      // Neutron mass in grams
    constexpr double kChargeUnit              = 4.803204673e-10;      // Elementary charge in statcoulombs
    constexpr double kElectronVolt            = 1.60218e-12;          // Electron volt in erg
    constexpr double kBohrRadius              = 5.2917720859e-9;      // Bohr radius in cm
    constexpr double kAvogadroConstant        = 6.0221e23;            // Avogadro constant in mol^-1
    constexpr double kGasConstantMol          = 8.3145e7;             // Gas constant in erg/mol K
    constexpr double kRadiationConstant       = 7.5646e-15;           // Radiation constant in erg cm^-3 K^-4
    constexpr double kAstronomicalUnit        = 1.495979e+13;         // Astronomical unit in cm
    constexpr double kParsec                  = 3.085677e18;          // Parsec in cm
    constexpr double kLightYear               = 9.460730473e17;       // Light year in cm
    constexpr double kSolarYear               = 3.1556925e7;          // Solar year in seconds
    constexpr double kSolarMass               = 1.9884e33;            // Solar mass in grams
    constexpr double kSolarRadius             = 6.955080e10;          // Solar radius in cm
    constexpr double kSolarLuminosity         = 3.828e33;             // Solar luminosity in erg/s
    constexpr double kMeanMolecularWeight     = 2.34;                 // Mean molecular weight (for a typical gas)
    constexpr double kGasMolecularMass        = kMeanMolecularWeight * kProtonMass; // Molecular mass of the gas in grams
    constexpr double kH2CrossSection          = 2.0e-15;              // H2 cross-section in cm^2
    constexpr double kInvMg                   = 1.0 / kGasMolecularMass; // Inverse of the gas molecular mass

    // Mathematical constants
    constexpr double kSqrt2Pi                 = 2.5066282746310002;    // sqrt(2 * pi)
    constexpr double kInv3Pi                  = 1.0 / (3.0 * M_PI);    // 1 / (3 * pi)
    constexpr double kInv4Pi                  = 1.0 / (4.0 * M_PI);    // 1 / (4 * pi)
    constexpr double kSqrt8Pi3                = 2.8944050182330705;    // sqrt(8 * pi / 3)
    constexpr double kPi8                     = M_PI / 8.0;            // pi / 8
    constexpr double kSqrt8Pi                 = 1.5957691216057308;    // sqrt(8 / pi)
    constexpr double kSqrt2overPi             = 0.7978845608028654;    // sqrt(2 / pi)
    constexpr double kInvSqrt2Pi              = 0.3989422804014327;    // 1 / sqrt(2 * pi) 
}

#endif /* CONSTANTS_HPP */