#ifndef NICOLE_DEFS_HPP_
#define NICOLE_DEFS_HPP_

// C++ standard headers
#include <cmath>
#include <string>
#include <vector>

// Nicole code specific headers
#include "nicole/constants.hpp"
#include "nicole/types.hpp"

namespace nicole {
    // Nicole code-specific constants
    constexpr std::size_t kNotFoundElement = 9999;              // Element not found identifier
    constexpr std::size_t kNotFoundSpecies = 9999;              // Species not found identifier
    constexpr Real kDustSurfaceSitesDensity = 1.0e15;           // Dust surface sites density [cm^-2]
    constexpr Real kWidthOfBarrier = 1.0e-8;                    // Width of dust surface/mantle reaction activation energy barrier [cm]
    constexpr Real kPolarizabilityH = 0.667;                    // Polarizability of Hydrogen (H)
    constexpr Real kPolarizabilityH2 = 0.804;                   // Polarizability of Hydrogen (H2)
    constexpr Real kPolarizabilityHe = 0.207;                   // Polarizability of Helium (He)
    constexpr Real kBindingDiffusionRatioOfDustSurface = 0.4;   // Binding diffusion ratio of dust surface
    constexpr Real kBindingDiffusionRatioOfDustMantle = 0.8;    // Binding diffusion ratio of dust mantle
    constexpr Real kVisualExtinction = 1.0e1;                   // Visual extinction coefficient
    constexpr Real kNumberOfActiveSurfaceLayer = 4.0;           // Number of active surface layers on dust
    constexpr Real kStickingProbabilityElectron = 0.6;          // Sticking probability for electrons
    constexpr Real kStickingProbabilityIon = 1.0;               // Sticking probability for ions
    constexpr Real kMinimumRateCoefficient = 1.0e-99;           // Minimum rate coefficient for reactions
    constexpr Real kMinimumSpeciesAbundance = 1.0e-99;          // Minimum species abundance for calculations
    constexpr Real kOdepackMinimumAbsoluteTolerance = 1.0e-99;  // Minimum absolute tolerance for ODE solver


    // Prefixes for dust surface and mantle species
    const std::string kDustSurfaceSpeciesPrefix = "s";           // Prefix for dust surface species
    const std::string kDustMantleSpeciesPrefix  = "m";           // Prefix for dust mantle species
    const std::string kDustSpeciesPrefix        = "D";           // Prefix for dust species


    // Specific species names
    const std::string kElectron = "e-";  // Electron
    const std::string kH2 = "H2";        // Hydrogen molecule (H2)
    const std::string kH  = "H";         // Hydrogen atom (H)
    const std::string kHe = "He";        // Helium atom (He)


    // Special species names
    const std::string kEmptyString = "";  // Empty string (used for certain special cases)
    const std::string kAsterisk = "*";    // Asterisk symbol (wildcard)
    const std::string kCR = "CR";         // Cosmic ray
    const std::string kCRP = "CRP";       // Cosmic ray proton
    const std::string kPhoton = "Photon"; // Photon
    const std::vector<std::string> SPECIAL_SPECIES_LIST = {
        kEmptyString, 
        kAsterisk,
        kCR, 
        kCRP, 
        kPhoton
    };


    /**
     * @struct EnvironmentParameters
     * @brief Structure to hold environmental parameters such as gas density and temperature.
     */
    struct EnvironmentParameters {
        Real gas_number_density;          // gas (hydrogen nuclei) number density [cm^-3]
        Real gas_mass_density;
        Real gas_temperature;             // gas temperature [K]
        Real cosmic_ray_ionization_rate;  // Cosmic ray ionization rate [s^-1]
        Real x_rays_ionization_rate;      // X-ray ionization rate [s^-1]
        Real visual_extinction;           // Visual extinction coefficient
        Real scaling_factor_uv_field;     // Scaling factor for UV field
        Real magnetic_field;              // Magnetic field strength [Gauss]
    };


    /**
     * @namespace reaction_type_id
     * @brief Namespace to define IDs for different reaction types.
     */
    namespace reaction_type_id {
        constexpr int kReactionDummy = 0; // Dummy reaction ID

        // Gas-phase reaction IDs
        constexpr int kGasPhase1  = 1;
        constexpr int kGasPhase2  = 2;
        constexpr int kGasPhase3  = 3;
        constexpr int kGasPhase4  = 4;
        constexpr int kGasPhase5  = 5;
        constexpr int kGasPhase6  = 6;
        constexpr int kGasPhase7  = 7;
        constexpr int kGasPhase8  = 8;
        constexpr int kGasPhase9  = 9;
        constexpr int kGasPhase10 = 10;
        constexpr int kGasPhase11 = 11;
        constexpr int kNumberOfGasPhase = 12;
        constexpr int kGasPhaseStart    = kGasPhase1;
        constexpr int kGasPhaseEnd      = kGasPhase11;

        // Dust-related reaction IDs
        constexpr int kDustAndChargedGasParticleCollison   = 15;
        constexpr int kDustCollision                       = 16;
        constexpr int kAccretionGasParticleOnDustSurfaces  = 17;
        constexpr int kThermalDesorptionOnDustSurfaces     = 18;
        constexpr int kCosmicRayDesorptionOnDustSurfaces   = 19;
        constexpr int kPhotoDesorptionByExternalUV         = 20;
        constexpr int kPhotoDesorptionByCRGeneratedUV      = 21;
        constexpr int kPhotoDissociationByCROnDustSurfaces = 22;
        constexpr int kPhotoDissociationByUVOnDustSurfaces = 23;
        constexpr int kDustSurfaceReaction                 = 24;
        constexpr int kDustMantleReaction                  = 25;
        constexpr int kDustSurfaceToMantleSwapping         = 26;
        constexpr int kDustMantleToSurfaceSwapping         = 27;

        constexpr int kNumberOfTypeID = 28;
    }


    /**
     * @fn SQR
     * @brief Computes the square of a number.
     * @param x The input value.
     * @return The square of the input value.
     */
    template <typename T>
#if __cplusplus >= 201402L  // C++14 or later uses constexpr
    constexpr 
#else                       // C++11 uses inline
    inline 
    #endif
    T SQR(T x) {
        return x * x;
    }


    /**
     * @fn CUB
     * @brief Computes the cube of a number.
     * @param x The input value.
     * @return The cube of the input value.
     */
    template <typename T>
#if __cplusplus >= 201402L  // C++14 or later uses constexpr
    constexpr 
#else                       // C++11 uses inline
    inline 
#endif
    T CUB(T x) {
        return x * x * x;
    }

    
    /**
     * @fn SIGN
     * @brief Returns the sign of a number.
     * @param x The input value.
     * @return 1 if x >= 0, -1 if x < 0.
     */
    template <typename T>
#if __cplusplus >= 201402L  // C++14 or later uses constexpr
    constexpr 
#else                       // C++11 uses inline
    inline 
#endif
    T SIGN(T x) {
        return (x >= static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(-1));
    }

}

#endif /* NICOLE_DEFS_HPP_ */
