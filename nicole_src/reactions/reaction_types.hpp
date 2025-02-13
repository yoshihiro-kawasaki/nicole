/**
 * @file reaction_types.hpp
 * @brief Defines parameter indices for different reaction types in the simulation.
 * 
 * This header file includes parameter index constants for various reaction types,
 * such as gas phase reactions, dust collision, accretion, desorption processes, etc.
 * 
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#ifndef REACTION_TYPES_HPP
#define REACTION_TYPES_HPP

#include "../nicole_defs.hpp"

namespace nicole
{
    /**
     * @namespace gas_phase_reaction_params
     * @brief Parameter indices for gas-phase reactions
     * 
     * This namespace contains constant indices for the parameters associated with 
     * gas-phase reactions, including rate coefficients and temperature limits.
     */
    namespace gas_phase_reaction_params 
    {
        constexpr std::size_t kAlpha = 0;   // Rate coefficient alpha
        constexpr std::size_t kBeta  = 1;   // Rate coefficient beta
        constexpr std::size_t kGamma = 2;   // Rate coefficient gamma
        constexpr std::size_t kTemperatureLowerLimit = 3;  // Lower temperature limit for reaction validity
        constexpr std::size_t kTemperatureUpperLimit = 4;  // Upper temperature limit for reaction validity
        constexpr std::size_t kFormulaID = 5; // Formula ID for the reaction
        constexpr std::size_t kID = 6;      // Reaction ID
        constexpr std::size_t kNumParams = 7; // Total number of parameters
    };

    /**
     * @namespace dust_and_charged_particle_collison_params
     * @brief Parameter indices for dust and charged particle collision
     * 
     * This namespace contains constant indices for the parameters related to 
     * the collision between dust particles and charged particles, such as 
     * gas mass, gas charge, dust charge, dust radius, and sticking probability.
     */
    namespace dust_and_charged_particle_collison_params 
    {
        constexpr std::size_t kGasMass = 0;     // Gas mass involved in the collision
        constexpr std::size_t kGasCharge = 1;   // Charge of the gas particle
        constexpr std::size_t kDustCharge = 2;  // Charge of the dust particle
        constexpr std::size_t kDustRadius = 3;  // Radius of the dust particle
        constexpr std::size_t kDustCrossSection = 4; // Cross-section of the dust particle
        constexpr std::size_t kStickingProbability = 5; // Probability of sticking after collision
        constexpr std::size_t kNumParams = 6;   // Total number of parameters
    };

    /**
     * @namespace dust_collision_params
     * @brief Parameter indices for dust-dust collisions
     * 
     * This namespace contains constant indices for the parameters related to 
     * collisions between dust particles, including radii, charges, and masses.
     */
    namespace dust_collision_params {
        constexpr std::size_t kDustRadius1 = 0; // Radius of dust particle 1
        constexpr std::size_t kDustRadius2 = 1; // Radius of dust particle 2
        constexpr std::size_t kDustCharge1 = 2; // Charge of dust particle 1
        constexpr std::size_t kDustCharge2 = 3; // Charge of dust particle 2
        constexpr std::size_t kDustMass1 = 4;   // Mass of dust particle 1
        constexpr std::size_t kDustMass2 = 5;   // Mass of dust particle 2
        constexpr std::size_t kNumParams = 6;   // Total number of parameters
    }

    /**
     * @namespace accretion_on_dusts_params
     * @brief Parameter index for accretion of gas particle on dust surfaces processes
     * 
     * This namespace contains the parameter for accretion of gas particle on dust surfaces processes.
     */
    namespace accretion_on_dusts_params 
    {
        constexpr std::size_t kNumParams = 1;   // Total number of parameters
    };

    /**
     * @namespace thermal_desorption_params
     * @brief Parameter indices for thermal desorption processes
     * 
     * This namespace contains constant indices for thermal desorption processes, 
     * including vibration frequency and binding energy on H2O ice.
     */
    namespace thermal_desorption_params 
    {
        constexpr std::size_t kVibrationFrequency = 0;  // Vibration frequency for thermal desorption
        constexpr std::size_t kBindingEnergyOnH2Oice = 1; // Binding energy on H2O ice for desorption
        constexpr std::size_t kNumParams = 2;  // Total number of parameters
    };

    /**
     * @namespace cosmic_ray_desorption_on_dusts_surfaces_params
     * @brief Parameter indices for cosmic ray desorption on dust surfaces
     * 
     * This namespace defines parameters for the desorption process caused by cosmic rays 
     * on dust surfaces, including vibration frequency and binding energy.
     */
    namespace cosmic_ray_desorption_on_dusts_surfaces_params 
    {
        constexpr std::size_t kVibrationFrequency = 0;  // Vibration frequency for cosmic ray desorption
        constexpr std::size_t kBindingEnergyOnH2Oice = 1; // Binding energy on H2O ice for desorption
        constexpr std::size_t kNumPramas = 2;  // Total number of parameters
    }

    /**
     * @namespace photo_desorption_by_external_UV_params
     * @brief Parameter index for photo desorption by external UV radiation
     * 
     * This namespace defines a parameter for photo desorption processes triggered 
     * by external UV radiation.
     */
    namespace photo_desorption_by_external_UV_params 
    {
        constexpr std::size_t kNumParams = 1;   // Total number of parameters
    };

    /**
     * @namespace photo_desorption_by_CR_generated_UV_params
     * @brief Parameter index for photo desorption by UV radiation generated by cosmic rays
     * 
     * This namespace defines a parameter for photo desorption processes triggered 
     * by UV radiation generated by cosmic rays.
     */
    namespace photo_desorption_by_CR_generated_UV_params 
    {
        constexpr std::size_t kNumParams = 1;   // Total number of parameters
    };

    /**
     * @namespace DustSurfaceReactionParams
     * @brief Parameter indices for dust surface reactions
     * 
     * This namespace contains constant indices for parameters in dust surface reactions, 
     * including vibration frequencies, diffusion barriers, activation energies, and probabilities.
     */
    namespace dust_surface_reaction_params 
    {
        constexpr std::size_t kVibrationFrequency1 = 0;  // Vibration frequency of species 1
        constexpr std::size_t kVibrationFrequency2 = 1;  // Vibration frequency of species 2
        constexpr std::size_t kDiffusionBarrier1 = 2;    // Diffusion barrier for species 1
        constexpr std::size_t kDiffusionBarrier2 = 3;    // Diffusion barrier for species 2
        constexpr std::size_t kActivationEnergy = 4;     // Activation energy for the reaction
        constexpr std::size_t kTunnelingEffectProbability = 5;  // Probability of tunneling effect
        constexpr std::size_t kChemicalDesorptionOnH2OProbability = 6;  // Probability of chemical desorption on H2O ice
        constexpr std::size_t kChemicalDesorptionOnSilicateProbability = 7;  // Probability of chemical desorption on silicates
        constexpr std::size_t kNumParams = 8;  // Total number of parameters
    };

    /**
     * @namespace dust_mantle_reation_params
     * @brief Parameter indices for dust mantle reactions
     * 
     * This namespace contains constant indices for parameters in dust mantle reactions,
     * similar to dust surface reactions, but for the mantle.
     */
    namespace dust_mantle_reaction_params 
    {
        constexpr std::size_t kVibrationFrequency1 = 0;  // Vibration frequency of species 1
        constexpr std::size_t kVibrationFrequency2 = 1;  // Vibration frequency of species 2
        constexpr std::size_t kDiffusionBarrier1 = 2;    // Diffusion barrier for species 1
        constexpr std::size_t kDiffusionBarrier2 = 3;    // Diffusion barrier for species 2
        constexpr std::size_t kActivationEnergy = 4;     // Activation energy for the reaction
        constexpr std::size_t kTunnelingEffectProbability = 5;  // Probability of tunneling effect
        constexpr std::size_t kNumParams = 6;  // Total number of parameters
    };

    /**
     * @namespace photo_dissociation_by_CR_on_dusts_params
     * @brief Parameter indices for photo-dissociation by cosmic rays on dust particles
     * 
     * This namespace defines parameters related to photo-dissociation caused by cosmic rays on dust particles.
     */
    namespace photo_dissociation_by_CR_on_dusts_params 
    {
        constexpr std::size_t kAlpha = 0;  // Rate coefficient alpha for photo-dissociation
        constexpr std::size_t kBeta = 1;   // Rate coefficient beta for photo-dissociation
        constexpr std::size_t kGamma = 2;  // Rate coefficient gamma for photo-dissociation
        constexpr std::size_t kNumParams = 3;  // Total number of parameters
    };

    /**
     * @namespace photo_dissociation_by_UV_on_dusts_params
     * @brief Parameter indices for photo-dissociation by UV radiation on dust particles
     * 
     * This namespace defines parameters related to photo-dissociation caused by UV radiation on dust.
     */
    namespace photo_dissociation_by_UV_on_dusts_params 
    {
        constexpr std::size_t kAlpha = 0;  // Rate coefficient alpha for photo-dissociation
        constexpr std::size_t kBeta = 1;   // Rate coefficient beta for photo-dissociation
        constexpr std::size_t kGamma = 2;  // Rate coefficient gamma for photo-dissociation
        constexpr std::size_t kNumParams = 3;  // Total number of parameters
    };

    /**
     * @namespace dust_surface_to_mantle_params
     * @brief Parameter indices for dust surface to mantle processes
     * 
     * This namespace contains parameters for processes where species move from the dust surface 
     * to the dust mantle, such as vibration frequency and binding energy.
     */
    namespace dust_surface_to_mantle_params 
    {
        constexpr std::size_t kVibrationFrequency = 0;  // Vibration frequency for the transition
        constexpr std::size_t kBindingEnergyOnH2Oice = 1;  // Binding energy on H2O ice
        constexpr std::size_t kNumPramas = 2;  // Total number of parameters
    };

    /**
     * @namespace dust_mantle_to_surface_params
     * @brief Parameter indices for dust mantle to surface processes
     * 
     * This namespace contains parameters for processes where species move from the dust mantle 
     * to the dust surface, including vibration frequency and binding energy.
     */
    namespace dust_mantle_to_surface_params 
    {
        constexpr std::size_t kVibrationFrequency = 0;  // Vibration frequency for the transition
        constexpr std::size_t kBindingEnergyOnH2Oice = 1;  // Binding energy on H2O ice
        constexpr std::size_t kNumPramas = 2;  // Total number of parameters
    };
}

#endif /* REACTION_TYPES_HPP */