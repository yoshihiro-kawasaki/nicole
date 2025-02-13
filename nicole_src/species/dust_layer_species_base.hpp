/**
 * @file dust_layer_species_base.hpp
 * @brief This file defines the DustLayerSpeciesBase class, which serves as the base class 
 * for species that are involved in dust layers, such as dust surface species and dust mantle species.
 * 
 * It includes properties and methods for calculating vibration frequencies and diffusion barriers
 * specific to different materials (e.g., H2O ice and silicate).
 * 
 * The class is a derived class of Species and provides functionality specific to dust layer species.
 */

#ifndef DUST_LAYER_SPECIES_BASE_HPP
#define DUST_LAYER_SPECIES_BASE_HPP

#include "../nicole_defs.hpp"
#include "species.hpp"

namespace nicole
{
    /**
     * @class DustLayerSpeciesBase
     * @brief Base class for species involved in dust layers (e.g., surface or mantle of dust).
     * 
     * This class extends the Species class and provides specific functionality for dust layer species,
     * including vibration frequencies and diffusion barriers on different materials (H2O ice and silicate).
     * It also has methods to calculate these properties.
     */
    class DustLayerSpeciesBase : public Species
    {
    protected:
        // Vibration frequencies for the species on H2O ice and silicate surfaces
        double vibration_frequency_on_H2O_ice_;
        double vibration_frequency_on_silicate_;

        // Diffusion barriers for the species on H2O ice and silicate surfaces
        double diffusion_barrier_on_H2O_ice_;
        double diffusion_barrier_on_silicate_;

        /**
         * @brief Calculates the vibration frequency of the species.
         * This is a placeholder for the actual implementation.
         * The function might be specialized in derived classes.
         */
        void CalculateVibrationFrequency();

        /**
         * @brief Calculates the diffusion barrier for the species.
         * This is a virtual function that can be overridden by derived classes.
         * The default behavior does nothing, but the derived class may provide its own implementation.
         */
        virtual void CalculateDiffusionBarrier() { };

    public:
        /**
         * @brief Constructor for DustLayerSpeciesBase
         * 
         * Initializes a dust layer species with the given properties.
         * 
         * @param index The index of the species.
         * @param name The name of the species.
         * @param charge The charge of the species.
         * @param element_composition The composition of elements in the species.
         * @param ptr_element_manager A pointer to the ElementManager managing elements.
         * @param type The type of the species.
         * @param binding_energy_on_H2O_ice Binding energy of the species on H2O ice.
         * @param binding_energy_on_silicate Binding energy of the species on silicate.
         */
        DustLayerSpeciesBase(
            std::size_t index,
            const std::string& name, 
            int charge, 
            const std::vector<std::size_t> element_composition, 
            ElementManager* ptr_element_manager,
            SpeciesType type,
            double binding_energy_on_H2O_ice,
            double binding_energy_on_silicate
        );

        // Getters for the vibration frequency and diffusion barrier on H2O ice and silicate
        inline double GetVibrationFrequencyOnH2Oice() const { return vibration_frequency_on_H2O_ice_; }
        inline double GetVibrationFrequencyOnSilicate() const { return vibration_frequency_on_silicate_; }
        inline double GetDiffusionBarrierOnH2Oice() const { return diffusion_barrier_on_H2O_ice_; }
        inline double GetDiffusionBarrierOnSilicate() const { return diffusion_barrier_on_silicate_; }

        // Setters for the vibration frequency and diffusion barrier on H2O ice and silicate
        inline void SetVibrationFrequencyOnH2Oice(const double vibration_frequency_on_H2O_ice) { 
            vibration_frequency_on_H2O_ice_ = vibration_frequency_on_H2O_ice;
        }
        inline void SetVibrationFrequencyOnSilicate(const double vibration_frequency_on_silicate) {
            vibration_frequency_on_silicate_ = vibration_frequency_on_silicate;
        }
        inline void SetDiffusionBarrierOnH2Oice(const double diffusion_barrier_on_H2O_ice) {
            diffusion_barrier_on_H2O_ice_ = diffusion_barrier_on_H2O_ice;
        }
        inline void SetDiffusionBarrierOnSilicate(const double diffusion_barrier_on_silicate) {
            diffusion_barrier_on_silicate_ = diffusion_barrier_on_silicate;
        }
    };
}

#endif /* DUST_LAYER_SPECIES_BASE_HPP */