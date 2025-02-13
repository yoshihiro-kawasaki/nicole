/**
 * @file dust_layer_species_base.cpp
 * @brief This file contains the implementation of the DustLayerSpeciesBase class, 
 * which is responsible for managing species involved in dust layers.
 */

#include "dust_layer_species_base.hpp"

namespace nicole
{
    /**
     * @brief Constructor for DustLayerSpeciesBase
     * 
     * Initializes a DustLayerSpeciesBase object with the specified properties such as 
     * index, name, charge, element composition, element manager, species type, 
     * binding energies on H2O ice and silicate, and calls the base class (Species) constructor.
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
    DustLayerSpeciesBase::DustLayerSpeciesBase(
        std::size_t index,
        const std::string& name, 
        int charge, 
        const std::vector<std::size_t> element_composition, 
        ElementManager* ptr_element_manager,
        SpeciesType type,
        double binding_energy_on_H2O_ice,
        double binding_energy_on_silicate
    ) : Species(
            index,
            name, 
            charge, 
            element_composition, 
            ptr_element_manager, 
            type,
            binding_energy_on_H2O_ice,
            binding_energy_on_silicate
        )
    { }

    /**
     * @brief Calculates the vibration frequency of the species on H2O ice and silicate surfaces.
     * 
     * This method calculates the vibration frequency using the following formula:
     * 
     * frequency = sqrt( (2 * density * kB * binding_energy) / (π^2 * mass) )
     * 
     * Where:
     * - density is a constant value representing the dust surface site density,
     * - kB is the Boltzmann constant,
     * - binding_energy is the energy of the species' binding to the surface (H2O ice or silicate),
     * - mass is the mass of the species.
     */
    void DustLayerSpeciesBase::CalculateVibrationFrequency() 
    {
        // Check if the mass is non-zero to avoid division by zero
        if (mass_ != 0.0) {
            // Calculate vibration frequency for H2O ice surface
            vibration_frequency_on_H2O_ice_ = std::sqrt(2.0 * kDustSurfaceSitesDensity 
                * constants::kBoltzmannConstant * binding_energy_on_H2O_ice_ 
                / (M_PI * M_PI * mass_));

            // Calculate vibration frequency for silicate surface
            vibration_frequency_on_silicate_ = std::sqrt(2.0 * kDustSurfaceSitesDensity
                * constants::kBoltzmannConstant * binding_energy_on_silicate_
                / (M_PI * M_PI * mass_));
        }
    }
}
