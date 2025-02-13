/**
 * @file dust_surface_species.cpp
 * 
 * Implementation of the DustSurfaceSpecies class. This file defines the methods of the DustSurfaceSpecies class,
 * including the constructor and the methods related to diffusion barrier calculations and vibration frequency calculations.
 */

#include "dust_surface_species.hpp"

namespace nicole
{
    /**
     * @brief Constructor for DustSurfaceSpecies.
     * 
     * Initializes a DustSurfaceSpecies object by setting up the index, name, charge, element composition, and
     * associated gas species properties. It also calculates the vibration frequency and the diffusion barrier for the species.
     * 
     * @param index The index of the dust surface species.
     * @param name The name of the dust surface species.
     * @param corresponding_gas_species The corresponding gas species related to this dust surface species.
     * @param binding_energy_on_H2O_ice Binding energy of the dust surface species on H2O ice.
     * @param binding_energy_on_silicate Binding energy of the dust surface species on silicate.
     */
    DustSurfaceSpecies::DustSurfaceSpecies(
        std::size_t index,
        const std::string& name,
        std::shared_ptr<GasSpecies> corresponding_gas_species,
        double binding_energy_on_H2O_ice,
        double binding_energy_on_silicate
    ) : DustLayerSpeciesBase( // Call to base class constructor (DustLayerSpeciesBase)
            index,
            name, 
            corresponding_gas_species->GetCharge(), 
            corresponding_gas_species->GetElementComposition(), 
            corresponding_gas_species->GetPtrElementManager(),
            SpeciesType::Surface,
            binding_energy_on_H2O_ice,
            binding_energy_on_silicate
        ),
        corresponding_gas_species_(corresponding_gas_species),
        corresponding_dust_mantle_species_(nullptr) // Dust mantle species is initialized to nullptr
    { 
        // Calculate vibration frequency for the dust surface species
        CalculateVibrationFrequency();

        // Calculate the diffusion barrier for the dust surface species based on the given binding energies
        CalculateDiffusionBarrier();
    }

    /**
     * @brief Calculate the diffusion barrier for the dust surface species.
     * 
     * This method calculates the diffusion barrier based on the binding energies of the species
     * on H2O ice and silicate. The diffusion barrier is used to estimate the activation energy required
     * for species diffusion on dust surfaces.
     * 
     * @ref 
     * Ruaud et al. 2016
     */
    void DustSurfaceSpecies::CalculateDiffusionBarrier() 
    {
        // If the species name is "sH", set diffusion barrier to specific values.
        if (name_ == "sH") {
            // Assigning diffusion barrier values for "sH" species on H2O ice and silicate
            diffusion_barrier_on_H2O_ice_ = 230.0;
            diffusion_barrier_on_silicate_ = 230.0;
        } else {
            // For other species, the diffusion barrier is calculated as a ratio of the binding energy
            diffusion_barrier_on_H2O_ice_ = kBindingDiffusionRatioOfDustSurface * binding_energy_on_H2O_ice_;
            diffusion_barrier_on_silicate_ = kBindingDiffusionRatioOfDustSurface * binding_energy_on_silicate_;
        }
    }

}