/**
 * @file dust_mantle_species.cpp
 * @brief Implementation of the DustMantleSpecies class.
 */

#include "dust_mantle_species.hpp"

namespace nicole
{
    /**
     * @brief Constructor for DustMantleSpecies.
     * 
     * Initializes a dust mantle species with its corresponding gas-phase 
     * and dust surface species. It also calculates the vibration frequency 
     * and diffusion barrier during initialization.
     * 
     * @param index The index of the species.
     * @param name The name of the species.
     * @param corresponding_gas_species Shared pointer to the corresponding gas-phase species.
     * @param corresponding_dust_surface_species Shared pointer to the corresponding dust surface species.
     */
    DustMantleSpecies::DustMantleSpecies(
        std::size_t index,
        const std::string& name,
        std::shared_ptr<GasSpecies> corresponding_gas_species,
        std::shared_ptr<DustSurfaceSpecies> corresponding_dust_surface_species
    ) : DustLayerSpeciesBase(
            index,
            name,
            corresponding_gas_species->GetCharge(),                      // Use charge from gas species
            corresponding_gas_species->GetElementComposition(),          // Use element composition from gas species
            corresponding_gas_species->GetPtrElementManager(),           // Use element manager from gas species
            SpeciesType::Mantle,                                         // Set species type as Mantle
            corresponding_dust_surface_species->GetBindingEnergyOnH2Oice(),  // Use binding energy on H2O ice from dust surface species
            corresponding_dust_surface_species->GetBindingEnergyOnBareSilicate() // Use binding energy on silicate from dust surface species
        ),
        corresponding_gas_species_(corresponding_gas_species),
        corresponding_dust_surface_species_(corresponding_dust_surface_species)
    { 
        // Calculate vibration frequency and diffusion barrier upon creation
        CalculateVibrationFrequency();
        CalculateDiffusionBarrier();
    }

    /**
     * @brief Calculates the diffusion barrier for species in the dust mantle.
     * 
     * The diffusion barrier is computed as a fraction of the binding energy.
     * The fraction is determined by the constant kBindingDiffusionRatioOfDustMantle.
     */
    void DustMantleSpecies::CalculateDiffusionBarrier() 
    {
        diffusion_barrier_on_H2O_ice_ = kBindingDiffusionRatioOfDustMantle * binding_energy_on_H2O_ice_;
        diffusion_barrier_on_silicate_ = kBindingDiffusionRatioOfDustMantle * binding_energy_on_silicate_;
    }
}
