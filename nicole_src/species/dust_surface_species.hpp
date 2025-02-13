/**
 * @file dust_surface_species.hpp
 * 
 * Header file defining the DustSurfaceSpecies class. This class represents a species that exists on the surface
 * of a dust particle. It inherits from the DustLayerSpeciesBase class and is used to manage properties and behavior
 * specific to dust species on a surface (e.g., diffusion barriers, corresponding gas species, and dust mantle species).
 * 
 * @date 2025-02-12
 */


#ifndef DUST_SURFACE_SPECIES_HPP
#define DUST_SURFACE_SPECIES_HPP

#include "../nicole_defs.hpp"
#include "species.hpp"
#include "gas_species.hpp"
#include "dust_layer_species_base.hpp"
#include "dust_mantle_species.hpp"

namespace nicole
{
    // Forward declaration of related classes for circular dependencies
    class GasSpecies;
    class DustMantleSpecies;
    
    /**
     * @class DustSurfaceSpecies
     * 
     * Represents the surface species on a dust particle. This class is a specialized form of DustLayerSpeciesBase,
     * and it manages the diffusion barrier as well as relationships with corresponding gas species and dust mantle species.
     * It allows for accessing and setting the corresponding gas species and dust mantle species.
     */
    class DustSurfaceSpecies
        : public DustLayerSpeciesBase
    {
    private:

        // Member variables to hold references to the corresponding gas species and dust mantle species
        std::weak_ptr<GasSpecies> corresponding_gas_species_;  // Weak pointer to avoid circular reference
        std::shared_ptr<DustMantleSpecies> corresponding_dust_mantle_species_;  // Shared pointer to dust mantle species

        // Private method to calculate the diffusion barrier
        void CalculateDiffusionBarrier();

    public:

        /**
         * @brief Constructor for DustSurfaceSpecies.
         * Initializes a DustSurfaceSpecies object with the given parameters.
         * 
         * @param index The index of the dust surface species.
         * @param name The name of the dust surface species.
         * @param corresponding_gas_species The corresponding gas species related to this dust surface species.
         * @param binding_energy_on_H2O_ice Binding energy of the dust surface species on H2O ice.
         * @param binding_energy_on_silicate Binding energy of the dust surface species on silicate.
         */
        DustSurfaceSpecies(
            std::size_t index,
            const std::string& name,
            std::shared_ptr<GasSpecies> corresponding_gas_species,
            double binding_energy_on_H2O_ice,
            double binding_energy_on_silicate
        );

        // Getter functions for corresponding gas species and dust mantle species

        /**
         * @brief Get the corresponding gas species.
         * @return A shared pointer to the corresponding gas species.
         */
        inline std::shared_ptr<GasSpecies> GetCorrespondingGasSpecies() const {
            return corresponding_gas_species_.lock();
        }

        /**
         * @brief Get the corresponding dust mantle species.
         * @return A shared pointer to the corresponding dust mantle species.
         */
        inline const std::shared_ptr<DustMantleSpecies>& GetCorrespondingDustMantleSpecies() const {
            return corresponding_dust_mantle_species_;
        }

        // Setter functions for corresponding gas species and dust mantle species

        /**
         * @brief Set the corresponding gas species.
         * @param gas_species A shared pointer to the corresponding gas species.
         */
        inline void SetCorrespondingGasSpecies(const std::shared_ptr<GasSpecies>& gas_species) {
            corresponding_gas_species_= gas_species;
        }

        /**
         * @brief Set the corresponding dust mantle species.
         * @param dust_mantle_species A shared pointer to the corresponding dust mantle species.
         */
        inline void SetCorrespondingDustMantleSpecies(const std::shared_ptr<DustMantleSpecies>& dust_mantle_species) {
            corresponding_dust_mantle_species_ = dust_mantle_species;
        }
        
    };
}

#endif /* #include DUST_SURFACE_SPECIES_HPP */