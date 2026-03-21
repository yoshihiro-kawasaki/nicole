/**
 * @file dust_mantle_species.hpp
 * @brief This file defines the DustMantleSpecies class, representing chemical species embedded in the dust mantle.
 */

#ifndef DUST_MANTLE_SPECIES_HPP
#define DUST_MANTLE_SPECIES_HPP

#include "species.hpp"
#include "gas_species.hpp"
#include "dust_surface_species.hpp"
#include "dust_layer_species_base.hpp"

namespace nicole
{
    // Forward declarations to avoid circular dependencies
    class GasSpecies;
    class DustSurfaceSpecies;

    /**
     * @class DustMantleSpecies
     * @brief Represents a species embedded in the dust mantle.
     * 
     * This class is derived from DustLayerSpeciesBase and represents a species that resides
     * within the ice mantle of dust grains. It maintains references to the corresponding 
     * gas-phase species and dust surface species.
     */
    class DustMantleSpecies
        : public DustLayerSpeciesBase
    {
    private:
        /// Weak pointer to the corresponding gas-phase species
        std::weak_ptr<GasSpecies> corresponding_gas_species_;
        
        /// Weak pointer to the corresponding dust surface species
        std::weak_ptr<DustSurfaceSpecies> corresponding_dust_surface_species_;

        /**
         * @brief Calculates the diffusion barrier for the dust mantle species.
         */
        void CalculateDiffusionBarrier();

    public:

        /**
         * @brief Constructor for DustMantleSpecies
         * 
         * @param index The index of the species.
         * @param name The name of the species.
         * @param corresponding_gas_species Shared pointer to the corresponding gas-phase species.
         * @param corresponding_dust_surface_species Shared pointer to the corresponding dust surface species.
         */
        DustMantleSpecies(
            std::size_t index,
            const std::string& name,
            std::shared_ptr<GasSpecies> corresponding_gas_species,
            std::shared_ptr<DustSurfaceSpecies> corresponding_dust_surface_species
        );

        ~DustMantleSpecies() { }

        /**
         * @brief Get the corresponding gas-phase species.
         * 
         * @return std::shared_ptr<GasSpecies> Shared pointer to the gas-phase species.
         */
        inline std::shared_ptr<GasSpecies> GetCorrespondingGasSpecies() const {
            return corresponding_gas_species_.lock();
        }

        /**
         * @brief Get the corresponding dust surface species.
         * 
         * @return std::shared_ptr<DustSurfaceSpecies> Shared pointer to the dust surface species.
         */
        inline std::shared_ptr<DustSurfaceSpecies> GetCorrespondingDustSurfaceSpecies() const {
            return corresponding_dust_surface_species_.lock();
        }

        /**
         * @brief Set the corresponding gas-phase species.
         * 
         * @param gas_species Shared pointer to the gas-phase species.
         */
        inline void SetCorrespondingGasSpecies(const std::shared_ptr<GasSpecies>& gas_species) {
            corresponding_gas_species_ = gas_species;
        }

        /**
         * @brief Set the corresponding dust surface species.
         * 
         * @param dust_surface_species Shared pointer to the dust surface species.
         */
        inline void SetCorrespondingDustSurfaceSpecies(const std::shared_ptr<DustSurfaceSpecies>& dust_surface_species) {
            corresponding_dust_surface_species_ = dust_surface_species;
        }
    };
}

#endif /* DUST_MANTLE_SPECIES_HPP */
