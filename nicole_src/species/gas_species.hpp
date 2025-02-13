/**
 * @file gas_species.hpp
 * @brief Definition of GasSpecies class, which represents gas-phase chemical species.
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#ifndef GAS_SPECIES_HPP
#define GAS_SPECIES_HPP

#include "../nicole_defs.hpp"
#include "species.hpp"
#include "dust_surface_species.hpp"
#include "dust_mantle_species.hpp"

namespace nicole
{
    // Forward declaration
    class DustSurfaceSpecies;
    class DustMantleSpecies;

    /**
     * @class GasSpecies
     * @brief A class representing gas-phase chemical species.
     * 
     * This class inherits from `Species` and represents chemical species in the gas phase.
     * Each gas species has corresponding surface and mantle species, which are handled using shared pointers.
     */
    class GasSpecies
        : public Species
    {
    private:

        std::shared_ptr<DustSurfaceSpecies> corresponding_dust_surface_species_;
        std::shared_ptr<DustMantleSpecies>  corresponding_dust_mantle_species_;

    public:

        /**
         * @brief Constructor for GasSpecies.
         * 
         * @param index Unique species index.
         * @param name Name of the gas species.
         * @param charge Electrical charge of the species.
         * @param element_composition Composition of elements in the species.
         * @param ptr_element_manager Pointer to the ElementManager.
         */
        GasSpecies(
            std::size_t index,
            const std::string& name, 
            int charge, 
            const std::vector<std::size_t>& element_composition, 
            ElementManager* ptr_element_manager
        );

        /**
         * @brief Get the corresponding dust surface species.
         * @return Shared pointer to the corresponding DustSurfaceSpecies.
         */
        const std::shared_ptr<DustSurfaceSpecies>& GetCorrespondingSurfaceSpecies() const {
            return corresponding_dust_surface_species_;
        }

        /**
         * @brief Get the corresponding dust mantle species.
         * @return Shared pointer to the corresponding DustMantleSpecies.
         */
        const std::shared_ptr<DustMantleSpecies>& GetCorrespondingMantleSpecies() const {
            return corresponding_dust_mantle_species_;
        }

        /**
         * @brief Set the corresponding dust surface species.
         * @param dust_surface_species Shared pointer to the DustSurfaceSpecies.
         */
        void SetCorrespondingSurfaceSpecies(const std::shared_ptr<DustSurfaceSpecies>& dust_surface_species) {
            corresponding_dust_surface_species_ = dust_surface_species;
        }

        /**
         * @brief Set the corresponding dust mantle species.
         * @param dust_mantle_species Shared pointer to the DustMantleSpecies.
         */
        void SetCorrespondingMantleSpecies(const std::shared_ptr<DustMantleSpecies>& dust_mantle_species) {
            corresponding_dust_mantle_species_ = dust_mantle_species;
        }
    };
    
} // namespace nicole


#endif /* GAS_SPECIES_HPP */