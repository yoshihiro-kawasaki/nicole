#ifndef GAS_SPECIES_HPP_
#define GAS_SPECIES_HPP_

#include "nicole/species/species.hpp"

namespace nicole {
    // Forward declaration
    class DustSurfaceSpecies;
    class DustMantleSpecies;

    class GasSpecies : public Species {
    public:
        GasSpecies(
            SpeciesID id,
            const std::string& name,
            int charge,
            const std::vector<std::size_t>& element_composition,
            ElementManager* ptr_element_manager
        );

        // Getter
        std::shared_ptr<DustSurfaceSpecies> GetCorrespondingSurfeceSpecies() const {
            return corresponding_dust_surface_species_;
        }
        std::shared_ptr<DustMantleSpecies> GetCorrespondingMantleSpecies() const {
            return corresponding_dust_mantle_species_;
        }

        // Setter
        void SetCorrespondingSurfeceSpecies(const std::shared_ptr<DustSurfaceSpecies>& dust_surface_species) {
            corresponding_dust_surface_species_ = dust_surface_species;
        }
        void SetCorrespondingMantleSpecies(const std::shared_ptr<DustMantleSpecies>& dust_mantle_species) {
            corresponding_dust_mantle_species_ = dust_mantle_species;
        }
    private:
        std::shared_ptr<DustSurfaceSpecies> corresponding_dust_surface_species_;
        std::shared_ptr<DustMantleSpecies>  corresponding_dust_mantle_species_;
    };
}

#endif /* GAS_SPECIES_HPP_ */
