#ifndef DUST_MANTLE_SPECIES_HPP_
#define DUST_MANTLE_SPECIES_HPP_

#include "nicole/nicole_defs.hpp"
#include "nicole/species/gas_species.hpp"
#include "nicole/species/dust_layer_species_base.hpp"
#include "nicole/species/dust_surface_species.hpp"
#include "nicole/species/species.hpp"

namespace nicole {
    // Forward declarations
    class GasSpecies;
    class DustSurfaceSpecies;

    class DustMantleSpecies : public DustLayerSpeciesBase {
    public:
        DustMantleSpecies(
            SpeciesID id,
            const std::string& name,
            std::shared_ptr<GasSpecies> corresponding_gas_species,
            std::shared_ptr<DustSurfaceSpecies> corresponding_dust_surface_species
        );

        // Getter
        inline std::shared_ptr<GasSpecies> GetCorrespondingGasSpecies() const {
            return corresponding_gas_species_.lock();
        }
        inline std::shared_ptr<DustSurfaceSpecies> GetCorrespondingDustSurfaceSpecies() const {
            return corresponding_dust_surface_species_.lock();
        }

        // Setter
        inline void SetCorrespondingGasSpecies(const std::shared_ptr<GasSpecies>& gas_species) {
            corresponding_gas_species_ = gas_species;
        }
        inline void SetCorrespondingDustSurfaceSpecies(const std::shared_ptr<DustSurfaceSpecies>& dust_surface_species) {
            corresponding_dust_surface_species_ = dust_surface_species;
        }

    private:
        void CalculateDiffusionBarrier();
        std::weak_ptr<GasSpecies> corresponding_gas_species_;
        std::weak_ptr<DustSurfaceSpecies> corresponding_dust_surface_species_;
    };
}

#endif /* DUST_MANTLE_SPECIES_HPP_ */
