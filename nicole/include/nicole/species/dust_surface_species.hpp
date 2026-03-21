#ifndef DUST_SURFACE_SPECIES_HPP_
#define DUST_SURFACE_SPECIES_HPP_

#include "nicole/nicole_defs.hpp"
#include "nicole/species/gas_species.hpp"
#include "nicole/species/dust_layer_species_base.hpp"
#include "nicole/species/dust_mantle_species.hpp"
#include "nicole/species/species.hpp"

namespace nicole {
    // Forward declaration
    class GasSpecies;
    class DustMantleSpecies;

    class DustSurfaceSpecies : public DustLayerSpeciesBase {
    public:
        DustSurfaceSpecies(
            SpeciesID id,
            const std::string& name,
            std::shared_ptr<GasSpecies> corresponding_gas_species,
            Real binding_energy_on_H2O_ice,
            Real binding_energy_on_silicate
        );

        // Getter
        inline std::shared_ptr<GasSpecies> GetCorrespondingGasSpecies() const {
            return corresponding_gas_species_.lock();
        }
        inline std::shared_ptr<DustMantleSpecies> GetCorrespondingDustMantleSpecies() const {
            return corresponding_dust_mantle_species_;
        }

        // Setter
        inline void SetCorrespondingGasSpecies(const std::shared_ptr<GasSpecies>& gas_species) {
            corresponding_gas_species_= gas_species;
        }
        inline void SetCorrespondingDustMantleSpecies(const std::shared_ptr<DustMantleSpecies>& dust_mantle_species) {
            corresponding_dust_mantle_species_ = dust_mantle_species;
        }

    private:
        void CalculateDiffusionBarrier();
        std::weak_ptr<GasSpecies> corresponding_gas_species_;
        std::shared_ptr<DustMantleSpecies> corresponding_dust_mantle_species_;
    };

}

#endif /* DUST_SURFACE_SPECIES_HPP_ */
