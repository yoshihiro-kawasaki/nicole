#include "nicole/species/dust_layer_species_base.hpp"

namespace nicole {
    DustLayerSpeciesBase::DustLayerSpeciesBase(
        SpeciesID id,
        const std::string& name, 
        int charge, 
        const std::vector<std::size_t> element_composition,
        SpeciesType type,
        ElementManager* ptr_element_manager,
        double binding_energy_on_H2O_ice,
        double binding_energy_on_silicate
    ) : Species(
            id,
            name, 
            charge, 
            element_composition,
            type, 
            ptr_element_manager,
            binding_energy_on_H2O_ice,
            binding_energy_on_silicate
        )
    { }


    void DustLayerSpeciesBase::CalculateVibrationFrequency() {
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
