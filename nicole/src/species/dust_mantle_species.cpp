#include "nicole/species/dust_mantle_species.hpp"

namespace nicole {
    DustMantleSpecies::DustMantleSpecies(
        SpeciesID id,
        const std::string& name,
        std::shared_ptr<GasSpecies> corresponding_gas_species,
        std::shared_ptr<DustSurfaceSpecies> corresponding_dust_surface_species
    ) : DustLayerSpeciesBase(
            id,
            name,
            corresponding_gas_species->GetCharge(),
            corresponding_gas_species->GetElementComposition(),
            SpeciesType::Mantle,
            corresponding_gas_species->GetPtrElementManager(),
            corresponding_dust_surface_species->GetBindingEnergyOnH2Oice(),
            corresponding_dust_surface_species->GetBindingEnergyOnBareSilicate()
        ),
        corresponding_gas_species_(corresponding_gas_species),
        corresponding_dust_surface_species_(corresponding_dust_surface_species)
    {
        CalculateVibrationFrequency();
        CalculateDiffusionBarrier();
    }


    void DustMantleSpecies::CalculateDiffusionBarrier() {
        diffusion_barrier_on_H2O_ice_ = kBindingDiffusionRatioOfDustMantle * binding_energy_on_H2O_ice_;
        diffusion_barrier_on_silicate_ = kBindingDiffusionRatioOfDustMantle * binding_energy_on_silicate_;
    }
}
