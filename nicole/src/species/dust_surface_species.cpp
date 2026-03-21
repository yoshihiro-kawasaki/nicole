#include "nicole/species/dust_surface_species.hpp"

namespace nicole {
    DustSurfaceSpecies::DustSurfaceSpecies(
        SpeciesID id,
        const std::string& name,
        std::shared_ptr<GasSpecies> corresponding_gas_species,
        Real binding_energy_on_H2O_ice,
        Real binding_energy_on_silicate
    ) : DustLayerSpeciesBase(
            id,
            name, 
            corresponding_gas_species->GetCharge(), 
            corresponding_gas_species->GetElementComposition(),
            SpeciesType::Surface,
            corresponding_gas_species->GetPtrElementManager(),
            binding_energy_on_H2O_ice,
            binding_energy_on_silicate
        ),
        corresponding_gas_species_(corresponding_gas_species),
        corresponding_dust_mantle_species_(nullptr)
    { 
        CalculateVibrationFrequency();
        CalculateDiffusionBarrier();
    }


    void DustSurfaceSpecies::CalculateDiffusionBarrier() 
    {
        if (name_ == "sH") {
            diffusion_barrier_on_H2O_ice_ = 230.0;
            diffusion_barrier_on_silicate_ = 230.0;
        } else {
            diffusion_barrier_on_H2O_ice_ = kBindingDiffusionRatioOfDustSurface * binding_energy_on_H2O_ice_;
            diffusion_barrier_on_silicate_ = kBindingDiffusionRatioOfDustSurface * binding_energy_on_silicate_;
        }
    }
}
