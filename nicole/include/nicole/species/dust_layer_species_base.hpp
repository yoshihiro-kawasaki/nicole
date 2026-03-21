#ifndef DUST_LAYER_SPECIES_BASE_HPP_
#define DUST_LAYER_SPECIES_BASE_HPP_

#include "nicole/nicole_defs.hpp"
#include "nicole/species/species.hpp"

namespace nicole {
    class DustLayerSpeciesBase : public Species {
    public:
        DustLayerSpeciesBase(
            SpeciesID id,
            const std::string& name, 
            int charge, 
            const std::vector<std::size_t> element_composition,
            SpeciesType type,
            ElementManager* ptr_element_manager,
            Real binding_energy_on_H2O_ice,
            Real binding_energy_on_silicate
        );

        // Getter
        inline Real GetVibrationFrequencyOnH2Oice() const { return vibration_frequency_on_H2O_ice_; }
        inline Real GetVibrationFrequencyOnSilicate() const { return vibration_frequency_on_silicate_; }
        inline Real GetDiffusionBarrierOnH2Oice() const { return diffusion_barrier_on_H2O_ice_; }
        inline Real GetDiffusionBarrierOnSilicate() const { return diffusion_barrier_on_silicate_; }

        // Setter
        inline void SetVibrationFrequencyOnH2Oice(const Real vibration_frequency_on_H2O_ice) { 
            vibration_frequency_on_H2O_ice_ = vibration_frequency_on_H2O_ice;
        }
        inline void SetVibrationFrequencyOnSilicate(const Real vibration_frequency_on_silicate) {
            vibration_frequency_on_silicate_ = vibration_frequency_on_silicate;
        }
        inline void SetDiffusionBarrierOnH2Oice(const Real diffusion_barrier_on_H2O_ice) {
            diffusion_barrier_on_H2O_ice_ = diffusion_barrier_on_H2O_ice;
        }
        inline void SetDiffusionBarrierOnSilicate(const Real diffusion_barrier_on_silicate) {
            diffusion_barrier_on_silicate_ = diffusion_barrier_on_silicate;
        }
    protected:
        void CalculateVibrationFrequency();
        virtual void CalculateDiffusionBarrier() { };
        Real vibration_frequency_on_H2O_ice_;
        Real vibration_frequency_on_silicate_;
        Real diffusion_barrier_on_H2O_ice_;
        Real diffusion_barrier_on_silicate_;
    };
}

#endif /* DUST_LAYER_SPECIES_BASE_HPP_ */
