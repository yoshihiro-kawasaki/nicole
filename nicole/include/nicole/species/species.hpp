#ifndef SPECIES_HPP_
#define SPECIES_HPP_

#include "nicole/element/element_manager.hpp"
#include "nicole/nicole_defs.hpp"

namespace nicole {
    enum class SpeciesType {
        Gas,
        Surface,
        Mantle,
        Dust
    };

    class Species {
    public:
        // Constructor for GasSpecies
        Species(
            SpeciesID id,
            const std::string& name,
            int charge,
            const std::vector<std::size_t>& element_composition,
            SpeciesType type,
            ElementManager* ptr_element_manager
        );

        // Constructor for DustSurfaceSpecies and DustMantleSpecies.
        Species(
            SpeciesID id,
            const std::string& name,
            int charge,
            const std::vector<std::size_t>& element_composition,
            SpeciesType type,
            ElementManager* ptr_element_manager,
            Real binding_energy_on_H2O_ice,
            Real binding_energy_on_silicate
        );

        // Constructor for DustSpecies.
        Species(
            SpeciesID id,
            const std::string& name,
            int charge,
            SpeciesType type
        );

        virtual ~Species() = 0;

        // Getter
        inline SpeciesID GetID() const { return id_; }
        inline const std::string& GetName() const { return name_; }
        inline Real GetMass() const { return mass_; }
        inline int GetCharge() const { return charge_; }
        inline const std::vector<std::size_t>& GetElementComposition() const {
            return element_composition_;
        }
        inline SpeciesType GetSpeciesType() const { return type_; }
        inline Real GetBindingEnergyOnH2Oice() const { return binding_energy_on_H2O_ice_; }
        inline Real GetBindingEnergyOnBareSilicate() const { return binding_energy_on_silicate_; }
        inline Real GetEnthalpyOfFormation() const { return enthalpy_of_formation_; }
        inline ElementManager* GetPtrElementManager() const { return ptr_element_manager_; }

        // Setter
        inline void SetID(const SpeciesID id) { id_ = id; }
        inline void SetBindingEnergyOnH2Oice(const Real binding_energy_on_H2O_ice) {
            binding_energy_on_H2O_ice_ = binding_energy_on_H2O_ice;
        }
        inline void SetBindingEnergyOnBareSilicate(const Real binding_energy_on_silicate) {
            binding_energy_on_silicate_ = binding_energy_on_silicate;
        }
        inline void SetEnthalpyOfFormation(const Real enthalpy_of_formation) {
            enthalpy_of_formation_ = enthalpy_of_formation; 
        }

        void CalculateMass();

    protected:
        SpeciesID id_;
        std::string name_;
        Real mass_;
        int charge_;
        std::vector<std::size_t> element_composition_;
        SpeciesType type_;
        Real binding_energy_on_H2O_ice_;
        Real binding_energy_on_silicate_;
        Real enthalpy_of_formation_;
        ElementManager* ptr_element_manager_;
    };
}

#endif /* SPECIES_HPP_ */
