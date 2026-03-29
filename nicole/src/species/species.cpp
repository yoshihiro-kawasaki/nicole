#include "nicole/species/species.hpp"

namespace nicole {
    // Constructor for GasSpecies
    Species::Species(
        SpeciesID id,
        const std::string& name,
        int charge,
        const std::vector<std::size_t>& element_composition,
        SpeciesType type,
        ElementManager* ptr_element_manager
    ) : id_(id),
        name_(name),
        mass_(0.0),
        charge_(charge),
        element_composition_(element_composition),
        type_(type),
        binding_energy_on_H2O_ice_(0.0),
        binding_energy_on_silicate_(0.0),
        enthalpy_of_formation_(0.0),
        ptr_element_manager_(ptr_element_manager)
    {
        CalculateMass();
    }

    
    // Constructor for DustSurfaceSpecies and DustMantleSpecies.
    Species::Species(
        SpeciesID id,
        const std::string& name,
        int charge,
        const std::vector<std::size_t>& element_composition,
        SpeciesType type,
        ElementManager* ptr_element_manager,
        Real binding_energy_on_H2O_ice,
        Real binding_energy_on_silicate
    ) : id_(id),
        name_(name),
        mass_(0.0),
        charge_(charge),
        element_composition_(element_composition),
        type_(type),
        binding_energy_on_H2O_ice_(binding_energy_on_H2O_ice),
        binding_energy_on_silicate_(binding_energy_on_silicate),
        enthalpy_of_formation_(0.0),
        ptr_element_manager_(ptr_element_manager)
    {
        CalculateMass();
    }


    // Constructor for DustSpecies.
    Species::Species(
        SpeciesID id,
        const std::string& name,
        int charge,
        SpeciesType type
    ) : id_(id),
        name_(name),
        charge_(charge),
        type_(type),
        ptr_element_manager_(nullptr)
    { }

    Species::~Species() { }


    void Species::CalculateMass() {
        if (!ptr_element_manager_) {
            mass_ = 0.0;
            return;
        }

        if (element_composition_.empty()) {
            mass_ = 0.0;
            return;
        }

        if (name_ == nicole::kElectron) {
            mass_ = constants::kElectronMass;
            return;
        }

        std::size_t number_of_elements = element_composition_.size();
        mass_ = 0.0;
        for (ElementID id = 0; id < number_of_elements; ++id) {
            mass_ += ptr_element_manager_->GetElementMass(id) * static_cast<Real>(element_composition_[id]);
        }
        mass_ *= constants::kProtonMass;
        return;
    }
}
