/**
 * @file species.cpp
 * @brief Implimentation of class species, which is the base class for different chemical species
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#include "species.hpp"

namespace nicole
{
    /**
     * @brief Constructor for gas-phase species.
     * @param index Species index
     * @param name Species name
     * @param charge Electronic charge
     * @param element_composition List of element indices composing the species
     * @param ptr_element_manager Pointer to ElementManager
     * @param type Species type
     */
    Species::Species(
        std::size_t index,
        const std::string& name,
        int charge, 
        const std::vector<std::size_t>& element_composition,
        ElementManager* ptr_element_manager,
        SpeciesType type
    ) : index_(index),
        name_(name), 
        charge_(charge),
        mass_(0.0),
        element_composition_(element_composition),
        ptr_element_manager_(ptr_element_manager),
        type_(type),
        binding_energy_on_H2O_ice_(0.0),
        binding_energy_on_silicate_(0.0),
        enthalpy_of_formation_(0.0)
    {
        // Check if ptr_element_manager_ is nullptr
        if (!ptr_element_manager_) {
            throw std::runtime_error("Error: ElementManager pointer is null.");
        }
        CalculateMass();
    }

    /**
     * @brief Constructor for dust surface and mantle species.
     * @param index Species index
     * @param name Species name
     * @param charge Electronic charge
     * @param element_composition List of element indices composing the species
     * @param ptr_element_manager Pointer to ElementManager
     * @param type Species type
     * @param binding_energy_on_H2O_ice Binding energy on H2O ice surface
     * @param binding_energy_on_silicate Binding energy on bare silicate surface
     */
    Species::Species(
        std::size_t index,
        const std::string& name,
        int charge, 
        const std::vector<std::size_t>& element_composition,
        ElementManager* ptr_element_manager,
        SpeciesType type,
        double binding_energy_on_H2O_ice,
        double binding_energy_on_silicate
    ) : index_(index),
        name_(name), 
        charge_(charge),
        mass_(0.0),
        element_composition_(element_composition),
        ptr_element_manager_(ptr_element_manager),
        type_(type),
        binding_energy_on_H2O_ice_(binding_energy_on_H2O_ice),
        binding_energy_on_silicate_(binding_energy_on_silicate),
        enthalpy_of_formation_(0.0)
    {
        // Check if ptr_element_manager_ is nullptr
        if (!ptr_element_manager_) {
            throw std::runtime_error("Error: ElementManager pointer is null.");
        }
        CalculateMass();
    }

    /**
     * @brief Constructor for dust species.
     * @param index Species index
     * @param name Species name
     * @param charge Electronic charge
     * @param type Species type
     */
    Species::Species(
        std::size_t index,
        const std::string& name,
        int charge,
        SpeciesType type
    ) : index_(index),
        name_(name),
        charge_(charge),
        ptr_element_manager_(nullptr),
        type_(type)
    { }

    /**
     * @brief Display species properties.
     */
    void Species::DisplayProperties() const 
    {
        const std::size_t num_setw = 10;
        std::cout << std::setw(num_setw) << "Index  : " << index_  << " "
                  << std::setw(num_setw) << "Name   : " << name_   << " "
                  << std::setw(num_setw) << "Mass   : " << mass_   << " "
                  << std::setw(num_setw) << "Charge : " << charge_ << std::endl;
    }

    /**
     * @brief Calculate species mass based on elemental composition.
     */
    void Species::CalculateMass()
    {
        if (!ptr_element_manager_) {
            std::cerr << "Error: ElementManager is not initialized!" << std::endl;
            throw std::runtime_error("Error: ElementManager is not initialized!");
            mass_ = 0.0;
            return;
        }

        if (element_composition_.empty()) {
            mass_ = 0.0;
            return;
        }

        size_t number_of_elements = element_composition_.size();
        mass_ = 0.0;
        for (std::size_t i = 0; i < number_of_elements; ++i) {
            mass_ += ptr_element_manager_->GetElementMass(i) * static_cast<double>(element_composition_[i]);
        }
        mass_ *= constants::kProtonMass;

        if (name_ == nicole::kElectron) {
            mass_ = constants::kElectronMass;
        }
    }
}