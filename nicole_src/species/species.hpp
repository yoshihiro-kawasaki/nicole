/**
 * @file species.hpp
 * @brief Base class for different chemical species
 * @date 2025/02/12
 * @author Y. Kawasaki
 */

#ifndef SPECIES_HPP
#define SPECIES_HPP

#include "../nicole_defs.hpp"
#include "../elements/element_manager.hpp"

namespace nicole
{
    /**
     * @enum SpeciesType
     * @brief Represents different types of species in Nicole code
     */
    enum class SpeciesType {
        Gas,
        Surface,
        Mantle,
        Dust
    };

    /**
     * @class Species
     * @brief Base class representing a chemical species.
     */
    class Species 
    {
    protected:

        // Unique index of the species
        std::size_t index_;
        // Species name
        std::string name_;
        // Species electronic charge
        int charge_;
        // Species charge in gram
        double mass_;
        // Number of each elements composing species
        std::vector<std::size_t> element_composition_;
        // Pointer to ElementManager
        ElementManager *ptr_element_manager_;
        // Species type
        SpeciesType type_;
        // Binding energy on dust surfaces of H2O ice [K]
        double binding_energy_on_H2O_ice_;
        // Binding energy on dust surfaces of (bare) silicate [K]
        double binding_energy_on_silicate_;

        double enthalpy_of_formation_;  // [K]

    public:

        /**
         * @brief Constructor for GasSpecies.
         * @param index Unique index of the species.
         * @param name Species name.
         * @param charge Species charge.
         * @param element_composition Element composition indices.
         * @param ptr_element_manager Pointer to ElementManager.
         * @param type Species type (Gas).
         */
        Species(
            std::size_t index,
            const std::string& name, 
            int charge, 
            const std::vector<std::size_t>& element_composition, 
            ElementManager* ptr_element_manager,
            SpeciesType type
        );

        /**
         * @brief Constructor for DustSurfaceSpecies and DustMantleSpecies.
         * @param index Unique index of the species.
         * @param name Species name.
         * @param charge Species charge.
         * @param element_composition Element composition indices.
         * @param ptr_element_manager Pointer to ElementManager.
         * @param type Species type (Surface or Mantle).
         * @param binding_energy_on_H2O_ice Binding energy on H2O ice [K].
         * @param binding_energy_on_silicate Binding energy on silicate [K].
         */
        Species(
            std::size_t index,
            const std::string& name, 
            int charge, 
            const std::vector<std::size_t>& element_composition, 
            ElementManager* ptr_element_manager,
            SpeciesType type,
            double binding_energy_on_H2O_ice,
            double binding_energy_on_silicate
        );

        /**
         * @brief Constructor for DustSpecies.
         * @param index Unique index of the species.
         * @param name Species name.
         * @param charge Species charge.
         * @param type Species type (Dust).
         */
        Species(
            std::size_t index,
            const std::string& name,
            int charge,
            SpeciesType type
        );

        virtual ~Species() = 0;

        // Get functions
        inline std::size_t GetIndex() const { return index_; }
        inline const std::string& GetName() const { return name_; }
        inline double GetMass() const { return mass_; }
        inline int GetCharge() const { return charge_; }
        inline SpeciesType GetSpeciesType() const { return type_; }
        inline const std::vector<std::size_t>& GetElementComposition() const {
            return element_composition_;
        };
        inline ElementManager* GetPtrElementManager() const {
            return ptr_element_manager_;
        }
        inline double GetBindingEnergyOnH2Oice() const { return binding_energy_on_H2O_ice_; }
        inline double GetBindingEnergyOnBareSilicate() const { return binding_energy_on_silicate_; }
        inline double GetEnthalpyOfFormation() const { return enthalpy_of_formation_; }

        // Set functions
        inline void SetIndex(const std::size_t index) { index_ = index; }
        inline void SetBindingEnergyOnH2Oice(const double binding_energy_on_H2O_ice) {
            binding_energy_on_H2O_ice_ = binding_energy_on_H2O_ice;
        }
        inline void SetBindingEnergyOnBareSilicate(const double binding_energy_on_silicate) {
            binding_energy_on_silicate_ = binding_energy_on_silicate;
        }
        inline void SetEnthalpyOfFormation(const double enthalpy_of_formation) {
            enthalpy_of_formation_ = enthalpy_of_formation; 
        }

        /**
         * @brief Calculate species mass based on element_composition_
         */
        void CalculateMass();

        /**
         * @brief display species properties
         */
        virtual void DisplayProperties() const;
    };

}

#endif /* SPECIES_HPP */