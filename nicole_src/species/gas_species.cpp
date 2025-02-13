/**
 * @file gas_species.cpp
 * @brief Implimentation of class GasSpecies
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#include "gas_species.hpp"

namespace nicole
{
    /**
     * @brief Constructor for GasSpecies.
     * 
     * Initializes a gas-phase species with its properties.
     *
     * @param index Unique species index.
     * @param name Name of the gas species.
     * @param charge Electrical charge of the species.
     * @param element_composition Composition of elements in the species.
     * @param ptr_element_manager Pointer to the ElementManager.
     */
    GasSpecies::GasSpecies(
        std::size_t index,
        const std::string& name, 
        int charge, 
        const std::vector<std::size_t>& element_composition, 
        ElementManager* ptr_element_manager
    ) : Species(
            index,
            name, 
            charge, 
            element_composition, 
            ptr_element_manager, 
            SpeciesType::Gas
        ),
        corresponding_dust_surface_species_(nullptr),
        corresponding_dust_mantle_species_(nullptr)
    { }
}