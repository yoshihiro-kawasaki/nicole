/**
 * @file reaction.cpp
 * @brief Implementation of the Reaction class
 * 
 * This file contains the implementation of the Reaction class, which handles
 * chemical reactions by storing relevant indices for reactants and products,
 * and their associated rate parameters.
 * 
 * @date 2025-02-12
 */

#include "reaction.hpp"

namespace nicole
{
    /**
     * @brief Constructor for the Reaction class
     * 
     * Initializes a Reaction object with the given indices for reactants, 
     * products, rate parameters, and type ID. The branching ratio is set to 1.0 by default.
     * 
     * @param index Reaction index
     * @param reactant_indices Indices of the reactants
     * @param product_indices Indices of the products
     * @param rate_parameters Rate parameters for the reaction
     * @param type_id Type ID of the reaction
     */
    Reaction::Reaction(
        std::size_t index,
        const std::vector<std::size_t>& reactant_indices,
        const std::vector<std::size_t>& product_indices,
        const std::vector<double>& rate_parameters,
        std::size_t type_id
    ) :
        index_(index),
        reactant_indices_(reactant_indices),
        product_indices_(product_indices),
        rate_parameters_(rate_parameters),
        type_id_(type_id),
        branching_ratio_(1.0) // // Default branching ratio
    { }

    /**
     * @brief Prints the information of the reaction
     * 
     * This function prints the reaction details, including the indices of reactants 
     * and products, as well as the type ID. It uses the provided species name list 
     * to display the names of the reactants and products.
     * 
     * @param species_name_list A list of species names corresponding to species indices
     */
    void Reaction::PrintInfo(const std::vector<std::string>& species_name_list) const 
    {
        std::string species_name; // Temporary variable to store species name

        std::cout << std::setw(6) << index_ << " ";
        // Loop through and print reactants
        for (const std::size_t index : reactant_indices_) {
            if (index == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[index];
            }
            std::cout << std::setw(11) << species_name << " ";
        }

        // Loop through and print products
        for (const std::size_t index : product_indices_) {
            if (index == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[index];
            }
            std::cout << std::setw(11) << species_name << " ";
        }
        std::cout << std::setw(4) << type_id_ << std::endl;
    }

    /**
     * @brief Writes the information of the reaction to a file
     * 
     * This function writes the reaction details to an output file, including the 
     * indices of reactants and products, and the reaction type ID. It uses the 
     * provided species name list to display the names of the reactants and products.
     * 
     * @param species_name_list A list of species names corresponding to species indices
     * @param file An output file stream to write the reaction data
     */
    void Reaction::WriteInfoToFile(const std::vector<std::string>& species_name_list, std::ofstream &file) const 
    {
        std::string species_name;

        file << std::setw(6) << index_ << " ";
        // Loop through and write reactants to file
        for (const std::size_t index : reactant_indices_) {
            if (index == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[index];
            }
            file << std::setw(11) << species_name << " ";
        }

        // Loop through and write products to file
        for (const std::size_t index : product_indices_) {
            if (index == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[index];
            }
            file << std::setw(11) << species_name << " ";
        }
        file << std::setw(4) << type_id_ << std::endl;
    }
}