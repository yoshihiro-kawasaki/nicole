#include <fstream>
#include <iomanip>
#include <iostream>

#include "nicole/reaction/reaction.hpp"

namespace nicole {
    Reaction::Reaction(
        ReactionID id,
        const std::vector<std::size_t>& reactant_indices,
        const std::vector<std::size_t>& product_indices,
        const std::vector<double>& rate_parameters,
        std::size_t type_id
    ) :
        id_(id),
        reactant_indices_(reactant_indices),
        product_indices_(product_indices),
        rate_parameters_(rate_parameters),
        type_id_(type_id),
        branching_ratio_(1.0)  // Default branching ratio
    { }


    void Reaction::PrintInfo(const std::vector<std::string>& species_name_list) const {
        std::string species_name; // Temporary variable to store species name

        std::cout << std::setw(6) << id_ << " ";
        for (const std::size_t index : reactant_indices_) {
            if (index == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[index];
            }
            std::cout << std::setw(11) << species_name << " ";
        }

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


    void Reaction::WriteInfoToFile(const std::vector<std::string>& species_name_list, std::ofstream &file) const {
        std::string species_name;

        file << std::setw(6) << id_ << " ";
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
