#include <fstream>
#include <iomanip>
#include <iostream>

#include "nicole/reaction/reaction.hpp"

namespace nicole {
    Reaction::Reaction(
        ReactionID id,
        const std::vector<SpeciesID>& reactant_ids,
        const std::vector<SpeciesID>& product_ids,
        const std::vector<double>& rate_parameters,
        std::size_t type_id
    ) :
        id_(id),
        reactant_ids_(reactant_ids),
        product_ids_(product_ids),
        rate_parameters_(rate_parameters),
        type_id_(type_id),
        branching_ratio_(1.0)  // Default branching ratio
    { }


    void Reaction::PrintInfo(const std::vector<std::string>& species_name_list) const {
        std::string species_name;

        std::cout << std::setw(6) << id_ << " ";
        for (const SpeciesID id : reactant_ids_) {
            if (id == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[id];
            }
            std::cout << std::setw(11) << species_name << " ";
        }

        for (const SpeciesID id : product_ids_) {
            if (id == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[id];
            }
            std::cout << std::setw(11) << species_name << " ";
        }
        std::cout << std::setw(4) << type_id_ << std::endl;
    }


    void Reaction::WriteInfoToFile(const std::vector<std::string>& species_name_list, std::ofstream &file) const {
        std::string species_name;

        file << std::setw(6) << id_ << " ";
        for (const SpeciesID id : reactant_ids_) {
            if (id == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[id];
            }
            file << std::setw(11) << species_name << " ";
        }

        for (const SpeciesID id : product_ids_) {
            if (id == kNotFoundSpecies) {
                species_name = kEmptyString;
            } else {
                species_name = species_name_list[id];
            }
            file << std::setw(11) << species_name << " ";
        }
        file << std::setw(4) << type_id_ << std::endl;
    }
}
