/**
 * @file reaction.hpp
 * @brief Header file for the Reaction class.
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#ifndef REACTION_HPP
#define REACTION_HPP

#include "../nicole_defs.hpp"

namespace nicole
{
    /**
     * @class Reaction
     * @brief A class representing a chemical reaction.
     * 
     * This class stores information related to a chemical reaction, including the indices of reactant and product species, 
     * rate parameters, reaction type, and branching ratio of the reaction.
     */
    class Reaction
    {
    protected:

        // Index of the reaction
        std::size_t index_;

        // Indices of the reactant species
        std::vector<std::size_t> reactant_indices_;

        // Indices of the product species
        std::vector<std::size_t> product_indices_;

        // Rate parameters for calculating reaction rate coefficient of the reaction
        std::vector<double> rate_parameters_;

        // Type ID of the reaction
        std::size_t type_id_;

        // Branching ratio of the reaction
        double branching_ratio_;
        
    public:

        /**
         * @brief Constructor for initializing a Reaction object.
         * 
         * @param index Index of the reaction.
         * @param reactant_indices Indices of the reactant species.
         * @param product_indices Indices of the product species.
         * @param rate_parameters Rate parameters for the reaction.
         * @param type_id Type ID of the reaction.
         */
        Reaction(
            std::size_t index,
            const std::vector<std::size_t>& reactant_indices,
            const std::vector<std::size_t>& product_indices,
            const std::vector<double>& rate_parameters,
            std::size_t type_id
        );

        /**
         * @brief Get the branching ratio of the reaction.
         * 
         * @return Branching ratio.
         */
        inline double GetBranchingRatio() const { return branching_ratio_; }

        /**
         * @brief Set the branching ratio of the reaction.
         * 
         * @param branching_ratio Branching ratio value to set.
         */
        inline void SetBranchingRatio(const double branching_ratio) { branching_ratio_ = branching_ratio; }

        /**
         * @brief Print information about the reaction.
         * 
         * This method prints detailed information about the reaction, including the names of the reactant and product species.
         * 
         * @param species_name_list List of species names.
         */
        void PrintInfo(const std::vector<std::string>& species_name_list) const;

        /**
         * @brief Write information about the reaction to a file.
         * 
         * This method writes detailed information about the reaction to the provided file, including the names of the reactant 
         * and product species.
         * 
         * @param species_name_list List of species names.
         * @param file Output file stream to write information to.
         */
        void WriteInfoToFile(const std::vector<std::string>& species_name_list, std::ofstream &file) const;

        friend class ReactionManager;
        friend class ReactionSimulator;
    };
}

#endif /* REACTION_HPP */