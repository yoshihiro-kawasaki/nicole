/**
 * @file reaction_manager.hpp
 * @brief Header file for the ReactionManager class, which manages chemical reactions.
 * @date 2025-02-12
 */

#ifndef REACTION_MANAGER_HPP
#define REACTION_MANAGER_HPP

#include "../nicole_defs.hpp"
#include "../species/species_manager.hpp"
#include "reaction.hpp"
#include "reaction_types.hpp"

namespace nicole
{
    /**
     * @class ReactionManager
     * @brief Manages chemical reactions, including reading reaction data, generating reaction lists,
     *        and calculating probabilities and branching ratios for different types of reactions.
     */
    class ReactionManager
    {
    public:

        /**
         * @brief Constructor that initializes the ReactionManager with the SpeciesManager and input configuration.
         * @param ptr_species_manager A pointer to the SpeciesManager, which manages chemical species.
         * @param input The input configuration object that provides simulation parameters.
         */
        ReactionManager(SpeciesManager *ptr_species_manager, InputConfig& input);

        /**
         * @brief Checks the consistency of the ReactionManager by verifying reaction data from the specified file.
         * @param filename The name of the file containing reaction data to verify.
         */
        void CheckReactionManager(const std::string& filename);

        friend class ReactionSimulator;
        
    protected:

        /**
         * @brief Retrieves the indices list of species based on a list of species names.
         * @param species_name_list A vector containing the names of species.
         * @return A vector containing the indices of the species corresponding to the names in the input list.
         */
        std::vector<std::size_t> GetSpeciesListIndices(const std::vector<std::string>& species_name_list);

        /**
         * @brief Retrieves the names of species based on a list of species indices.
         * @param indices A vector containing the indices of species.
         * @return A vector containing the names of the species corresponding to the indices in the input list.
         */
        std::vector<std::string> GetSpeciesNameList(const std::vector<std::size_t>& indices);

        /**
         * @brief Checks whether all reactants and products are present in the species name list.
         * @param reactants A vector containing the names of reactant species.
         * @param products A vector containing the names of product species.
         * @return Returns true if all species in reactants and products are found in the species list or special species list, false otherwise.
         */
        bool AreReactionSpeciesInSpeciesNameList(const std::vector<std::string>& reactants, const std::vector<std::string>& products);

        /**
         * @brief Reads the gas-phase reaction data from a file and stores the relevant information.
         * @param filename The name of the file containing the gas-phase reaction data.
         */
        void ReadGasPhaseReactionFile(std::string filename);

        /**
         * @brief Splits a gas-phase reaction line into individual components.
         * @param line The line from the gas-phase reaction file to be split.
         * @return A vector of strings containing the individual components of the reaction.
         */
        std::vector<std::string> SplitGasReactionLine(const std::string& line);

        /**
         * @brief Reads the dust surface reaction data from a file and stores the relevant information.
         * @param filename The name of the file containing the dust surface reaction data.
         */
        void ReadDustSurfaceReactionFile(std::string filename);

        /**
         * @brief Splits a line from the dust surface reaction file into individual components.
         * @param line A single line from the dust surface reaction file containing the reaction information.
         * @return A vector of strings, where each string represents a different component of the reaction:
         *         - Reactants (3 items)
         *         - Products (5 items)
         *         - Rate parameters (branching ratio and uncertainty)
         */
        std::vector<std::string> SplitDustSurfaceReactionLine(const std::string& line);

        /**
         * @brief Reads the surface activation energy file and updates the activation energy and tunneling effect probability for reactions.
         * @param filename The name of the file containing the activation energy data for surface reactions.
         */
        void ReadSurfaceActivationEnergyFile(std::string filename);

        /**
         * @brief Calculates the tunneling effect probability for a reaction.
         * @param mass1 The mass of the first reactant.
         * @param mass2 The mass of the second reactant.
         * @param activation_energy The activation energy for the reaction.
         * @return The calculated tunneling effect probability.
         */
        double CalculateTunnellingEffectProbabilty(const double mass1, const double mass2, const double activation_energy);

        /**
         * @brief Splits a line from the surface activation energy file into its components.
         * @param line The line from the surface activation energy file.
         * @return A vector containing the split components (reactants, products, and rate parameters).
         */
        std::vector<std::string> SplitSurfaceActivationEnergyLine(const std::string& line);

        // Generate reaction list related to dust particles
        void GenerateReactionListForDustAndChargedParticleCollision();
        void GenerateReactionListForDustCollision();
        void GenerateReactionListForNeutralSpeciesAccretionOnDustSurfaces();
        void GenerateReactionListForThermalDesorptionOnDustSurfaces();
        void GenerateReactionListForCosmicRayDesorptionOnDustSurfaces();
        void GenerateReactionListForPhotoDesorptionByExternalUV();
        void GenerateReactionListForPhotoDesorptionByCosmicRayGeneratedUV();
        void GenerateReactionListForPhotoDissociationInducedByCRsOnDustSurfaces();
        void GenerateReactionListForPhotoDissociationByExternalUVOnDustSurfaces();
        void GenerateReactionListForDustMantleReaction();
        void GenerateReactionListForDustSurfaceToMantleSwapping();
        void GenerateReactionListForDustMantleToSurfaceSwapping();

        /**
         * @brief Calculates the branching ratio for each reaction.
         */
        void CalculateReactionBranchingRatio();

        /**
         * @brief Calculates the chemical desorption probabilities.
         */
        void CalculateChemicalDesorptionProbabilities();

        /**
         * @brief Counts the number of reactions of each type.
         */
        void CountNumberOfEachTypeReactions();

        /**
         * @brief Sets the start and end indices for each reaction type ID.
         */
        void SetReactionTypeIdStartAndEnd();

        /**
         * @brief Sets the reactions involving specific species.
         */
        void SetReactionsInvolvedWithSpecies();

        // Pointer to the SpeciesManager for managing species.
        SpeciesManager *ptr_species_manager_;

        // List of chemical reactions.
        std::vector<std::shared_ptr<Reaction> > reaction_list_;

        // number of reactions
        std::array<std::size_t, reaction_type_id::kNumberOfTypeID> number_of_each_type_reactions_ = {0};
        std::size_t total_number_of_reactions_;
        std::size_t total_number_of_gas_phase_reactions_;

        // Start and end indices for each reaction type.
        std::array<std::size_t, reaction_type_id::kNumberOfTypeID> reaction_type_id_start_ = {0};
        std::array<std::size_t, reaction_type_id::kNumberOfTypeID> reaction_type_id_end_   = {0};

        // Number of reactions each species is involved in.
        std::vector<std::size_t> number_of_reactions_involved_with_species_;

        // List of reaction indices each species is involved in.
        std::vector<std::vector<std::size_t> > reaction_index_list_involved_with_species_;

        // Flags
        bool is_dust_collision_;
        bool is_dust_surface_reaction_;
        bool is_chemical_desorption_;
        bool is_H2_desorption_;
        bool is_three_phase_reaction_;
        bool is_H2_self_shielding_;
        bool is_CO_self_shielding_;

    };
}

#endif /* REACTION_MANAGER_HPP */