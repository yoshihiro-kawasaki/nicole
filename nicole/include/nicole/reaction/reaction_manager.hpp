#ifndef REACTION_MANAGER_HPP_
#define REACTION_MANAGER_HPP_

#include <array>
#include <string>
#include <vector>

#include "nicole/nicole_defs.hpp"
#include "nicole/reaction/reaction.hpp"
#include "nicole/reaction/reaction_types.hpp"
#include "nicole/species/species_manager.hpp"

namespace nicole {
    class ReactionManager {
    public:
        ReactionManager(SpeciesManager *ptr_species_manager, InputConfig& input);

        void CheckReactionManager(const std::string& filename);

        friend class ReactionSimulator;
        
    protected:

        std::vector<SpeciesID> GetSpeciesIDList(const std::vector<std::string>& species_name_list) const;
        std::vector<std::string> GetSpeciesNameList(const std::vector<SpeciesID>& ids) const;

        bool AreReactionSpeciesInSpeciesNameList(const std::vector<std::string>& reactants, const std::vector<std::string>& products) const;
        
        std::vector<std::string> SplitGasReactionLine(const std::string& line) const;
        void ReadGasPhaseReactionFile(const std::string& file_path);

        std::vector<std::string> SplitDustSurfaceReactionLine(const std::string& line) const;
        void ReadDustSurfaceReactionFile(const std::string& file_path);

        std::vector<std::string> SplitSurfaceActivationEnergyLine(const std::string& line) const;
        Real CalculateTunnellingEffectProbabilty(const Real mass1, const Real mass2, const Real activation_energy) const;
        void ReadSurfaceActivationEnergyFile(const std::string& file_path);

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

        void CalculateReactionBranchingRatio();
        void CalculateChemicalDesorptionProbabilities();
        void CountNumberOfEachTypeReactions();
        void SetReactionTypeIdStartAndEnd();
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

#endif /* REACTION_MANAGER_HPP_ */
