#ifndef REACTION_HPP_
#define REACTION_HPP_

#include "nicole/nicole_defs.hpp"

namespace nicole {
    class Reaction {    
    public:
        Reaction(
            ReactionID id,
            const std::vector<SpeciesID>& reactant_ids,
            const std::vector<SpeciesID>& product_ids,
            const std::vector<Real>& rate_parameters,
            std::size_t type_id
        );

        inline Real GetBranchingRatio() const { return branching_ratio_; }
        inline void SetBranchingRatio(const Real branching_ratio) { branching_ratio_ = branching_ratio; }

        void PrintInfo(const std::vector<std::string>& species_name_list) const;
        void WriteInfoToFile(const std::vector<std::string>& species_name_list, std::ofstream &file) const;

        friend class ReactionManager;
        friend class ReactionSimulator;

    protected:
        ReactionID id_;
        std::vector<SpeciesID> reactant_ids_;
        std::vector<SpeciesID> product_ids_;
        std::vector<Real> rate_parameters_;
        std::size_t type_id_;
        Real branching_ratio_;
    };
}

#endif /* REACTION_HPP_ */