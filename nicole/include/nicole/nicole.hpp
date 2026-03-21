/**
 * @file nicole.hpp
 * @brief NICOLE: A framework for computing chemical reactions and non-ideal MHD resistivity in astrophysics.
 * 
 * This header provides access to all major modules of the NICOLE framework.
 * 
 * @date 2025-02-11
 * @author Y. Kawasaki
 */

#ifndef NICOLE_HPP
#define NICOLE_HPP

// Manages chemical elements
#include "elements/element_manager.hpp"

// Manages chemical species
#include "species/species_manager.hpp"

// Manages chemical reactions
#include "reactions/reaction_manager.hpp"

// Runs chemical reaction simulations
#include "reaction_simulator/reaction_simulator.hpp"

// Computes non-ideal MHD resistivity
#include "mhd_resistivity/non_ideal_mhd_effect.hpp"

// namespace nicole {

//     class Nicole
//     {
//     public:
//         Nicole(InputConfig& input);
//         Nicole(InputConfig& input, const std::vector<std::string>& user_gas_species_list);

//     private:

//         ElementManager element_manager_;
//         SpeciesManager species_manager_;
//         ReactionManager reaction_manager_;
//         ReactionSimulator reaction_simulator_;
//         NonIdealMHDeffect non_ideal_mhd_effect_;

//         bool is_check_;
//     };
// }

#endif /* NICOLE_HPP */