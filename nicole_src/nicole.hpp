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

#endif /* NICOLE_HPP */