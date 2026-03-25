/**
 * @file self_shielding_factor.hpp
 * @brief Provides shielding factor data for H2 and CO molecules.
 * 
 * This file contains precomputed shielding factor tables for molecular hydrogen (H2) 
 * and carbon monoxide (CO), based on the study by Lee et al. (1996).
 * 
 * @ref Lee et al. (1996), "Self-shielding and CO photodissociation in molecular clouds", 
 *      Astronomy & Astrophysics, 311, 690.
 *      https://ui.adsabs.harvard.edu/abs/1996A%26A...311..690L/abstract
 */

#ifndef SELF_SHIELDING_FACTOR_HPP_
#define SELF_SHIELDING_FACTOR_HPP_

#include "nicole/nicole_defs.hpp"

namespace nicole {
    // Lee et al. (1996) - Shielding factors for H2 and CO
    // Table 10 and Table 11 from the reference paper

    /**
     * @brief Number of tabulated shielding factors for H2 self-shielding.
     */
    const std::size_t kNumberOfH2ShieldingFactors = 105;

    /**
     * @brief Precomputed H2 column density values (N(H2)) used for shielding calculations.
     */
    extern const Real kH2ColumnDensity[kNumberOfH2ShieldingFactors];

    /**
     * @brief Shielding factors θ[N(H2)] for H2 self-shielding.
     */
    extern const Real kH2ShieldingFactors[kNumberOfH2ShieldingFactors];

    /**
     * @brief Number of tabulated shielding factors for CO self-shielding.
     */
    const std::size_t kNumberOfCOShieldingFactors = 52;

    /**
     * @brief Precomputed CO column density values (N(CO)) used for shielding calculations.
     */
    extern const Real kCOColumnDensity[kNumberOfCOShieldingFactors];

    /**
     * @brief Shielding factors θ[N(CO)] for CO self-shielding.
     */
    extern const Real kCOShieldingFactors[kNumberOfCOShieldingFactors];

    /**
     * @brief Number of tabulated H2 shielding factors used for CO shielding calculations.
     * 
     * This is different from H2 self-shielding because CO is shielded by both CO and H2.
     */
    const std::size_t kNumberOfH2ShieldingFactorsForCOShielding = 43;

    /**
     * @brief Precomputed H2 column density values (N(H2)) used for CO shielding calculations.
     */
    extern const Real kH2ColumnDensityForCOShielding[kNumberOfH2ShieldingFactorsForCOShielding];

    /**
     * @brief Shielding factors θ2[N(H2)] for H2 contribution to CO shielding.
     */
    extern const Real kH2ShieldingFactorsForCOShielding[kNumberOfH2ShieldingFactorsForCOShielding];

    /**
     * @brief Number of tabulated visual extinction (Aν) factors used for CO shielding calculations.
     */
    const std::size_t kNumberOfVisualExtinctionFactorsForCOShielding = 43;

    /**
     * @brief Precomputed visual extinction (Aν) values used for CO shielding calculations.
     */
    extern const Real kVisualExtinctionForCOShielding[kNumberOfVisualExtinctionFactorsForCOShielding];

    /**
     * @brief Shielding factors θ3(Aν) for visual extinction contribution to CO shielding.
     */
    extern const Real kVisualExtinctionFactorsForCOShielding[kNumberOfVisualExtinctionFactorsForCOShielding];

} // namespace nicole

#endif /* SELF_SHIELDING_FACTOR_HPP_ */
