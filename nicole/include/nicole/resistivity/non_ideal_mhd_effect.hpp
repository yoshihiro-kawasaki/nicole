/**
 * @file non_ideal_mhd_effect.hpp
 * @brief Calculate non-ideal MHD resistivities and conductivities.
 * 
 * This class computes various MHD resistivities and conductivities, such as 
 * Ohmic, Hall, and Pedersen conductivities, along with the corresponding 
 * resistivities for non-ideal MHD effects.
 * 
 * @ref 
 * https://ui.adsabs.harvard.edu/abs/2002ApJ...573..199N/abstract
 * https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.5588K/abstract
 */
#ifndef NON_IDEAL_MHD_EFFECT_HPP_
#define NON_IDEAL_MHD_EFFECT_HPP_

#include "nicole/nicole_defs.hpp"
#include "nicole/species/species_manager.hpp"

namespace nicole {
    /**
     * @struct Conductivity
     * @brief Stores conductivity information for each species and total conductivity.
     */
    struct Conductivity {
        std::vector<Real> species_;   // Conductivity of each species
        Real total_;                  // Total conductivity (sum of all species conductivities)
    };

    /**
     * @struct MagneticResistivity
     * @brief Stores resistivity data for Ohmic, Hall, and Ambipolar effects.
     */
    struct MagneticResistivity {
        Real ohmic_;      // Ohmic resistivity [cm^2 s^-1]
        Real hall_;       // Hall resistivity [cm^2 s^-1]
        Real ambipolar_;  // Ambipolar resistivity [cm^ s^-1]
    };

    /**
     * @class NonIdealMHDeffect
     * @brief Computes non-ideal MHD resistivities and conductivities.
     * 
     * This class calculates resistivities (Ohmic, Hall, and Ambipolar) and 
     * conductivities (Ohmic, Hall, and Pedersen) for a plasma, taking into account 
     * various species and environmental parameters.
     */
    class NonIdealMHDeffect {
    public:

        /**
         * @brief Constructor for NonIdealMHDeffect class.
         * 
         * Initializes the non-ideal MHD effect calculations with species manager 
         * and environmental parameters.
         * 
         * @param ptr_species_manager Pointer to the SpeciesManager that manages species data.
         * @param ptr_environment_conditions Pointer to the environmental parameters.
         */
        NonIdealMHDeffect(SpeciesManager *ptr_species_manager, EnvironmentParameters* ptr_environment_conditions);

        /**
         * @brief Calculate the Hall parameters for each species based on their abundances.
         * 
         * The Hall parameter is the ratio of the drag force acting on the species due to neutral particles 
         * (such as H2, H, He) to the Lorentz force. 
         * - If the Hall parameter is greater than 1, the Lorentz force dominates.
         * - If the Hall parameter is less than 1, the drag force due to neutrals dominates.
         * 
         * @param species_abundances The abundances of each species.
         */
        void CalculateHallParameters(const Real *species_abundances);

        /**
         * @brief Calculates the conductivities for Ohmic, Hall, and Pedersen effects.
         * 
         * This function computes the conductivities for the species and sums them to get 
         * the total conductivity for each type (Ohmic, Hall, and Pedersen).
         * 
         * @param species_abundances Array of species abundances.
         */
        void CalculateConductivities(const Real *species_abundances);

        /**
         * @brief Calculates the resistivities for Ohmic, Hall, and Ambipolar effects.
         * 
         * This function computes the resistivity for Ohmic dispation, the Hall effect, and 
         * Ambipolar diffusion.
         */
        void CalculateResistivities();

        /**
         * @brief Retrieves the total conductivities.
         * 
         * @param ohmic_conductivity The computed Ohmic conductivity.
         * @param hall_conductivity The computed Hall conductivity.
         * @param pedersen_conductivity The computed Pedersen conductivity.
         */
        void GetTotalConductivites(Real &ohmic_conductivity, Real &hall_conductivity, Real &pedersen_conductivity) const {
            ohmic_conductivity    = ohmic_conductivity_.total_;
            hall_conductivity     = hall_conductivity_.total_;
            pedersen_conductivity = pedersen_conductivity_.total_;
            return;
        }

        /**
         * @brief Get the Ohmic conductivity of a species by its name.
         * @param species_name The name of the species.
         * @return The Ohmic conductivity of the species, or 0.0 if the species is not found.
         */
        Real GetSpeciesOhmicConductivityByName(const std::string& species_name) const;

        /**
         * @brief Get the Ohmic conductivity of a species by its id.
         * @param index The index of the species.
         * @return The Ohmic conductivity of the species, or 0.0 if the index is invalid.
         */
        Real GetSpeciesOhmicConductivityByID(const std::size_t id) const;

        /**
         * @brief Get the Hall conductivity of a species by its name.
         * @param species_name The name of the species.
         * @return The Hall conductivity of the species, or 0.0 if the species is not found.
         */
        Real GetSpeciesHallConductivityByName(const std::string& species_name) const;

        /**
         * @brief Get the Hall conductivity of a species by its id.
         * @param index The index of the species.
         * @return The Hall conductivity of the species, or 0.0 if the index is invalid.
         */
        Real GetSpeciesHallConductivityByID(const std::size_t id) const;

        /**
         * @brief Get the Pedersen conductivity of a species by its name.
         * @param species_name The name of the species.
         * @return The Pedersen conductivity of the species, or 0.0 if the species is not found.
         */
        Real GetSpeciesPedersenConductivityByName(const std::string& species_name) const;

        /**
         * @brief Get the Pedersen conductivity of a species by its id.
         * @param index The index of the species.
         * @return The Pedersen conductivity of the species, or 0.0 if the index is invalid.
         */
        Real GetSpeciesPedersenConductivityByID(const std::size_t id) const;

        /**
         * @brief Get the resistivity values.
         */
        void GetResistivite(Real &ohmic_resistivity, Real &hall_resistivity, Real &ambipolar_resistivity) const {
            ohmic_resistivity     = resistivity_.ohmic_;
            hall_resistivity      = resistivity_.hall_;
            ambipolar_resistivity = resistivity_.ambipolar_;
            return;
        }

    private:
        SpeciesManager* ptr_species_manager_;
        EnvironmentParameters* ptr_environment_parameters_;
        
        std::vector<Real> hall_parameters_;
        
        Conductivity ohmic_conductivity_;
        Conductivity hall_conductivity_;
        Conductivity pedersen_conductivity_;
        
        MagneticResistivity resistivity_;

        // specific sepcies id
        std::size_t id_H_;
        std::size_t id_H2_;
        std::size_t id_He_;

        // specific species masses
        Real mass_H_;
        Real mass_H2_;
        Real mass_He_;
    };
}

#endif /* NON_IDEAL_MHD_EFFECT_HPP_ */
