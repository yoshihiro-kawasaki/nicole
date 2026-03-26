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

#ifndef NON_IDEAL_MHD_EFFECT_HPP
#define NON_IDEAL_MHD_EFFECT_HPP

#include "nicole/nicole_defs.hpp"
#include "nicole/species/species_manager.hpp"

namespace nicole {
    /**
     * @struct Conductivity
     * @brief Stores conductivity information for each species and total conductivity.
     */
    struct Conductivity {
        std::vector<double> species_;   // Conductivity of each species
        double total_;                  // Total conductivity (sum of all species conductivities)
    };

    /**
     * @struct MagneticResistivity
     * @brief Stores resistivity data for Ohmic, Hall, and Ambipolar effects.
     */
    struct MagneticResistivity {
        double ohmic_;      // Ohmic resistivity [cm^2 s^-1]
        double hall_;       // Hall resistivity [cm^2 s^-1]
        double ambipolar_;  // Ambipolar resistivity [cm^ s^-1]
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
    private:
        // Pointer to SpeciesManager
        SpeciesManager* ptr_species_manager_;
        
        // Pointer to environmental parameters
        EnvironmentParameters* ptr_environment_parameters_;
        
        // Hall parameters for species
        std::vector<double> hall_parameters_;
        
        // conductivity data
        Conductivity ohmic_conductivity_;     // Ohmic conductivity data
        Conductivity hall_conductivity_;      // Hall conductivity data
        Conductivity pedersen_conductivity_;  // Pedersen conductivity data;
        
        // Magnetic resistivity data
        MagneticResistivity resistivity_;

        // specific sepcies index
        std::size_t index_H_;
        std::size_t index_H2_;
        std::size_t index_He_;

        // specific species masses
        double mass_H_;
        double mass_H2_;
        double mass_He_;

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
         * @brief Calculates Hall parameters based on species abundances.
         * 
         * @param species_abundances Array of species abundances.
         */
        void CalculateHallParameters(const double *species_abundances);

        /**
         * @brief Calculates the conductivities for Ohmic, Hall, and Pedersen effects.
         * 
         * This function computes the conductivities for the species and sums them to get 
         * the total conductivity for each type (Ohmic, Hall, and Pedersen).
         * 
         * @param species_abundances Array of species abundances.
         */
        void CalculateConductivities(const double *species_abundances);

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
        void GetTotalConductivites(double &ohmic_conductivity, double &hall_conductivity, double &pedersen_conductivity) {
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
        double GetSpeciesOhmicConductivityByName(const std::string& species_name) const;

        /**
         * @brief Get the Ohmic conductivity of a species by its index.
         * @param index The index of the species.
         * @return The Ohmic conductivity of the species, or 0.0 if the index is invalid.
         */
        double GetSpeciesOhmicConductivityByIndex(const std::size_t index) const;

        /**
         * @brief Get the Hall conductivity of a species by its name.
         * @param species_name The name of the species.
         * @return The Hall conductivity of the species, or 0.0 if the species is not found.
         */
        double GetSpeciesHallConductivityByName(const std::string& species_name) const;

        /**
         * @brief Get the Hall conductivity of a species by its index.
         * @param index The index of the species.
         * @return The Hall conductivity of the species, or 0.0 if the index is invalid.
         */
        double GetSpeciesHallConductivityByIndex(const std::size_t index) const;

        /**
         * @brief Get the Pedersen conductivity of a species by its name.
         * @param species_name The name of the species.
         * @return The Pedersen conductivity of the species, or 0.0 if the species is not found.
         */
        double GetSpeciesPedersenConductivityByName(const std::string& species_name) const;

        /**
         * @brief Get the Pedersen conductivity of a species by its index.
         * @param index The index of the species.
         * @return The Pedersen conductivity of the species, or 0.0 if the index is invalid.
         */
        double GetSpeciesPedersenConductivityByIndex(const std::size_t index) const;

        /**
         * @brief Retrieves the resistivity values.
         * 
         * @param ohmic_resistivity The computed Ohmic resistivity.
         * @param hall_resistivity The computed Hall resistivity.
         * @param ambipolar_resistivity The computed Ambipolar resistivity.
         */
        void GetResistivite(double &ohmic_resistivity, double &hall_resistivity, double &ambipolar_resistivity) {
            ohmic_resistivity     = resistivity_.ohmic_;
            hall_resistivity      = resistivity_.hall_;
            ambipolar_resistivity = resistivity_.ambipolar_;
            return;
        }
    };
}

#endif /* NON_IDEAL_MHD_EFFECT_HPP */