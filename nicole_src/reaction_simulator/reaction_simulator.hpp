/**
 * @file reaction_simulator.hpp
 * @brief Header file for the ReactionSimulator class, handling the integration and simulation of chemical reactions.
 * @date 2025-02-12
 */

#ifndef REACTION_SIMULATOR_HPP
#define REACTION_SIMULATOR_HPP

#include "../nicole_defs.hpp"
#include "../species/species_manager.hpp"
#include "../reactions/reaction_manager.hpp"
#include "../shielding/self_shielding_factor.hpp"
#include "../odepack_cpp/odepack.hpp"

namespace nicole
{
    /**
     * @struct LsodesParameters
     * @brief Parameters for the LSODES solver.
     */
    struct LsodesParameters
    {
        int liw;
        int lrw;
        double *rwork;
        int *iwork;
        int itol;
        double *rtol;
        double *atol;
        int itask;
        int iopt;
        int mf;
        int istate;
        bool is_allocate_work_arrays;
    };

    /**
     * @struct LsodeParameters
     * @brief Parameters for the LSODE solver (alternative to LSODES).
     */
    struct LsodeParameters
    {
        int liw;
        int lrw;
        double *rwork;
        int *iwork;
        bool is_allocate_work_arrays;
        int itol;
        double *rtol;
        double *atol;
        int itask;
        int iopt;
        int mf;
        int istate;
    };

    /**
     * @class ReactionSimulator
     * @brief A class responsible for simulating and integrating chemical reactions.
     * 
     * This class is used to calculate reaction rates and integrate the abundances of species.
     */
    class ReactionSimulator
    {
    public:

        /**
         * @brief Constructor for initializing the ReactionSimulator.
         * 
         * @param ptr_species_manager Pointer to the species manager.
         * @param ptr_reaction_manager Pointer to the reaction manager.
         * @param ptr_environment_parameters Pointer to environment parameters.
         * @param input Configuration input for the simulation.
         */
        ReactionSimulator(
            SpeciesManager* ptr_species_manager,
            ReactionManager* ptr_reaction_manager,
            EnvironmentParameters* ptr_environment_parameters,
            InputConfig& input
        );

        /**
         * @brief Destructor for the ReactionSimulator.
         */
        ~ReactionSimulator();

        /**
         * @brief Calculate the rate coefficients for reactions.
         * 
         * This function calculate the rate coefficient independent on species abundances.
         */
        void CalculateRateCoefficient();

        /**
         * @brief Check the reaction rate coefficients and write results to file.
         * 
         * @param filename The name of the file to write the results to.
         */
        void CheckReactionRateCoefficient(const std::string filename);

        /**
         * @brief Set the initial abundances for each species.
         * 
         * @param species_abundances Array of species abundances.
         */
        void SetInitialSpeciesAbundances(double* species_abundances);

        /**
         * @brief Integrate species abundances over time using LSODE or LSODES.
         * 
         * @param t Current time.
         * @param tout Time to integrate to.
         * @param species_abundances Array to store the integrated species abundances.
         * @return true if integration was successful, false otherwise.
         */
        bool Integrate(double &t, const double tout, double *species_abundances);

        /**
         * @brief Integrate species abundances over time and write output to a file.
         * 
         * @param t Current time.
         * @param tout Time to integrate to.
         * @param species_abundances Array to store the integrated species abundances.
         * @param file Output file stream to write results.
         * @return true if integration was successful, false otherwise.
         */
        bool Integrate(double &t, const double tout, double *species_abundances, std::ofstream& file);

        /**
         * @brief Check the calculation result by evaluating total charge, dust number density, and other quantities.
         * 
         * @param species_abundances Array of species abundances to check.
         */
        void CheckCalculationResult(const double *species_abundances) const;

    private:

        // Pointer to the species manager
        SpeciesManager* ptr_species_manager_;

        // Pointer to the reaction manager
        ReactionManager* ptr_reaction_manager_;

        // Pointer to the environment parameters
        EnvironmentParameters *ptr_environment_parameters_;

        // Initial abundances of species
        std::vector<double> initial_species_abundances_;

        // List of reaction rate coefficients
        std::vector<double> reaction_rate_coefficient_;

        // Flag for LSODE integrator, if true use LSODE, use LDODES otherwise.
        bool is_lsode_integrator_;

        // LSODE parameters
        LsodeParameters lsode_parameters_;

        // LSODES parameters
        LsodesParameters lsodes_parameters_;

        // Tolerance values for numerical integration
        double relative_tolerance_;
        double absolute_tolerance_;

        // Variables for dust surface/mantle related reactions
        double total_desorption_rate_;
        double total_accretion_rate_;
        double total_abundances_of_dust_surface_species_;
        double total_abundances_of_dust_mantle_species_;
        double number_of_surface_layers_;
        double number_of_mantle_layers_;
        double number_of_total_layers_;
        double coverage_of_H2O_on_dust_surface_;
        double coverage_of_silicate_on_dust_surface_;

        /**
         * @brief Read species abundances from a file.
         * 
         * @param filename The name of the file containing species abundances data.
         */
        void ReadAbundancesFile(const std::string& filename);

        // Reaction rate coefficient calculation functions
        void CalculateGasPhaseReactionRateCoefficient();
        using CalculateGasReactionRateFunction = double(ReactionSimulator::*)(const std::shared_ptr<Reaction> reaction, const double temperature);
        double CalculateGasPhaseModifiedArrhenius(const std::shared_ptr<Reaction> reaction, const double temperature);
        double CalculateGasPhaseIonpol1(const std::shared_ptr<Reaction> reaction, const double temperature);
        double CalculateGasPhaseIonpol2(const std::shared_ptr<Reaction> reaction, const double temperature);

        // Dust related reaction rate coefficient calculation functions
        void CalculateDustAndChargedParticleCollisionRateCoefficient();
        void CalculateDustCollisionRateCoefficient();
        void CalculateNeutralSpeciesAccretionOnDustSurfacesRateCoefficient();
        void CalculateThermalDesorptionOnDustSurfacesRateCoefficient();
        void CalculateCosmicRayDesorptionOnDustSurfacesRateCoefficient();
        void CalculatePhotoDesorptionByExternalUVRateCoefficient();
        void CalculatePhotoDesorptionByCosmicRayGeneratedUVRateCoefficient();
        void CalculateDustSurfaceReactionRateCoefficient();
        void CalculateDustMantleReactionRateCoefficient();
        void CalculateSpeciesAbundancesDependentRateCoefficient(const double *species_abundances);
        void CalculateDustSurfaceAndMantleSwappingRateCoefficient(const double *species_abundances);

        /**
         * @brief Check if the Jacobian is sparse based on a given threshold.
         * 
         * @param threshold Threshold for sparsity check.
         * @return true if the Jacobian is sparse, false otherwise.
         */
        bool IsSparseJacobian(const double threshold);

        // Allocate and set work arrays for LSODE or LSODES
        void AllocateAndSetLsodeWorkArrays();
        void AllocateAndSetLsodesWorkArrays();

        // Reset LSODE or LSODES work arrays.
        void ResetLsodeWorkArrays();
        void ResetLsodesWorkArrays();

        /**
         * @brief Ordinary Differential Equation (ODE) function for LSODE or LSODES integration.
         * 
         * @param neq Number of equations.
         * @param t Current time.
         * @param y State vector of species abundances.
         * @param ydot Derivative of the state vector.
         * @param user_data User-defined data passed to the function.
         */
        static void OrdinaryDifferentialEquation(const int neq, const double t, const double *y, double *ydot, void *user_data);

        /**
         * @brief Jacobian function for LSODES integration.
         * 
         * @param neq Number of equations.
         * @param t Current time.
         * @param y State vector of species abundances.
         * @param j Index of the Jacobian to compute.
         * @param ian, jan Index arrays for the Jacobian.
         * @param pdj Jacobian matrix values.
         * @param user_data User-defined data passed to the function.
         */
        static void JacobianJth(const int neq, const double t, const double *y, const int j, int *ian, int *jan, double *pdj, void *user_data);
        
        /**
         * @brief Full Jacobian matrix function for LSODE integration.
         * 
         * @param neq Number of equations.
         * @param t Current time.
         * @param y State vector of species abundances.
         * @param ml, mu Number of subdiagonal and superdiagonal elements.
         * @param pd Jacobian matrix.
         * @param nrowpd Number of rows in the Jacobian.
         * @param user_data User-defined data passed to the function.
         */
        static void Jacobian(const int neq, const double t, const double *y, const int ml, const int mu, double *pd, const int nrowpd, void *user_data);
    };
    
} // namespace nicole

#endif /* REACTION_SIMULATOR_HPP */