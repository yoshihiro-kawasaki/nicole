#ifndef REACTION_SIMULATOR_HPP_
#define REACTION_SIMULATOR_HPP_

#include <vector>

#include "nicole/nicole_defs.hpp"
#include "nicole/species/species_manager.hpp"
#include "nicole/reaction/reaction_manager.hpp"
#include "odepack_cpp/odepack.hpp"

namespace nicole {
    /**
     * @struct LsodesParameters
     * @brief Parameters for the LSODES solver.
     */
    struct LsodesParameters {
        int liw;
        int lrw;
        std::vector<Real> rwork;
        std::vector<int> iwork;
        int itol;
        std::vector<Real> rtol;
        std::vector<Real> atol;
        int itask;
        int iopt;
        int mf;
        int istate;
        bool is_allocate_arrays;
    };

    /**
     * @struct LsodeParameters
     * @brief Parameters for the LSODE solver (alternative to LSODES).
     */
    struct LsodeParameters {
        int liw;
        int lrw;
        std::vector<Real> rwork;
        std::vector<int> iwork;
        int itol;
        std::vector<Real> rtol;
        std::vector<Real> atol;
        int itask;
        int iopt;
        int mf;
        int istate;
        bool is_allocate_arrays;
    };

    /**
     * @class ReactionSimulator
     * @brief A class responsible for simulating and integrating chemical reactions.
     * 
     * This class is used to calculate reaction rates and integrate the abundances of species.
     */
    class ReactionSimulator {
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

        ~ReactionSimulator();

        void CheckCalculationResult(const Real *species_abundances) const;
        void CheckReactionRateCoefficient(const std::string filename);

        /**
         * @brief Set the initial abundances for each species.
         * 
         * @param species_abundances Array of species abundances.
         */
        void SetInitialSpeciesAbundances(Real* species_abundances);

        /**
         * @brief Integrate species abundances over time using LSODE or LSODES.
         * 
         * @param t Current time.
         * @param tout Time to integrate to.
         * @param species_abundances Array to store the integrated species abundances.
         * @return true if integration was successful, false otherwise.
         */
        bool Integrate(Real &t, const Real tout, Real *species_abundances);

        /**
         * @brief Integrate species abundances over time and write output to a file.
         * 
         * @param t Current time.
         * @param tout Time to integrate to.
         * @param species_abundances Array to store the integrated species abundances.
         * @param file Output file stream to write results.
         * @return true if integration was successful, false otherwise.
         */
        bool Integrate(Real &t, const Real tout, Real *species_abundances, std::ofstream& file);

        /**
         * @brief Calculate the rate coefficients for reactions.
         * 
         * This function calculate the rate coefficient independent on species abundances.
         */
        void CalculateRateCoefficient();

    private:
        void ReadAbundancesFile(const std::string& filename);

        // Gas phase reaction rate coefficient calculation functions
        void CalculateGasPhaseReactionRateCoefficient();
        using CalculateGasReactionRateFunction = Real(ReactionSimulator::*)(const std::shared_ptr<Reaction> reaction, const Real temperature);
        Real CalculateGasPhaseModifiedArrhenius(const std::shared_ptr<Reaction> reaction, const Real temperature);
        Real CalculateGasPhaseIonpol1(const std::shared_ptr<Reaction> reaction, const Real temperature);
        Real CalculateGasPhaseIonpol2(const std::shared_ptr<Reaction> reaction, const Real temperature);

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
        void CalculateSpeciesAbundancesDependentRateCoefficient(const Real *species_abundances);
        void CalculateDustSurfaceAndMantleSwappingRateCoefficient(const Real *species_abundances);

        /**
         * @brief Check if the Jacobian is sparse based on a given threshold.
         * 
         * @param threshold Threshold for sparsity check.
         * @return true if the Jacobian is sparse, false otherwise.
         */
        bool IsSparseJacobian(const Real threshold);

        // Allocate and set work arrays for LSODE or LSODES
        void AllocateAndSetLsodeArrays();
        void AllocateAndSetLsodesArrays();

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
        static void OrdinaryDifferentialEquation(int neq, Real t, Real *y, Real *ydot, void *user_data);

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
        static void JacobianJth(int neq, Real t, Real *y, int j, int *ian, int *jan, Real *pdj, void *user_data);
        
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
        static void Jacobian(int neq, Real t, Real *y, int ml, int mu, Real *pd, int nrowpd, void *user_data);

        SpeciesManager* ptr_species_manager_;
        ReactionManager* ptr_reaction_manager_;
        EnvironmentParameters *ptr_environment_parameters_;

        std::vector<Real> initial_species_abundances_;
        std::vector<Real> reaction_rate_coefficient_;

        bool is_lsode_integrator_;  // Flag for LSODE integrator, if true use LSODE, use LDODES otherwise.
        LsodeParameters lsode_parameters_;
        LsodesParameters lsodes_parameters_;

        Real relative_tolerance_;
        Real absolute_tolerance_;

        // Variables for dust surface/mantle related reactions
        Real total_desorption_rate_;
        Real total_accretion_rate_;
        Real total_abundances_of_dust_surface_species_;
        Real total_abundances_of_dust_mantle_species_;
        Real number_of_surface_layers_;
        Real number_of_mantle_layers_;
        Real number_of_total_layers_;
        Real coverage_of_H2O_on_dust_surface_;
        Real coverage_of_silicate_on_dust_surface_;
    };
} // namespace nicole

#endif /* REACTION_SIMULATOR_HPP */
