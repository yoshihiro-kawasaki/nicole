/**
 * @file non_ideal_mhd_effect.cpp
 * @brief Implementation of class NonIdealMHDeffect, calculating non-ideal MHD resistivities and conductivities
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#include "non_ideal_mhd_effect.hpp"

namespace nicole
{
    /**
     * @brief Constructs the NonIdealMHDeffect class with pointers to the species manager and environment parameters
     * 
     * @param ptr_species_manager Pointer to the species manager.
     * @param ptr_environment_conditions Pointer to the environment parameters.
     * @throws std::runtime_error If any pointer is null.
     */
    NonIdealMHDeffect::NonIdealMHDeffect(
        SpeciesManager *ptr_species_manager, 
        EnvironmentParameters* ptr_environment_conditions
    ) : ptr_species_manager_(ptr_species_manager),
        ptr_environment_parameters_(ptr_environment_conditions)
    {
        if (ptr_species_manager_ == nullptr) {
            std::cerr << "Error: ptr_species_manager_ is nullptr" << std::endl;
            throw std::runtime_error("ptr_species_manager_ is nullptr");
        }

        if (ptr_environment_parameters_ == nullptr) {
            std::cerr << "Error: ptr_environment_parameters_ is nullptr" << std::endl;
            throw std::runtime_error("ptr_environment_parameters_ is nullptr");
        }

        // allocate arrays
        int number_of_species = ptr_species_manager_->number_of_total_species_;
        hall_parameters_.resize(number_of_species);
        ohmic_conductivity_.species_.resize(number_of_species);
        hall_conductivity_.species_.resize(number_of_species);
        pedersen_conductivity_.species_.resize(number_of_species);

        // Set species-specific indices
        index_H_  = ptr_species_manager_->index_H_;
        index_H2_ = ptr_species_manager_->index_H2_;
        index_He_ = ptr_species_manager_->index_He_;

        // Set species-specific masses
        mass_H_   = ptr_species_manager_->GetSpeciesMass(index_H_);
        mass_H2_  = ptr_species_manager_->GetSpeciesMass(index_H2_);
        mass_He_  = ptr_species_manager_->GetSpeciesMass(index_He_);

    }

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
    void NonIdealMHDeffect::CalculateHallParameters(const double *species_abundances)
    {
        const std::size_t number_of_species = ptr_species_manager_->number_of_total_species_;
        const double gas_number_density = ptr_environment_parameters_->gas_number_density;

        // Calculate mass density of H, H2, He based on species abundances.
        const double rho_H  = mass_H_  * species_abundances[index_H_] * gas_number_density;
        const double rho_H2 = mass_H2_ * species_abundances[index_H2_] * gas_number_density;
        const double rho_He = mass_He_ * species_abundances[index_He_] * gas_number_density;

        // Get environment temperature and calculate auxiliary values
        const double T = ptr_environment_parameters_->gas_temperature;
        const double sqrtT = std::sqrt(T);
        const double logT = std::log10(T);

        // Constants for collision cross-sections
        const double sigmavH_coef = 2.81e-9 * std::sqrt(kPolarizabilityH);
        const double sigmavH2_coef = 2.81e-9 * std::sqrt(kPolarizabilityH2);
        const double sigmavHe_coef = 2.81e-9 * std::sqrt(kPolarizabilityHe);
        const double delta_g = 1.3;
        const double sigmavg_coe = delta_g * std::sqrt(128.0 * constants::kBoltzmannConstant * T / (9.0 * M_PI));

        // Magnetic field strength from environment parameters
        const double magnetic_field = ptr_environment_parameters_->magnetic_field;

        // Iterate over species and calculate Hall parameters
        for (std::size_t ispe = 0; ispe < number_of_species; ++ispe) {
            const auto &species = ptr_species_manager_->species_list_[ispe];

            // Skip neutral species (charge == 0)
            if (species->GetCharge() == 0.0) continue;

            const std::string& species_name = species->GetName();
            const double species_mass = species->GetMass();
            const double species_charge = species->GetCharge();
            
            // Calculate Hall parameters based on species name and charge
            if (species_name == "e-") { // Electron - collisions with H, H2, He

                // sigmav from Pinto & Galli (2008) table1
                // electron - H2
                double sigmavH2 = 1.0e-9 * sqrtT * (0.535 + 0.203*logT - 0.136*logT*logT + 0.050*logT*logT*logT);
                double tauH2inv = (sigmavH2 / (constants::kElectronMass + mass_H2_)) * rho_H2;
                // electron - He
                double sigmavHe = 1.0e-9 * 0.428 * sqrtT;
                double tauHeinv = (sigmavHe / (constants::kElectronMass + mass_He_)) * rho_He;
                // electron - H
                double sigmavH = 1.0e-9 * sqrtT * (2.841 + 0.093*logT + 0.245*logT*logT - 0.089*logT*logT*logT);
                double tauHinv = (sigmavH / (constants::kElectronMass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                double tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                double omega_cyc = species_charge * constants::kChargeUnit * magnetic_field / (constants::kElectronMass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species_name == "H3+") { // H3+ - collisions with H, H2, He

                // H3+ - H2 sigmav from Pinto & Galli (2008) table1
                double sigmavH2 = 1.0e-9 * (2.693 - 1.238*logT + 0.664*logT*logT - 0.089*logT*logT*logT);
                double tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // H3+ - He sigmav from Pinto & Galli (2008) Appendix (A.5)
                double reduced_mass = species_mass * mass_He_ / (species_mass + mass_He_);
                double sigmavHe = sigmavHe_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // H3+ - H He sigmav from Pinto & Galli (2008) Appendix (A.5)
                reduced_mass = species_mass * mass_H_ / (species_mass + mass_H_);
                double sigmavH = sigmavH_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                double tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                double omega_cyc = species_charge * constants::kChargeUnit * magnetic_field / (species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species_name == "H+") { // H+ - collisions with H, H2, He

                // sigmav from Pinto & Galli (2008) table1
                // H+ - H2 
                double sigmavH2 = 1.0e-9 * (1.003 + 0.050*logT + 0.136*logT*logT - 0.014*logT*logT*logT);
                double tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // H+ - He
                double sigmavHe = 1.0e-9 * (1.424 + 7.438e-6*T - 6.734e-9*T*T);
                double tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // H+ - H 
                double sigmavH = 1.0e-9 * 0.649*std::pow(T, 0.375);
                double tauHinv = (sigmavH / (species_charge + mass_H_)) * rho_H;
                // tau, cyclotron freq
                double tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                double omega_cyc = species_charge * constants::kChargeUnit * magnetic_field / (species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species_name == "HCO+") { // HCO+ - collisions with H, H2, He

                // HCO+ - H2 sigmav from Pinto & Galli (2008) table1
                double sigmavH2 = 1.0e-9 * sqrtT * (1.476 - 1.409*logT + 0.555*logT*logT - 0.0775*logT*logT*logT);
                double tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // HCO+ - He sigmav from Pinto & Galli (2008) Appendix (A.5)
                double reduced_mass = species_mass * mass_He_ / (species_mass + mass_He_);
                double sigmavHe = sigmavHe_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // HCO+ - H sigmav from Pinto & Galli (2008) Appendix (A.5)
                reduced_mass = species_mass * mass_H_ / (species_mass + mass_H_);
                double sigmavH = sigmavH_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                double tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                double omega_cyc = species_charge * constants::kChargeUnit * magnetic_field /(species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species->GetSpeciesType() == SpeciesType::Dust) { // Dust species - collisions with H, H2, He

                // Pinto & Galli (2008) eq.(25)
                double cross_section = std::dynamic_pointer_cast<DustSpecies>(species)->GetCrossSection();
                // collision with H2
                double sigmavH2 = cross_section * sigmavg_coe / std::sqrt(mass_H2_);
                double tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // collision with He
                double sigmavHe = cross_section * sigmavg_coe / std::sqrt(mass_He_);
                double tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // collision with H
                double sigmavH = cross_section * sigmavg_coe / std::sqrt(mass_H_);
                double tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                double tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                double omega_cyc = species_charge * constants::kChargeUnit * magnetic_field /(species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else { // Other ions

                if (species_charge < 0.0) continue;

                // other ions sigmav from Pinto & Galli (2008) Appendix (A.5)
                // collision with H2
                double reduced_mass = species_mass * mass_H2_ / (species_mass + mass_H2_);
                double sigmavH2 = sigmavH2_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // collision with He
                reduced_mass = species_mass * mass_He_ / (species_mass + mass_He_);
                double sigmavHe = sigmavHe_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // collision with H
                reduced_mass = species_mass * mass_H_ / (species_mass + mass_H_);
                double sigmavH = sigmavH_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                double tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                double tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                double omega_cyc = species_charge * constants::kChargeUnit * magnetic_field /(species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            }

        } // end for ispe

        return;
    }

    /**
     * @brief Calculate conductivities (Ohmic, Hall, and Pedersen) based on species abundances.
     * 
     * @param species_abundances Abundance of each species.
     */
    void NonIdealMHDeffect::CalculateConductivities(const double *species_abundances)
    {
        // Initialize conductivities to zero
        ohmic_conductivity_.total_    = 0.0;
        hall_conductivity_.total_     = 0.0;
        pedersen_conductivity_.total_ = 0.0;

        const double number_density = ptr_environment_parameters_->gas_number_density;
        const double magnetic_field = ptr_environment_parameters_->magnetic_field;
        const double ecn_B = constants::kChargeUnit * constants::kSpeedOfLight * number_density / magnetic_field;
        const std::size_t number_of_species = ptr_species_manager_->number_of_total_species_;

        // Calculate individual species conductivities
        for (std::size_t ispe = 0; ispe < number_of_species; ++ispe) {

            const double species_charge = static_cast<double>(ptr_species_manager_->GetSpeciesCharge(ispe));

            // Skip neutral species
            if (species_charge == 0.0) continue;

            ohmic_conductivity_.species_[ispe] = ecn_B * species_abundances[ispe] * species_charge * hall_parameters_[ispe];
            ohmic_conductivity_.total_ += ohmic_conductivity_.species_[ispe];

            hall_conductivity_.species_[ispe] = - ecn_B * species_abundances[ispe] * species_charge * SQR(hall_parameters_[ispe]) / (1.0 + SQR(hall_parameters_[ispe]));
            hall_conductivity_.total_ += hall_conductivity_.species_[ispe];

            pedersen_conductivity_.species_[ispe] = ecn_B * species_abundances[ispe] * species_charge * hall_parameters_[ispe] / (1.0 + SQR(hall_parameters_[ispe]));
            pedersen_conductivity_.total_ += pedersen_conductivity_.species_[ispe];
        }
    }

    /**
     * @brief Calculates the resistivities for Ohmic, Hall, and Ambipolar effects.
     * 
     * This function computes the resistivity for Ohmic dispation, the Hall effect, and 
     * Ambipolar diffusion.
     */
    void NonIdealMHDeffect::CalculateResistivities()
    {
        const double sigma_perp2 = hall_conductivity_.total_ * hall_conductivity_.total_ + pedersen_conductivity_.total_ * pedersen_conductivity_.total_;
        const double coef_res = (constants::kSpeedOfLight * constants::kSpeedOfLight / (4.0 * M_PI));
        resistivity_.ohmic_     = coef_res / ohmic_conductivity_.total_;
        resistivity_.hall_      = coef_res * hall_conductivity_.total_ / sigma_perp2;
        resistivity_.ambipolar_ = coef_res * pedersen_conductivity_.total_ / sigma_perp2 - resistivity_.ohmic_;
    }
}