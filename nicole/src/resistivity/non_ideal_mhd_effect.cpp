#include <iostream>

#include "nicole/resistivity/non_ideal_mhd_effect.hpp"

namespace nicole {
    NonIdealMHDeffect::NonIdealMHDeffect(
        SpeciesManager *ptr_species_manager, 
        EnvironmentParameters* ptr_environment_conditions
    ) : ptr_species_manager_(ptr_species_manager),
        ptr_environment_parameters_(ptr_environment_conditions)
    {
        // allocate arrays
        int number_of_species = ptr_species_manager_->total_number_of_species_;
        hall_parameters_.resize(number_of_species);
        ohmic_conductivity_.species_.resize(number_of_species);
        hall_conductivity_.species_.resize(number_of_species);
        pedersen_conductivity_.species_.resize(number_of_species);

        // Set species-specific indices
        id_H_  = ptr_species_manager_->id_H_;
        id_H2_ = ptr_species_manager_->id_H2_;
        id_He_ = ptr_species_manager_->id_He_;

        // Set species-specific masses
        mass_H_   = ptr_species_manager_->GetSpeciesMass(id_H_);
        mass_H2_  = ptr_species_manager_->GetSpeciesMass(id_H2_);
        mass_He_  = ptr_species_manager_->GetSpeciesMass(id_He_);
    }


    void NonIdealMHDeffect::CalculateHallParameters(const Real *species_abundances, const std::size_t number_of_species) {
        if (number_of_species != ptr_species_manager_->total_number_of_species_) {
            std::cout << "Warning: 引数のnumber_of_speciesの値が計算の化学種の数とは異なっています。" << std::endl;
            return;
        }
        const Real gas_number_density = ptr_environment_parameters_->gas_number_density;

        // Calculate mass density of H, H2, He based on species abundances.
        const Real rho_H  = mass_H_  * species_abundances[id_H_] * gas_number_density;
        const Real rho_H2 = mass_H2_ * species_abundances[id_H2_] * gas_number_density;
        const Real rho_He = mass_He_ * species_abundances[id_He_] * gas_number_density;

        // Get environment temperature and calculate auxiliary values
        const Real T = ptr_environment_parameters_->gas_temperature;
        const Real sqrtT = std::sqrt(T);
        const Real logT = std::log10(T);

        // Constants for collision cross-sections
        const Real sigmavH_coef = 2.81e-9 * std::sqrt(kPolarizabilityH);
        const Real sigmavH2_coef = 2.81e-9 * std::sqrt(kPolarizabilityH2);
        const Real sigmavHe_coef = 2.81e-9 * std::sqrt(kPolarizabilityHe);
        const Real delta_g = 1.3;
        const Real sigmavg_coe = delta_g * std::sqrt(128.0 * constants::kBoltzmannConstant * T / (9.0 * M_PI));

        // Magnetic field strength from environment parameters
        const Real magnetic_field = ptr_environment_parameters_->magnetic_field;

        // Calculate Hall parameters
        for (std::size_t ispe = 0; ispe < number_of_species; ++ispe) {
            const auto &species = ptr_species_manager_->species_list_[ispe];

            if (species->GetCharge() == 0.0) continue;  // Skip neutral species (charge == 0)

            const std::string& species_name = species->GetName();
            const Real species_mass = species->GetMass();
            const Real species_charge = species->GetCharge();
            
            // Calculate Hall parameters based on species name and charge
            if (species_name == "e-") { // Electron - collisions with H, H2, He

                // sigmav from Pinto & Galli (2008) table1
                // electron - H2
                Real sigmavH2 = 1.0e-9 * sqrtT * (0.535 + 0.203*logT - 0.136*logT*logT + 0.050*logT*logT*logT);
                Real tauH2inv = (sigmavH2 / (constants::kElectronMass + mass_H2_)) * rho_H2;
                // electron - He
                Real sigmavHe = 1.0e-9 * 0.428 * sqrtT;
                Real tauHeinv = (sigmavHe / (constants::kElectronMass + mass_He_)) * rho_He;
                // electron - H
                Real sigmavH = 1.0e-9 * sqrtT * (2.841 + 0.093*logT + 0.245*logT*logT - 0.089*logT*logT*logT);
                Real tauHinv = (sigmavH / (constants::kElectronMass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                Real tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                Real omega_cyc = species_charge * constants::kChargeUnit * magnetic_field / (constants::kElectronMass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species_name == "H3+") { // H3+ - collisions with H, H2, He

                // H3+ - H2 sigmav from Pinto & Galli (2008) table1
                Real sigmavH2 = 1.0e-9 * (2.693 - 1.238*logT + 0.664*logT*logT - 0.089*logT*logT*logT);
                Real tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // H3+ - He sigmav from Pinto & Galli (2008) Appendix (A.5)
                Real reduced_mass = species_mass * mass_He_ / (species_mass + mass_He_);
                Real sigmavHe = sigmavHe_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // H3+ - H He sigmav from Pinto & Galli (2008) Appendix (A.5)
                reduced_mass = species_mass * mass_H_ / (species_mass + mass_H_);
                Real sigmavH = sigmavH_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                Real tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                Real omega_cyc = species_charge * constants::kChargeUnit * magnetic_field / (species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species_name == "H+") { // H+ - collisions with H, H2, He

                // sigmav from Pinto & Galli (2008) table1
                // H+ - H2 
                Real sigmavH2 = 1.0e-9 * (1.003 + 0.050*logT + 0.136*logT*logT - 0.014*logT*logT*logT);
                Real tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // H+ - He
                Real sigmavHe = 1.0e-9 * (1.424 + 7.438e-6*T - 6.734e-9*T*T);
                Real tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // H+ - H 
                Real sigmavH = 1.0e-9 * 0.649*std::pow(T, 0.375);
                Real tauHinv = (sigmavH / (species_charge + mass_H_)) * rho_H;
                // tau, cyclotron freq
                Real tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                Real omega_cyc = species_charge * constants::kChargeUnit * magnetic_field / (species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species_name == "HCO+") { // HCO+ - collisions with H, H2, He

                // HCO+ - H2 sigmav from Pinto & Galli (2008) table1
                Real sigmavH2 = 1.0e-9 * sqrtT * (1.476 - 1.409*logT + 0.555*logT*logT - 0.0775*logT*logT*logT);
                Real tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // HCO+ - He sigmav from Pinto & Galli (2008) Appendix (A.5)
                Real reduced_mass = species_mass * mass_He_ / (species_mass + mass_He_);
                Real sigmavHe = sigmavHe_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // HCO+ - H sigmav from Pinto & Galli (2008) Appendix (A.5)
                reduced_mass = species_mass * mass_H_ / (species_mass + mass_H_);
                Real sigmavH = sigmavH_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                Real tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                Real omega_cyc = species_charge * constants::kChargeUnit * magnetic_field /(species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else if (species->GetSpeciesType() == SpeciesType::Dust) { // Dust species - collisions with H, H2, He

                // Pinto & Galli (2008) eq.(25)
                Real cross_section = std::dynamic_pointer_cast<DustSpecies>(species)->GetCrossSection();
                // collision with H2
                Real sigmavH2 = cross_section * sigmavg_coe / std::sqrt(mass_H2_);
                Real tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // collision with He
                Real sigmavHe = cross_section * sigmavg_coe / std::sqrt(mass_He_);
                Real tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // collision with H
                Real sigmavH = cross_section * sigmavg_coe / std::sqrt(mass_H_);
                Real tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                Real tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                Real omega_cyc = species_charge * constants::kChargeUnit * magnetic_field /(species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            } else { // Other ions

                if (species_charge < 0.0) continue;

                // other ions sigmav from Pinto & Galli (2008) Appendix (A.5)
                // collision with H2
                Real reduced_mass = species_mass * mass_H2_ / (species_mass + mass_H2_);
                Real sigmavH2 = sigmavH2_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauH2inv = (sigmavH2 / (species_mass + mass_H2_)) * rho_H2;
                // collision with He
                reduced_mass = species_mass * mass_He_ / (species_mass + mass_He_);
                Real sigmavHe = sigmavHe_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauHeinv = (sigmavHe / (species_mass + mass_He_)) * rho_He;
                // collision with H
                reduced_mass = species_mass * mass_H_ / (species_mass + mass_H_);
                Real sigmavH = sigmavH_coef * std::sqrt(species_charge) / std::sqrt(reduced_mass / mass_H_);
                Real tauHinv = (sigmavH / (species_mass + mass_H_)) * rho_H;
                // tau, cyclotron freq
                Real tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
                Real omega_cyc = species_charge * constants::kChargeUnit * magnetic_field /(species_mass * constants::kSpeedOfLight);
                // hall parameter
                hall_parameters_[ispe] = tau * omega_cyc;

            }
        }

        return;
    }


    void NonIdealMHDeffect::CalculateHallParameters(const std::vector<Real>& species_abundances) {
        CalculateHallParameters(species_abundances.data(), species_abundances.size());
        return;
    }


    void NonIdealMHDeffect::CalculateConductivities(const Real *species_abundances, const std::size_t number_of_species) {
        if (number_of_species != ptr_species_manager_->total_number_of_species_) {
            std::cout << "Warning: 引数のnumber_of_speciesの値が計算の化学種の数とは異なっています。" << std::endl;
            return;
        }

        // Initialize conductivities to zero
        ohmic_conductivity_.total_    = 0.0;
        hall_conductivity_.total_     = 0.0;
        pedersen_conductivity_.total_ = 0.0;

        const Real number_density = ptr_environment_parameters_->gas_number_density;
        const Real magnetic_field = ptr_environment_parameters_->magnetic_field;
        const Real ecn_B = constants::kChargeUnit * constants::kSpeedOfLight * number_density / magnetic_field;

        // Calculate individual species conductivities
        for (std::size_t ispe = 0; ispe < number_of_species; ++ispe) {
            const Real species_charge = static_cast<Real>(ptr_species_manager_->GetSpeciesCharge(ispe));

            if (species_charge == 0.0) continue;  // Skip neutral species

            ohmic_conductivity_.species_[ispe] = ecn_B * species_abundances[ispe] * species_charge * hall_parameters_[ispe];
            ohmic_conductivity_.total_ += ohmic_conductivity_.species_[ispe];

            hall_conductivity_.species_[ispe] = - ecn_B * species_abundances[ispe] * species_charge * SQR(hall_parameters_[ispe]) / (1.0 + SQR(hall_parameters_[ispe]));
            hall_conductivity_.total_ += hall_conductivity_.species_[ispe];

            pedersen_conductivity_.species_[ispe] = ecn_B * species_abundances[ispe] * species_charge * hall_parameters_[ispe] / (1.0 + SQR(hall_parameters_[ispe]));
            pedersen_conductivity_.total_ += pedersen_conductivity_.species_[ispe];
        }
    }


    void NonIdealMHDeffect::CalculateConductivities(const std::vector<Real>& species_abundances) {
        CalculateConductivities(species_abundances.data(), species_abundances.size());
        return;
    }


    void NonIdealMHDeffect::CalculateResistivities() {
        const Real sigma_perp2 = hall_conductivity_.total_ * hall_conductivity_.total_ + pedersen_conductivity_.total_ * pedersen_conductivity_.total_;
        const Real coef_res = (constants::kSpeedOfLight * constants::kSpeedOfLight / (4.0 * M_PI));
        resistivity_.ohmic_     = coef_res / ohmic_conductivity_.total_;
        resistivity_.hall_      = coef_res * hall_conductivity_.total_ / sigma_perp2;
        resistivity_.ambipolar_ = coef_res * pedersen_conductivity_.total_ / sigma_perp2 - resistivity_.ohmic_;
    }


    Real NonIdealMHDeffect::GetSpeciesOhmicConductivityByName(const std::string& species_name) const {
        if (ohmic_conductivity_.species_.empty()) return 0.0;
        std::size_t id = ptr_species_manager_->FindSpeciesID(species_name);
        if (id >= ptr_species_manager_->total_number_of_species_) return 0.0;
        return ohmic_conductivity_.species_[id];
    }


    Real NonIdealMHDeffect::GetSpeciesOhmicConductivityByID(const std::size_t id) const {
        if (id >= ptr_species_manager_->total_number_of_species_) return 0.0;
        if (ohmic_conductivity_.species_.empty()) return 0.0;
        return ohmic_conductivity_.species_[id];
    }


    Real NonIdealMHDeffect::GetSpeciesHallConductivityByName(const std::string& species_name) const {
        if (hall_conductivity_.species_.empty()) return 0.0;
        std::size_t id = ptr_species_manager_->FindSpeciesID(species_name);
        if (id >= ptr_species_manager_->total_number_of_species_) return 0.0;
        return ohmic_conductivity_.species_[id];
    }


    Real NonIdealMHDeffect::GetSpeciesHallConductivityByID(const std::size_t id) const {
        if (id >= ptr_species_manager_->total_number_of_species_) return 0.0;
        if (hall_conductivity_.species_.empty()) return 0.0;
        return hall_conductivity_.species_[id];
    }


    Real NonIdealMHDeffect::GetSpeciesPedersenConductivityByName(const std::string& species_name) const {
        if (pedersen_conductivity_.species_.empty()) return 0.0;
        std::size_t id = ptr_species_manager_->FindSpeciesID(species_name);
        if (id >= ptr_species_manager_->total_number_of_species_) return 0.0;
        return pedersen_conductivity_.species_[id];
    }


    Real NonIdealMHDeffect::GetSpeciesPedersenConductivityByID(const std::size_t id) const {
        if (id >= ptr_species_manager_->total_number_of_species_) return 0.0;
        if (pedersen_conductivity_.species_.empty()) return 0.0;
        return pedersen_conductivity_.species_[id];
    }
}
