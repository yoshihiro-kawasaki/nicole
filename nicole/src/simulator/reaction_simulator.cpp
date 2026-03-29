#include <fstream>

#include "odepack_cpp/odepack.hpp"

#include "nicole/shielding/self_shielding_factor.hpp"
#include "nicole/simulator/reaction_simulator.hpp"
#include "nicole/utils/string_utils.hpp"
#include "nicole/utils/vector_utils.hpp"

namespace nicole {
    ReactionSimulator::ReactionSimulator(
        SpeciesManager *ptr_species_manager, 
        ReactionManager*ptr_reaction_manager, 
        EnvironmentParameters *ptr_environment_parameters,
        InputConfig& config
    ) : ptr_species_manager_(ptr_species_manager),
        ptr_reaction_manager_(ptr_reaction_manager),
        ptr_environment_parameters_(ptr_environment_parameters),
        relative_tolerance_(config.GetReal("reltol")),
        absolute_tolerance_(config.GetReal("abstol"))
    {
        if (relative_tolerance_ <= 0.0) {
            std::cout << "Warning: relative tolerance <= 0.0, set relative tolerance=1e-4." << std::endl;
            relative_tolerance_ = 1e-4;
        }

        if (absolute_tolerance_ <= 0.0) {
            std::cout << "Warning: absolute_tolerance <= 0.0, set absolute_tolerance=1e-15." << std::endl;
            absolute_tolerance_ = 1e-15;
        }

        ReadAbundancesFile(config.GetString("abundances_file"));

        reaction_rate_coefficient_.resize(ptr_reaction_manager->total_number_of_reactions_);

        // Initialize LSODE/LSODES parameters
        Real threshold_specre = 0.15;
        bool is_sparce = IsSparseJacobian(threshold_specre);
        if (is_sparce) {
            AllocateAndSetLsodesArrays();
            is_lsode_integrator_ = false;
            std::cout << "Using LSODES integrator (sparse Jacobian)" << std::endl;
        } else {
            AllocateAndSetLsodeArrays();
            is_lsode_integrator_ = true;
            std::cout << "Using LSODE integrator" << std::endl;
        }
    }


    ReactionSimulator::~ReactionSimulator() { }


    void ReactionSimulator::CheckCalculationResult(const Real *species_abundances, const std::size_t number_of_species) const {
        if (number_of_species != ptr_species_manager_->total_number_of_species_) {
            std::cout << "Warning: 引数のnumber_of_speciesの値が計算の化学種の数とは異なっています。" << std::endl;
            return;
        }
        // const std::size_t number_of_total_species = ptr_species_manager_->total_number_of_species_;
        const Real gas_number_density = ptr_environment_parameters_->gas_number_density;

        Real total_charge = 0.0;
        Real cation_density = 0.0;
        Real anion_density = 0.0;
        Real total_dust_number_density_result = 0.0;
        Real mean_dust_charge = 0.0;

        // Compute charge
        for (std::size_t i = 0; i < number_of_species; ++i) {
            const auto& species = ptr_species_manager_->species_list_[i];
            Real charge = static_cast<Real>(species->GetCharge());

            total_charge += charge * species_abundances[i];
            if (charge > 0.0) {
                cation_density += charge * species_abundances[i];
            } else if (charge < 0.0) {
                anion_density += charge * species_abundances[i];
            }

            if (species->GetSpeciesType() == SpeciesType::Dust) {
                total_dust_number_density_result += species_abundances[i];
                mean_dust_charge += charge * species_abundances[i];
            }
        }

        // Calculate mean dust charge.
        if (total_dust_number_density_result == 0.0) {
            mean_dust_charge = 0.0; // Avoid division by zero if no dust particles exist
        } else {
            mean_dust_charge /= total_dust_number_density_result;
        }

        // Compute the error in dust number density
        total_dust_number_density_result *= gas_number_density;
        Real dust_number_density = ptr_species_manager_->dust_species_model_parameters_.GetDustTotalAbundance() * gas_number_density;
        Real err_dust = std::abs((total_dust_number_density_result - dust_number_density) / dust_number_density);

        // Print the results of the calculation
        for (std::size_t index = 0; index < number_of_species; ++index) {
            std::cout << "x[" << std::setw(12) << ptr_species_manager_->GetSpeciesName(index) << "] = " 
                      << std::scientific << std::setw(15) << species_abundances[index] << std::endl;
        }
        std::cout << "initial total dust number density = " << std::setw(15) << dust_number_density << std::endl;
        std::cout << "finish total dust number density  = " << std::setw(15) << total_dust_number_density_result << std::endl;
        std::cout << "error dust number density         = " << std::setw(15) << err_dust << std::endl;
        std::cout << "mean dust charge                  = " << std::setw(15) << mean_dust_charge << std::endl;
        std::cout << "total charge                      = " << std::setw(15) << total_charge << std::endl;
        std::cout << "total cation number density       = " << std::setw(15) << cation_density << std::endl;
        std::cout << "total antion number density       = " << std::setw(15) << anion_density << std::endl;
        std::cout << "electron number density           = " << std::setw(15) << species_abundances[ptr_species_manager_->id_electron_] << std::endl;
    }


    void ReactionSimulator::CheckCalculationResult(const std::vector<Real>& species_abundances) const {
        CheckCalculationResult(species_abundances.data(), species_abundances.size());
        return;
    }


    void ReactionSimulator::CheckReactionRateCoefficient(const std::string filename) {
        std::ofstream file(filename, std::ios::out | std::ios::trunc);
        if (!file.is_open()) {
            return;
        }

        const int n = ptr_species_manager_->total_number_of_species_;
        std::vector<Real> y(n), ydot(n);
        for (int i = 0; i < n; ++i) {
            y[i] = 1.0e-5;
            ydot[i] = 0.0;
        }
        SetInitialSpeciesAbundances(y);
        CalculateRateCoefficient();
        OrdinaryDifferentialEquation(n, 0.0, y.data(), ydot.data(), this);

        std::size_t number_of_total_reaction = reaction_rate_coefficient_.size();
        file << std::scientific;
        for (std::size_t ireac = 0; ireac < number_of_total_reaction; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            file << reaction->type_id_ << " " << reaction_rate_coefficient_[ireac] << std::endl;
        }

        file.close();
    }


    void ReactionSimulator::SetInitialSpeciesAbundances(Real* species_abundances, const std::size_t number_of_species) {
        if (number_of_species != ptr_species_manager_->total_number_of_species_) {
            std::cout << "Warning: 引数のnumber_of_speciesの値が計算の化学種の数とは異なっています。" << std::endl;
            return;
        }

        // Set initialize species abundances. See also ReadAbundancesFile.
        for (std::size_t id = 0; id < number_of_species; ++id) {
            species_abundances[id] = initial_species_abundances_[id];
        }

        // Calculate the electron abundance based on ionization equilibrium.
        Real electron_abundance = 0.0;
        for (std::size_t id = 0; id < number_of_species; ++id) {
            const auto& species = ptr_species_manager_->species_list_[id];
            const int charge = species->GetCharge();

            // Sum the charge contributions from all charged species to calculate the electron abundance
            if (id != ptr_species_manager_->id_electron_ && charge != 0.0) {
                electron_abundance += (charge * species_abundances[id]);
            }
        }

        // Set initial dust species abundances
        if (ptr_species_manager_->dust_species_model_parameters_.IsSizeDistributionModel()) {
            // If using a size distribution model, set the abundances for each dust size bin
            for (std::size_t id = 0; id < number_of_species; ++id) {
                const auto& species = std::dynamic_pointer_cast<DustSpecies>(ptr_species_manager_->species_list_[id]);
                if (!species) continue;
                const int charge = species->GetCharge();
                const std::size_t bin_number = species->GetBinNumber();

                // Set abundances for electrically neutral dust species
                if (species->GetSpeciesType() == SpeciesType::Dust && charge == 0.0) {
                    species_abundances[id] = ptr_species_manager_->dust_species_model_parameters_.GetDustAbundancesForBin(bin_number);
                }
            }
        } else {
            // If using a single size model, set the total dust abundance
            for (std::size_t id = 0; id < number_of_species; ++id) {
                const auto& species = ptr_species_manager_->species_list_[id];
                const int charge = species->GetCharge();
                if (species->GetSpeciesType() == SpeciesType::Dust && charge == 0.0) {
                    species_abundances[id] = ptr_species_manager_->dust_species_model_parameters_.GetDustTotalAbundance();
                }
            }
        }

        // Set the electron abundance
        species_abundances[ptr_species_manager_->id_electron_] = electron_abundance;

        for (std::size_t id = 0; id < number_of_species; ++id) {
            if (species_abundances[id] < kMinimumSpeciesAbundance) {
                species_abundances[id] = kMinimumSpeciesAbundance;
            }
        }
    }


    void ReactionSimulator::SetInitialSpeciesAbundances(std::vector<Real>& species_abundances) {
        SetInitialSpeciesAbundances(species_abundances.data(), species_abundances.size());
        return;
    }


    void ReactionSimulator::ReadAbundancesFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::in);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        std::size_t number_of_species = ptr_species_manager_->total_number_of_species_;
        initial_species_abundances_.resize(number_of_species, 0.0);

        std::string line;
        int line_number = 0;
        while (std::getline(file, line)) {
            line_number++;

            // Skip empty lines or comment lines
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;

            std::vector<std::string> split_result = string_utils::Split(line, ' ', true);
            
            std::string species_name = split_result[0];

            std::size_t id = ptr_species_manager_->FindSpeciesID(species_name);
            if (id == kNotFoundSpecies) continue;

            try {
                initial_species_abundances_[id] = std::stod(split_result[2]);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid abundance value for species " << species_name 
                        << " at line " << line_number << ": " << split_result[2] << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Error: Abundance value out of range for species " << species_name
                        << " at line " << line_number << ": " << split_result[2] << std::endl;
            }
        }
    }


    /**
     * @brief Calculate the gas phase reaction rate coefficients for all relevant reactions.
     * 
     * This function calculates the reaction rate coefficients for gas phase reactions based on the temperature
     * and the specific reaction rate formula used for each reaction (Modified Arrhenius, ionpol1, ionpol2).
     * 
     * This calculation method(code) is based on the code Nahoon_kida.uva.2014.
     * 
     * @ref
     * https://kida.astrochem-tools.org/codes.html
     * https://kida.astrochem-tools.org/help.html
     * https://ui.adsabs.harvard.edu/abs/2015ApJS..217...20W/abstract
     * https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..63W/abstract
     */
    void ReactionSimulator::CalculateGasPhaseReactionRateCoefficient() {
        const Real temperature = ptr_environment_parameters_->gas_temperature;

        // Arrays for storing reaction indices and distance from temperature limits.
        std::vector<std::size_t> indice(10);
        std::vector<Real> distmin(10), distmax(10);
        for (std::size_t i = 0; i < 10; ++i) {
            indice[i]  = 0;
            distmin[i] = 9999.0;
            distmax[i] = 9999.0;
        }
        std::size_t j = 0;

        const std::size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kGasPhase4];
        const std::size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kGasPhase8];

        // Function pointer for reaction rate calculation
        CalculateGasReactionRateFunction calc_rate_func = nullptr;

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];

            // Select the appropriate rate calculation function based on the formula ID
            const int formula_id = static_cast<int>(reaction->rate_parameters_[gas_phase_reaction_params::kFormulaID]);
            if (formula_id == 3) { // Modified Arrhenius
                calc_rate_func = &ReactionSimulator::CalculateGasPhaseModifiedArrhenius;
            } else if (formula_id == 4) { // ionpol1
                calc_rate_func = &ReactionSimulator::CalculateGasPhaseIonpol1;
            } else if (formula_id == 5) { // ionpol2
                calc_rate_func = &ReactionSimulator::CalculateGasPhaseIonpol2;
            }

            // Check temperature limits and calculate rate coefficient
            const Real temperature_lower_limit = reaction->rate_parameters_[gas_phase_reaction_params::kTemperatureLowerLimit];
            const Real temperature_upper_limit = reaction->rate_parameters_[gas_phase_reaction_params::kTemperatureUpperLimit];
            if (temperature < temperature_lower_limit) {
                reaction_rate_coefficient_[ireac] = (this->*calc_rate_func)(reaction, temperature_lower_limit);
            } else if (temperature > temperature_upper_limit) {
                reaction_rate_coefficient_[ireac] = (this->*calc_rate_func)(reaction, temperature_upper_limit);
            } else {
                reaction_rate_coefficient_[ireac] = (this->*calc_rate_func)(reaction, temperature);
            }

            // Handle reactions with multiple rate coefficients
            int gas_phase_id = static_cast<int>(reaction->rate_parameters_[gas_phase_reaction_params::kID]);
            int gas_phase_id_next;
            if (ireac + 1 <= ireac_end) {
                gas_phase_id_next = static_cast<int>(ptr_reaction_manager_->reaction_list_[ireac+1]->rate_parameters_[gas_phase_reaction_params::kID]);
            } else {
                break;
            }

            // Check for the presence of several rate coefficients present in the network for the same reaction
            if (gas_phase_id == gas_phase_id_next) {
                indice[j]  = ireac;
                distmin[j] = temperature_lower_limit - temperature;
                distmax[j] = temperature - temperature_upper_limit;
                j++;
            }

            // If no further reactions with the same ID, process the stored reactions
            if (gas_phase_id != gas_phase_id_next && j != 0) {
                indice[j]  = ireac;
                distmin[j] = temperature_lower_limit - temperature;
                distmax[j] = temperature - temperature_upper_limit;

                // Set the rate coefficients for the stored reactions
                for (std::size_t k = 0; k <= j; ++k) {
                    std::size_t n = indice[k];
                    if (temperature < ptr_reaction_manager_->reaction_list_[n]->rate_parameters_[gas_phase_reaction_params::kTemperatureLowerLimit]) {
                        reaction_rate_coefficient_[n] = 0.0;
                    }
                    if (temperature > ptr_reaction_manager_->reaction_list_[n]->rate_parameters_[gas_phase_reaction_params::kTemperatureUpperLimit]) {
                        reaction_rate_coefficient_[n] = 0.0;
                    }
                }

                // Set the rate coefficient for the reaction closest to the limits
                if (vector_utils::GetMaxValue(reaction_rate_coefficient_, indice, 0, j) < 1.0e-99) {
                    if (vector_utils::GetMinAbsValue(distmin) <  vector_utils::GetMinAbsValue(distmax)) {
                        std::size_t n = indice[vector_utils::GetMinAbsValueIndex(distmin)];
                        reaction_rate_coefficient_[n] 
                            = (this->*calc_rate_func)(
                                ptr_reaction_manager_->reaction_list_[n], 
                                ptr_reaction_manager_->reaction_list_[n]->rate_parameters_[gas_phase_reaction_params::kTemperatureLowerLimit]
                            );

                    } else {
                        std::size_t n = indice[vector_utils::GetMinValueIndex(distmax)];
                        reaction_rate_coefficient_[n] 
                            = (this->*calc_rate_func)(
                                ptr_reaction_manager_->reaction_list_[n], 
                                ptr_reaction_manager_->reaction_list_[n]->rate_parameters_[gas_phase_reaction_params::kTemperatureUpperLimit]
                            );
                    }
                } // end if

                // Reset and prepare for the next set of reactions
                j = 0;
                for (std::size_t i = 0; i < 10; ++i) {
                    indice[i] = 0;
                    distmin[i] = 9999.0;
                    distmax[i] = 9999.0;
                }
            } // end if
        } // end for ireac
    }


    /**
     * @brief Calculate the reaction rate coefficient using the Modified Arrhenius formula.
     * 
     * @param reaction The reaction object containing rate parameters.
     * @param temperature The temperature at which to calculate the rate coefficient.
     * @return The calculated reaction rate coefficient.
     * 
     * @ref
     * https://kida.astrochem-tools.org/help.html
     * https://ui.adsabs.harvard.edu/abs/2015ApJS..217...20W/abstract
     * https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..63W/abstract
     */
    Real ReactionSimulator::CalculateGasPhaseModifiedArrhenius(const std::shared_ptr<Reaction> reaction, const Real temperature) {
        Real alpha = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
        Real beta  = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
        Real gamma = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];
        return alpha * std::pow(temperature / 300.0, beta) * std::exp(- gamma / temperature);
    }

    /**
     * @brief Calculate the reaction rate coefficient using the ionpol1 formula.
     * 
     * @param reaction The reaction object containing rate parameters.
     * @param temperature The temperature at which to calculate the rate coefficient.
     * @return The calculated reaction rate coefficient.
     * 
     * @ref
     * https://kida.astrochem-tools.org/help.html
     * https://ui.adsabs.harvard.edu/abs/2015ApJS..217...20W/abstract
     * https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..63W/abstract
     * https://kida.astrochem-tools.org/uploads/documents/ionpol_notice.pdf
     */
    Real ReactionSimulator::CalculateGasPhaseIonpol1(const std::shared_ptr<Reaction> reaction, const Real temperature) {
        Real alpha = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
        Real beta  = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
        Real gamma = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];
        return alpha * beta * (0.62 + 0.4767 * gamma * std::sqrt(300.0 / temperature));
    }

    /**
     * @brief Calculate the reaction rate coefficient using the ionpol2 formula.
     * 
     * @param reaction The reaction object containing rate parameters.
     * @param temperature The temperature at which to calculate the rate coefficient.
     * @return The calculated reaction rate coefficient.
     * 
     * @ref
     * https://kida.astrochem-tools.org/help.html
     * https://ui.adsabs.harvard.edu/abs/2015ApJS..217...20W/abstract
     * https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..63W/abstract
     * https://kida.astrochem-tools.org/uploads/documents/ionpol_notice.pdf
     */
    Real ReactionSimulator::CalculateGasPhaseIonpol2(const std::shared_ptr<Reaction> reaction, const Real temperature) {
        Real alpha = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
        Real beta  = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
        Real gamma = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];
        return alpha * beta * (1.0 + 0.0967 * gamma * std::sqrt(300.0 / temperature) + gamma * gamma * 300.0 / (10.526 * temperature));
    }

    /**
     * @brief Calculates the collision rate coefficient between dust particles and charged gas-phase particles.
     * 
     * This function implements the rate coefficient calculation based on the model by Draine & Sutin (1987).
     *
     * @note This function assumes that all required reaction parameters are properly initialized before calling.
     * The function uses the following reaction model to compute the collision rate coefficient.
     *
     * @ref https://ui.adsabs.harvard.edu/abs/1987ApJ...320..803D/abstract
     */
    void ReactionSimulator::CalculateDustAndChargedParticleCollisionRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kDustAndChargedGasParticleCollison] == 0) return;

        const Real temperature = ptr_environment_parameters_->gas_temperature;
        const Real thermal_velocity_coefficient = std::sqrt(8.0 * constants::kBoltzmannConstant * temperature / (M_PI));
        const Real tau_coefficient = constants::kBoltzmannConstant * temperature / (SQR(constants::kChargeUnit));
        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kDustAndChargedGasParticleCollison];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kDustAndChargedGasParticleCollison];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const Real gas_mass = reaction->rate_parameters_[dust_and_charged_particle_collison_params::kGasMass];
            const Real gas_charge = reaction->rate_parameters_[dust_and_charged_particle_collison_params::kGasCharge];
            const Real dust_size = reaction->rate_parameters_[dust_and_charged_particle_collison_params::kDustRadius];
            const Real dust_charge = reaction->rate_parameters_[dust_and_charged_particle_collison_params::kDustCharge];
            const Real dust_cross_section = reaction->rate_parameters_[dust_and_charged_particle_collison_params::kDustCrossSection];
            const Real sticking_probability = reaction->rate_parameters_[dust_and_charged_particle_collison_params::kStickingProbability];

            // Calculate the collision rate coefficient based on the charge interaction model
            Real thermal_velocity = thermal_velocity_coefficient / std::sqrt(gas_mass);
            Real tau = tau_coefficient * dust_size / SQR(gas_charge);
            Real nu = dust_charge / gas_charge;
            Real theta_nu = (nu > 0.0 ? nu / (1.0 + 1.0 / std::sqrt(nu)) : 0.0);

            if (nu == 0.0) {
                reaction_rate_coefficient_[ireac] = sticking_probability * thermal_velocity * dust_cross_section
                    * (1.0 + std::sqrt(M_PI / (2.0 * tau)));
            } else if (nu < 0.0) {
                reaction_rate_coefficient_[ireac] = sticking_probability * thermal_velocity * dust_cross_section
                    * (1.0 - nu / tau) * (1.0 + std::sqrt(2.0 / (tau - 2.0 * nu)));
            } else {
                reaction_rate_coefficient_[ireac] = sticking_probability * thermal_velocity * dust_cross_section
                    * SQR((1.0 + 1.0 / std::sqrt(4.0 * tau + 3.0 * nu))) * std::exp(-theta_nu / tau);
            }
        }
    }


    /**
     * @brief Calculates the collision rate between dust particles
     * 
     * This function calculates the collision rate coefficient for dust particles using a model based on 
     * Umebayashi and Nakano (1990).
     * 
     * @ref 
     * https://ui.adsabs.harvard.edu/abs/1990MNRAS.243..103U/abstract
     */
    void ReactionSimulator::CalculateDustCollisionRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kDustCollision] == 0) return;

        const Real gas_temperature = ptr_environment_parameters_->gas_temperature;
        const Real thermal_velocity_coef = std::sqrt(8.0 * constants::kBoltzmannConstant * gas_temperature / M_PI);
        const std::size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kDustCollision];
        const std::size_t ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kDustCollision];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const Real dust_radius1 = reaction->rate_parameters_[dust_collision_params::kDustRadius1];
            const Real dust_radius2 = reaction->rate_parameters_[dust_collision_params::kDustRadius2];
            const Real dust_charge1 = reaction->rate_parameters_[dust_collision_params::kDustCharge1];
            const Real dust_charge2 = reaction->rate_parameters_[dust_collision_params::kDustCharge2];
            const Real dust_mass1 = reaction->rate_parameters_[dust_collision_params::kDustMass1];
            const Real dust_mass2 = reaction->rate_parameters_[dust_collision_params::kDustMass2];

            // Calculate the collision rate coefficient based on the formula
            const Real reduced_mass = dust_mass1 * dust_mass2 / (dust_mass1 + dust_mass2);
            const Real thermal_velocity = thermal_velocity_coef / std::sqrt(reduced_mass);
            reaction_rate_coefficient_[ireac]
                = M_PI * SQR(dust_radius1 + dust_radius2) * thermal_velocity 
                * (1.0 - dust_charge1 * dust_charge2 / ((dust_radius1 + dust_radius2) * constants::kBoltzmannConstant * gas_temperature));
        }
    }


    /**
     * @brief Calculates the accretion rate coefficient of neutral gas particles on dust surfaces.
     * 
     * This function calculates the rate coefficient for the accretion of neutral gas-phase species
     * onto dust surfaces. The rate is determined by the sticking probability, the thermal velocity of 
     * the gas species, and the total dust cross-section. The sticking probability is computed using 
     * temperature-dependent functions for H and H2, and the coverage of H2O and silicate on dust surfaces. 
     * For other species, it is assumed that the sticking probability is unity.
     * 
     * @ref
     * Chaabouni et al. 2012
     */
    void ReactionSimulator::CalculateNeutralSpeciesAccretionOnDustSurfacesRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kAccretionGasParticleOnDustSurfaces] == 0) return;

        const Real gas_number_density = ptr_environment_parameters_->gas_number_density;
        const Real temperature = ptr_environment_parameters_->gas_temperature;
        const Real thermal_velocity_coefficient = std::sqrt(8.0 * constants::kBoltzmannConstant * temperature / M_PI);
        const Real dust_total_cross_section = ptr_species_manager_->dust_species_model_parameters_.GetDustTotalCrossSectionPernH() * gas_number_density;
        const std::size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kAccretionGasParticleOnDustSurfaces];
        const std::size_t ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kAccretionGasParticleOnDustSurfaces];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const Real gas_mass = ptr_species_manager_->GetSpeciesMass(idx_r1);
            const Real thermal_velocity = thermal_velocity_coefficient / std::sqrt(gas_mass);

            // Sticking probability
            Real sticking_probability = 1.0;
            Real si_ice, si_bare;

            // Calculate sticking probability for H species (Chaabouni et al. 2012)
            if (idx_r1 == ptr_species_manager_->id_H_) {
                // Ice surface
                si_ice = (1.0 + 2.5 * temperature / 52.0) / std::pow(1.0 + temperature / 52.0, 2.5);
                // Silicate surface
                si_bare = (1.0 + 2.5 * temperature / 25.0) / std::pow(1.0 + temperature / 25.0, 2.5);
                // Overall sticking probability
                sticking_probability = coverage_of_H2O_on_dust_surface_ * si_ice + coverage_of_silicate_on_dust_surface_ * si_bare;
            }
            // Calculate sticking probability for H2 species (Chaabouni et al. 2012)
            else if (idx_r1 == ptr_species_manager_->id_H2_) {
                // Ice surface
                si_ice = (1.0 + 2.5 * temperature / 87.0) / std::pow(1.0 + temperature / 87.0, 2.5);
                // Silicate surface
                si_bare = (1.0 + 2.5 * temperature / 56.0) / std::pow(1.0 + temperature / 56.0, 2.5);
                // Overall sticking probability
                sticking_probability = coverage_of_H2O_on_dust_surface_ * si_ice + coverage_of_silicate_on_dust_surface_ * si_bare;
            }

            // Calculate the reaction rate coefficient
            reaction_rate_coefficient_[ireac] = sticking_probability * thermal_velocity * dust_total_cross_section;
        }
    }


    /**
     * @brief Calculates the thermal desorption rate coefficient on dust surfaces.
     * 
     * This function calculates the rate coefficient for the thermal desorption of species from dust 
     * surfaces. The rate coefficient is determined by the vibrational frequency of the species and 
     * its binding energy on the dust surface. If the rate coefficient is below a certain threshold 
     * (`kMinimumRateCoefficient`), it is set to zero.
     */
    void ReactionSimulator::CalculateThermalDesorptionOnDustSurfacesRateCoefficient() {
        // If there are no reactions of the type 'ThermalDesorptionOnDustSurfaces', return early
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kThermalDesorptionOnDustSurfaces] == 0) return;

        const Real dust_temperature = ptr_environment_parameters_->gas_temperature;
        const Real invT = 1.0 / dust_temperature;
        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kThermalDesorptionOnDustSurfaces];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kThermalDesorptionOnDustSurfaces];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];

            // Calculate the rate coefficient using the Arrhenius-like equation
            const Real vibrational_frequency = reaction->rate_parameters_[thermal_desorption_params::kVibrationFrequency];
            const Real binding_energy = reaction->rate_parameters_[thermal_desorption_params::kBindingEnergyOnH2Oice];
            reaction_rate_coefficient_[ireac] = vibrational_frequency * std::exp(-binding_energy * invT);
            
            // If the rate coefficient is lower than the minimum threshold, set it to zero
            if (reaction_rate_coefficient_[ireac] < kMinimumRateCoefficient) {
                reaction_rate_coefficient_[ireac] = 0.0;
            }
        }
    }


    /**
     * @brief Calculates the cosmic ray desorption rate coefficient on dust surfaces.
     * 
     * This function calculates the rate coefficient for cosmic ray desorption of species from 
     * dust surfaces. 
     */
    void ReactionSimulator::CalculateCosmicRayDesorptionOnDustSurfacesRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kCosmicRayDesorptionOnDustSurfaces] == 0) return;

        // Constants for the calculation, assuming the dust temperature is 70K
        const Real invT = 1.0 / 70.0;
        const Real f70K = 3.16e-19;

        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kCosmicRayDesorptionOnDustSurfaces];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kCosmicRayDesorptionOnDustSurfaces];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            
            // Calculate the rate coefficient based on the formula for cosmic ray desorption
            const Real vibrational_frequency = reaction->rate_parameters_[cosmic_ray_desorption_on_dusts_surfaces_params::kVibrationFrequency];
            const Real binding_energy_H2Oice = reaction->rate_parameters_[cosmic_ray_desorption_on_dusts_surfaces_params::kBindingEnergyOnH2Oice];
            reaction_rate_coefficient_[ireac] = f70K * vibrational_frequency * std::exp(-binding_energy_H2Oice * invT);
        }
    }


    /**
     * @brief 
     */
    void ReactionSimulator::CalculatePhotoDesorptionByExternalUVRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kPhotoDesorptionByExternalUV] == 0) return;

        const Real S_UV = 1.0;
        const Real F_UV = 1.0e8;
        const Real Ypd = 1.0e-4;
        const Real c = F_UV * S_UV * Ypd / (4.0 * kDustSurfaceSitesDensity);

        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kPhotoDesorptionByExternalUV];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kPhotoDesorptionByExternalUV];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            reaction_rate_coefficient_[ireac] = c * std::exp(-2.0 * kVisualExtinction);
            if (number_of_total_layers_ >= kNumberOfActiveSurfaceLayer) {
                reaction_rate_coefficient_[ireac] *= (kNumberOfActiveSurfaceLayer / number_of_total_layers_);
            }
        }
    }


    /**
     * @brief
     */
    void ReactionSimulator::CalculatePhotoDesorptionByCosmicRayGeneratedUVRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kPhotoDesorptionByCRGeneratedUV] == 0) return;

        const Real S_UV_CR = 1.0;
        const Real F_UV_CR = 1.0e4; 
        const Real Ypd = 1.0e-4;
        const Real c = F_UV_CR * S_UV_CR * Ypd / (4.0 * kDustSurfaceSitesDensity);
        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kPhotoDesorptionByCRGeneratedUV];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kPhotoDesorptionByCRGeneratedUV];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            reaction_rate_coefficient_[ireac] = c;
            if (number_of_total_layers_ >= kNumberOfActiveSurfaceLayer) {
                reaction_rate_coefficient_[ireac] *= (kNumberOfActiveSurfaceLayer / number_of_total_layers_);
            }
        }
    }


    /**
     * @brief
     */
    void ReactionSimulator::CalculateDustSurfaceReactionRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kDustSurfaceReaction] == 0) return;

        const Real temperature = ptr_environment_parameters_->gas_temperature;
        const Real invT = 1.0 / temperature;
        const Real gas_number_density = ptr_environment_parameters_->gas_number_density;
        const Real total_dust_abundances = ptr_species_manager_->dust_species_model_parameters_.GetDustTotalAbundance();
        const Real number_of_sites_per_a_dust = ptr_species_manager_->dust_species_model_parameters_.GetNumberOfSitesPerDust();
        const Real inv_number_of_dust_surface_sites = 1.0 / (number_of_sites_per_a_dust * total_dust_abundances * gas_number_density);

        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kDustSurfaceReaction];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kDustSurfaceReaction];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const Real vibrational_frequency1 = reaction->rate_parameters_[dust_surface_reaction_params::kVibrationFrequency1];
            const Real vibrational_frequency2 = reaction->rate_parameters_[dust_surface_reaction_params::kVibrationFrequency2];
            const Real diffusion_barrier1 = reaction->rate_parameters_[dust_surface_reaction_params::kDiffusionBarrier1];
            const Real diffusion_barrier2 = reaction->rate_parameters_[dust_surface_reaction_params::kDiffusionBarrier2];

            const Real khop1 = vibrational_frequency1 * std::exp(-diffusion_barrier1 * invT);
            const Real khop2 = vibrational_frequency2 * std::exp(-diffusion_barrier2 * invT);

            Real kappa = 1.0;
            const Real activation_energy = reaction->rate_parameters_[dust_surface_reaction_params::kActivationEnergy];
            const Real tunneling_effect_probability = reaction->rate_parameters_[dust_surface_reaction_params::kTunnelingEffectProbability];
            if (activation_energy > 0.0) {
                Real act_Tinv = activation_energy * invT;
                if (act_Tinv > tunneling_effect_probability) act_Tinv = tunneling_effect_probability;
                Real kappa0 = std::exp(-act_Tinv);
                Real nu_max = std::max(vibrational_frequency1, vibrational_frequency2);
                kappa = nu_max * kappa0 / (nu_max * kappa0 + khop1 + khop2);
            }

            const Real pcd_silicate = reaction->rate_parameters_[dust_surface_reaction_params::kChemicalDesorptionOnSilicateProbability];
            const Real pcd_H2Oice = reaction->rate_parameters_[dust_surface_reaction_params::kChemicalDesorptionOnH2OProbability];
            const Real pcd = coverage_of_silicate_on_dust_surface_ * pcd_silicate + coverage_of_H2O_on_dust_surface_ * pcd_H2Oice;

            reaction_rate_coefficient_[ireac] = pcd * kappa * (khop1 + khop2) * inv_number_of_dust_surface_sites;
        }
    }


    /**
     * @brief
     */
    void ReactionSimulator::CalculateDustMantleReactionRateCoefficient() {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kDustMantleReaction] == 0) return;

        const Real temperature = ptr_environment_parameters_->gas_temperature;
        const Real invT = 1.0 / temperature;
        const Real gas_number_density = ptr_environment_parameters_->gas_number_density;
        const Real total_dust_abundances = ptr_species_manager_->dust_species_model_parameters_.GetDustTotalAbundance();
        const Real number_of_sites_per_a_dust = ptr_species_manager_->dust_species_model_parameters_.GetNumberOfSitesPerDust();
        const Real inv_number_of_dust_surface_sites = 1.0 / (number_of_sites_per_a_dust * total_dust_abundances * gas_number_density);

        const size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kDustMantleReaction];
        const size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kDustMantleReaction];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const Real vibrational_frequency1 = reaction->rate_parameters_[dust_mantle_reaction_params::kVibrationFrequency1];
            const Real vibrational_frequency2 = reaction->rate_parameters_[dust_mantle_reaction_params::kVibrationFrequency2];
            const Real diffusion_barrier1 = reaction->rate_parameters_[dust_mantle_reaction_params::kDiffusionBarrier1];
            const Real diffusion_barrier2 = reaction->rate_parameters_[dust_mantle_reaction_params::kDiffusionBarrier2];

            const Real khop1 = vibrational_frequency1 * std::exp(-diffusion_barrier1 * invT);
            const Real khop2 = vibrational_frequency2 * std::exp(-diffusion_barrier2 * invT);

            Real kappa = 1.0;
            const Real activation_energy = reaction->rate_parameters_[dust_mantle_reaction_params::kActivationEnergy];
            const Real tunneling_effect_probability = reaction->rate_parameters_[dust_mantle_reaction_params::kTunnelingEffectProbability];
            if (activation_energy > 0.0) {
                Real act_Tinv = activation_energy * invT;
                if (act_Tinv > tunneling_effect_probability) act_Tinv = tunneling_effect_probability;
                Real kappa0 = std::exp(-act_Tinv);
                Real nu_max = std::max(vibrational_frequency1, vibrational_frequency2);
                kappa = nu_max * kappa0 / (nu_max * kappa0 + khop1 + khop2);
            }

            reaction_rate_coefficient_[ireac] = reaction->branching_ratio_ * kappa * (khop1 + khop2) * inv_number_of_dust_surface_sites;

            if (number_of_mantle_layers_ > 1.0) {
                reaction_rate_coefficient_[ireac] /= number_of_mantle_layers_;
            }
        }
    }

    void ReactionSimulator:: CalculateRateCoefficient() {
        std::size_t ireac_start, ireac_end;
        const Real cosmic_ray_ionization_rate = ptr_environment_parameters_->cosmic_ray_ionization_rate;
        const Real x_rays_ionization_rate = ptr_environment_parameters_->x_rays_ionization_rate;

        // type1 Dissociation or ionization of species due to direct collision with cosmic-ray particles.
        ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kGasPhase1];
        ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kGasPhase1];
        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            reaction_rate_coefficient_[ireac] = reaction->rate_parameters_[0] 
                * (cosmic_ray_ionization_rate + x_rays_ionization_rate);
        }

        // type3 Dissociation or ionization of neutral species by UV photons with a standard interstellar UV field.
        // const Real visual_extinction  = ptr_environment_parameters_->visual_extinction;
        // const Real scaling_factor_uv_field = ptr_environment_parameters_->scaling_factor_uv_field;
        // ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kGasPhase3];
        // ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kGasPhase3];
        // for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
        //     const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
        //     reaction_rate_coefficient_[ireac] = reaction->rate_parameters_[0] 
        //         * std::exp(-reaction->rate_parameters_[2] * visual_extinction) * scaling_factor_uv_field;
        // }

        // type4-8 Bimolecular reactions includes all chemical reactions between two species.
        CalculateGasPhaseReactionRateCoefficient();

        if (ptr_species_manager_->is_dust_species_) {
            CalculateDustAndChargedParticleCollisionRateCoefficient();
            CalculateDustCollisionRateCoefficient();
        }

        if (ptr_reaction_manager_->is_dust_surface_reaction_) {

            CalculateThermalDesorptionOnDustSurfacesRateCoefficient();
            CalculateCosmicRayDesorptionOnDustSurfacesRateCoefficient();

            if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kPhotoDissociationByUVOnDustSurfaces] > 0) {
                const Real visual_extinction = ptr_environment_parameters_->visual_extinction;
                const Real scaling_factor_uv_field = ptr_environment_parameters_->scaling_factor_uv_field;

                ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kPhotoDissociationByUVOnDustSurfaces];
                ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kPhotoDissociationByUVOnDustSurfaces];

                for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
                    const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
                    reaction_rate_coefficient_[ireac] = reaction->rate_parameters_[0] 
                        * std::exp(-reaction->rate_parameters_[2] * visual_extinction) * scaling_factor_uv_field;
                }
            }
        }
    }


    void ReactionSimulator::CalculateSpeciesAbundancesDependentRateCoefficient(const Real *species_abundances) {
        std::size_t ireac_start, ireac_end;

        // gas phase type2 Dissociation or ionization of species due to UV photons emitted following H2 excitation.
        const Real cosmic_ray_ionization_rate = ptr_environment_parameters_->cosmic_ray_ionization_rate;
        const Real x_rays_ionization_rate = ptr_environment_parameters_->x_rays_ionization_rate;
        const Real x_H2 = species_abundances[ptr_species_manager_->id_H2_];

        ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kGasPhase2];
        ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kGasPhase2];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            // (1/(1-omega) = 2, omega = 0.5)
            reaction_rate_coefficient_[ireac] = reaction->rate_parameters_[0] 
                * (cosmic_ray_ionization_rate + x_rays_ionization_rate) * 2.0 * x_H2;
        }

        // type3 Dissociation or ionization of neutral species by UV photons with a standard interstellar UV field.
        const Real visual_extinction  = ptr_environment_parameters_->visual_extinction;
        const Real scaling_factor_uv_field = ptr_environment_parameters_->scaling_factor_uv_field;

        Real H2_column_density = 0.0, CO_column_density = 0.0;
        if (ptr_species_manager_->id_H2_ != kNotFoundSpecies) {
            H2_column_density = species_abundances[ptr_species_manager_->id_H2_];
        }
        if (ptr_species_manager_->id_CO_ != kNotFoundSpecies) {
            CO_column_density = species_abundances[ptr_species_manager_->id_CO_];
        }

        ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kGasPhase3];
        ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kGasPhase3];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {
            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            reaction_rate_coefficient_[ireac] = reaction->rate_parameters_[0] 
                * std::exp(-reaction->rate_parameters_[2] * visual_extinction) * scaling_factor_uv_field;

            // H2 self-shielding
            if (ptr_reaction_manager_->is_H2_self_shielding_ && reaction->reactant_ids_[0] == ptr_species_manager_->id_H2_) {
                Real thetaH2 = 1.0;
                if (H2_column_density <= kH2ColumnDensity[kNumberOfH2ShieldingFactors-1]) {
                    // Linear extrapolation of the shielding factors
                    for (std::size_t ncol = 0; ncol < kNumberOfH2ShieldingFactors-1; ++ncol) {
                        if (kH2ColumnDensity[ncol] <= H2_column_density && H2_column_density < kH2ColumnDensity[ncol+1]) {
                            thetaH2 = kH2ShieldingFactors[ncol]
                                + (H2_column_density - kH2ColumnDensity[ncol])
                                * (kH2ShieldingFactors[ncol+1] - kH2ShieldingFactors[ncol]) 
                                / (kH2ColumnDensity[ncol+1] - kH2ColumnDensity[ncol]);
                            break;
                        }
                    }
                } else {
                    thetaH2 = kH2ShieldingFactors[kNumberOfH2ShieldingFactors-1];
                }
                // pH20 = 2.54e-11 Lee et al. (1996) Appendix
                reaction_rate_coefficient_[ireac] = 2.54e-11 * thetaH2 * scaling_factor_uv_field;
            }

            // CO self-shielding
            if (ptr_reaction_manager_->is_CO_self_shielding_ && reaction->reactant_ids_[0] == ptr_species_manager_->id_CO_) {

                Real thetaCO = 1.0;
                if (CO_column_density <= kCOColumnDensity[kNumberOfCOShieldingFactors-1]) {
                    for (std::size_t ncol = 0; ncol < kNumberOfCOShieldingFactors-1; ++ncol) {
                        if (kCOColumnDensity[ncol] <= CO_column_density && CO_column_density < kCOColumnDensity[ncol+1]) {
                            thetaCO = kCOShieldingFactors[ncol]
                                + (CO_column_density - kCOColumnDensity[ncol])
                                * (kCOShieldingFactors[ncol+1] - kCOShieldingFactors[ncol])
                                / (kCOColumnDensity[ncol+1] - kCOColumnDensity[ncol]);
                            break;
                        }
                    }
                } else {
                    thetaCO = kCOShieldingFactors[kNumberOfCOShieldingFactors-1];
                }

                Real thetaH2 = 1.0;
                if (H2_column_density <= kH2ColumnDensityForCOShielding[kNumberOfH2ShieldingFactorsForCOShielding-1]) {
                    for (std::size_t ncol = 0; ncol < kNumberOfH2ShieldingFactorsForCOShielding-1; ++ncol) {
                        if (kH2ColumnDensityForCOShielding[ncol] <= H2_column_density && H2_column_density < kH2ColumnDensityForCOShielding[ncol+1]) {
                            thetaH2 = kH2ShieldingFactorsForCOShielding[ncol]
                                + (H2_column_density - kH2ColumnDensityForCOShielding[ncol])
                                * (kH2ShieldingFactorsForCOShielding[ncol+1] - kH2ShieldingFactors[ncol])
                                / (kH2ColumnDensityForCOShielding[ncol+1] - kH2ColumnDensityForCOShielding[ncol]);
                            break;
                        }
                    }
                } else {
                    thetaH2 = kH2ShieldingFactorsForCOShielding[kNumberOfH2ShieldingFactorsForCOShielding-1];
                }

                Real thetaAv = 1.0;
                if (visual_extinction <= kVisualExtinctionForCOShielding[kNumberOfVisualExtinctionFactorsForCOShielding-1]) {
                    for (std::size_t ncol = 0; ncol < kNumberOfVisualExtinctionFactorsForCOShielding-1; ++ncol) {
                        thetaAv = kVisualExtinctionFactorsForCOShielding[ncol]
                            + (visual_extinction - kVisualExtinctionForCOShielding[ncol])
                            * (kVisualExtinctionFactorsForCOShielding[ncol+1] - kVisualExtinctionFactorsForCOShielding[ncol])
                            / (kVisualExtinctionForCOShielding[ncol+1] - kVisualExtinctionForCOShielding[ncol]);
                    }
                } else {
                    thetaAv = kVisualExtinctionForCOShielding[kNumberOfVisualExtinctionFactorsForCOShielding-1];
                }

                // pCO0 = 1.03e-10 Lee et al. (1996) Appendix
                reaction_rate_coefficient_[ireac] = 1.03e-10 * thetaCO * thetaH2 * thetaAv * scaling_factor_uv_field;
            }
        }


        if (!ptr_reaction_manager_->is_dust_surface_reaction_) return;

        // ダスト表面種とダストマントル種の全存在量を計算
        total_abundances_of_dust_surface_species_ = 0.0;
        total_abundances_of_dust_mantle_species_  = 0.0;
        const size_t number_of_gas_species = ptr_species_manager_->number_of_gas_species_;
        const size_t number_of_total_species = ptr_species_manager_->total_number_of_species_;
        for (std::size_t ispe = number_of_gas_species; ispe < number_of_total_species; ++ispe) {
            if (ptr_species_manager_->IsDustSurfaceSpecies(ispe)) {
                total_abundances_of_dust_surface_species_ += species_abundances[ispe];
            }
            if (ptr_species_manager_->IsDustMantleSpecies(ispe)) {
                total_abundances_of_dust_mantle_species_  += species_abundances[ispe];
            }
        }

        // ダストの表面層とマントル層の数の計算
        const Real total_dust_abundances = ptr_species_manager_->dust_species_model_parameters_.GetDustTotalAbundance();
        const Real number_of_sites_per_a_dust = ptr_species_manager_->dust_species_model_parameters_.GetNumberOfSitesPerDust();
        const Real dsite = 1.0 / (number_of_sites_per_a_dust * total_dust_abundances);
        number_of_surface_layers_ = total_abundances_of_dust_surface_species_ * dsite;
        number_of_mantle_layers_ = total_abundances_of_dust_mantle_species_  * dsite;
        number_of_total_layers_ = number_of_surface_layers_ + number_of_mantle_layers_;

        // ダスト表面のH2Oとsilicateの被覆率の計算
        const Real x_sH2O = species_abundances[ptr_species_manager_->id_sH2O_];
        coverage_of_H2O_on_dust_surface_ = x_sH2O * dsite;
        coverage_of_H2O_on_dust_surface_ = std::min(coverage_of_H2O_on_dust_surface_, 1.0);
        coverage_of_silicate_on_dust_surface_ = 1.0 - coverage_of_H2O_on_dust_surface_;

        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kPhotoDissociationByCROnDustSurfaces] > 0) {
            ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kPhotoDissociationByCROnDustSurfaces];
            ireac_end = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kPhotoDissociationByCROnDustSurfaces];
            for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ ireac) {
                const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
                // (1/(1-omega) = 2, omega = 0.5)
                reaction_rate_coefficient_[ireac] = reaction->rate_parameters_[photo_dissociation_by_CR_on_dusts_params::kAlpha] 
                    * (cosmic_ray_ionization_rate + x_rays_ionization_rate) * 2.0 * x_H2;
            }
        }

        CalculateNeutralSpeciesAccretionOnDustSurfacesRateCoefficient();
        CalculatePhotoDesorptionByExternalUVRateCoefficient();
        CalculatePhotoDesorptionByCosmicRayGeneratedUVRateCoefficient();
        CalculateDustSurfaceReactionRateCoefficient();

        if (ptr_reaction_manager_->is_three_phase_reaction_) {
            CalculateDustMantleReactionRateCoefficient();
        }

        return;
    }


    void ReactionSimulator::CalculateDustSurfaceAndMantleSwappingRateCoefficient(const Real *species_abundances) {
        if (ptr_reaction_manager_->number_of_each_type_reactions_[reaction_type_id::kDustSurfaceToMantleSwapping] == 0) return;

        const Real temperature = ptr_environment_parameters_->gas_temperature;
        const Real invT = 1.0 / temperature;

        Real swap_mantle_to_surface = 0.0;
        Real sum_swap_mantle_to_surface = 0.0;

        // calculate mantle to surface swpping rates
        Real alpha_loss = 0.0, r_loss = 0.0;
        if (total_abundances_of_dust_surface_species_ > 0.0) {
            alpha_loss = total_abundances_of_dust_mantle_species_ / total_abundances_of_dust_surface_species_;
            if (alpha_loss >= 1.0) alpha_loss = 1.0;
        }
        if (total_abundances_of_dust_mantle_species_ > 0.0) {
            r_loss = - alpha_loss * total_desorption_rate_ / total_abundances_of_dust_mantle_species_;
        }

        size_t ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kDustMantleToSurfaceSwapping];
        size_t ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kDustMantleToSurfaceSwapping];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {

            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const Real abund_r1 = species_abundances[idx_r1];
            const Real vibrational_frequency = reaction->rate_parameters_[dust_mantle_to_surface_params::kVibrationFrequency];
            const Real binding_energy_on_H2Oice = reaction->rate_parameters_[dust_mantle_to_surface_params::kBindingEnergyOnH2Oice];

            if (abund_r1 > kMinimumSpeciesAbundance) {
                swap_mantle_to_surface = vibrational_frequency * std::exp(-binding_energy_on_H2Oice * invT);
                if (number_of_mantle_layers_ >= 1.0) swap_mantle_to_surface /= number_of_mantle_layers_;
                sum_swap_mantle_to_surface += swap_mantle_to_surface * abund_r1;
            } else {
                swap_mantle_to_surface = 0.0;
            }

            reaction_rate_coefficient_[ireac] = r_loss + swap_mantle_to_surface;
        }

        // surface to mantle
        Real swap_surface_to_mantle = 0.0;
        Real alpha_gain = number_of_surface_layers_ / kNumberOfActiveSurfaceLayer;
        Real r_gain = 0.0;
        if (total_abundances_of_dust_surface_species_ > 0.0) {
            r_gain = alpha_gain * total_accretion_rate_ / total_abundances_of_dust_surface_species_;
        }

        ireac_start = ptr_reaction_manager_->reaction_type_id_start_[reaction_type_id::kDustSurfaceToMantleSwapping];
        ireac_end   = ptr_reaction_manager_->reaction_type_id_end_[reaction_type_id::kDustSurfaceToMantleSwapping];

        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {

            const auto& reaction = ptr_reaction_manager_->reaction_list_[ireac];
            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const Real abund_r1 = species_abundances[idx_r1];

            if (abund_r1 > kMinimumSpeciesAbundance && total_abundances_of_dust_surface_species_ > 0.0) {
                swap_surface_to_mantle = sum_swap_mantle_to_surface / total_abundances_of_dust_surface_species_;
            } else {
                swap_surface_to_mantle = 0.0;
            }

            reaction_rate_coefficient_[ireac] = r_gain + swap_surface_to_mantle;
        }
    }


    bool ReactionSimulator::IsSparseJacobian(const Real threshold) {
        const std::size_t number_of_species = ptr_species_manager_->total_number_of_species_;
        const std::size_t nn = number_of_species * number_of_species;
        Real *y = new Real[number_of_species];
        Real *ydot = new Real[number_of_species];
        Real *jac = new Real[nn];

        for (std::size_t i = 0; i < number_of_species; ++i) {
            y[i]    = 1.0e-5;
            ydot[i] = 0.0;
        }

        for (std::size_t i = 0; i < nn; ++i) {
            jac[i] = 0.0;
        }

        CalculateRateCoefficient();
        OrdinaryDifferentialEquation(number_of_species, 0.0, y, ydot, this);

        Jacobian(number_of_species, 0.0, y, 0, 0, jac, number_of_species, this);
        int count_nonzero = 0;
        for (std::size_t i = 0; i < nn; ++i) {
            if (jac[i] != 0.0) count_nonzero++;
        }

        delete[] y;
        delete[] ydot;
        delete[] jac;

        Real r = static_cast<Real>(count_nonzero) / static_cast<Real>(nn);

        if (r < threshold) {
            return true; // is sparce
        } else {
            return false;  // is not sparce
        }
    }


    void ReactionSimulator::AllocateAndSetLsodeArrays() {
        const std::size_t number_of_species = ptr_species_manager_->total_number_of_species_;
        lsode_parameters_.liw = 20 + number_of_species;
        lsode_parameters_.lrw = 22 + 9*number_of_species + number_of_species*number_of_species;

        lsode_parameters_.rtol.resize(number_of_species);
        lsode_parameters_.atol.resize(number_of_species);
        lsode_parameters_.rwork.resize(lsode_parameters_.lrw);
        lsode_parameters_.iwork.resize(lsode_parameters_.liw);
        lsode_parameters_.is_allocate_arrays = true;

        for (std::size_t i = 0; i < number_of_species; ++i) {
            lsode_parameters_.rtol[i] = relative_tolerance_;
            lsode_parameters_.atol[i] = absolute_tolerance_;
        }
    }


    void ReactionSimulator::AllocateAndSetLsodesArrays() {
        const std::size_t number_of_species = ptr_species_manager_->total_number_of_species_;
        std::vector<Real> y(number_of_species);
        std::vector<Real> ydot(number_of_species);
        std::vector<Real> pdj(number_of_species);
        std::vector<int> ian(number_of_species);
        std::vector<int> jan(number_of_species);
        Real t = 0.0;
        int max_nonzeros = 0;

        for (std::size_t i = 0; i < number_of_species; ++i) {
            y[i]    = 1.0e-5;
            ydot[i] = 0.0;
            ian[i]  = 0;
            jan[i]  = 0;
            pdj[i]  = 0;
        }

        CalculateRateCoefficient();
        OrdinaryDifferentialEquation(number_of_species, t, y.data(), ydot.data(), this);

        for (std::size_t j = 1; j <= number_of_species; ++j) {
            JacobianJth(number_of_species, t, y.data(), j, ian.data(), jan.data(), pdj.data(), this);

            int count_nonzeros = 0;
            for (std::size_t i = 0; i < number_of_species; ++i) {
                if (pdj[i] != 0.0) {
                    count_nonzeros++;
                }
            }

            if (count_nonzeros > max_nonzeros) {
                max_nonzeros = count_nonzeros;
            }
        }

        lsodes_parameters_.liw = 30;
        lsodes_parameters_.lrw = 20 + 3 * max_nonzeros * number_of_species + 16 * number_of_species;

        lsodes_parameters_.rtol.resize(number_of_species);
        lsodes_parameters_.atol.resize(number_of_species);
        lsodes_parameters_.rwork.resize(lsodes_parameters_.lrw);
        lsodes_parameters_.iwork.resize(lsodes_parameters_.liw);
        lsodes_parameters_.is_allocate_arrays = true;

        for (std::size_t i = 0; i < number_of_species; ++i) {
            lsodes_parameters_.rtol[i] = relative_tolerance_;
            lsodes_parameters_.atol[i] = absolute_tolerance_;
        }

        return;
    }


    void ReactionSimulator::ResetLsodeWorkArrays() {
        if (!lsode_parameters_.is_allocate_arrays) return;

        for (int i = 0; i < lsode_parameters_.liw; ++i) {
            lsode_parameters_.iwork[i] = 0;
        }

        for (int i = 0; i < lsode_parameters_.lrw; ++i) {
            lsode_parameters_.rwork[i] = 0.0;
        }

        lsode_parameters_.rwork[5] = 3.154e14;
        lsode_parameters_.iwork[5] = 3000;
    }


    void ReactionSimulator::ResetLsodesWorkArrays() {
        if (!lsodes_parameters_.is_allocate_arrays) return;

        for (int i = 0; i < lsodes_parameters_.liw; ++i) {
            lsodes_parameters_.iwork[i] = 0;
        }

        for (int i = 0; i < lsodes_parameters_.lrw; ++i) {
            lsodes_parameters_.rwork[i] = 0.0;
        }

        lsodes_parameters_.rwork[5] = 3.154e14;
        lsodes_parameters_.iwork[5] = 3000;
    }


    bool ReactionSimulator::Integrate(Real &t, const Real tout, Real *species_abundance, const std::size_t number_of_species) {
        // const size_t number_of_species = ptr_species_manager_->total_number_of_species_;
        if (number_of_species != ptr_species_manager_->total_number_of_species_) {
            std::cout << "Warning: 引数のnumber_of_speciesの値が計算の化学種の数とは異なっています。" << std::endl;
            return false;
        }

        if (is_lsode_integrator_) {
            lsode_parameters_.itol = 2;
            lsode_parameters_.iopt = 1;
            lsode_parameters_.itask = 1;
            lsode_parameters_.mf = 21;
            lsode_parameters_.istate = 1;

            for (std::size_t i = 0; i < number_of_species; ++i) {
                lsode_parameters_.atol[i] = std::max(kOdepackMinimumAbsoluteTolerance, absolute_tolerance_ * species_abundance[i]);
            }

            ResetLsodeWorkArrays();

            odepack_cpp::Odepack odepack;
            odepack.DLSODE(
                OrdinaryDifferentialEquation, 
                number_of_species,
                species_abundance, 
                t, 
                tout, 
                lsode_parameters_.itol,
                lsode_parameters_.rtol.data(), 
                lsode_parameters_.atol.data(), 
                lsode_parameters_.itask, 
                lsode_parameters_.istate, 
                lsode_parameters_.iopt, 
                lsode_parameters_.rwork.data(), 
                lsode_parameters_.lrw, 
                lsode_parameters_.iwork.data(), 
                lsode_parameters_.liw,
                Jacobian, 
                lsode_parameters_.mf, 
                this
            );

            for (std::size_t i = 0; i < number_of_species; ++i) {
                if (species_abundance[i] < kMinimumSpeciesAbundance) {
                    species_abundance[i] = kMinimumSpeciesAbundance;
                }
            }

            if (lsode_parameters_.istate < 2) {
                std::cerr << "Error: LSODE istate = " << lsode_parameters_.istate << std::endl;
                return false;
            }
            return true;

        } else {
            lsodes_parameters_.itol = 2;
            lsodes_parameters_.iopt = 1;
            lsodes_parameters_.itask = 1;
            lsodes_parameters_.mf = 121;
            lsodes_parameters_.istate = 1;

            for (std::size_t i = 0; i < number_of_species; ++i) {
                lsodes_parameters_.atol[i] = std::max(kOdepackMinimumAbsoluteTolerance, absolute_tolerance_ * species_abundance[i]);
            }

            ResetLsodesWorkArrays();
            
            odepack_cpp::Odepack odepack;
            odepack.DLSODES(
                OrdinaryDifferentialEquation, 
                number_of_species,
                species_abundance, 
                t, 
                tout, 
                lsodes_parameters_.itol,
                lsodes_parameters_.rtol.data(), 
                lsodes_parameters_.atol.data(), 
                lsodes_parameters_.itask, 
                lsodes_parameters_.istate, 
                lsodes_parameters_.iopt, 
                lsodes_parameters_.rwork.data(), 
                lsodes_parameters_.lrw, 
                lsodes_parameters_.iwork.data(), 
                lsodes_parameters_.liw,
                JacobianJth, 
                lsodes_parameters_.mf, 
                this
            );

            for (std::size_t i = 0; i < number_of_species; ++i) {
                if (species_abundance[i] < kMinimumSpeciesAbundance) {
                    species_abundance[i] = kMinimumSpeciesAbundance;
                }
            }

            if (lsodes_parameters_.istate < 2) {
                std::cerr << "Error: DLSODES istate = " << lsodes_parameters_.istate << std::endl;
                return false;
            }
            return true;
        }
    }


    bool ReactionSimulator::Integrate(Real &t, const Real tout, std::vector<Real>& species_abundance) {
        return Integrate(t, tout, species_abundance.data(), species_abundance.size());
    }


    bool ReactionSimulator::Integrate(Real &t, const Real tout, Real *species_abundance, const std::size_t number_of_species, std::ofstream& file) {
        if (number_of_species != ptr_species_manager_->total_number_of_species_) {
            std::cout << "Warning: 引数のnumber_of_speciesの値が計算の化学種の数とは異なっています。" << std::endl;
            return false;
        }

        if (!file.is_open()) {
            return false;
        }

        if (is_lsode_integrator_) {
            lsode_parameters_.itol = 2;
            lsode_parameters_.iopt = 1;
            lsode_parameters_.itask = 1;
            lsode_parameters_.mf = 21;
            lsode_parameters_.istate = 1;

            for (std::size_t i = 0; i < number_of_species; ++i) {
                lsode_parameters_.atol[i] = std::max(kOdepackMinimumAbsoluteTolerance, absolute_tolerance_ * species_abundance[i]);
            }

            ResetLsodeWorkArrays();

            odepack_cpp::Odepack odepack;
            odepack.DLSODE(
                OrdinaryDifferentialEquation, 
                number_of_species,
                species_abundance, 
                t, 
                tout, 
                lsode_parameters_.itol,
                lsode_parameters_.rtol.data(), 
                lsode_parameters_.atol.data(), 
                lsode_parameters_.itask, 
                lsode_parameters_.istate, 
                lsode_parameters_.iopt, 
                lsode_parameters_.rwork.data(), 
                lsode_parameters_.lrw, 
                lsode_parameters_.iwork.data(), 
                lsode_parameters_.liw,
                Jacobian, 
                lsode_parameters_.mf, 
                this
            );

            file << std::scientific << std::setw(12) << t << " ";
            for (std::size_t i = 0; i < number_of_species; ++i) {
                file << std::setw(12) << species_abundance[i] << " ";
            }
            file << std::endl;

            if (lsode_parameters_.istate < 0) {
                std::cerr << "Error: DLSODE istate = " << lsode_parameters_.istate << std::endl;
                return false;
            }
            return true;

        } else {
            lsodes_parameters_.itol = 2;
            lsodes_parameters_.iopt = 1;
            lsodes_parameters_.itask = 1;
            lsodes_parameters_.mf = 121;
            lsodes_parameters_.istate = 1;

            for (std::size_t i = 0; i < number_of_species; ++i) {
                lsodes_parameters_.atol[i] = std::max(kOdepackMinimumAbsoluteTolerance, absolute_tolerance_ * species_abundance[i]);
            }

            ResetLsodesWorkArrays();

            odepack_cpp::Odepack odepack;
            odepack.DLSODES(
                OrdinaryDifferentialEquation, 
                number_of_species,
                species_abundance, 
                t, 
                tout,
                lsodes_parameters_.itol,
                lsodes_parameters_.rtol.data(), 
                lsodes_parameters_.atol.data(), 
                lsodes_parameters_.itask, 
                lsodes_parameters_.istate, 
                lsodes_parameters_.iopt, 
                lsodes_parameters_.rwork.data(), 
                lsodes_parameters_.lrw, 
                lsodes_parameters_.iwork.data(), 
                lsodes_parameters_.liw,
                JacobianJth, 
                lsodes_parameters_.mf, 
                this
            );

            file << std::scientific << std::setw(12) << t << " ";
            for (std::size_t i = 0; i < number_of_species; ++i) {
                if (species_abundance[i] < kMinimumSpeciesAbundance) {
                    species_abundance[i] = kMinimumSpeciesAbundance;
                }
                file << std::setw(12) << species_abundance[i] << " ";
            }
            file << std::endl;

            if (lsodes_parameters_.istate < 0) {
                std::cerr << "Error: DLSODES istate = " << lsodes_parameters_.istate << std::endl;
                return false;
            }
            return true;
        }
    }


    bool ReactionSimulator::Integrate(Real &t, const Real tout, std::vector<Real>& species_abundance, std::ofstream& file) {
        return Integrate(t, tout, species_abundance.data(), species_abundance.size(), file);
    }


    void ReactionSimulator::OrdinaryDifferentialEquation(int neq, Real t, Real *y, Real *ydot, void *user_data) {
        ReactionSimulator* ptr_recation_driver  = (ReactionSimulator*)(user_data);
        ReactionManager*   ptr_reaction_manager = ptr_recation_driver->ptr_reaction_manager_;
        SpeciesManager*    ptr_species_manager  = ptr_recation_driver->ptr_species_manager_;
        EnvironmentParameters*  ptr_environment_parameters = ptr_recation_driver->ptr_environment_parameters_;

        // 環境パラメータの取得
        const std::size_t number_of_total_reactions = ptr_reaction_manager->total_number_of_reactions_;
        const Real number_density = ptr_environment_parameters->gas_number_density;

        // 変化率を格納する配列の初期化
        Real ydot_loss[neq], ydot_gain[neq];
        for (int i = 0; i < neq; ++i) {
            ydot[i]      = 0.0;
            ydot_loss[i] = 0.0;
            ydot_gain[i] = 0.0;
        }

        // 存在量に依存する反応率係数の計算
        ptr_recation_driver->CalculateSpeciesAbundancesDependentRateCoefficient(y);

        // 反応ごとのループして変化率を計算
        for (std::size_t ireac = 0; ireac < number_of_total_reactions; ++ireac) {
            const auto& reaction = ptr_reaction_manager->reaction_list_[ireac];

            // 特定の反応をスキップ
            const size_t type_id = reaction->type_id_;
            if (type_id == reaction_type_id::kDustSurfaceToMantleSwapping || 
                type_id == reaction_type_id::kDustMantleToSurfaceSwapping ) {
                continue;
            }

            // 反応係数の取得
            const Real rate_coef = ptr_recation_driver->reaction_rate_coefficient_[ireac];
            if (rate_coef <= kMinimumRateCoefficient) continue;

            // 反応物と生成物のインデックスを取得
            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const std::size_t idx_r2 = reaction->reactant_ids_[1];
            const std::size_t idx_r3 = reaction->reactant_ids_[2];
            const std::size_t idx_p1 = reaction->product_ids_[0];
            const std::size_t idx_p2 = reaction->product_ids_[1];
            const std::size_t idx_p3 = reaction->product_ids_[2];
            const std::size_t idx_p4 = reaction->product_ids_[3];
            const std::size_t idx_p5 = reaction->product_ids_[4];

            // 反応率の計算
            Real rate;
            if (idx_r2 == kNotFoundSpecies) {
                rate = rate_coef * y[idx_r1];
            } else {
                if (idx_r3 == kNotFoundSpecies) {
                    rate = rate_coef * y[idx_r1] * y[idx_r2] * number_density;
                } else {
                    rate = rate_coef * y[idx_r1] * y[idx_r2] * y[idx_r3] * number_density * number_density;
                }
            }

            // 変化率の計算
            ydot[idx_r1] -= rate;
            if (idx_r2 != kNotFoundSpecies) ydot[idx_r2] -= rate;
            if (idx_r3 != kNotFoundSpecies) ydot[idx_r3] -= rate;

            ydot[idx_p1] += rate;
            if (idx_p2 != kNotFoundSpecies) ydot[idx_p2] += rate;
            if (idx_p3 != kNotFoundSpecies) ydot[idx_p3] += rate;
            if (idx_p4 != kNotFoundSpecies) ydot[idx_p4] += rate;
            if (idx_p5 != kNotFoundSpecies) ydot[idx_p5] += rate;

            // ロスとゲインの更新 (三相反応の計算に利用)
            if (ptr_reaction_manager->is_three_phase_reaction_) {
                ydot_loss[idx_r1] -= rate;
                if (idx_r2 != kNotFoundSpecies) ydot_loss[idx_r2] -= rate;
                if (idx_r3 != kNotFoundSpecies) ydot_loss[idx_r3] -= rate;

                ydot_gain[idx_p1] += rate;
                if (idx_p2 != kNotFoundSpecies) ydot_gain[idx_p2] += rate;
                if (idx_p3 != kNotFoundSpecies) ydot_gain[idx_p3] += rate;
                if (idx_p4 != kNotFoundSpecies) ydot_gain[idx_p4] += rate;
                if (idx_p5 != kNotFoundSpecies) ydot_gain[idx_p5] += rate;
            }
        }

        // 三相反応がある場合
        if (ptr_reaction_manager->is_three_phase_reaction_) {

            Real total_ydot_gain = 0.0;
            Real total_ydot_loss = 0.0;

            const std::size_t number_of_gas_species   = ptr_species_manager->number_of_gas_species_;
            const std::size_t number_of_total_species = ptr_species_manager->total_number_of_species_;

            // ダスト表面種に関して、各々のゲインとロスの和の計算
            for (std::size_t ispe = number_of_gas_species; ispe < number_of_total_species; ++ispe) {
                if (ptr_species_manager->IsDustSurfaceSpecies(ispe)) {
                    total_ydot_gain += ydot_gain[ispe];
                    total_ydot_loss += ydot_loss[ispe];
                }
            }

            // 反応率係数の計算
            ptr_recation_driver->total_desorption_rate_ = total_ydot_loss;
            ptr_recation_driver->total_accretion_rate_  = total_ydot_gain;
            ptr_recation_driver->CalculateDustSurfaceAndMantleSwappingRateCoefficient(y);

            // ダスト表面->ダストマントルの変化率の計算
            std::size_t ireac_start = ptr_reaction_manager->reaction_type_id_start_[reaction_type_id::kDustSurfaceToMantleSwapping];
            std::size_t ireac_end   = ptr_reaction_manager->reaction_type_id_end_[reaction_type_id::kDustSurfaceToMantleSwapping];

            for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {

                const auto& reaction = ptr_reaction_manager->reaction_list_[ireac];
                const Real rate_coef = ptr_recation_driver->reaction_rate_coefficient_[ireac];
                
                const std::size_t idx_r1 = reaction->reactant_ids_[0];
                const std::size_t idx_p1 = reaction->product_ids_[0];

                Real rate = rate_coef * y[idx_r1];
                ydot[idx_r1] -= rate;
                ydot[idx_p1] += rate;
            }

            // ダストマントル種->ダスト表面種の変化率の計算
            ireac_start = ptr_reaction_manager->reaction_type_id_start_[reaction_type_id::kDustMantleToSurfaceSwapping];
            ireac_end   = ptr_reaction_manager->reaction_type_id_end_[reaction_type_id::kDustMantleToSurfaceSwapping];

            for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {

                const auto& reaction = ptr_reaction_manager->reaction_list_[ireac];
                Real rate = ptr_recation_driver->reaction_rate_coefficient_[ireac];
                
                const std::size_t idx_r1 = reaction->reactant_ids_[0];
                const std::size_t idx_p1 = reaction->product_ids_[0];

                rate = rate * y[idx_r1];
                ydot[idx_r1] -= rate;
                ydot[idx_p1] += rate;
            }
        }
    }


    void ReactionSimulator::JacobianJth(int neq, Real t, Real *y, int j, int *ian, int *jan, Real *pdj, void *user_data) {
        const std::size_t idx_j = static_cast<std::size_t>(j-1);

        // 必要なポインタをuser_dataから取得
        ReactionSimulator* ptr_recation_driver  = (ReactionSimulator*)(user_data);
        ReactionManager*   ptr_reaction_manager = ptr_recation_driver->ptr_reaction_manager_;
        EnvironmentParameters*  ptr_environment_parameters = ptr_recation_driver->ptr_environment_parameters_;

        // 環境パラメータの取得
        const Real number_density = ptr_environment_parameters->gas_number_density;

        //　初期化
        for (std::size_t i = 0; i < neq; ++i) pdj[i] = 0.0;

        // ループ
        const std::size_t number_of_reaction_species = ptr_reaction_manager->number_of_reactions_involved_with_species_[idx_j];
        for (std::size_t i = 0; i < number_of_reaction_species; ++i) {

            const std::size_t idx_reaction = ptr_reaction_manager->reaction_index_list_involved_with_species_[idx_j][i];
            const auto& reaction = ptr_reaction_manager->reaction_list_[idx_reaction];
            const Real rate_coef = ptr_recation_driver->reaction_rate_coefficient_[idx_reaction];
            // if (rate_coef == 0.0) continue;
            if (rate_coef <= kMinimumRateCoefficient) continue;

            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const std::size_t idx_r2 = reaction->reactant_ids_[1];
            const std::size_t idx_r3 = reaction->reactant_ids_[2];
            const std::size_t idx_p1 = reaction->product_ids_[0];
            const std::size_t idx_p2 = reaction->product_ids_[1];
            const std::size_t idx_p3 = reaction->product_ids_[2];
            const std::size_t idx_p4 = reaction->product_ids_[3];
            const std::size_t idx_p5 = reaction->product_ids_[4];

            if (idx_r2 == kNotFoundSpecies) {

                if (idx_r1 == idx_j) {
                    const Real rate = rate_coef;
                    pdj[idx_r1] -= rate;
                    pdj[idx_p1] += rate;
                    if (idx_p2 != kNotFoundSpecies) pdj[idx_p2] += rate;
                    if (idx_p3 != kNotFoundSpecies) pdj[idx_p3] += rate;
                    if (idx_p4 != kNotFoundSpecies) pdj[idx_p4] += rate;
                    if (idx_p5 != kNotFoundSpecies) pdj[idx_p5] += rate;
                }

            } else if (idx_r3 == kNotFoundSpecies) {

                if (idx_r1 == idx_j) {
                    const Real rate = rate_coef * y[idx_r2] * number_density;
                    pdj[idx_r1] -= rate;
                    pdj[idx_r2] -= rate;
                    pdj[idx_p1] += rate;
                    if (idx_p2 != kNotFoundSpecies) pdj[idx_p2] += rate;
                    if (idx_p3 != kNotFoundSpecies) pdj[idx_p3] += rate;
                    if (idx_p4 != kNotFoundSpecies) pdj[idx_p4] += rate;
                    if (idx_p5 != kNotFoundSpecies) pdj[idx_p5] += rate;
                }

                if (idx_r2 == idx_j) {
                    const Real rate = rate_coef * y[idx_r1] * number_density;
                    pdj[idx_r1] -= rate;
                    pdj[idx_r2] -= rate;
                    pdj[idx_p1] += rate;
                    if (idx_p2 != kNotFoundSpecies) pdj[idx_p2] += rate;
                    if (idx_p3 != kNotFoundSpecies) pdj[idx_p3] += rate;
                    if (idx_p4 != kNotFoundSpecies) pdj[idx_p4] += rate;
                    if (idx_p5 != kNotFoundSpecies) pdj[idx_p5] += rate;
                }

            } else {

                if (idx_r1 == idx_j) {
                    const Real rate = rate_coef * y[idx_r2] * y[idx_r3] * number_density * number_density;
                    pdj[idx_r1] -= rate;
                    pdj[idx_r2] -= rate;
                    pdj[idx_r3] -= rate;
                    pdj[idx_p1] += rate;
                    if (idx_p2 != kNotFoundSpecies) pdj[idx_p2] += rate;
                    if (idx_p3 != kNotFoundSpecies) pdj[idx_p3] += rate;
                    if (idx_p4 != kNotFoundSpecies) pdj[idx_p4] += rate;
                    if (idx_p5 != kNotFoundSpecies) pdj[idx_p5] += rate;
                }

                if (idx_r2 == idx_j) {
                    const Real rate = rate_coef * y[idx_r1] * y[idx_r3] * number_density * number_density;
                    pdj[idx_r1] -= rate;
                    pdj[idx_r2] -= rate;
                    pdj[idx_r3] -= rate;
                    pdj[idx_p1] += rate;
                    if (idx_p2 != kNotFoundSpecies) pdj[idx_p2] += rate;
                    if (idx_p3 != kNotFoundSpecies) pdj[idx_p3] += rate;
                    if (idx_p4 != kNotFoundSpecies) pdj[idx_p4] += rate;
                    if (idx_p5 != kNotFoundSpecies) pdj[idx_p5] += rate;
                }

                if (idx_r3 == idx_j) {
                    const Real rate = rate_coef * y[idx_r1] * y[idx_r2] * number_density * number_density;
                    pdj[idx_r1] -= rate;
                    pdj[idx_r2] -= rate;
                    pdj[idx_r3] -= rate;
                    pdj[idx_p1] += rate;
                    if (idx_p2 != kNotFoundSpecies) pdj[idx_p2] += rate;
                    if (idx_p3 != kNotFoundSpecies) pdj[idx_p3] += rate;
                    if (idx_p4 != kNotFoundSpecies) pdj[idx_p4] += rate;
                    if (idx_p5 != kNotFoundSpecies) pdj[idx_p5] += rate;
                }
            }
        }
    }


    void ReactionSimulator::Jacobian(int neq, Real t, Real *y, int ml, int mu, Real *pd, int nrowpd, void *user_data) {
    #ifndef JAC
    #define JAC(i, j) ARRAY2D(pd, neq, i+1, j+1)
    #endif

        // 必要なポインタをuser_dataから取得
        ReactionSimulator* ptr_recation_driver  = (ReactionSimulator*)(user_data);
        ReactionManager*   ptr_reaction_manager = ptr_recation_driver->ptr_reaction_manager_;
        EnvironmentParameters*  ptr_environment_parameters = ptr_recation_driver->ptr_environment_parameters_;

        // 環境パラメータの取得
        const Real number_density = ptr_environment_parameters->gas_number_density;

        std::size_t number_of_total_reaction = ptr_reaction_manager->total_number_of_reactions_;

        for (std::size_t ireac = 0; ireac < number_of_total_reaction; ++ireac) {

            const auto& reaction = ptr_reaction_manager->reaction_list_[ireac];
            const Real rate_coef = ptr_recation_driver->reaction_rate_coefficient_[ireac];
            if (rate_coef <= kMinimumRateCoefficient) continue;

            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const std::size_t idx_r2 = reaction->reactant_ids_[1];
            const std::size_t idx_r3 = reaction->reactant_ids_[2];
            
            const std::size_t idx_p1 = reaction->product_ids_[0];
            const std::size_t idx_p2 = reaction->product_ids_[1];
            const std::size_t idx_p3 = reaction->product_ids_[2];
            const std::size_t idx_p4 = reaction->product_ids_[3];
            const std::size_t idx_p5 = reaction->product_ids_[4];

            if (idx_r2 == kNotFoundSpecies) {
                const Real rate = rate_coef;
                JAC(idx_r1, idx_r1) -= rate;
                JAC(idx_p1, idx_r1) += rate;
                if (idx_p2 != kNotFoundSpecies) JAC(idx_p2, idx_r1) += rate;
                if (idx_p3 != kNotFoundSpecies) JAC(idx_p3, idx_r1) += rate;
                if (idx_p4 != kNotFoundSpecies) JAC(idx_p4, idx_r1) += rate;
                if (idx_p5 != kNotFoundSpecies) JAC(idx_p5, idx_r1) += rate;
                continue;
            }

            if (idx_r3 == kNotFoundSpecies) {

                Real rate = rate_coef * y[idx_r2] * number_density;
                JAC(idx_r1, idx_r1) -= rate;
                JAC(idx_r2, idx_r1) -= rate;
                JAC(idx_p1, idx_r1) += rate;
                if (idx_p2 != kNotFoundSpecies) JAC(idx_p2, idx_r1) += rate;
                if (idx_p3 != kNotFoundSpecies) JAC(idx_p3, idx_r1) += rate;
                if (idx_p4 != kNotFoundSpecies) JAC(idx_p4, idx_r1) += rate;
                if (idx_p5 != kNotFoundSpecies) JAC(idx_p5, idx_r1) += rate;

                rate = rate_coef * y[idx_r1] * number_density;
                JAC(idx_r1, idx_r2) -= rate;
                JAC(idx_r2, idx_r2) -= rate;
                JAC(idx_p1, idx_r2) += rate;
                if (idx_p2 != kNotFoundSpecies) JAC(idx_p2, idx_r2) += rate;
                if (idx_p3 != kNotFoundSpecies) JAC(idx_p3, idx_r2) += rate;
                if (idx_p4 != kNotFoundSpecies) JAC(idx_p4, idx_r2) += rate;
                if (idx_p5 != kNotFoundSpecies) JAC(idx_p5, idx_r2) += rate;

                continue;
            }

            if (idx_r3 != kNotFoundSpecies) {
                Real rate = rate_coef * y[idx_r2] * y[idx_r3] * number_density * number_density;
                JAC(idx_r1, idx_r1) -= rate;
                JAC(idx_r2, idx_r1) -= rate;
                JAC(idx_r3, idx_r1) -= rate;
                JAC(idx_p1, idx_r1) += rate;
                if (idx_p2 != kNotFoundSpecies) JAC(idx_p2, idx_r1) += rate;
                if (idx_p3 != kNotFoundSpecies) JAC(idx_p3, idx_r1) += rate;
                if (idx_p4 != kNotFoundSpecies) JAC(idx_p4, idx_r1) += rate;
                if (idx_p5 != kNotFoundSpecies) JAC(idx_p5, idx_r1) += rate;

                rate = rate_coef * y[idx_r1] * y[idx_r3] * number_density * number_density;
                JAC(idx_r1, idx_r2) -= rate;
                JAC(idx_r2, idx_r2) -= rate;
                JAC(idx_r3, idx_r2) -= rate;
                JAC(idx_p1, idx_r2) += rate;
                if (idx_p2 != kNotFoundSpecies) JAC(idx_p2, idx_r2) += rate;
                if (idx_p3 != kNotFoundSpecies) JAC(idx_p3, idx_r2) += rate;
                if (idx_p4 != kNotFoundSpecies) JAC(idx_p4, idx_r2) += rate;
                if (idx_p5 != kNotFoundSpecies) JAC(idx_p5, idx_r2) += rate;

                rate = rate_coef * y[idx_r1] * y[idx_r2] * number_density * number_density;
                JAC(idx_r1, idx_r3) -= rate;
                JAC(idx_r2, idx_r3) -= rate;
                JAC(idx_r3, idx_r3) -= rate;
                JAC(idx_p1, idx_r3) += rate;
                if (idx_p2 != kNotFoundSpecies) JAC(idx_p2, idx_r3) += rate;
                if (idx_p3 != kNotFoundSpecies) JAC(idx_p3, idx_r3) += rate;
                if (idx_p4 != kNotFoundSpecies) JAC(idx_p4, idx_r3) += rate;
                if (idx_p5 != kNotFoundSpecies) JAC(idx_p5, idx_r3) += rate;

                continue;
            }
        }

        return;
    #ifdef JAC
    #undef JAC
    #endif
    }
}
