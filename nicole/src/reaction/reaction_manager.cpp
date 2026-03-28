#include <fstream>
#include <iostream>
#include <memory>

#include "nicole/reaction/reaction_manager.hpp"
#include "nicole/utils/string_utils.hpp"
#include "nicole/utils/vector_utils.hpp"

namespace nicole {
    ReactionManager::ReactionManager(
        SpeciesManager *ptr_species_manager, 
        InputConfig& input
    ) : ptr_species_manager_(ptr_species_manager),
        is_dust_collision_(input.GetBool("is_dust_collision")),
        is_dust_surface_reaction_(input.GetBool("is_dust_surface_reaction")),
        is_chemical_desorption_(input.GetBool("is_chemical_desorption")),
        is_H2_desorption_(input.GetBool("is_H2_desorption")),
        is_three_phase_reaction_(input.GetBool("is_three_phase_reaction")),
        is_H2_self_shielding_(false),
        is_CO_self_shielding_(false)
    {
        ReadGasPhaseReactionFile(input.GetString("gas_reaction_file"));

        if (ptr_species_manager_->is_dust_species_) {
            GenerateReactionListForDustAndChargedParticleCollision();
            if (is_dust_collision_) {
                GenerateReactionListForDustCollision();
            }
        }

        if (is_dust_surface_reaction_) {
            // Generate reactions for various processes on dust surfaces
            GenerateReactionListForNeutralSpeciesAccretionOnDustSurfaces();
            GenerateReactionListForThermalDesorptionOnDustSurfaces();
            GenerateReactionListForCosmicRayDesorptionOnDustSurfaces();
            GenerateReactionListForPhotoDesorptionByExternalUV();
            GenerateReactionListForPhotoDesorptionByCosmicRayGeneratedUV();
            GenerateReactionListForPhotoDissociationInducedByCRsOnDustSurfaces();
            GenerateReactionListForPhotoDissociationByExternalUVOnDustSurfaces();

            ReadDustSurfaceReactionFile(input.GetString("dust_surface_reaction_file"));

            if (is_three_phase_reaction_) {
                GenerateReactionListForDustMantleReaction();
                GenerateReactionListForDustSurfaceToMantleSwapping();
                GenerateReactionListForDustMantleToSurfaceSwapping();
            }

            ReadSurfaceActivationEnergyFile(input.GetString("activation_energy_file"));
        }

        total_number_of_reactions_ = reaction_list_.size();
        CountNumberOfEachTypeReactions();
        SetReactionTypeIdStartAndEnd();
        SetReactionsInvolvedWithSpecies();

        CalculateReactionBranchingRatio();
        if (is_chemical_desorption_) CalculateChemicalDesorptionProbabilities();
    }


    void ReactionManager::CheckReactionManager(const std::string& filename) {
        std::ofstream file(filename, std::ios::out | std::ios::trunc);
        if (!file.is_open()) {
            return;
        }

        for (const auto& reaction : reaction_list_) {

            reaction->WriteInfoToFile(ptr_species_manager_->species_name_list_, file);

            // if (reaction->type_id_ == reaction_type_id::kDustSurfaceReaction) {
            //     reaction->WriteInfoToFile(ptr_species_manager_->species_name_list_, file);
            // }

            // if (reaction->type_id_ == reaction_type_id::kDustAndChargedGasParticleCollison) {
            //     reaction->WriteInfoToFile(ptr_species_manager_->species_name_list_, file);
            // }

            // if (reaction->type_id_ == reaction_type_id::kDustCollision) {
            //     reaction->WriteInfoToFile(ptr_species_manager_->species_name_list_, file);
            // }

            // if (reaction->type_id_ == reaction_type_id::kPhotoDissociationByCROnDustSurfaces) {
            //     reaction->WriteInfoToFile(ptr_species_manager_->species_name_list_, file);
            // }

            // if (reaction->type_id_ == reaction_type_id::kDustCollision) {
            //     reaction->WriteInfoToFile(ptr_species_manager_->species_name_list_, file);
            // }

        }

        file.close();

        std::cout << "number of gas phase species = " << total_number_of_gas_phase_reactions_ << std::endl;
        std::cout << "number of total reactions   = " << total_number_of_reactions_ << std::endl;
    }


    std::vector<SpeciesID> ReactionManager::GetSpeciesIDList(const std::vector<std::string>& species_name_list) const {
        // Initialize the result vector with the same size as the input species list
        std::vector<SpeciesID> ids(species_name_list.size());

        std::transform(
            species_name_list.begin(), 
            species_name_list.end(), 
            ids.begin(),
            [&](const std::string& species) {
                return ptr_species_manager_->FindSpeciesID(species);
            }
        );

        return ids;
    }


    std::vector<std::string> ReactionManager::GetSpeciesNameList(const std::vector<SpeciesID>& ids) const {
        // Initialize the result vector with the same size as the input indices list
        std::vector<std::string> species_name_list(ids.size());

        std::transform(
            ids.begin(), 
            ids.end(), 
            species_name_list.begin(),
            [&](std::size_t index) {
                return (index == kNotFoundSpecies) ? kEmptyString : ptr_species_manager_->GetSpeciesName(index); 
            }
        );

        return species_name_list;
    }


    bool ReactionManager::AreReactionSpeciesInSpeciesNameList(const std::vector<std::string>& reactants, const std::vector<std::string>& products) const {       
        // Check if all reactants are present in the species name list or special species list
        for (const auto& reactant : reactants) {
            if (!string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant) && 
                !string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) {
                return false;  // If a reactant is not found, return false
            }
        }

        // Check if all products are present in the species name list or special species list
        for (const auto& product : products) {
            if (!string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product) &&
                !string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) {
                return false;  // If a product is not found, return false
            }
        }

        // If all reactants and products are found, return true
        return true;
    }


    std::vector<std::string> ReactionManager::SplitGasReactionLine(const std::string& line) const {
        // Extract reactant names (trim leading and trailing spaces)
        std::string reactant1 = string_utils::Trim(line.substr( 0, 11));
        std::string reactant2 = string_utils::Trim(line.substr(11, 11));
        std::string reactant3 = string_utils::Trim(line.substr(22, 11));

        // Extract product names (trim leading and trailing spaces)
        std::string product1  = string_utils::Trim(line.substr(34, 11));
        std::string product2  = string_utils::Trim(line.substr(45, 11));
        std::string product3  = string_utils::Trim(line.substr(56, 11));
        std::string product4  = string_utils::Trim(line.substr(67, 11));
        std::string product5  = string_utils::Trim(line.substr(78, 11));

        // Extract rate parameters (trim spaces)
        std::string rate_param1 = string_utils::Trim(line.substr(90, 10));   // alpha (reaction rate coefficient)
        std::string rate_param2 = string_utils::Trim(line.substr(101, 10));  // beta (reaction rate coefficient)
        std::string rate_param3 = string_utils::Trim(line.substr(112, 10));  // gamma (reaction rate coefficient)
        std::string rate_param4 = string_utils::Trim(line.substr(123, 9));   // uncertainty factor on rate coefficient
        std::string rate_param5 = string_utils::Trim(line.substr(132, 9));   // temperature dependence of uncertainty factor

        // Extract additional information (trim spaces)
        std::string logn        = string_utils::Trim(line.substr(141, 5));  // Type of uncertainty (e.g., 'logn')
        std::string type_id     = string_utils::Trim(line.substr(146, 4));  // Type ID of the gas-phase reaction
        std::string tmin        = string_utils::Trim(line.substr(150, 8));  // Minimum temperature for validity
        std::string tmax        = string_utils::Trim(line.substr(158, 4));  // Maximum temperature for validity
        std::string formula_id  = string_utils::Trim(line.substr(164, 1));  // Formula ID
        std::string reac_id     = string_utils::Trim(line.substr(165, 6));  // Reaction ID

        // Return all extracted components as a vector of strings
        return std::vector<std::string>{
            reactant1,
            reactant2,
            reactant3,
            product1,
            product2,
            product3,
            product4,
            product5,
            rate_param1,
            rate_param2,
            rate_param3,
            rate_param4,
            rate_param5,
            logn,
            type_id,
            tmin,
            tmax,
            formula_id,
            reac_id
        };
    }

    
    void ReactionManager::ReadGasPhaseReactionFile(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << file_path << std::endl;
            return;
        }

        // Temporary list to store gas phase reactions
        std::vector<std::shared_ptr<Reaction> > local_gas_reaction_list;

        // Read the file line by line
        std::string line;
        int line_number = 0;
        while (std::getline(file, line)) {
            line_number++;

            // Skip empty lines or comment lines (starting with '!' or '#')
            if (line.empty() || line[0] == '!' || line[1] == '#') continue;
            
            std::vector<std::string> split_result = SplitGasReactionLine(line);
            if (split_result.size() < 19) {
                std::cerr << "Warning: Invalid format in line: " << line << " : line number = " << line_number << std::endl;
                continue;
            }

            // Extract reactants and products
            std::vector<std::string> reactants = {
                split_result[0], 
                split_result[1], 
                split_result[2]
            };
            
            std::vector<std::string> products = {
                split_result[3], 
                split_result[4], 
                split_result[5],
                split_result[6], 
                split_result[7]
            };

            if (!AreReactionSpeciesInSpeciesNameList(reactants, products)) continue;
            
            std::vector<SpeciesID> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<SpeciesID> product_ids  = GetSpeciesIDList(products);

            // Extract reaction parameters 
            const Real alpha = std::stod(split_result[8]);
            const Real beta  = std::stod(split_result[9]);
            const Real gamma = std::stod(split_result[10]);
            std::size_t type_id = std::stoul(split_result[14]);
            const Real temperature_lower_limit = std::stod(split_result[15]);
            const Real temperature_upper_limit = std::stod(split_result[16]);
            const std::size_t gas_phase_formula_id = std::stoul(split_result[17]);
            const std::size_t gas_phase_id = std::stoul(split_result[18]);

            std::vector<Real> rate_parameters(gas_phase_reaction_params::kNumParams);
            rate_parameters[gas_phase_reaction_params::kAlpha] = alpha;
            rate_parameters[gas_phase_reaction_params::kBeta]  = beta;
            rate_parameters[gas_phase_reaction_params::kGamma] = gamma;
            rate_parameters[gas_phase_reaction_params::kTemperatureLowerLimit] = temperature_lower_limit;
            rate_parameters[gas_phase_reaction_params::kTemperatureUpperLimit] = temperature_upper_limit;
            rate_parameters[gas_phase_reaction_params::kFormulaID] = gas_phase_formula_id;
            rate_parameters[gas_phase_reaction_params::kID] = gas_phase_id;

            // Create a new Reaction object and add it to the temporary list
            const std::size_t dummy_index = 0;
            local_gas_reaction_list.emplace_back(std::make_shared<Reaction>(
                dummy_index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        }

        // Reorder the reaction list by reaction type ID and update internal data structures
        ReactionID id = 0;
        total_number_of_gas_phase_reactions_ = 0;
        for (std::size_t itype = 0; itype < reaction_type_id::kNumberOfTypeID; ++itype) {
            for (const auto& reaction : local_gas_reaction_list) {
                if (reaction->type_id_ == itype) {
                    id++;
                    reaction->id_ = id;
                    reaction_list_.emplace_back(reaction);
                    total_number_of_gas_phase_reactions_++;
                }
            }
        }
    }


    std::vector<std::string> ReactionManager::SplitDustSurfaceReactionLine(const std::string& line) const {
        // Trim and extract reactants
        std::string reactant1 = string_utils::Trim(line.substr(0, 11));
        std::string reactant2 = string_utils::Trim(line.substr(11, 11));
        std::string reactant3 = string_utils::Trim(line.substr(22, 11));

        // Trim and extract products
        std::string product1  = string_utils::Trim(line.substr(34, 11));
        std::string product2  = string_utils::Trim(line.substr(45, 11));
        std::string product3  = string_utils::Trim(line.substr(56, 11));
        std::string product4  = string_utils::Trim(line.substr(67, 11));
        std::string product5  = string_utils::Trim(line.substr(78, 11));

        // Trim and extract rate parameters (branching ratio and uncertainty)
        std::string rate_param1 = string_utils::Trim(line.substr(91, 4)); // branching ratio
        std::string rate_param2 = string_utils::Trim(line.substr(96, 8)); // uncertainty on the branching ratio

        // Return all extracted values as a vector of strings
        return std::vector<std::string>{
            reactant1,
            reactant2,
            reactant3,
            product1,
            product2,
            product3,
            product4,
            product5,
            rate_param1,
            rate_param2,
        };
    }


    void ReactionManager::ReadDustSurfaceReactionFile(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << file_path << std::endl;
            return;
        }

        const std::size_t type_id = reaction_type_id::kDustSurfaceReaction;
        std::size_t index = reaction_list_.size();

        std::string line;
        std::vector<std::string> previous_reactants;
        std::vector<std::string> previous_products;
        while (std::getline(file, line)) {

            // Skip empty lines or comment lines
            if (line.empty() || line[0] == '!' || line[1] == '#') continue;

            std::vector<std::string> split_result = SplitDustSurfaceReactionLine(line);

            // Store the reactants
            std::vector<std::string> reactants = {
                split_result[0], 
                split_result[1], 
                split_result[2]
            };

            // If any reactant exists in the gas species list, add 's' prefix
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(ptr_species_manager_->gas_species_name_list_, reactant)) {
                    reactant = kDustSurfaceSpeciesPrefix + reactant;
                }
            }
            
            // Store the products
            std::vector<std::string> products = {
                split_result[3], 
                split_result[4], 
                split_result[5],
                split_result[6], 
                split_result[7]
            };

            // If any product exists in the gas species list, add 's' prefix
            for (auto& product : products) {
                if (string_utils::IsInStringVector(ptr_species_manager_->gas_species_name_list_, product)) {
                    product = kDustSurfaceSpeciesPrefix + product;
                }
            }

            if (!AreReactionSpeciesInSpeciesNameList(reactants, products)) continue;

            // Skip if the current reaction is identical to the previous one to avoid duplicates
            if (reactants == previous_reactants && products == previous_products) continue;

            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Get DustSurfaceSpecies Object for reactants
            std::shared_ptr<DustSurfaceSpecies> dust_surface_species1 = ptr_species_manager_->FindDustSurfaceSpeciesByName(reactants[0]);
            std::shared_ptr<DustSurfaceSpecies> dust_surface_species2 = ptr_species_manager_->FindDustSurfaceSpeciesByName(reactants[1]);

            // Reaction parameters
            std::vector<Real> rate_parameters(dust_surface_reaction_params::kNumParams);
            rate_parameters[dust_surface_reaction_params::kVibrationFrequency1] = dust_surface_species1->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kVibrationFrequency2] = dust_surface_species2->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kDiffusionBarrier1] = dust_surface_species1->GetDiffusionBarrierOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kDiffusionBarrier2] = dust_surface_species2->GetDiffusionBarrierOnH2Oice();

            // Store the reaction data
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
            
            // Update previous reactants and products
            previous_reactants = reactants;
            previous_products = products;

            // Chemical desorption reactions
            if (is_chemical_desorption_) {
                if (products[1] == "" && products[2] == "") {
                    std::vector<std::string> products_cd = products;
                    // sA + sB → C
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );
                }

                if (products[1] != "" && products[2] == "") {
                    std::vector<std::string> products_cd;

                    // sA + sB → C + sD
                    products_cd     = products;
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB → sC + D
                    products_cd     = products;
                    products_cd[1]  = products_cd[1].substr(1); // sD -> D
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB → C + D
                    products_cd     = products;
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    products_cd[1]  = products_cd[1].substr(1); // sD -> D
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );
                }

                if (products[1] != "" && products[2] != "") {
                    std::vector<std::string> products_cd;

                    // sA + sB → sC + sD + E
                    products_cd     = products;
                    products_cd[2]  = products_cd[2].substr(1); // sE -> E
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB → sC + D + sE
                    products_cd     = products;
                    products_cd[1]  = products_cd[1].substr(1); // sD -> D
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB → sC + D +  E
                    products_cd     = products;
                    products_cd[1]  = products_cd[1].substr(1); // sD -> D
                    products_cd[2]  = products_cd[2].substr(1); // sE -> E
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB →  C + sD + sE
                    products_cd     = products;
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB →  C + sD +  E
                    products_cd     = products;
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    products_cd[2]  = products_cd[2].substr(1); // sE -> E
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB →  C +  D + sE
                    products_cd     = products;
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    products_cd[1]  = products_cd[1].substr(1); // sD -> D
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                    // sA + sB →  C +  D + E
                    products_cd     = products;
                    products_cd[0]  = products_cd[0].substr(1); // sC -> C
                    products_cd[1]  = products_cd[1].substr(1); // sD -> D
                    products_cd[1]  = products_cd[1].substr(1); // sE -> E
                    product_ids = GetSpeciesIDList(products_cd);
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );
                }
            } // End chemical desorption
        } // End while loop

        if (is_H2_desorption_) {
            // H2 desorption reaction on dust surface
            // sH2 + sH2 → sH2 + H2

            std::vector<std::string> reactants = {
                "sH2",
                "sH2",
                ""
            };

            std::vector<std::string> products = {
                "sH2",
                "H2",
                "",
                "",
                ""
            };

            if (!AreReactionSpeciesInSpeciesNameList(reactants, products)) return;

            // Get the indices of reactants and products in the species list
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Get DustSurfaceSpecies object for reactants
            std::shared_ptr<DustSurfaceSpecies> dust_surface_species1 = ptr_species_manager_->FindDustSurfaceSpeciesByName(reactants[0]);
            std::shared_ptr<DustSurfaceSpecies> dust_surface_species2 = ptr_species_manager_->FindDustSurfaceSpeciesByName(reactants[1]);

            // Set up reaction parameters
            std::vector<Real> rate_parameters(dust_surface_reaction_params::kNumParams);
            rate_parameters[dust_surface_reaction_params::kVibrationFrequency1] = dust_surface_species1->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kVibrationFrequency2] = dust_surface_species2->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kDiffusionBarrier1] = dust_surface_species1->GetDiffusionBarrierOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kDiffusionBarrier2] = dust_surface_species2->GetDiffusionBarrierOnH2Oice();

            // Store the reaction data
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        } // End H2 desorption
    }


    std::vector<std::string> ReactionManager::SplitSurfaceActivationEnergyLine(const std::string& line) const {
        std::string reactant1 = string_utils::Trim(line.substr( 0, 11));
        std::string reactant2 = string_utils::Trim(line.substr(11, 11));
        std::string reactant3 = string_utils::Trim(line.substr(22, 11));

        std::string product1  = string_utils::Trim(line.substr(34, 11));
        std::string product2  = string_utils::Trim(line.substr(45, 11));
        std::string product3  = string_utils::Trim(line.substr(56, 11));
        std::string product4  = string_utils::Trim(line.substr(67, 11));
        std::string product5  = string_utils::Trim(line.substr(78, 11));

        std::string rate_param1 = line.substr( 91, 11); // pre-exponential
        std::string rate_param2 = line.substr(102, 11); // activation barrier
        std::string rate_param3 = line.substr(113,  9); // uncertainty on the activation barrier

        return std::vector<std::string>{
            reactant1,
            reactant2,
            reactant3,
            product1,
            product2,
            product3,
            product4,
            product5,
            rate_param2
        };
    }



    Real ReactionManager::CalculateTunnellingEffectProbabilty(const Real mass1, const Real mass2, const Real activation_energy) const {
        const Real reduced_mass = mass1 * mass2 / (mass1 + mass2);
        // Calculate and return the tunneling effect probability using the WKB approximation
        return 2.0 * kWidthOfBarrier / constants::kDiracConstant *
            std::sqrt(2.0 * reduced_mass * constants::kBoltzmannConstant * activation_energy);
    }


    void ReactionManager::ReadSurfaceActivationEnergyFile(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << file_path << std::endl;
            return;
        }

        std::string line;
        std::vector<std::string> previous_reactants;
        std::vector<std::string> previous_products;
        while (std::getline(file, line)) {
            std::vector<std::string> split_result = SplitSurfaceActivationEnergyLine(line);
            if (split_result.size() < 9) {
                std::cerr << "Error: Insufficient data in line: " << line << std::endl;
                continue;
            }

            // Store reactants
            std::vector<std::string> reactants = {
                split_result[0],
                split_result[1],
                split_result[2]
            };

            // Check if reactants exist in the species list or SPECIAL_SPECIES_LIST
            bool skip = false;
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) continue;
                reactant = kDustSurfaceSpeciesPrefix + reactant;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) continue;
                skip = true;
            }
            if (skip) continue;

            // Store products
            std::vector<std::string> products = {
                split_result[3],
                split_result[4],
                split_result[5],
                split_result[6],
                split_result[7]
            };

            // Check if products exist in the species list
            skip = false;
            for (auto& product : products) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) continue;
                product = kDustSurfaceSpeciesPrefix + product;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) continue;
                skip = true;
            }
            if (skip) continue;

            // Prevent duplicate reactions by comparing reactants and products
            if (reactants == previous_reactants && products == previous_products) continue;

            // Get indices for reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids  = GetSpeciesIDList(products);

            // Get activation energy
            const Real activation_energy = std::stod(split_result[8]);

            // Search for the reaction in the dust surface reaction list
            for (auto& reaction : reaction_list_) {
                if (reaction->type_id_ != reaction_type_id::kDustSurfaceReaction) continue;

                // Check if reactant and product indices match
                if (reactant_ids == reaction->reactant_ids_ && product_ids  == reaction->product_ids_) {
                    // Update mass and rate parameters for the reaction
                    const Real mass1 = ptr_species_manager_->species_list_[reactant_ids[0]]->GetMass();
                    const Real mass2 = ptr_species_manager_->species_list_[reactant_ids[1]]->GetMass();
                    reaction->rate_parameters_[dust_surface_reaction_params::kActivationEnergy] = activation_energy;
                    reaction->rate_parameters_[dust_surface_reaction_params::kTunnelingEffectProbability]
                        = CalculateTunnellingEffectProbabilty(mass1, mass2, activation_energy);
                    
                    break; // Exit loop once the reaction is found
                }
            }

            // Update previous reactants and products for next iteration
            previous_reactants = reactants;
            previous_products  = products;
        }

        if (!is_three_phase_reaction_) return;

        // Reset the file pointer to the beginning
        file.clear();
        file.seekg(0);

        while (std::getline(file, line)) {
            std::vector<std::string> split_result = SplitSurfaceActivationEnergyLine(line);
            if (split_result.size() < 9) {
                std::cerr << "Error: Insufficient data in line: " << line << std::endl;
                continue;
            }

            // Store reactants
            std::vector<std::string> reactants = {
                split_result[0],
                split_result[1],
                split_result[2]
            };

            // Check if reactants exist in the species list
            bool skip = false;
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) continue;
                reactant = kDustMantleSpeciesPrefix + reactant;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) continue;
                skip = true;
            }
            if (skip) continue;

            // Store products
            std::vector<std::string> products = {
                split_result[3],
                split_result[4],
                split_result[5],
                split_result[6],
                split_result[7]
            };

            // Check if products exist in the species list
            skip = false;
            for (auto& product : products) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) continue;
                product = kDustMantleSpeciesPrefix + product;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) continue;
                skip = true;
            }
            if (skip) continue;

            // Prevent duplicate reactions by comparing reactants and products
            if (reactants == previous_reactants && products == previous_products) continue;

            // Get indices for reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids  = GetSpeciesIDList(products);

            // Get activation energy
            const Real activation_energy = std::stod(split_result[8]);

            // Search for the reaction in the dust mantle reaction list
            for (auto& reaction : reaction_list_) {
                if (reaction->type_id_ == reaction_type_id::kDustMantleReaction) continue;

                // Check if reactant and product indices match
                if (reactant_ids == reaction->reactant_ids_ && product_ids  == reaction->product_ids_) {
                    // Update mass and rate parameters for the reaction
                    const Real mass1 = ptr_species_manager_->species_list_[reactant_ids[0]]->GetMass();
                    const Real mass2 = ptr_species_manager_->species_list_[reactant_ids[1]]->GetMass();
                    reaction->rate_parameters_[dust_mantle_reaction_params::kActivationEnergy] = activation_energy;
                    reaction->rate_parameters_[dust_mantle_reaction_params::kTunnelingEffectProbability] 
                        = CalculateTunnellingEffectProbabilty(mass1, mass2, activation_energy);
                    
                    break;
                }
            }

            // Update previous reactants and products for next iteration
            previous_reactants = reactants;
            previous_products  = products;
        }
    }

    
    /**
     * @brief Generate a list of reactions for dust particle and charged gas particle collisions.
     *
     * This function iterates through each charged gas species and each dust species to generate
     * reaction channels based on the charge of the gas particle. It handles:
     * - Electron interactions: Electron + Dust(x) -> Dust(x-1)
     * - Ion interactions: M+ + Dust(x) -> M + Dust(x+1) or M- + Dust(x) -> M + Dust(x-1)
     * - Handling special cases where there are no corresponding neutral particles.
     *      For H3+, there is no corresponding neutral species (H3). In this case, look for gas-phase reactions 
     *      where dust replaces electrons, and substitute dust for electrons. In a gas-phase reaction H3+ + e- -> H2 + H, 
     *      the reaction would become: H3+ + Dust(x) -> H2 + H + Dust(x+1).
     *
     * It also ensures that dust particles with charges beyond the specified limits are excluded.
     * 
     * @ref
     * https://ui.adsabs.harvard.edu/abs/2006A%26A...445..205I/abstract
     * https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.5588K/abstract
     */
    void ReactionManager::GenerateReactionListForDustAndChargedParticleCollision() {
        const std::size_t type_id = reaction_type_id::kDustAndChargedGasParticleCollison;
        std::size_t index = reaction_list_.size();

        // Get the maximum and minimum charge numbers for dust particles
        const int max_dust_charge_number = ptr_species_manager_->dust_species_model_parameters_.GetMaxDustChargeNumber();
        const int min_dust_charge_number = ptr_species_manager_->dust_species_model_parameters_.GetMinDustChargeNumber();

        // Iterate over each gas species to process reactions with dust species
        for (auto const& gas_species : ptr_species_manager_->gas_species_list_) {

            // Skip neutral gas particles
            if (gas_species->GetCharge() == 0) continue;

            // Get gas species name and index
            const std::string& gas_species_name = gas_species->GetName();
            const std::size_t gas_species_index = ptr_species_manager_->FindSpeciesID(gas_species_name);

            // Handle electron interactions (e- + Dust(x) -> DUST(x-1))
            if (gas_species_name == kElectron) {
                // Iterate over each dust species
                for (auto const& dust_species : ptr_species_manager_->dust_species_list_) {
                    // Get the bin number and charge of the dust species
                    const std::size_t bin_number = dust_species->GetBinNumber();
                    const int dust_charge = dust_species->GetCharge();

                    // Skip dust species if its charge is beyond the minimum allowed charge
                    if (dust_charge == min_dust_charge_number) continue;

                    // Get the dust species with a charge that is one less than the current dust species' charge
                    const int dust_charge_m1 = dust_charge - 1;
                    const auto& dust_species_m1 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number, dust_charge_m1);
                    if (!dust_species_m1) continue;

                    // Set up the list of reactants for the reaction: gas species and dust species
                    std::vector<std::string> reactants = {
                        gas_species_name,
                        dust_species->GetName(),
                        ""
                    };

                    // Set up the list of products for the reaction: reduced dust species
                    std::vector<std::string> products = {
                        dust_species_m1->GetName(),
                        "",
                        "",
                        "",
                        ""
                    };
                    
                    // Get the indices of the reactants and products from their names
                    std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
                    std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

                    // Set up rate parameters for the reaction
                    std::vector<Real> rate_parameters(dust_and_charged_particle_collison_params::kNumParams);
                    rate_parameters[dust_and_charged_particle_collison_params::kGasMass] = gas_species->GetMass();
                    rate_parameters[dust_and_charged_particle_collison_params::kGasCharge] = gas_species->GetCharge();
                    rate_parameters[dust_and_charged_particle_collison_params::kDustCharge] = dust_species->GetCharge();
                    rate_parameters[dust_and_charged_particle_collison_params::kDustRadius] = dust_species->GetRadius();
                    rate_parameters[dust_and_charged_particle_collison_params::kDustCrossSection] = dust_species->GetCrossSection();
                    rate_parameters[dust_and_charged_particle_collison_params::kStickingProbability] = kStickingProbabilityElectron;

                    // Store the reaction in the reaction list
                    reaction_list_.emplace_back(std::make_shared<Reaction>(
                        ++index,
                        reactant_ids,
                        product_ids,
                        rate_parameters,
                        type_id
                        )
                    );

                }

                continue; // next gas species
            } // end if gas_species_name == "e-"

            // Handle reactions with other charged particles
            // Get neutral gas species name corresponding the current charged gas particle (ex. H+ and H)
            const std::string neutral_gas_species_name = gas_species_name.substr(0, gas_species_name.length() - 1);

            // Check if the corresponding neutral gas species exists in the list
            if (string_utils::IsInStringVector(ptr_species_manager_->gas_species_name_list_, neutral_gas_species_name)) {
                // If it exists, branch based on whether the gas species is positively or negatively charged
                if (gas_species->GetCharge() > 0) {
                    // For positive charges: M+ + Dust(x) -> M + Dust(x+1), charge exchange

                    // Iterate over each dust species
                    for (const auto& dust_species : ptr_species_manager_->dust_species_list_) {
                        // Get dust species information (bin number and charge)
                        const std::size_t bin_number = dust_species->GetBinNumber();
                        const int dust_charge = dust_species->GetCharge();

                        // Skip if dust charge is at the maximum positive value (dust cannot have charge greater than this in this code)
                        if (dust_charge == max_dust_charge_number) continue;

                        // Get the dust species with one charge unit greater than the current dust species
                        const int dust_charge_p1 = dust_charge + 1;
                        const auto& dust_species_p1 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number, dust_charge_p1);
                        if (!dust_species_p1) continue;

                        // Set up the reactants (gas species and dust species)
                        std::vector<std::string> reactants = {
                            gas_species_name,
                            dust_species->GetName(),
                            ""
                        };

                        // Set up the products (neutral gas species and dust species with increased charge)
                        std::vector<std::string> products = {
                            neutral_gas_species_name,
                            dust_species_p1->GetName(),
                            "",
                            "",
                            ""
                        };
                        
                        // Get the indices for the reactants and products
                        std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
                        std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

                        // Set up rate parameters for the reaction
                        std::vector<Real> rate_parameters(dust_and_charged_particle_collison_params::kNumParams);
                        rate_parameters[dust_and_charged_particle_collison_params::kGasMass] = gas_species->GetMass();
                        rate_parameters[dust_and_charged_particle_collison_params::kGasCharge] = gas_species->GetCharge();
                        rate_parameters[dust_and_charged_particle_collison_params::kDustCharge] = dust_species->GetCharge();
                        rate_parameters[dust_and_charged_particle_collison_params::kDustRadius] = dust_species->GetRadius();
                        rate_parameters[dust_and_charged_particle_collison_params::kDustCrossSection] = dust_species->GetCrossSection();
                        rate_parameters[dust_and_charged_particle_collison_params::kStickingProbability] = kStickingProbabilityIon;

                        // Store the reaction data
                        reaction_list_.emplace_back(std::make_shared<Reaction>(
                            ++index,
                            reactant_ids,
                            product_ids,
                            rate_parameters,
                            type_id
                            )
                        );
                    } // End loop over dust species (positive charge case)
                } else {
                    // For negative charges: M- + Dust(x) -> M + Dust(x-1), charge exchange

                    // Iterate over each dust species
                    for (const auto& dust_species : ptr_species_manager_->dust_species_list_) {
                        // Get dust species information (bin number and charge)
                        const std::size_t bin_number = dust_species->GetBinNumber();
                        const int dust_charge = dust_species->GetCharge();

                        // Skip if dust charge is at the minimum value (dust cannot have charge less than this in this code)
                        if (dust_charge == min_dust_charge_number) continue;

                        // Get the dust species with one charge unit less than the current dust species
                        const int dust_charge_m1 = dust_charge - 1;
                        const auto& dust_species_m1 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number, dust_charge_m1);
                        if (!dust_species_m1) continue;

                        // Set up the reactants (gas species and dust species)
                        std::vector<std::string> reactants = {
                            gas_species_name,
                            dust_species->GetName(),
                            ""
                        };

                        // Set up the products (neutral gas species and dust species with decreased charge)
                        std::vector<std::string> products = {
                            neutral_gas_species_name,
                            dust_species_m1->GetName(),
                            "",
                            "",
                            ""
                        };
                        
                        // Get the indices for the reactants and products
                        std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
                        std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

                        // Set up rate parameters for the reaction
                        std::vector<Real> rate_parameters(dust_and_charged_particle_collison_params::kNumParams);
                        rate_parameters[dust_and_charged_particle_collison_params::kGasMass] = gas_species->GetMass();
                        rate_parameters[dust_and_charged_particle_collison_params::kGasCharge] = gas_species->GetCharge();
                        rate_parameters[dust_and_charged_particle_collison_params::kDustCharge] = dust_species->GetCharge();
                        rate_parameters[dust_and_charged_particle_collison_params::kDustRadius] = dust_species->GetRadius();
                        rate_parameters[dust_and_charged_particle_collison_params::kDustCrossSection] = dust_species->GetCrossSection();
                        rate_parameters[dust_and_charged_particle_collison_params::kStickingProbability] = kStickingProbabilityIon;

                        // Store the reaction data
                        reaction_list_.emplace_back(std::make_shared<Reaction>(
                            ++index,
                            reactant_ids,
                            product_ids,
                            rate_parameters,
                            type_id
                            )
                        );

                    } // End loop over dust species (negative charge case)
                } // End charge check
            } else {
                // If the charged particle does not exist in the list of neutral particles 
                // (e.g., for H3+, there is no H3 in the considered species list)
                
                // Branch based on the charge of the charged particle
                if (gas_species->GetCharge() > 0) {
                    // Look for gas-phase reactions where dust replaces electrons, and substitute dust for electrons.
                    // For example, in a gas-phase reaction H3+ + e- -> H2 + H, the reaction would become:
                    // H3+ + Dust(x) -> H2 + H + Dust(x+1)

                    // Search for reactions where the charged particle interacts with electrons
                    for (std::size_t ireac = 0; ireac < total_number_of_gas_phase_reactions_; ++ireac) {

                        const auto reaction = reaction_list_[ireac];
                        if (!reaction) continue;

                        // Check the reactants
                        bool skip = true;
                        // Check if there is a reaction between the charged particle and electrons
                        if ((reaction->reactant_ids_[0] == gas_species_index && reaction->reactant_ids_[1] == ptr_species_manager_->id_electron_) ||
                            (reaction->reactant_ids_[1] == gas_species_index && reaction->reactant_ids_[0] == ptr_species_manager_->id_electron_)) {
                            skip = false;
                        }
                        if (skip) continue; // Skip if the reaction doesn't exist

                        // If it exists, loop over dust species to create the corresponding reaction
                        for (const auto& dust_species : ptr_species_manager_->dust_species_list_) {

                            // Get the information of the dust species interacting with the charged particle
                            const std::size_t bin_number = dust_species->GetBinNumber();
                            const int dust_charge = dust_species->GetCharge();

                            // If the dust charge is the maximum (positive maximum), we cannot have dust with a charge greater than that.
                            if (dust_charge == max_dust_charge_number) continue;

                            // Get the dust species with a charge one unit higher than the current dust species within the same bin。
                            const int dust_charge_p1 = dust_charge + 1;
                            const auto& dust_species_p1 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number, dust_charge_p1);
                            if (!dust_species_p1) continue;

                            // Get the indices of the reactants
                            std::size_t number_of_reactant = reaction->reactant_ids_.size();
                            std::vector<std::size_t> reactant_ids(number_of_reactant, kNotFoundSpecies);
                            // Replace the electron in the reactants with the dust species
                            for (std::size_t index = 0; index < number_of_reactant; ++index) {
                                if (reaction->reactant_ids_[index] == ptr_species_manager_->id_electron_) {
                                    reactant_ids[index] = ptr_species_manager_->FindSpeciesID(dust_species->GetName());
                                } else {
                                    reactant_ids[index] = reaction->reactant_ids_[index];
                                }
                            }

                            // Get the indices of the products
                            std::size_t number_of_product = reaction->product_ids_.size();
                            std::vector<std::size_t> product_ids(number_of_product, kNotFoundSpecies);
                            // Add the dust species after the collision
                            for (std::size_t index = 0; index < number_of_product; ++index) {
                                if (reaction->product_ids_[index] == kNotFoundSpecies) {
                                    product_ids[index] = ptr_species_manager_->FindSpeciesID(dust_species_p1->GetName());
                                    break;
                                } else {
                                    product_ids[index] = reaction->product_ids_[index];
                                }
                            }

                            // Set the values related to the reaction
                            std::vector<Real> rate_parameters(dust_and_charged_particle_collison_params::kNumParams);
                            rate_parameters[dust_and_charged_particle_collison_params::kGasMass] = gas_species->GetMass();
                            rate_parameters[dust_and_charged_particle_collison_params::kGasCharge] = gas_species->GetCharge();
                            rate_parameters[dust_and_charged_particle_collison_params::kDustCharge] = dust_species->GetCharge();
                            rate_parameters[dust_and_charged_particle_collison_params::kDustRadius] = dust_species->GetRadius();
                            rate_parameters[dust_and_charged_particle_collison_params::kDustCrossSection] = dust_species->GetCrossSection();
                            rate_parameters[dust_and_charged_particle_collison_params::kStickingProbability] = kStickingProbabilityIon;

                            // Store the reaction data
                            reaction_list_.emplace_back(std::make_shared<Reaction>(
                                ++index,
                                reactant_ids,
                                product_ids,
                                rate_parameters,
                                type_id
                                )
                            );

                        } // end for dust_species

                    } // end for ireac

                }

                // Anions are not considered in this case
                // (There is no processing for anions here)

            } // end else (gas_species does not have a corresponding neutral species)
        } // end for gas_species
    }

    /**
     * @brief Generate a reaction list for dust collisions
     * 
     * This method generates reactions when dust particles collide. It ignores dust coagulation (dust growth) after collisions.
     * During collisions, if dust particles have different charge signs, they will exchange charge and neutralize. 
     * Example reactions:
     *      Dust(+)  + Dust(-)  -> Dust(0)  + Dust(0)
     *      Dust(2+) + Dust(-)  -> Dust(+)  + Dust(0)
     *      Dust(+)  + Dust(2-) -> Dust(0)  + Dust(-)
     *      Dust(3+) + Dust(-)  -> Dust(2+) + Dust(0)
     * 
     * @ref
     * https://ui.adsabs.harvard.edu/abs/1990MNRAS.243..103U/abstract
     * https://ui.adsabs.harvard.edu/abs/2006A%26A...445..205I/abstract
     * https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.5588K/abstract
     */
    void ReactionManager::GenerateReactionListForDustCollision() {
        std::size_t index = reaction_list_.size();
        std::size_t type_id = reaction_type_id::kDustCollision;

        // Get dust information
        const std::size_t number_of_dust_bins = ptr_species_manager_->dust_species_model_parameters_.GetNumberOfDustBins();
        const int max_dust_charge_number = ptr_species_manager_->dust_species_model_parameters_.GetMaxDustChargeNumber();
        const int min_dust_charge_number = ptr_species_manager_->dust_species_model_parameters_.GetMinDustChargeNumber();

        // Loop through all combinations of dust particles
        for (std::size_t bin_number1 = 1; bin_number1 <= number_of_dust_bins; ++bin_number1) { // Collision dust 1

            // Loop through positive dust
            for (int dust_charge1 = 1; dust_charge1 <= max_dust_charge_number; ++dust_charge1) {

                const auto& reactant_dust_species1 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number1, dust_charge1);
                if (!reactant_dust_species1) continue;

                for (std::size_t bin_number2 = 1; bin_number2 <= number_of_dust_bins; ++bin_number2) { // Collision dust 2

                    // Loop through negative dust
                    for (int dust_charge2 = -1; dust_charge2 >= min_dust_charge_number; --dust_charge2) { // negativeダストでloop

                        const auto& reactant_dust_species2 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number2, dust_charge2);
                        if (!reactant_dust_species2) continue;

                        std::vector<std::string> reactants = {
                            reactant_dust_species1->GetName(),
                            reactant_dust_species2->GetName(),
                            ""
                        };

                        int net_charge = reactant_dust_species1->GetCharge() + reactant_dust_species2->GetCharge();
                        int product_dust_charge1, product_dust_charge2;
                        // Calculate the charges of the product dust species based on net charge
                        if (net_charge < 0) {
                            product_dust_charge1 = 0;
                            product_dust_charge2 = net_charge;
                        } else if (net_charge == 0) {
                            product_dust_charge1 = 0;
                            product_dust_charge2 = 0;
                        } else {
                            product_dust_charge1 = net_charge;
                            product_dust_charge2 = 0;
                        }

                        const auto& product_dust_species1 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number1, product_dust_charge1);
                        const auto& product_dust_species2 = ptr_species_manager_->FindDustSpeciesByBinNumberAndCharge(bin_number2, product_dust_charge2);

                        std::vector<std::string> products = {
                            product_dust_species1->GetName(),
                            product_dust_species2->GetName(),
                            "",
                            "",
                            ""
                        };

                        std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
                        std::vector<std::size_t> product_ids  = GetSpeciesIDList(products);

                        std::vector<Real> rate_parameters(dust_collision_params::kNumParams);
                        rate_parameters[dust_collision_params::kDustRadius1] = reactant_dust_species1->GetRadius();
                        rate_parameters[dust_collision_params::kDustRadius2] = reactant_dust_species2->GetRadius();
                        rate_parameters[dust_collision_params::kDustCharge1] = static_cast<Real>(reactant_dust_species1->GetCharge()) * constants::kChargeUnit;
                        rate_parameters[dust_collision_params::kDustCharge2] = static_cast<Real>(reactant_dust_species2->GetCharge()) * constants::kChargeUnit;
                        rate_parameters[dust_collision_params::kDustMass1] = reactant_dust_species1->GetMass();
                        rate_parameters[dust_collision_params::kDustMass2] = reactant_dust_species2->GetMass();

                        // Store the reaction
                        reaction_list_.emplace_back(std::make_shared<Reaction>(
                            ++index,
                            reactant_ids,
                            product_ids,
                            rate_parameters,
                            type_id
                            )
                        );

                    } // end for charge 2
                } // end for bin_number 2

            } // end for charge 2
        } // end for bin_number 2
    }

    /**
     * @brief Generate a reaction list for neutral gas species accretion onto dust surfaces
     * 
     * The reaction represents the accretion of neutral gas species onto dust surfaces.
     * Example:
     *      M -> sM
     */
    void ReactionManager::GenerateReactionListForNeutralSpeciesAccretionOnDustSurfaces() {
        std::size_t index = reaction_list_.size();
        const std::size_t type_id = reaction_type_id::kAccretionGasParticleOnDustSurfaces;

        // Loop over all dust surface species
        for (const auto& dust_surface_species : ptr_species_manager_->dust_surface_species_list_) {

            // Get the corresponding gas species for the dust surface species
            const auto& corresponding_gas_species = dust_surface_species->GetCorrespondingGasSpecies();

            // Set of reactants
            std::vector<std::string> reactants = {
                corresponding_gas_species->GetName(),
                "",
                ""
            };

            // Set of products
            std::vector<std::string> products = {
                dust_surface_species->GetName(),
                "",
                "",
                "",
                ""
            };
                    
            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Reaction parameters
            std::vector<Real> rate_parameters(accretion_on_dusts_params::kNumParams); 

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));

        } // end for dust_surface_species
    }

    /**
     * @brief Generate a reaction list for thermal desorption on dust surfaces
     * 
     * The reaction represents thermal desorption of species from dust surfaces.
     * Example:
     *      sM -> M
     */
    void ReactionManager::GenerateReactionListForThermalDesorptionOnDustSurfaces() {
        std::size_t index = reaction_list_.size();
        const std::size_t type_id = reaction_type_id::kThermalDesorptionOnDustSurfaces;

        // Loop over all dust surface species
        for (const auto& dust_surface_species : ptr_species_manager_->dust_surface_species_list_) {

            // Get the corresponding gas species for the dust surface species
            const auto& corresponding_gas_species = dust_surface_species->GetCorrespondingGasSpecies();

            // Set of reactants
            std::vector<std::string> reactants = {
                dust_surface_species->GetName(),
                "",
                ""
            };

            // Set of products
            std::vector<std::string> products = {
                corresponding_gas_species->GetName(),
                "",
                "",
                "",
                ""
            };
                    
            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Reaction parameters
            std::vector<Real> rate_parameters(thermal_desorption_params::kNumParams);
            rate_parameters[thermal_desorption_params::kVibrationFrequency] = dust_surface_species->GetVibrationFrequencyOnH2Oice();
            rate_parameters[thermal_desorption_params::kBindingEnergyOnH2Oice] = dust_surface_species->GetBindingEnergyOnH2Oice();

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));

        } // end for dust_surface_species
    }

    /**
     * @brief Generate a reaction list for cosmic ray desorption on dust surfaces
     * 
     * This reaction represents the desorption of species from dust surfaces due to cosmic rays.
     * Example:
     *      sM -> M
     */
    void ReactionManager::GenerateReactionListForCosmicRayDesorptionOnDustSurfaces() {
        std::size_t index = reaction_list_.size();
        const std::size_t type_id = reaction_type_id::kCosmicRayDesorptionOnDustSurfaces;

        // Loop over all dust surface species
        for (const auto& dust_surface_species : ptr_species_manager_->dust_surface_species_list_) {

            // Get the corresponding gas species for the dust surface species
            const auto& corresponding_gas_species = dust_surface_species->GetCorrespondingGasSpecies();

            // Set of reactants
            std::vector<std::string> reactants = {
                dust_surface_species->GetName(),
                "",
                ""
            };

            // Set of products
            std::vector<std::string> products = {
                corresponding_gas_species->GetName(),
                "",
                "",
                "",
                ""
            };
                    
            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Reaction parameters
            std::vector<Real> rate_parameters(cosmic_ray_desorption_on_dusts_surfaces_params::kNumPramas);
            rate_parameters[cosmic_ray_desorption_on_dusts_surfaces_params::kVibrationFrequency] = dust_surface_species->GetVibrationFrequencyOnH2Oice();
            rate_parameters[cosmic_ray_desorption_on_dusts_surfaces_params::kBindingEnergyOnH2Oice] = dust_surface_species->GetBindingEnergyOnH2Oice();

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));

        } // end for dust_surface_species
    }

    /**
     * @brief Generate a reaction list for photo desorption by external UV
     * 
     * This reaction represents the photo desorption of species from dust surfaces by external UV radiation.
     * Example:
     *      sM -> M
     */
    void ReactionManager::GenerateReactionListForPhotoDesorptionByExternalUV() {
        std::size_t index = reaction_list_.size();
        const std::size_t type_id = reaction_type_id::kPhotoDesorptionByExternalUV;

        // Loop over all dust surface species
        for (const auto& dust_surface_species : ptr_species_manager_->dust_surface_species_list_) {

            // Get the corresponding gas species for the dust surface species
            const auto& corresponding_gas_species = dust_surface_species->GetCorrespondingGasSpecies();

            // Set of reactants
            std::vector<std::string> reactants = {
                dust_surface_species->GetName(),
                "",
                ""
            };

            // Set of products
            std::vector<std::string> products = {
                corresponding_gas_species->GetName(),
                "",
                "",
                "",
                ""
            };
                
            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Reaction parameters
            std::vector<Real> rate_parameters(photo_desorption_by_external_UV_params::kNumParams);

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));

        } // end for dust_surface_species
    }

    /**
     * @brief Generate a reaction list for photo desorption by cosmic ray generated UV
     * 
     * This reaction represents the photo desorption of species from dust surfaces induced by UV radiation
     * generated by cosmic rays.
     * Example:
     *      sM -> M
     */
    void ReactionManager::GenerateReactionListForPhotoDesorptionByCosmicRayGeneratedUV() {
        std::size_t index = reaction_list_.size();
        const std::size_t type_id = reaction_type_id::kPhotoDesorptionByCRGeneratedUV;

        // Loop over all dust surface species
        for (const auto& dust_surface_species : ptr_species_manager_->dust_surface_species_list_) {

            // Get the corresponding gas species for the dust surface species
            const auto& corresponding_gas_species = dust_surface_species->GetCorrespondingGasSpecies();

            // Set of reactants
            std::vector<std::string> reactants = {
                dust_surface_species->GetName(),
                "",
                ""
            };

            // Set of products
            std::vector<std::string> products = {
                corresponding_gas_species->GetName(),
                "",
                "",
                "",
                ""
            };
                
            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Reaction parameters
            std::vector<Real> rate_parameters(photo_desorption_by_CR_generated_UV_params::kNumParams);

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));

        } // end for dust_surface_species
    }

    /**
     * @brief Generate a reaction list for photo dissociation induced by cosmic rays on dust surfaces.
     * 
     * This reaction represents the dissociation of species on dust surfaces caused by photoionization 
     * or photodissociation induced by cosmic rays.
     * The reaction list generated here is the same as the gas-phase photo dissociation reaction induced 
     * by cosmic rays, but with the reactants and products replaced by the corresponding dust surface species 
     * of the gas-phase species.
     */
    void ReactionManager::GenerateReactionListForPhotoDissociationInducedByCRsOnDustSurfaces() {
        std::size_t index = reaction_list_.size();
        const std::size_t type_id = reaction_type_id::kPhotoDissociationByCROnDustSurfaces;

        // Loop through the reaction list to generate reactions for dust surface species
        for (const auto& reaction : reaction_list_) {

            // Only consider reactions of type kGasPhase2
            if (reaction->type_id_ != reaction_type_id::kGasPhase2) continue;

            // Get the list of reactants based on the indices
            std::vector<std::string> reactants = GetSpeciesNameList(reaction->reactant_ids_);

            bool skip = false;
            // Check each reactant, if it's not special, prepend dust surface species prefix
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) continue;
                reactant = kDustSurfaceSpeciesPrefix + reactant;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) continue;
                skip = true;
            }
            if (skip) continue;

            // Get the list of products based on the indices
            std::vector<std::string> products = GetSpeciesNameList(reaction->product_ids_);

            skip = false;
            // Check each product, if it's not special, prepend dust surface species prefix
            for (auto& product : products) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) continue;
                product = kDustSurfaceSpeciesPrefix + product;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) continue;
                skip = true;
            }
            if (skip) continue;

            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Set the reaction rate parameters
            std::vector<Real> rate_parameters(photo_dissociation_by_CR_on_dusts_params::kNumParams);
            rate_parameters[photo_dissociation_by_CR_on_dusts_params::kAlpha] = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
            rate_parameters[photo_dissociation_by_CR_on_dusts_params::kBeta] = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
            rate_parameters[photo_dissociation_by_CR_on_dusts_params::kGamma] = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));

        } // end for reactions

        // If three-phase reaction is not enabled, return
        if (!is_three_phase_reaction_) return;

        // If three-phase reactions are enabled, generate reactions for dust mantle species
        for (const auto& reaction : reaction_list_) {

            // Only consider reactions of type kGasPhase2
            if (reaction->type_id_ != reaction_type_id::kGasPhase2) continue;

            // Get the list of reactants based on the indices
            std::vector<std::string> reactants = GetSpeciesNameList(reaction->reactant_ids_);

            bool skip = false;
            // Check each reactant, if it's not special, prepend dust mantle species prefix
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) continue;
                reactant = kDustMantleSpeciesPrefix + reactant;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) continue;
                skip = true;
            }
            if (skip) continue;

            // Get the list of products based on the indices
            std::vector<std::string> products = GetSpeciesNameList(reaction->product_ids_);

            skip = false;
            // Check each product, if it's not special, prepend dust mantle species prefix
            for (auto& product : products) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) continue;
                product = kDustMantleSpeciesPrefix + product;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) continue;
                skip = true;
            }
            if (skip) continue;

            // Get the indices of the reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Set the reaction rate parameters
            std::vector<Real> rate_parameters(photo_dissociation_by_CR_on_dusts_params::kNumParams);
            rate_parameters[photo_dissociation_by_CR_on_dusts_params::kAlpha] = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
            rate_parameters[photo_dissociation_by_CR_on_dusts_params::kBeta] = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
            rate_parameters[photo_dissociation_by_CR_on_dusts_params::kGamma] = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];

            // Store the reaction in the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
            ));
                    
        } // end for reactions
    }

    /**
     * @brief Generate a reaction list for photo dissociation by external UV on dust surfaces.
     * 
     * This reaction represents the dissociation of species on dust surfaces caused by photodissociation 
     * induced by external UV radiation. The reaction list generated here mirrors the gas-phase 
     * photo dissociation reaction induced by external UV, but with the reactants and products replaced 
     * by the corresponding dust surface species of the gas-phase species.
     */
    void ReactionManager::GenerateReactionListForPhotoDissociationByExternalUVOnDustSurfaces() {
        // Initialize the index for new reactions
        std::size_t index = reaction_list_.size();
        // Define the reaction type for photo dissociation by external UV on dust surfaces
        std::size_t type_id = reaction_type_id::kPhotoDissociationByUVOnDustSurfaces;

        // Loop through each reaction in the existing reaction list
        for (const auto& reaction : reaction_list_) {

            // Skip reactions that are not of type kGasPhase3 (gas phase photo dissociation)
            if (reaction->type_id_ != reaction_type_id::kGasPhase3) continue;

            // Retrieve the list of reactants (species involved in the reaction)
            std::vector<std::string> reactants = GetSpeciesNameList(reaction->reactant_ids_);

            // Check if reactants can be converted to corresponding dust surface species
            bool skip = false;
            for (auto& reactant : reactants) {
                // Skip special species (species listed in SPECIAL_SPECIES_LIST)
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) continue;
                // Add the dust surface prefix to the reactant species name
                reactant = kDustSurfaceSpeciesPrefix + reactant;
                // If the dust surface species does not exist, skip this reaction
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) continue;
                skip = true;
            }
            if (skip) continue;

            // Retrieve the list of products (species formed by the reaction)
            std::vector<std::string> products = GetSpeciesNameList(reaction->product_ids_);

            // Check if products can be converted to corresponding dust surface species
            skip = false;
            for (auto& product : products) {
                // Skip special species
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) continue;
                // Add the dust surface prefix to the product species name
                product = kDustSurfaceSpeciesPrefix + product;
                // If the dust surface species does not exist, skip this reaction
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) continue;
                skip = true;
            }
            if (skip) continue;

            // Get indices for the reactants and products in the species list
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Set up rate parameters for the reaction
            std::vector<Real> rate_parameters(photo_dissociation_by_UV_on_dusts_params::kNumParams);
            // Assign values from the original gas phase reaction parameters
            rate_parameters[photo_dissociation_by_UV_on_dusts_params::kAlpha] = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
            rate_parameters[photo_dissociation_by_UV_on_dusts_params::kBeta] = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
            rate_parameters[photo_dissociation_by_UV_on_dusts_params::kGamma] = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];

            // Add the newly created reaction to the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        } // End loop over reactions

        // If three-phase reactions are not enabled, skip further processing
        if (!is_three_phase_reaction_) return;

        // Process reactions for dust mantle species in a similar way as for dust surface species
        for (const auto& reaction : reaction_list_) {

            // Skip reactions that are not of type kGasPhase3
            if (reaction->type_id_ != reaction_type_id::kGasPhase3) continue;

            // Retrieve the list of reactants for this reaction
            std::vector<std::string> reactants = GetSpeciesNameList(reaction->reactant_ids_);

            // Check if reactants can be converted to corresponding dust mantle species
            bool skip = false;
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, reactant)) continue;
                reactant = kDustMantleSpeciesPrefix + reactant;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) continue;
                skip = true;
            }
            if (skip) continue;

            // Retrieve the list of products for this reaction
            std::vector<std::string> products = GetSpeciesNameList(reaction->product_ids_);

            // Check if products can be converted to corresponding dust mantle species
            skip = false;
            for (auto& product : products) {
                if (string_utils::IsInStringVector(SPECIAL_SPECIES_LIST, product)) continue;
                product = kDustMantleSpeciesPrefix + product;
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) continue;
                skip = true;
            }
            if (skip) continue;

            // Get indices for the dust mantle species reactants and products
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Set up rate parameters for the reaction for dust mantle species
            std::vector<Real> rate_parameters(photo_dissociation_by_UV_on_dusts_params::kNumParams);
            rate_parameters[photo_dissociation_by_UV_on_dusts_params::kAlpha] = reaction->rate_parameters_[gas_phase_reaction_params::kAlpha];
            rate_parameters[photo_dissociation_by_UV_on_dusts_params::kBeta] = reaction->rate_parameters_[gas_phase_reaction_params::kBeta];
            rate_parameters[photo_dissociation_by_UV_on_dusts_params::kGamma] = reaction->rate_parameters_[gas_phase_reaction_params::kGamma];

            // Add the reaction for dust mantle species to the list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        } // End loop for dust mantle species reactions
    }

    /**
     * @brief Generate reaction list for dust mantle reactions.
     * 
     * This reaction list mirrors the reaction list for dust surface reactions, 
     * but with the reactants and products replaced by the corresponding dust mantle species 
     * of the dust surface species. Additionally, chemical desorption does not occur in the dust mantle.
     */
    void ReactionManager::GenerateReactionListForDustMantleReaction() {
        // Initialize the index for new reactions and the reaction type ID for dust mantle reactions
        std::size_t index = reaction_list_.size();
        std::size_t type_id = reaction_type_id::kDustMantleReaction;

        // Loop through each reaction in the existing reaction list
        for (const auto& reaction : reaction_list_) {

            // Skip reactions that are not of type kDustSurfaceReaction (dust surface reactions)
            if (reaction->type_id_ != reaction_type_id::kDustSurfaceReaction) continue;

            // Retrieve the list of reactants and products for the current reaction
            std::vector<std::string> reactants = GetSpeciesNameList(reaction->reactant_ids_);
            std::vector<std::string> products = GetSpeciesNameList(reaction->product_ids_);

            // Check if chemical desorption is enabled
            if (is_chemical_desorption_) {
                // In the dust mantle, chemical desorption does not occur. Therefore, check 
                // if any product is a gas-phase species (which would be generated by desorption).
                bool skip = false;
                for (const auto& product : products) {
                    // If a product is a gas-phase species, skip this reaction
                    if (string_utils::IsInStringVector(ptr_species_manager_->gas_species_name_list_, product)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue; // Skip this reaction if a gas-phase product is found
            }

            // Convert dust surface species names to corresponding dust mantle species
            for (auto& reactant : reactants) {
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, reactant)) {
                    // Replace the species name with the dust mantle version
                    reactant = kDustMantleSpeciesPrefix + reactant.substr(1); 
                }
            }

            // Similarly, convert product species names to corresponding dust mantle species
            for (auto& product : products) {
                if (string_utils::IsInStringVector(ptr_species_manager_->species_name_list_, product)) {
                    // Replace the species name with the dust mantle version
                    product = kDustMantleSpeciesPrefix + product.substr(1); 
                }
            }

            // Get the indices for the reactants and products in the species list
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids = GetSpeciesIDList(products);

            // Retrieve the corresponding dust mantle species objects for the first two reactants
            std::shared_ptr<DustMantleSpecies> dust_mantle_species1 = ptr_species_manager_->FindDustMantleSpeciesByName(reactants[0]);
            std::shared_ptr<DustMantleSpecies> dust_mantle_species2 = ptr_species_manager_->FindDustMantleSpeciesByName(reactants[1]);

            // Set up rate parameters specific to dust mantle reactions
            std::vector<Real> rate_parameters(dust_mantle_reaction_params::kNumParams);

            // Assign specific parameters from the dust mantle species, such as vibration frequencies and diffusion barriers
            rate_parameters[dust_surface_reaction_params::kVibrationFrequency1] = dust_mantle_species1->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kVibrationFrequency2] = dust_mantle_species2->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kDiffusionBarrier1] = dust_mantle_species1->GetDiffusionBarrierOnH2Oice();
            rate_parameters[dust_surface_reaction_params::kDiffusionBarrier2] = dust_mantle_species2->GetDiffusionBarrierOnH2Oice();

            // Add the newly created dust mantle reaction to the reaction list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        } // End loop over dust surface reactions
    }

    /**
     * @brief Generate reaction list for dust surface to mantle swapping.
     * 
     * This reaction represents the swapping of species from the dust surface layer to the dust mantle layer.
     * The reaction occurs from surface species (sM) to mantle species (mM).
     */
    void ReactionManager::GenerateReactionListForDustSurfaceToMantleSwapping() {
        // Initialize the index for new reactions and the reaction type ID for dust surface to mantle swapping
        std::size_t index = reaction_list_.size();
        std::size_t type_id = reaction_type_id::kDustSurfaceToMantleSwapping;

        // Loop through each species in the dust surface species list
        for (const auto& dust_surface_species : ptr_species_manager_->dust_surface_species_list_) {

            // Retrieve the corresponding dust mantle species for the current dust surface species
            const auto& corresponding_dust_mantle_species  = dust_surface_species->GetCorrespondingDustMantleSpecies();

            // Set up the reactants for this swapping reaction: the dust surface species
            std::vector<std::string> reactants = {
                dust_surface_species->GetName(),
                "",
                ""
            };

            // Set up the products for this swapping reaction: the corresponding dust mantle species
            std::vector<std::string> products = {
                corresponding_dust_mantle_species->GetName(),
                "",
                "",
                "",
                ""
            };

            // Get the indices of the reactants and products in the species list
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids  = GetSpeciesIDList(products);

            // Define rate parameters for the swapping reaction (currently empty)
            std::vector<Real> rate_parameters(dust_surface_to_mantle_params::kNumPramas);

            // Add the new reaction for dust surface to mantle swapping to the reaction list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        } // End loop over dust surface species list
    }

    /**
     * @brief Generate reaction list for dust mantle to surface swapping.
     * 
     * This reaction represents the swapping of species from the dust mantle layer to the dust surface layer.
     * The reaction occurs from mantle species (mM) to surface species (sM).
     */
    void ReactionManager::GenerateReactionListForDustMantleToSurfaceSwapping() {
        // Initialize the index for new reactions and the reaction type ID for dust mantle to surface swapping
        std::size_t index = reaction_list_.size();
        std::size_t type_id = reaction_type_id::kDustMantleToSurfaceSwapping;

        // Loop through each species in the dust mantle species list
        for (const auto& dust_mantle_species : ptr_species_manager_->dust_mantle_species_list_) {

            // Retrieve the corresponding dust surface species for the current dust mantle species
            const auto& corresponding_dust_surface_species = dust_mantle_species->GetCorrespondingDustSurfaceSpecies();

            // Set up the reactants for this swapping reaction: the dust mantle species
            std::vector<std::string> reactants = {
                dust_mantle_species->GetName(),
                "",
                ""
            };

            // Set up the products for this swapping reaction: the corresponding dust surface species
            std::vector<std::string> products = {
                corresponding_dust_surface_species->GetName(),
                "",
                "",
                "",
                ""
            };

            // Get the indices of the reactants and products in the species list
            std::vector<std::size_t> reactant_ids = GetSpeciesIDList(reactants);
            std::vector<std::size_t> product_ids  = GetSpeciesIDList(products);

            // Define rate parameters for the swapping reaction, which include vibration frequency and binding energy
            std::vector<Real> rate_parameters(dust_mantle_to_surface_params::kNumPramas);
            rate_parameters[dust_mantle_to_surface_params::kVibrationFrequency] = dust_mantle_species->GetVibrationFrequencyOnH2Oice();
            rate_parameters[dust_mantle_to_surface_params::kBindingEnergyOnH2Oice] = dust_mantle_species->GetBindingEnergyOnH2Oice();

            // Add the new reaction for dust mantle to surface swapping to the reaction list
            reaction_list_.emplace_back(std::make_shared<Reaction>(
                ++index,
                reactant_ids,
                product_ids,
                rate_parameters,
                type_id
                )
            );
        } // End loop over dust mantle species list
    }

    /**
     * @brief Calculate reaction branching ratio.
     * 
     * This function calculates the branching ratio for different types of reactions, 
     * including dust and charged gas particle collisions, dust surface reactions, 
     * and dust mantle reactions. The branching ratio is estimated either 
     * numerically or through laboratory experiments. However, such estimations 
     * can be difficult to make. In this code, the branching ratio is calculated 
     * by counting the number of reactions with identical reactants and dividing the count.
     * This means that the branching ratio of the reactions with identical reactants 
     * are assumed to have the same probability.
     * 
     * @ref
     * https://ui.adsabs.harvard.edu/abs/2006A%26A...445..205I/abstract
     * https://ui.adsabs.harvard.edu/abs/2021MNRAS.504.5588K/abstract
     */
    void ReactionManager::CalculateReactionBranchingRatio() {
        // Count the branching ratio for dust and charged particle collisions
        const size_t ireac_start = reaction_type_id_start_[reaction_type_id::kDustAndChargedGasParticleCollison];
        const size_t ireac_end = reaction_type_id_end_[reaction_type_id::kDustAndChargedGasParticleCollison];

        // Loop through reactions of dust and charged gas particle collisions
        for (std::size_t ireac = ireac_start; ireac <= ireac_end; ++ireac) {

            const auto& ith_reaction = reaction_list_[ireac];
            int number_of_branch = 0;
            std::size_t idx_ith_r1 = ith_reaction->reactant_ids_[0];
            std::size_t idx_ith_r2 = ith_reaction->reactant_ids_[1];

            // Count the number of reactions with matching reactants
            for (std::size_t jreac = ireac_start; jreac <= ireac_end; ++jreac) {

                const auto& jth_reaction = reaction_list_[jreac];
                std::size_t idx_jth_r1 = jth_reaction->reactant_ids_[0];
                std::size_t idx_jth_r2 = jth_reaction->reactant_ids_[1];

                // Check for matching reactants
                if ((idx_ith_r1 == idx_jth_r1 && idx_ith_r2 == idx_jth_r2) ||
                    (idx_ith_r1 == idx_jth_r2 && idx_ith_r2 == idx_jth_r1)) {
                    number_of_branch++;
                }
            }

            // Calculate the branching ratio based on the number of branches
            if (number_of_branch > 0) {
                ith_reaction->branching_ratio_ = 1.0 / static_cast<Real>(number_of_branch);
            } else {
                ith_reaction->branching_ratio_ = 0.0;
            }
        }

        // Count the branching ratio for dust surface reactions and dust mantle reactions
        for (const auto& ith_reaction : reaction_list_) {
            
            // Skip reactions that are neither dust surface nor dust mantle reactions
            if (ith_reaction->type_id_ != reaction_type_id::kDustSurfaceReaction &&
                ith_reaction->type_id_ != reaction_type_id::kDustMantleReaction) {
                continue;
            }

            ith_reaction->branching_ratio_ = 1.0;

            int number_of_branch = 0; // Initialize the number of branches
            std::size_t idx_ith_r1 = ith_reaction->reactant_ids_[0];
            std::size_t idx_ith_r2 = ith_reaction->reactant_ids_[1];

            // Count the branching reactions
            for (const auto& jth_reaction : reaction_list_) {

                // Skip reactions that are neither dust surface nor dust mantle reactions
                if (jth_reaction->type_id_ != reaction_type_id::kDustSurfaceReaction &&
                    jth_reaction->type_id_ != reaction_type_id::kDustMantleReaction) {
                    continue;
                }

                std::size_t idx_jth_r1 = jth_reaction->reactant_ids_[0];
                std::size_t idx_jth_r2 = jth_reaction->reactant_ids_[1];

                // Check for matching reactants
                if ((idx_ith_r1 == idx_jth_r1 && idx_ith_r2 == idx_jth_r2) || 
                    (idx_ith_r1 == idx_jth_r2 && idx_ith_r2 == idx_jth_r1)) {

                    std::size_t idx_jth_p1 = jth_reaction->product_ids_[0];
                    std::size_t idx_jth_p2 = jth_reaction->product_ids_[1];
                    std::size_t idx_jth_p3 = jth_reaction->product_ids_[2];

                    // Count branching reactions with surface or mantle products
                    if (idx_jth_p2 == kNotFoundSpecies && idx_jth_p3 == kNotFoundSpecies) {
                        if (ptr_species_manager_->IsDustSurfaceOrDustMantleSpecies(idx_jth_p1)) {
                            number_of_branch++;
                        }
                    } else if (idx_jth_p3 == kNotFoundSpecies) {
                        if (ptr_species_manager_->IsDustSurfaceOrDustMantleSpecies(idx_jth_p1) &&
                            ptr_species_manager_->IsDustSurfaceOrDustMantleSpecies(idx_jth_p2)) {
                            number_of_branch++;
                        }
                    } else {
                        if (ptr_species_manager_->IsDustSurfaceOrDustMantleSpecies(idx_jth_p1) &&
                            ptr_species_manager_->IsDustSurfaceOrDustMantleSpecies(idx_jth_p2) &&
                            ptr_species_manager_->IsDustSurfaceOrDustMantleSpecies(idx_jth_p3)) {
                            number_of_branch++;
                        }
                    }

                } // End of matching reactants check
            } // End of loop for jth_reaction

            // Calculate the branching ratio for dust surface and dust mantle reactions
            if (number_of_branch > 0) {
                ith_reaction->branching_ratio_ = 1.0 / static_cast<Real>(number_of_branch);
                // For reactions with identical reactants, halve the branching ratio
                if (idx_ith_r1 == idx_ith_r2) ith_reaction->branching_ratio_ *= 0.5;
            } else {
                ith_reaction->branching_ratio_ = 0.0;
            }

        } // End of loop for ith_reaction
    }

    /**
     * @brief Calculate chemical desorption probabilities
     * 
     * This function calculates the chemical desorption probabilities for dust surface reactions. 
     * The desorption occurs when reactions, such as:
     *  sA + sB → sC + sD
     *  sA + sB → C + D
     * exhibit exothermicity, where the released energy contributes to desorption. 
     * The exothermicity is evaluated by the enthalpy of formation of the reactants and products. 
     * In this implementation, chemical desorption probabilities are specifically calculated 
     * for H2O ice dust surfaces and bare silicate dust surfaces.
     */
    void ReactionManager::CalculateChemicalDesorptionProbabilities() {
        for (auto& reaction : reaction_list_) {

            // Skip if the reaction is not a dust surface reaction
            if (reaction->type_id_ != reaction_type_id::kDustSurfaceReaction) continue;

            // Retrieve the indices of reactants and products
            const std::size_t idx_r1 = reaction->reactant_ids_[0];
            const std::size_t idx_r2 = reaction->reactant_ids_[1];
            const std::size_t idx_p1 = reaction->product_ids_[0];
            const std::size_t idx_p2 = reaction->product_ids_[1];
            const std::size_t idx_p3 = reaction->product_ids_[2];

            // Retrieve the enthalpy of formation for reactants and products
            const Real enthalpy_r1 = ptr_species_manager_->GetSpeciesEnthalpyOfFormation(idx_r1);
            const Real enthalpy_r2 = ptr_species_manager_->GetSpeciesEnthalpyOfFormation(idx_r2);
            const Real enthalpy_p1 = ptr_species_manager_->GetSpeciesEnthalpyOfFormation(idx_p1);
            const Real enthalpy_p2 = (idx_p2 != kNotFoundSpecies ? ptr_species_manager_->GetSpeciesEnthalpyOfFormation(idx_p2) : 0.0);
            const Real enthalpy_p3 = (idx_p3 != kNotFoundSpecies ? ptr_species_manager_->GetSpeciesEnthalpyOfFormation(idx_p3) : 0.0);

            // Calculate the reaction enthalpy (kJ/mol)
            const Real reaction_enthalpy = (enthalpy_p1 + enthalpy_p2 + enthalpy_p3) - (enthalpy_r1 + enthalpy_r2);

            // Calculate exothermicity (convert from kJ to K using Boltzmann constant and Avogadro's number)
            Real exothermicity = -reaction_enthalpy; // Make exothermicity positive for exothermic reactions
            exothermicity = exothermicity * 1.0e3 / (constants::kBoltzmannConstant * 1.0e-7);
            exothermicity = exothermicity / constants::kAvogadroConstant;

            // If enthalpy is not valid, set exothermicity to -1 (invalid)
            if (enthalpy_r1 <= -9999.0 || enthalpy_r2 <= -9999.0 || enthalpy_p1 <= -9999.0) exothermicity = -1.0;
            if (idx_p2 != kNotFoundSpecies && enthalpy_p2 <= -9999.0) exothermicity = -1.0;
            if (idx_p3 != kNotFoundSpecies && enthalpy_p3 <= -9999.0) exothermicity = -1.0;

            // Initialize desorption probabilities for bare silicate and H2O ice surfaces
            Real probability_on_bare_silicate = 1.0;
            Real probability_on_H2O_ice = 1.0;

            if (exothermicity >= 0.0) { // Exothermic reactions

                // Calculate desorption rates for products
                int number_of_atoms_ = 0;
                for (std::size_t idx_p : reaction->product_ids_) {
                    if (idx_p == kNotFoundSpecies) continue;
                    const std::vector<std::size_t>& element_composition = ptr_species_manager_->GetSpeciesElementComposition(idx_p);
                    for (std::size_t nelm : element_composition) {
                        number_of_atoms_ += nelm;
                    }
                }
                number_of_atoms_ *= 3;

                for (std::size_t idx_p : reaction->product_ids_) {
                    if (idx_p == kNotFoundSpecies) continue;

                    // Calculate the desorption probability for bare silicate and H2O ice
                    Real eps_m = SQR(120.0 * constants::kProtonMass - ptr_species_manager_->GetSpeciesMass(idx_p)) 
                        / SQR(120.0 * constants::kProtonMass + ptr_species_manager_->GetSpeciesMass(idx_p));

                    Real pcb_bare = std::exp(-ptr_species_manager_->GetSpeciesBindingEnergyOnH2Oice(idx_p)
                        / (eps_m * exothermicity / static_cast<Real>(number_of_atoms_)));

                    Real pcb_ice = 0.1 * pcb_bare;

                    if (ptr_species_manager_->IsDustSurfaceSpecies(idx_p)) {
                        // For species that do not desorb chemically
                        probability_on_bare_silicate *= (1.0 - pcb_bare);
                        probability_on_H2O_ice *= (1.0 - pcb_ice);
                    } else {
                        // For species that desorb chemically
                        probability_on_bare_silicate *= pcb_bare;
                        probability_on_H2O_ice *= pcb_ice;
                    }
                }

                // Special cases for specific species (Cazaux et al. 2016)
                // Handle known species with specific desorption probabilities
                if ((ptr_species_manager_->GetSpeciesName(idx_r1) == "sO" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sH") ||
                    (ptr_species_manager_->GetSpeciesName(idx_r1) == "sH" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sO")) {
                    Real pcb_ice = 0.25;
                    probability_on_H2O_ice = (ptr_species_manager_->IsDustSurfaceSpecies(idx_p1)) ? 1.0 - pcb_ice : pcb_ice;
                }

                if ((ptr_species_manager_->GetSpeciesName(idx_r1) == "sOH" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sH") ||
                    (ptr_species_manager_->GetSpeciesName(idx_r1) == "sH" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sOH")) {
                    Real pcb_ice = 0.30;
                    probability_on_H2O_ice = (ptr_species_manager_->IsDustSurfaceSpecies(idx_p1)) ? 1.0 - pcb_ice : pcb_ice;
                }

                if (ptr_species_manager_->GetSpeciesName(idx_r1) == "sN" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sN") {
                    Real pcb_ice = 0.50;
                    probability_on_H2O_ice = (ptr_species_manager_->IsDustSurfaceSpecies(idx_p1)) ? 1.0 - pcb_ice : pcb_ice;
                }

                if (ptr_species_manager_->GetSpeciesName(idx_r1) == "sH" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sH") {
                    Real pcb_ice = 1.0, pcb_bare = 1.0;
                    probability_on_H2O_ice = (ptr_species_manager_->IsDustSurfaceSpecies(idx_p1)) ? 1.0 - pcb_ice : pcb_ice;
                    probability_on_bare_silicate = (ptr_species_manager_->IsDustSurfaceSpecies(idx_p1)) ? 1.0 - pcb_bare : pcb_bare;
                }

                // Special case for sH2 + sH2 -> sH2 + H2
                if (ptr_species_manager_->GetSpeciesName(idx_r1) == "sH2" && ptr_species_manager_->GetSpeciesName(idx_r2) == "sH2") {
                    reaction->branching_ratio_ = 1.0;
                    probability_on_H2O_ice = 1.0;
                    probability_on_bare_silicate = 1.0;
                }

            } else { // Endothermic reactions

                // No desorption occurs for endothermic reactions
                probability_on_bare_silicate = 0.0;
                probability_on_H2O_ice = 0.0;
            }

            // Update reaction rate parameters for chemical desorption probabilities
            reaction->rate_parameters_[dust_surface_reaction_params::kChemicalDesorptionOnSilicateProbability] = reaction->branching_ratio_ * probability_on_bare_silicate;
            reaction->rate_parameters_[dust_surface_reaction_params::kChemicalDesorptionOnH2OProbability] = reaction->branching_ratio_ * probability_on_H2O_ice;
        }

        return;
    }

    /**
     * @brief Set reaction type id start and end
     * 
     * This function sets the start and end indices for each reaction type in the reaction list.
     * The reaction list contains several different reaction types, and identifying the index range for each type helps with later calculations and processing.
     * For example, if different calculations need to be performed for different reaction types, this information allows the processing to be streamlined by focusing on the relevant indices for each type.
     * 
     * - The reactions are stored in `reaction_list_`, and each reaction is classified by an identifier called `type_id_`.
     * - The arrays `reaction_type_id_start_` and `reaction_type_id_end_` store the start and end indices for each reaction type.
     * 
     * The reaction list is grouped by reaction type. By identifying the index ranges for each reaction type, the processing for each type can be optimized later.
     */
    void ReactionManager::SetReactionTypeIdStartAndEnd() {
        // Set the start index for the first reaction type
        std::size_t type_id_start = reaction_list_[0]->type_id_;
        reaction_type_id_start_[type_id_start] = 0;

        // Iterate through the reaction list and set the start index for a new reaction type
        // Whenever the reaction type changes, also set the end index for the previous reaction type
        for (std::size_t ireac = 1; ireac < total_number_of_reactions_; ++ireac) {
            // When the reaction type changes from the previous one
            if (reaction_list_[ireac]->type_id_ != reaction_list_[ireac - 1]->type_id_) {
                // Set the start index for the new reaction type
                reaction_type_id_start_[reaction_list_[ireac]->type_id_] = ireac;
                // Set the end index for the previous reaction type
                reaction_type_id_end_[reaction_list_[ireac - 1]->type_id_] = ireac - 1;
            }
        }

        // Set the end index for the last reaction type
        std::size_t type_id_end = reaction_list_[total_number_of_reactions_ - 1]->type_id_;
        reaction_type_id_end_[type_id_end] = total_number_of_reactions_ - 1;
    }

    /**
    * @brief Count the number of reactions for each type
    * 
    * This function iterates through all reactions in `reaction_list_` and counts how many reactions 
    * there are for each reaction type. The results are stored in the `number_of_each_type_reactions_` array, 
    * where the index corresponds to a specific reaction type (`type_id_`), and the value at that index 
    * represents the total number of reactions of that type.
    */
    void ReactionManager::CountNumberOfEachTypeReactions() {
        // Iterate through all reactions in the reaction list
        for(std::size_t ireac = 0; ireac < total_number_of_reactions_; ++ireac) {
            // Get the type ID for the current reaction
            std::size_t itype = reaction_list_[ireac]->type_id_;
            // Increment the count for the corresponding reaction type
            number_of_each_type_reactions_[itype]++;
        }
    }
    
    /**
     * @brief Set reactions involved with species
     * 
     * This function identifies which reactions each species is involved in and stores this information in `reaction_index_list_involved_with_species_`.
     * Additionally, it stores the number of reactions each species is involved in in `number_of_reactions_involved_with_species_`.
     */
    void ReactionManager::SetReactionsInvolvedWithSpecies() {
        // Get the total number of species
        std::size_t number_of_total_species = ptr_species_manager_->total_number_of_species_;
        
        // Initialize a 2D array to hold the indices of reactions involving each species
        std::vector<std::vector<std::size_t>> use_species_for_reactions(number_of_total_species, std::vector<std::size_t>(total_number_of_reactions_));
        
        // Resize the array to store the number of reactions each species is involved in
        number_of_reactions_involved_with_species_.resize(number_of_total_species);

        // Loop through each species and count how many reactions it is involved in
        for (std::size_t ispe = 0; ispe < number_of_total_species; ++ispe) {

            // Initialize the count for the number of reactions involving the current species
            std::size_t count = 0;

            // Loop through all reactions and check if the species is involved in any of them
            for (std::size_t ireac = 0; ireac < total_number_of_reactions_; ++ireac) {
                
                const auto& reaction = reaction_list_[ireac];
                std::size_t idx_r1 = reaction->reactant_ids_[0];
                std::size_t idx_r2 = reaction->reactant_ids_[1];
                std::size_t idx_r3 = reaction->reactant_ids_[2];

                // If the species is part of the reaction (as a reactant), record the reaction index
                if (idx_r1 == ispe || idx_r2 == ispe || idx_r3 == ispe) {
                    use_species_for_reactions[ispe][count] = ireac;
                    count++;
                }
            }

            // Store the number of reactions the current species is involved in
            number_of_reactions_involved_with_species_[ispe] = count;
        }

        // Get the maximum number of reactions any species is involved in
        std::size_t max_nj = vector_utils::GetMaxValue(number_of_reactions_involved_with_species_);

        // Resize the array to hold the indices of reactions involved with each species
        reaction_index_list_involved_with_species_.resize(number_of_total_species, std::vector<std::size_t>(max_nj));

        // Store the indices of the reactions involved with each species
        for (std::size_t ispe = 0; ispe < number_of_total_species; ++ispe) {
            for (std::size_t j = 0; j < max_nj; ++j) {
                reaction_index_list_involved_with_species_[ispe][j] = use_species_for_reactions[ispe][j];
            }
        }
    }
}
