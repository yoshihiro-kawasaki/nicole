#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

#include "nicole/species/species_manager.hpp"
#include "nicole/utils/string_utils.hpp"

namespace nicole {
    SpeciesManager::SpeciesManager(
        ElementManager *ptr_element_manager, 
        InputConfig& config
    ) : ptr_element_manager_(ptr_element_manager),
        dust_species_model_parameters_(config),
        is_dust_species_(config.GetBool("is_dust")),
        is_dust_surface_reaction_(config.GetBool("is_dust_surface_reaction")),
        is_three_phase_reaction_(config.GetBool("is_three_phase_reaction"))
    {
        // Check ptr_element_manager
        if (ptr_element_manager_ == nullptr) {
            std::cerr << "Error: ptr_element_manager is nullptr" << std::endl;
            throw std::runtime_error("ElementManager pointer is null");
        }

        // Set up
        SetUpGasSpecies(config);
        if (is_dust_species_) SetUpDustRelatedSpecies(config);
        SetUpSpeciesNameList();
        SetUpSpecialSpeciesIndex();

        // Set the total number of species and the number of species in each category
        number_of_gas_species_ = gas_species_list_.size();
        number_of_dust_surface_species_ = dust_surface_species_list_.size();
        number_of_dust_mantle_species_ = dust_mantle_species_list_.size();
        number_of_dust_species_ = dust_species_list_.size();
        total_number_of_species_ = species_list_.size();
    }


    SpeciesManager::SpeciesManager(
        ElementManager *ptr_element_manager, 
        InputConfig& config, 
        const std::vector<std::string> user_gas_species_list
    ) : ptr_element_manager_(ptr_element_manager),
        dust_species_model_parameters_(config),
        is_dust_species_(config.GetBool("is_dust")),
        is_dust_surface_reaction_(config.GetBool("is_dust_surface_reaction")),
        is_three_phase_reaction_(config.GetBool("is_three_phase_reaction"))
    {
        if (ptr_element_manager_ == nullptr) {
            std::cerr << "Error: ptr_element_manager is nullptr" << std::endl;
            throw std::runtime_error("ElementManager pointer is null");
        }

        // Set up
        SetUpGasSpecies(config, user_gas_species_list);
        if (is_dust_species_) SetUpDustRelatedSpecies(config);
        SetUpSpeciesNameList();
        SetUpSpecialSpeciesIndex();

        // Set the total number of species and the number of species in each category
        number_of_gas_species_ = gas_species_list_.size();
        number_of_dust_surface_species_ = dust_surface_species_list_.size();
        number_of_dust_mantle_species_ = dust_mantle_species_list_.size();
        number_of_dust_species_ = dust_species_list_.size();
        total_number_of_species_ = species_list_.size();
    }


    void SpeciesManager::CheckSpeciesManager(const std::string& filename) const {
        std::ofstream file(filename, std::ios::out | std::ios::trunc);
        if (!file.is_open()) {
            return;
        }

        file << std::scientific;
        for (const auto& species : gas_species_list_) {
            file << std::setw(4)  << species->GetID() << " "
                 << std::setw(11) << species->GetName() << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass() << std::endl;
        }

        for (const auto& species : dust_surface_species_list_) {
            file << std::setw(4)  << species->GetID() << " "
                 << std::setw(11) << species->GetName() << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass() << std::endl;
        }

        for (const auto& species : dust_mantle_species_list_) {
            file << std::setw(4)  << species->GetID() << " "
                 << std::setw(11) << species->GetName() << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass() << std::endl;
        }

        for (const auto& species : dust_species_list_) {
            file << std::setw(4)  << species->GetID()     << " "
                 << std::setw(11) << species->GetName()   << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass()   << " "
                 << std::setw(11) << species->GetRadius() << " "
                 << std::setw(11) << species->GetCrossSection() << std::endl;
        }

        file.close();

        std::cout << "number of gas species          = " << number_of_gas_species_ << std::endl;
        std::cout << "number of dust surface species = " << number_of_dust_surface_species_ << std::endl;
        std::cout << "number of dust mantle species  = " << number_of_dust_mantle_species_ << std::endl;
        std::cout << "number of dust species         = " << number_of_dust_species_ << std::endl;
        std::cout << "total number of species        = " << total_number_of_species_ << std::endl;
    }


    std::shared_ptr<Species> SpeciesManager::GetSpecies(const SpeciesID id) const {
        if (id >= species_list_.size()) return nullptr;
        return species_list_[id];
    }


    std::string SpeciesManager::GetSpeciesName(const SpeciesID id) const {
        if (id >= species_list_.size()) return kEmptyString;
        return species_list_[id]->GetName();
    }


    int SpeciesManager::GetSpeciesCharge(const SpeciesID id) const {
        if (id >= species_list_.size()) return 0;
        return species_list_[id]->GetCharge();
    }

    
    Real SpeciesManager::GetSpeciesMass(const SpeciesID id) const {
        if (id >= species_list_.size()) return 0.0;
        return species_list_[id]->GetMass();
    }


    const std::vector<std::size_t>& SpeciesManager::GetSpeciesElementComposition(const SpeciesID id) const {
        // Default element composition (for dust species)
        static const std::vector<std::size_t> default_composition(ptr_element_manager_->GetNumberOfElements());

        // Check if the index is valid and within the bounds of the species list
        if (id >= species_list_.size()) return default_composition;
        
        // For dust species
        if (species_list_[id]->GetSpeciesType() == SpeciesType::Dust) return default_composition;

        // For gas species and dust-surface species, dust-mantle species
        return species_list_[id]->GetElementComposition();
    }


    Real SpeciesManager::GetSpeciesBindingEnergyOnH2Oice(const SpeciesID id) const {
        // If the index is out of bounds, return binding energy=0.0.
        if (id >= species_list_.size()) return 0.0;

        // Determine the species type and return the appropriate binding energy
        switch (species_list_[id]->GetSpeciesType()) {
            case SpeciesType::Gas:
            case SpeciesType::Surface:
            case SpeciesType::Mantle:
                return species_list_[id]->GetBindingEnergyOnH2Oice();
            default:
                return 0.0; // Dust species type does not have binding energy on H2O ice
        }
    }

    Real SpeciesManager::GetSpeciesBindingEnergyOnSilicate(const SpeciesID id) const {
        // If the index is out of bounds, return binding energy=0.0.
        if (id >= species_list_.size()) return 0.0;

        // Determine the species type and return the appropriate binding energy
        switch (species_list_[id]->GetSpeciesType()) {
            case SpeciesType::Gas:
            case SpeciesType::Surface:
            case SpeciesType::Mantle:
                return species_list_[id]->GetBindingEnergyOnBareSilicate();
            default:
                return 0.0; // Dust species type does not have binding energy on bare silicate
        }
    }


    Real SpeciesManager::GetSpeciesEnthalpyOfFormation(const SpeciesID id) const {
        // Check if the index is out of bounds, return enthalpy=0.0
        if (id >= species_list_.size()) return 0.0;

        // Determine the species type and return the appropriate enthalpy of formation
        switch (species_list_[id]->GetSpeciesType()) {
            case SpeciesType::Gas:
            case SpeciesType::Surface:
            case SpeciesType::Mantle:
                return species_list_[id]->GetEnthalpyOfFormation();
            default:
                return 0.0; // Dust species type does not have an enthalpy of formation value
        }
    }


    std::shared_ptr<Species> SpeciesManager::FindSpeciesByName(const std::string& name) const {
        auto it = std::find_if(
            species_list_.begin(), 
            species_list_.end(),
            [&name](const auto& species) { return species->GetName() == name; }
        );
        return (it != species_list_.end()) ? *it : nullptr;
    }


    std::shared_ptr<GasSpecies> SpeciesManager::FindGasSpeciesByName(const std::string& species_name) {
        auto it = std::find_if(
            gas_species_list_.begin(), 
            gas_species_list_.end(),
            [&species_name](const auto& gas_species) { return gas_species->GetName() == species_name; }
        );
        return (it != gas_species_list_.end()) ? *it : nullptr;
    }


    std::shared_ptr<DustSurfaceSpecies> SpeciesManager::FindDustSurfaceSpeciesByName(const std::string& species_name) {
        auto it = std::find_if(
            dust_surface_species_list_.begin(),
            dust_surface_species_list_.end(),
            [&species_name](const auto& dust_surface_species) { return dust_surface_species->GetName() == species_name; }
        );
        return (it != dust_surface_species_list_.end() ? *it : nullptr);
    }


    std::shared_ptr<DustMantleSpecies> SpeciesManager::FindDustMantleSpeciesByName(const std::string& species_name) {
        auto it = std::find_if(
            dust_mantle_species_list_.begin(),
            dust_mantle_species_list_.end(),
            [&species_name](const auto& dust_mantle_species) { return dust_mantle_species->GetName() == species_name; }
        );
        return (it != dust_mantle_species_list_.end() ? *it : nullptr);
    }


    std::shared_ptr<DustSpecies> SpeciesManager::FindDustSpeciesByBinNumberAndCharge(
        const std::size_t bin_number, 
        const int charge
    ) const {
        auto it = std::find_if(
            dust_species_list_.begin(),
            dust_species_list_.end(),
            [bin_number, charge](const std::shared_ptr<DustSpecies>& dust_species) {
                return (dust_species->GetBinNumber() == bin_number && dust_species->GetCharge() == charge);
            }
        );
        return (it != dust_species_list_.end()) ? *it : nullptr;
    }

    
    SpeciesID SpeciesManager::FindSpeciesID(const std::string& species_name) const {
        auto it = std::find(species_name_list_.begin(), species_name_list_.end(), species_name);
        return (it != species_name_list_.end()) ? std::distance(species_name_list_.begin(), it) : kNotFoundSpecies;
    }

    
    bool SpeciesManager::IsChemicalSpeciesByName(const std::string& species_name) const {
        return std::any_of(
            species_list_.begin(), 
            species_list_.end(),
            [&species_name](const auto& species) { return species->GetName() == species_name; }
        );
    }

    
    bool SpeciesManager::IsDustSurfaceSpecies(const SpeciesID id) const {
        if (id == kNotFoundSpecies || id >= species_list_.size()) return false;
        return (species_list_[id]->GetSpeciesType() == SpeciesType::Surface);
    }


    bool SpeciesManager::IsDustMantleSpecies(const SpeciesID id) const {
        if (id == kNotFoundSpecies || id >= species_list_.size()) return false;
        return (species_list_[id]->GetSpeciesType() == SpeciesType::Mantle);
    }


    bool SpeciesManager::IsDustSurfaceOrDustMantleSpecies(const SpeciesID id) const {
        if (id == kNotFoundSpecies || id >= species_list_.size()) return false;
        return ((species_list_[id]->GetSpeciesType() == SpeciesType::Surface) ||
                (species_list_[id]->GetSpeciesType() == SpeciesType::Mantle));
    }


    void SpeciesManager::ReadGasSpeciesFile(const std::string& file_path, std::vector<std::shared_ptr<GasSpecies>>& gas_species_list) {
        std::ifstream file(file_path, std::ios::in);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << file_path << std::endl;
            return;
        }

        std::size_t number_of_elements = ptr_element_manager_->GetNumberOfElements();
        std::string line;
        std::size_t dummy_index = 0;
        int line_number = 0;

        while (std::getline(file, line)) {
            line_number++;

            // Skip empty lines and comments
            if (line.empty() || line[0] == '#' || line[0] == '!') continue; 

            // Split the line into tokens
            std::vector<std::string> split_result = string_utils::Split(line, ' ', true);

            // Validate format: Expecting (name, charge, element composition)
            if (split_result.size() != number_of_elements + 2) {
                std::cerr << "Warning: Incorrect format in line " << line_number << ": " << line << std::endl;
                std::cout << split_result.size() << " " << number_of_elements << std::endl;
                for (const auto item : split_result) {
                    std::cout << "'" << item << "'" << " ";
                }
                std::cout << std::endl;
                continue;
            }

            // Extract species name
            const std::string& gas_species_name = split_result[0];

            // Skip species that start with "GRAIN"
            if (gas_species_name.rfind("GRAIN", 0) == 0) continue;

            // Extract species charge
            int gas_species_charge;
            try {
                gas_species_charge = std::stoi(split_result[1]);
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid charge format in line " << line_number << ": " << line << std::endl;
                continue;
            }

            // Extract element composition
            std::vector<std::size_t> gas_species_element_composition(number_of_elements);
            try {
                for (std::size_t i = 0; i < number_of_elements; ++i) {
                    gas_species_element_composition[i] = std::stoi(split_result[i + 2]);
                }
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid element composition in line " << line_number << ": " << line << std::endl;
                continue;
            }

            // Create GasSpecies object and add to list
            auto gas_species = std::make_shared<GasSpecies>(
                dummy_index,
                gas_species_name,
                gas_species_charge,
                gas_species_element_composition,
                ptr_element_manager_
            );

            if (!gas_species) {
                std::cerr << "Error: Failed to create GasSpecies: " << gas_species_name << std::endl;
                continue;
            }

            gas_species_list.emplace_back(gas_species);
        }

        file.close();
    }


    void SpeciesManager::ReadBindingEnergyFile(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << file_path << std::endl;
            return;
        }

        std::size_t index = species_list_.size();
        std::string line;
        int line_number = 0;

        while (std::getline(file, line)) {
            line_number++;

            // Skip empty lines and comments
            if (line.empty() || line[0] == '#') continue;

            // Split the line into tokens
            std::vector<std::string> split_result = string_utils::Split(line, ' ', true);

            // Ensure the expected number of tokens exist
            if (split_result.size() < 4) {
                std::cerr << "Warning: Incorrect format in line " << line_number << ": " << line << std::endl;
                continue;
            }

            // Extract species name
            const std::string species_name = split_result[0];

            // Generate dust surface species name
            const std::string dust_surface_species_name = kDustSurfaceSpeciesPrefix + species_name;

            // Extract binding energies
            Real binding_energy_on_H2O_ice = 0.0;
            Real binding_energy_on_bare_silicate = 0.0;
            try {
                binding_energy_on_H2O_ice = std::stod(split_result[2]);
                binding_energy_on_bare_silicate = std::stod(split_result[3]);
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid binding energy format in line " << line_number
                        << " : " << line << " (" << e.what() << ")" << std::endl;
                continue;
            }

            // Find corresponding gas-phase species
            const auto corresponding_gas_species = FindGasSpeciesByName(species_name);
            if (!corresponding_gas_species) {
                // std::cerr << "Warning: No corresponding gas species for " << species_name << " in line " << line_number << std::endl;
                continue;
            }

            // Set binding energies for the gas species
            corresponding_gas_species->SetBindingEnergyOnH2Oice(binding_energy_on_H2O_ice);
            corresponding_gas_species->SetBindingEnergyOnBareSilicate(binding_energy_on_bare_silicate);

            // Create dust surface species
            auto dust_surface_species = std::make_shared<DustSurfaceSpecies>(
                index++,
                dust_surface_species_name,
                corresponding_gas_species,
                binding_energy_on_H2O_ice,
                binding_energy_on_bare_silicate
            );

            // Store the created dust surface species
            species_list_.emplace_back(dust_surface_species);
            dust_surface_species_list_.emplace_back(dust_surface_species);

            // Link the gas species to its corresponding surface species
            corresponding_gas_species->SetCorrespondingSurfaceSpecies(dust_surface_species);

        }

        file.close();
    }


    void SpeciesManager::ReadEnthalpyOfFormationFile(const std::string& file_path) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << file_path << std::endl;
            return;
        }

        std::string line;
        int line_number = 0;
        
        while (std::getline(file, line)) {
            line_number++;

            // Skip empty lines and comments
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;

            // Split the line into tokens
            std::vector<std::string> split_result = string_utils::Split(line, ' ', true);

            // Validate format
            if (split_result.size() < 2) {
                std::cerr << "Warning: Invalid format at line " << line_number << " in " << file_path << std::endl;
                continue;
            }

            const std::string& species_name = split_result[0];
            Real species_enthalpy = 0.0;

            try {
                species_enthalpy = std::stod(split_result[1]); // Convert string to Real
            } catch (const std::invalid_argument& e) {
                std::cerr << "Warning: Invalid enthalpy value at line " << line_number << " in " << file_path << std::endl;
                continue;
            }

            // Retrieve the corresponding gas species
            const auto& gas_species = FindGasSpeciesByName(species_name);
            if (!gas_species) {
                // std::cerr << "Warning: No corresponding gas species for '" << species_name << "' at line " << line_number << std::endl;
                continue;
            }

            // Set enthalpy for gas species
            gas_species->SetEnthalpyOfFormation(species_enthalpy);

            // Retrieve and set enthalpy for corresponding dust surface species
            if (const auto& dust_surface_species = gas_species->GetCorrespondingSurfaceSpecies()) {
                dust_surface_species->SetEnthalpyOfFormation(species_enthalpy);
            } else {
                std::cerr << "Warning: No corresponding dust surface species for '" << species_name << "' at line " << line_number << std::endl;
            }

            // If three-phase reactions are enabled, set enthalpy for dust mantle species
            if (is_three_phase_reaction_) {
                if (const auto& dust_mantle_species = gas_species->GetCorrespondingMantleSpecies()) {
                    dust_mantle_species->SetEnthalpyOfFormation(species_enthalpy);
                } else {
                    std::cerr << "Warning: No corresponding dust mantle species for '" << species_name << "' at line " << line_number << std::endl;
                }
            }
        }

        file.close();
    }


    void SpeciesManager::SetUpGasSpecies(InputConfig &config) {
        const std::string file_path = config.GetString("species_file");

        // Temporary list to store gas species read from the file
        std::vector<std::shared_ptr<GasSpecies>> local_gas_species_list;

        // Read gas species from file into local_gas_species_list
        ReadGasSpeciesFile(file_path, local_gas_species_list);

        // Determine the starting id for the gas species
        SpeciesID id = 0;
        if (!species_list_.empty()) id = species_list_.size();

        // Add each gas species to the main species list and gas species list
        for (const auto& gas_species : local_gas_species_list) {
            gas_species->SetID(id++);
            species_list_.emplace_back(gas_species);
            gas_species_list_.emplace_back(gas_species);
        }
    }


    void SpeciesManager::SetUpGasSpecies(InputConfig& config, const std::vector<std::string>& user_gas_species_list) {
        const std::string file_path = config.GetString("species_file");

        // Temporary list to store gas species read from the file
        std::vector<std::shared_ptr<GasSpecies>> local_gas_species_list;

         // Read gas species from file into local_gas_species_list
        ReadGasSpeciesFile(file_path, local_gas_species_list);

        // Add gas species that match the user-specified list
        SpeciesID id = 0;
        if (!species_list_.empty()) id = species_list_.size();

        for (const auto& gas_species : local_gas_species_list) {
            const std::string& gas_species_name = gas_species->GetName();
            if (std::find(user_gas_species_list.begin(), user_gas_species_list.end(), gas_species_name) != user_gas_species_list.end()) {
                gas_species->SetID(id++);
                species_list_.emplace_back(gas_species);
                gas_species_list_.emplace_back(gas_species);
            }
        }

        // Identify species requested by the user but missing in the database
        std::vector<std::string> no_gas_species_list;
        for (const auto& gas_species_name : user_gas_species_list) {
            auto it = std::find_if(
                local_gas_species_list.begin(),
                local_gas_species_list.end(),
                [&gas_species_name](const std::shared_ptr<GasSpecies>& gas_species) {
                    return (gas_species->GetName() == gas_species_name);
                }
            );

            if (it == local_gas_species_list.end()) {
                no_gas_species_list.emplace_back(gas_species_name);
            }
        }

        // Report missing species
        if (!no_gas_species_list.empty()) {
            std::cout << "No species found in the database:" << std::endl;
            for (const auto& species_name : no_gas_species_list) {
                std::cout << std::setw(12) << species_name << std::endl;
            }
        }
    }

    
    void SpeciesManager::GenerateDustSpecies() {
        // Retrieve dust model parameters
        const std::size_t number_of_dust_bins = dust_species_model_parameters_.GetNumberOfDustBins();
        const int max_dust_charge_number = dust_species_model_parameters_.GetMaxDustChargeNumber();
        const int min_dust_charge_number = dust_species_model_parameters_.GetMinDustChargeNumber();
        
        SpeciesID id = species_list_.size(); // Start index from the current species list size

        // Loop through each dust bin
        for (std::size_t bin_number = 1; bin_number <= number_of_dust_bins; ++bin_number) {
            // Loop through each possible charge state
            for (int dust_charge_number = min_dust_charge_number; dust_charge_number <= max_dust_charge_number; ++dust_charge_number) {

                // Generate dust species name based on bin number and charge
                std::string dust_species_name = kDustSpeciesPrefix + "s" + std::to_string(bin_number) + 
                                                "c" + ((dust_charge_number >= 0) ? "+" : "") 
                                                + std::to_string(dust_charge_number);

                // Create dust species
                const std::shared_ptr<DustSpecies> dust_species = std::make_shared<DustSpecies>(
                    id++,
                    dust_species_name,
                    dust_charge_number,
                    bin_number,
                    &dust_species_model_parameters_
                );

                // Store dust species in the lists
                species_list_.emplace_back(dust_species);
                dust_species_list_.emplace_back(dust_species);
            }
        }
    }


    void SpeciesManager::SetUpDustMantleSpecies() {
        if (dust_surface_species_list_.empty() || !is_three_phase_reaction_) return;

        SpeciesID id = species_list_.size();

        for (const auto& dust_surface_species : dust_surface_species_list_) {
            // Retrieve the corresponding gas species
            const auto& gas_species = dust_surface_species->GetCorrespondingGasSpecies();
            if (!gas_species) {
                std::cerr << "Warning: Dust surface species has no corresponding gas species." << std::endl;
                continue;
            }

            // Generate the dust mantle species name
            const std::string dust_mantle_species_name = kDustMantleSpeciesPrefix + gas_species->GetName();

            // Create the dust mantle species
            const auto dust_mantle_species = std::make_shared<DustMantleSpecies>(
                id++,
                dust_mantle_species_name,
                gas_species,
                dust_surface_species
            );

            // Add the new mantle species to the lists
            species_list_.emplace_back(dust_mantle_species);
            dust_mantle_species_list_.emplace_back(dust_mantle_species);

            // Set corresponding relationships
            gas_species->SetCorrespondingMantleSpecies(dust_mantle_species);
            dust_surface_species->SetCorrespondingDustMantleSpecies(dust_mantle_species);
        }
    }


    void SpeciesManager::UpdateMantleSpeciesDiffusionBarriers() {
        if (dust_mantle_species_list_.empty() || !is_three_phase_reaction_) return;

        // List of dust mantle species to exclude from diffusion barrier update
        const std::vector<std::string> excluded_mantle_species = {"mH", "mH2", "mC", "mN", "mO"};

        // Get the diffusion barrier value for "mH2O"
        Real diffusion_barrier_H2O = 0.0;
        for (const auto& dust_mantle_species : dust_mantle_species_list_) {
            if (dust_mantle_species->GetName() == "mH2O") {
                diffusion_barrier_H2O = dust_mantle_species->GetDiffusionBarrierOnH2Oice();
                break;
            }
        }

        // Update the diffusion barrier for dust mantle species
        for (const auto& dust_mantle_species : dust_mantle_species_list_) {
            if (dust_mantle_species->GetBindingEnergyOnH2Oice() != 0.0 &&
                dust_mantle_species->GetDiffusionBarrierOnH2Oice() < diffusion_barrier_H2O) {
                if (!string_utils::IsInStringVector(excluded_mantle_species, dust_mantle_species->GetName())) {
                    dust_mantle_species->SetDiffusionBarrierOnH2Oice(diffusion_barrier_H2O);
                }
            }
        }
    }


    void SpeciesManager::SetUpDustRelatedSpecies(InputConfig& input) {
        GenerateDustSpecies();

        if (is_dust_surface_reaction_) {
            // If dust surface reactions are enabled, read the binding energy file to set up dust surface species
            ReadBindingEnergyFile(input.GetString("binding_energy_file"));

            // If three-phase reactions (gas-phase, dust surface, and dust mantle reactions) are enabled,
            // set up dust mantle species and update diffusion barriers for the mantle species.
            if (is_three_phase_reaction_) {
                SetUpDustMantleSpecies();
                UpdateMantleSpeciesDiffusionBarriers();
            }

            // Read the enthalpy of formation file. This file provides the enthalpy data used for calculating
            // chemical desorption reactions on the dust surface.
            ReadEnthalpyOfFormationFile(input.GetString("enthalpy_file"));
        }
    }


    void SpeciesManager::SetUpSpecialSpeciesIndex() {
        // Set the total number of species from the species list
        total_number_of_species_ = species_list_.size();

        // Initialize the indices for the special species to kNotFoundSpecies to signify that they are not yet assigned
        index_electron_ = index_H_ = index_H2_ = index_He_ = index_CO_ = kNotFoundSpecies;
        index_sH_ = index_sH2_ = index_mH_ = index_mH2_ = index_sH2O_ = kNotFoundSpecies;

        // Create a map that associates species names with their corresponding index variable
        std::unordered_map<std::string, std::size_t*> special_species_map = {
            {kElectron, &index_electron_},
            {"H", &index_H_},
            {"H2", &index_H2_},
            {"He", &index_He_},
            {"CO", &index_CO_},
            {"sH", &index_sH_},
            {"sH2", &index_sH2_},
            {"mH", &index_mH_},
            {"mH2", &index_mH2_},
            {"sH2O", &index_sH2O_}
        };

        // Iterate through the species list and check for special species names
        for (std::size_t ispe = 0; ispe < total_number_of_species_; ++ispe) {
            const std::string& species_name = species_list_[ispe]->GetName();
            auto it = special_species_map.find(species_name);
            if (it != special_species_map.end()) {
                *(it->second) = ispe;  // Set the index for the special species
            }
        }
    }


    void SpeciesManager::SetUpSpeciesNameList() {
        if (species_list_.empty()) return;
        for (const auto& species : species_list_) {
            species_name_list_.emplace_back(species->GetName());
            if (species->GetSpeciesType() == SpeciesType::Gas) {
                gas_species_name_list_.emplace_back(species->GetName());
            }
        }
    }
}
