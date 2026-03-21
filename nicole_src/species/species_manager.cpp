/**
 * @file species_manager.cpp
 * @date 2025-02-12
 * @author kawasaki
 * 
 * This file contains the implementation of the SpeciesManager class, which manages various species, including gas species, dust surface species,
 * dust mantle species, and dust species. The class provides functions for setting up species from input configurations, managing their properties,
 * and interacting with other components like the ElementManager.
 */

#include "species_manager.hpp"

namespace nicole
{
    /**
     * @brief Constructor for SpeciesManager.
     * 
     * Initializes a SpeciesManager object, setting up species based on the provided element manager and configuration.
     * 
     * @param ptr_element_manager Pointer to the ElementManager for handling elements.
     * @param input Configuration input.
     */
    SpeciesManager::SpeciesManager(
        ElementManager *ptr_element_manager, 
        InputConfig& input
    ) : ptr_element_manager_(ptr_element_manager),
        dust_species_model_parameters_(input),
        is_dust_species_(input.GetBool("is_dust")),
        is_dust_surface_reaction_(input.GetBool("is_dust_surface_reaction")),
        is_three_phase_reaction_(input.GetBool("is_three_phase_reaction"))
    {
        // Check if ptr_element_manager_ is nullptr
        // If the pointer is null, output an error message and throw an exception
        if (ptr_element_manager_ == nullptr) {
            std::cerr << "Error: ptr_element_manager is nullptr" << std::endl;
            throw std::runtime_error("ElementManager pointer is null");
        }

        // Read the gas species file specified in the input configuration
        SetUpGasSpecies(input.GetString("species_file"));

        // If dust species are enabled (as per input configuration), set up dust-related species
        if (is_dust_species_) SetUpDustRelatedSpecies(input);

        // Populate the species_name_list_ with names of all species
        // If a species is a gas species, also add it to the gas_species_name_list_
        for (const auto& species : species_list_) {
            species_name_list_.emplace_back(species->GetName());
            if (species->GetSpeciesType() == SpeciesType::Gas) {
                gas_species_name_list_.emplace_back(species->GetName());
            }
        }

        // Set up special species indices (like electron, H, H2, etc.)
        SetSpecialSpeciesIndex();

        // Set the total number of species and the number of species in each category
        number_of_gas_species_ = gas_species_list_.size();
        number_of_dust_surface_species_ = dust_surface_species_list_.size();
        number_of_dust_mantle_species_ = dust_mantle_species_list_.size();
        number_of_dust_species_ = dust_species_list_.size();
        total_number_of_species_ = species_list_.size();
    }

    /**
     * @brief Constructor for SpeciesManager with a user-provided gas species list.
     * 
     * This version of the constructor allows the user to provide a list of gas species names.
     * 
     * @param ptr_element_manager Pointer to the ElementManager for handling elements.
     * @param input Configuration input.
     * @param user_gas_species_list List of gas species names.
     */
    SpeciesManager::SpeciesManager(
        ElementManager *ptr_element_manager, 
        InputConfig& input, 
        const std::vector<std::string> user_gas_species_list
    ) : ptr_element_manager_(ptr_element_manager),
        dust_species_model_parameters_(input),
        is_dust_species_(input.GetBool("is_dust")),
        is_dust_surface_reaction_(input.GetBool("is_dust_surface_reaction")),
        is_three_phase_reaction_(input.GetBool("is_three_phase_reaction"))
    {
        // Read the gas species file, using both the species file and the user-provided gas species list
        SetUpGasSpecies(input.GetString("species_file"), user_gas_species_list);

        // If dust species are enabled, set up dust-related species
        if (is_dust_species_) SetUpDustRelatedSpecies(input);

        // Populate the species_name_list_ with names of all species
        // If a species is a gas species, also add it to the gas_species_name_list_
        for (const auto& species : species_list_) {
            species_name_list_.emplace_back(species->GetName());
            if (species->GetSpeciesType() == SpeciesType::Gas) {
                gas_species_name_list_.emplace_back(species->GetName());
            }
        }

        // Set up special species indices (like electron, H, H2, etc.)
        SetSpecialSpeciesIndex();

        // Set the total number of species and the number of species in each category
        number_of_gas_species_ = gas_species_list_.size();
        number_of_dust_surface_species_ = dust_surface_species_list_.size();
        number_of_dust_mantle_species_ = dust_mantle_species_list_.size();
        number_of_dust_species_ = dust_species_list_.size();
        total_number_of_species_ = species_list_.size();
    }

    /**
     * @brief Set up the dust-related species from the input configuration.
     * 
     * This function configures dust species, dust surface species, and dust mantle species based on the provided input configuration.
     * If dust surface reactions are enabled, it also loads the relevant binding energy and enthalpy files for further calculations.
     * If three-phase reactions are enabled, it also sets up dust mantle species and updates their diffusion barriers.
     * 
     * @param input Configuration input that provides the necessary files and parameters.
     */
    void SpeciesManager::SetUpDustRelatedSpecies(InputConfig& input)
    {
        // Set up dust species.
        GenerateDustSpecies();

        // Check if dust surface reactions are enabled (as per input configuration)
        if (is_dust_surface_reaction_) {

            // If dust surface reactions are enabled, read the binding energy file to set up dust surface species
            ReadBindingEnergyFile(input.GetString("binding_energy_file"));

            // If three-phase reactions (gas-phase, dust surface, and dust mantle reactions) are enabled,
            // set up dust mantle species and update diffusion barriers for the mantle species.
            if (is_three_phase_reaction_) {
                // Set up dust mantle species (e.g., species that are part of the dust mantle layer)
                SetDustMantleSpecies();

                // Update the diffusion barriers for dust mantle species based on the setup.
                UpdateMantleSpeciesDiffusionBarriers();
            }

            // Read the enthalpy of formation file. This file provides the enthalpy data used for calculating
            // chemical desorption reactions on the dust surface.
            ReadEnthalpyOfFormationFile(input.GetString("enthalpy_file"));
        }
    }

    /**
     * @brief Check the setup of species in the SpeciesManager.
     * 
     * This function checks if all species have been correctly initialized and sets up a validation report.
     * 
     * @param filename Name of the file used for checking.
     */
    void SpeciesManager::CheckSpeciesManager(const std::string& filename) const 
    {
        std::ofstream file(filename, std::ios::out | std::ios::trunc);
        if (!file.is_open()) {
            return;
        }

        file << std::scientific;
        for (const auto& species : gas_species_list_) {
            file << std::setw(4)  << species->GetIndex() << " "
                 << std::setw(11) << species->GetName() << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass() << std::endl;
        }

        for (const auto& species : dust_surface_species_list_) {
            file << std::setw(4)  << species->GetIndex() << " "
                 << std::setw(11) << species->GetName() << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass() << std::endl;
        }

        for (const auto& species : dust_mantle_species_list_) {
            file << std::setw(4)  << species->GetIndex() << " "
                 << std::setw(11) << species->GetName() << " "
                 << std::setw(4)  << species->GetCharge() << " "
                 << std::setw(11) << species->GetMass() << std::endl;
        }

        for (const auto& species : dust_species_list_) {
            file << std::setw(4)  << species->GetIndex()     << " "
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

    /**
     * @brief Get a species by index.
     * 
     * This function retrieves a species from the species list using the specified index.
     * It checks if the index is valid (i.e., within the bounds of the species list) and returns the corresponding species.
     * 
     * @param index The index of the species to retrieve from the species list.
     * @return A shared pointer to the species object, or nullptr if the index is out of bounds.
     */
    std::shared_ptr<Species> SpeciesManager::GetSpecies(const std::size_t index) const 
    {
        // Check if the index is within valid bounds of the species list
        if (index < 0 || index >= species_list_.size()) return nullptr;
        
        // Return the species object at the given index
        return species_list_[index];
    }

    /**
     * @brief Get the name of a species by index.
     * 
     * This function retrieves the name of a species by its index in the species list.
     * It checks if the index is valid and returns the name of the species at the specified index.
     * 
     * @param index The index of the species in the species list.
     * @return The name of the species, or an empty string if the index is out of bounds.
     */
    std::string SpeciesManager::GetSpeciesName(const std::size_t index) const 
    {
        // Check if the index is valid and within the bounds of the species list
        if (index >= species_list_.size()) return kEmptyString;
        
        // Return the name of the species at the given index
        return species_list_[index]->GetName();
    }

    /**
     * @brief Get the charge of a species by index.
     * 
     * This function retrieves the charge of a species by its index in the species list.
     * It checks if the index is valid and returns the charge of the species at the specified index.
     * 
     * @param index The index of the species in the species list.
     * @return The charge of the species, or 0 if the index is out of bounds.
     */
    int SpeciesManager::GetSpeciesCharge(const std::size_t index) const 
    {
        // Check if the index is valid and within the bounds of the species list
        if (index < 0 || index >= species_list_.size()) return 0;
        
        // Return the charge of the species at the given index
        return species_list_[index]->GetCharge();
    }

    /**
     * @brief Get the mass of a species by index.
     * 
     * This function retrieves the mass of a species at the specified index in the species list.
     * It checks if the index is valid and returns the mass of the species at the given index.
     * 
     * @param index The index of the species in the species list.
     * @return The mass of the species, or 0.0 if the index is out of bounds.
     */
    double SpeciesManager::GetSpeciesMass(const std::size_t index) const 
    {
        // Check if the index is within valid bounds of the species list
        if (index < 0 || index >= species_list_.size()) return 0.0;
        
        // Return the mass of the species at the given index
        return species_list_[index]->GetMass();
    }

    /**
     * @brief Get the element composition of a species by index.
     * 
     * This function retrieves the element composition of a species at the specified index.
     * It checks if the index is valid and returns the element composition of the species.
     * For dust species, a default composition is returned.
     * 
     * @param index The index of the species in the species list.
     * @return A reference to the vector containing the element composition, or a default composition if the species is of type Dust.
     */
    const std::vector<std::size_t>& SpeciesManager::GetSpeciesElementComposition(const std::size_t index) const 
    {
        // Default element composition (for dust species)
        static const std::vector<std::size_t> default_composition(ptr_element_manager_->GetNumberOfElements());

        // Check if the index is valid and within the bounds of the species list
        if (index < 0 || index >= species_list_.size()) return default_composition;

        // If the species type is not Dust, return its element composition
        if (species_list_[index]->GetSpeciesType() != SpeciesType::Dust) {
            return species_list_[index]->GetElementComposition();
        }

        // For dust species, return the default element composition
        return default_composition;
    }

    /**
     * @brief Get the binding energy of a species on H2O ice.
     * 
     * This function retrieves the binding energy of a species when adsorbed on H2O ice.
     * Gas-phase species, dust surface species, and dust mantle species have binding energy values.
     * Dust species types will return a default value of 0.0.
     * 
     * @param index The index of the species in the species list.
     * @return The binding energy on H2O ice, or 0.0 if the index is invalid or the species type is unsupported.
     */
    double SpeciesManager::GetSpeciesBindingEnergyOnH2Oice(const std::size_t index) const 
    {
        // Check if the index is out of bounds
        if (index >= species_list_.size()) return 0.0;

        // Determine the species type and return the appropriate binding energy
        switch (species_list_[index]->GetSpeciesType()) {
            case SpeciesType::Gas:
            case SpeciesType::Surface:
            case SpeciesType::Mantle:
                return species_list_[index]->GetBindingEnergyOnH2Oice();
            default:
                return 0.0; // Dust species type does not have binding energy on H2O ice
        }
    }

    /**
     * @brief Get the binding energy of a species on bare silicate.
     * 
     * This function retrieves the binding energy of a species when adsorbed on bare silicate.
     * Gas-phase species, dust surface species, and dust mantle species have binding energy values.
     * Dust species types will return a default value of 0.0.
     * 
     * @param index The index of the species in the species list.
     * @return The binding energy on bare silicate, or 0.0 if the index is invalid or the species type is unsupported.
     */
    double SpeciesManager::GetSpeciesBindingEnergyOnSilicate(const std::size_t index) const 
    {
        // Check if the index is out of bounds
        if (index >= species_list_.size()) return 0.0;

        // Determine the species type and return the appropriate binding energy
        switch (species_list_[index]->GetSpeciesType()) {
            case SpeciesType::Gas:
            case SpeciesType::Surface:
            case SpeciesType::Mantle:
                return species_list_[index]->GetBindingEnergyOnBareSilicate();
            default:
                return 0.0; // Dust species type does not have binding energy on bare silicate
        }
    }

    /**
     * @brief Get the enthalpy of formation of a species.
     * 
     * This function retrieves the enthalpy of formation of a species, which represents
     * the energy change when a species is formed from its constituent elements.
     * 
     * Gas-phase species, dust surface species, and dust mantle species have defined enthalpy values.
     * Dust species types will return a default value of 0.0.
     * 
     * @param index The index of the species in the species list.
     * @return The enthalpy of formation in K, or 0.0 if the index is invalid or the species type does 
     * not have enthalpy information.
     */
    double SpeciesManager::GetSpeciesEnthalpyOfFormation(const std::size_t index) const 
    {
        // Check if the index is out of bounds
        if (index >= species_list_.size()) return 0.0;

        // Determine the species type and return the appropriate enthalpy of formation
        switch (species_list_[index]->GetSpeciesType()) {
            case SpeciesType::Gas:
            case SpeciesType::Surface:
            case SpeciesType::Mantle:
                return species_list_[index]->GetEnthalpyOfFormation();
            default:
                return 0.0; // Dust species type does not have an enthalpy of formation value
        }
    }

    /**
     * @brief Find a species by its name.
     * 
     * This function searches for a species in the species list using its name.
     * If the species is found, a shared pointer to the species object is returned.
     * If not found, it returns nullptr.
     * 
     * @param name The name of the species to search for.
     * @return A shared pointer to the species if found, otherwise nullptr.
     */
    std::shared_ptr<Species> SpeciesManager::FindSpeciesByName(const std::string& name) const 
    {
        // Use std::find_if to search for the species by name
        auto it = std::find_if(
            species_list_.begin(), 
            species_list_.end(),
            [&name](const auto& species) { return species->GetName() == name; }
        );

        // If found, return the shared pointer to the species; otherwise, return nullptr
        return (it != species_list_.end()) ? *it : nullptr;
    }

    /**
     * @brief Find a gas-phase species by its name.
     * 
     * This function searches for a gas species in the gas species list using its name.
     * If the species is found, a shared pointer to the GasSpecies object is returned.
     * If not found, it returns nullptr.
     * 
     * @param species_name The name of the gas species to search for.
     * @return A shared pointer to the GasSpecies if found, otherwise nullptr.
     */
    std::shared_ptr<GasSpecies> SpeciesManager::FindGasSpeciesByName(const std::string& species_name) 
    {
        // Use std::find_if to search for the gas species by name
        auto it = std::find_if(
            gas_species_list_.begin(), 
            gas_species_list_.end(),
            [&species_name](const auto& gas_species) { return gas_species->GetName() == species_name; }
        );

        // If found, return the shared pointer to the gas species; otherwise, return nullptr
        return (it != gas_species_list_.end()) ? *it : nullptr;
    }

    /**
     * @brief Find a dust surface species by its name.
     * 
     * This function searches for a dust surface species in the list of dust surface species.
     * If the species is found, a shared pointer to the DustSurfaceSpecies object is returned.
     * If not found, it returns nullptr.
     * 
     * @param species_name The name of the dust surface species to search for.
     * @return A shared pointer to the DustSurfaceSpecies if found, otherwise nullptr.
     */
    std::shared_ptr<DustSurfaceSpecies> SpeciesManager::FindDustSurfaceSpeciesByName(const std::string& species_name)
    {
        // Use std::find_if to search for the dust surface species by name
        auto it = std::find_if(
            dust_surface_species_list_.begin(),
            dust_surface_species_list_.end(),
            [&species_name](const auto& dust_surface_species) { return dust_surface_species->GetName() == species_name; }
        );

        // If found, return the shared pointer to the dust surface species; otherwise, return nullptr
        return (it != dust_surface_species_list_.end() ? *it : nullptr);
    }

    /**
     * @brief Find a dust mantle species by its name.
     * 
     * This function searches for a dust mantle species in the list of dust mantle species.
     * If the species is found, a shared pointer to the DustMantleSpecies object is returned.
     * If not found, it returns nullptr.
     * 
     * @param species_name The name of the dust mantle species to search for.
     * @return A shared pointer to the DustMantleSpecies if found, otherwise nullptr.
     */
    std::shared_ptr<DustMantleSpecies> SpeciesManager::FindDustMantleSpeciesByName(const std::string& species_name)
    {
        // Use std::find_if to search for the dust mantle species by name
        auto it = std::find_if(
            dust_mantle_species_list_.begin(),
            dust_mantle_species_list_.end(),
            [&species_name](const auto& dust_mantle_species) { return dust_mantle_species->GetName() == species_name; }
        );

        // If found, return the shared pointer to the dust mantle species; otherwise, return nullptr
        return (it != dust_mantle_species_list_.end() ? *it : nullptr);
    }

    /**
     * @brief Find a dust species by its bin number and charge.
     * 
     * This function searches for a dust species in the dust species list that matches 
     * both the given bin number and charge. If a matching species is found, a shared 
     * pointer to the DustSpecies object is returned. Otherwise, it returns nullptr.
     * 
     * @param bin_number The bin number of the dust species.
     * @param charge The charge of the dust species.
     * @return A shared pointer to the DustSpecies if found, otherwise nullptr.
     */
    std::shared_ptr<DustSpecies> SpeciesManager::FindDustSpeciesByBinNumberAndCharge(
        const std::size_t bin_number, 
        const int charge
    ) const 
    {
        // Use std::find_if to search for a dust species with the specified bin number and charge
        auto it = std::find_if(
            dust_species_list_.begin(),
            dust_species_list_.end(),
            [bin_number, charge](const std::shared_ptr<DustSpecies>& dust_species) {
                return (dust_species->GetBinNumber() == bin_number && dust_species->GetCharge() == charge);
            }
        );

        // If found, return the shared pointer to the dust species; otherwise, return nullptr
        return (it != dust_species_list_.end()) ? *it : nullptr;
    }

    /**
     * @brief Find the index of a species by its name.
     * 
     * This function searches for a species in the species name list and returns its index.
     * If the species is not found, it returns a predefined constant `kNotFoundSpecies`.
     * 
     * @param species_name The name of the species to search for.
     * @return The index of the species if found; otherwise, `kNotFoundSpecies`.
     */
    std::size_t SpeciesManager::FindSpeciesIndex(const std::string& species_name) const 
    {
        // Search for the species name in the list
        auto it = std::find(species_name_list_.begin(), species_name_list_.end(), species_name);
        
        // If found, return the index of the species
        if (it != species_name_list_.end()) {
            return std::distance(species_name_list_.begin(), it);
        }

        // If not found, return a predefined constant indicating failure
        return kNotFoundSpecies;
    }

    /**
     * @brief Check if a species with a given name exists in the species list.
     * 
     * This function searches for a species with the specified name in the species list.
     * It returns true if the species is found, otherwise false.
     * 
     * @param species_name The name of the species to search for.
     * @return True if the species exists, false otherwise.
     */
    bool SpeciesManager::IsChemicalSpeciesByName(const std::string& species_name) const 
    {
        // Use std::any_of to check if any species in the list has the given name
        return std::any_of(
            species_list_.begin(), 
            species_list_.end(),
            [&species_name](const auto& species) { return species->GetName() == species_name; }
        );
    }

    /**
     * @brief Check if a species at a given index is a dust surface species.
     * 
     * This function verifies whether the species at the specified index 
     * belongs to the dust surface species.
     * 
     * @param index The index of the species in the list.
     * @return True if the species is a dust surface species, false otherwise.
     */
    bool SpeciesManager::IsDustSurfaceSpecies(const std::size_t index) const 
    {
        // Check if index is out of range
        if (index == kNotFoundSpecies || index >= species_list_.size()) return false;

        // Check if the species type is "Surface"
        return (species_list_[index]->GetSpeciesType() == SpeciesType::Surface);
    }

    /**
     * @brief Check if a species at a given index is a mdust antle species.
     * 
     * This function determines whether the species at the specified index 
     * belongs to the dust surface species.
     * 
     * @param index The index of the species in the list.
     * @return True if the species is a dust mantle species, false otherwise.
     */
    bool SpeciesManager::IsDustMantleSpecies(const std::size_t index) const 
    {
        // Check if index is out of range
        if (index == kNotFoundSpecies || index >= species_list_.size()) return false;
        
        // Check if the species type is "Mantle"
        return (species_list_[index]->GetSpeciesType() == SpeciesType::Mantle);
    }

    /**
     * @brief Check if a species at a given index is either a dust surface or a dust mantle species.
     * 
     * This function verifies whether the species at the specified index belongs 
     * to either the "dust surface species" or "dust mantle species" category.
     * 
     * @param index The index of the species in the list.
     * @return True if the species is a dust surface or dust mantle species, false otherwise.
     */
    bool SpeciesManager::IsDustSurfaceOrDustMantleSpecies(const std::size_t index) const 
    {
        // Check if the index is out of range
        if (index == kNotFoundSpecies || index >= species_list_.size()) return false;
        
        // Check if the species type is either "Surface" or "Mantle"
        return ((species_list_[index]->GetSpeciesType() == SpeciesType::Surface) ||
                (species_list_[index]->GetSpeciesType() == SpeciesType::Mantle));
    }

    /**
     * @brief Reads gas species from a file and adds them to the species list.
     * 
     * This function reads gas species data from the specified file and 
     * adds them to both the species list and the gas species list. 
     * 
     * @param filename The name of the file containing gas species data.
     */
    void SpeciesManager::SetUpGasSpecies(const std::string& filename)
    {
        // Temporary list to store gas species read from the file
        std::vector<std::shared_ptr<GasSpecies>> local_gas_species_list;

        // Read gas species from file into local_gas_species_list
        ReadGasSpeciesFile(filename, local_gas_species_list);

        // Determine the starting index for the gas species
        std::size_t index = 0;
        if (!species_list_.empty()) index = species_list_.size();

        // Add each gas species to the main species list and gas species list
        for (const auto& gas_species : local_gas_species_list) {
            gas_species->SetIndex(index++);  // Assign unique index
            species_list_.emplace_back(gas_species);  // Add to the species list
            gas_species_list_.emplace_back(gas_species);  // Add to the gas species list
        }
    }

    /**
     * @brief Reads gas species from a file and adds only those specified by the user.
     * 
     * This function reads gas species data from the specified file and filters 
     * them based on a user-provided gas species list. Only gas species that match the 
     * `user_gas_species_list` are added to the species list.
     * 
     * Additionally, it checks if any species in `user_gas_species_list` are missing 
     * from the database and reports them.
     * 
     * @param filename The name of the file containing gas species data.
     * @param user_gas_species_list A list of species names specified by the user.
     */
    void SpeciesManager::SetUpGasSpecies(const std::string& filename, const std::vector<std::string>& user_gas_species_list)
    {
        // Temporary list to store gas species read from the file
        std::vector<std::shared_ptr<GasSpecies>> local_gas_species_list;
        ReadGasSpeciesFile(filename, local_gas_species_list);

        // Add gas species that match the user-specified list
        std::size_t index = 0;
        if (!species_list_.empty()) index = species_list_.size();

        for (const auto& gas_species : local_gas_species_list) {
            const std::string& gas_species_name = gas_species->GetName();

            // Check if the gas species is in the user-provided list
            if (std::find(user_gas_species_list.begin(), user_gas_species_list.end(), gas_species_name) 
                != user_gas_species_list.end()) {
                gas_species->SetIndex(index++);  // Assign unique index
                species_list_.emplace_back(gas_species);  // Add to the general species list
                gas_species_list_.emplace_back(gas_species);  // Add to the gas species list
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

    /**
     * @brief Reads gas-phase species from a file and stores them in a list.
     * 
     * This function reads gas species data from a file, extracts species name, charge, 
     * and element composition, and stores the valid species in `gas_species_list`.
     * 
     * @param filename The file containing gas species data.
     * @param gas_species_list The list where the extracted gas species will be stored.
     */
    void SpeciesManager::ReadGasSpeciesFile(const std::string& filename, std::vector<std::shared_ptr<GasSpecies>>& gas_species_list) 
    {
        // Open file
        std::ifstream file(filename, std::ios::in);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        // Check if ElementManager is initialized
        if (!ptr_element_manager_) {
            std::cerr << "Error: ptr_element_manager_ is nullptr." << std::endl;
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
                continue;
            }

            // Extract species name
            const std::string& gas_species_name = split_result[0];

            // Skip species that start with "GRAIN" (likely dust species)
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

    /**
     * @brief Reads a binding energy file and creates dust surface species.
     * 
     * This function reads binding energy data from a file, extracts the binding energies 
     * for H2O ice and bare silicate surfaces, associates them with existing gas-phase species,
     * and creates corresponding dust surface species.
     * 
     * @param filename The file containing binding energy data.
     */
    void SpeciesManager::ReadBindingEnergyFile(const std::string& filename) 
    {
        // Open file
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
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
            double binding_energy_on_H2O_ice = 0.0;
            double binding_energy_on_bare_silicate = 0.0;
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

        } // end while

        file.close();
    }

    /**
     * @brief Creates and registers dust mantle species based on existing dust surface species.
     * 
     * If three-phase reactions are enabled, this function creates a corresponding 
     * dust mantle species for each dust surface species and stores them in the species lists.
     */
    void SpeciesManager::SetDustMantleSpecies() 
    {
        // If no dust surface species exist or three-phase reactions are disabled, return early.
        if (dust_surface_species_list_.empty() || !is_three_phase_reaction_) return;

        std::size_t index = species_list_.size();

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
                index++,
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

    /**
     * @brief Reads the enthalpy of formation data from a file and assigns it to the corresponding species.
     * 
     * This function reads enthalpy values from a file and sets them for gas-phase species, 
     * and if applicable, their corresponding dust surface and dust mantle species.
     * 
     * @param filename The file containing the enthalpy of formation information.
     */
    void SpeciesManager::ReadEnthalpyOfFormationFile(const std::string& filename) 
    {
        // Open the file
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
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
                std::cerr << "Warning: Invalid format at line " << line_number << " in " << filename << std::endl;
                continue;
            }

            const std::string& species_name = split_result[0];
            double species_enthalpy = 0.0;

            try {
                species_enthalpy = std::stod(split_result[1]); // Convert string to double
            } catch (const std::invalid_argument& e) {
                std::cerr << "Warning: Invalid enthalpy value at line " << line_number << " in " << filename << std::endl;
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

    /**
     * @brief Generates dust species based on bin number and charge state.
     * 
     * This function creates dust species for each dust bin and for all charge states 
     * within the specified charge range. The generated species are stored in `species_list_` 
     * and `dust_species_list_`.
     */
    void SpeciesManager::GenerateDustSpecies() 
    {
        // Retrieve dust model parameters
        const std::size_t number_of_dust_bins = dust_species_model_parameters_.GetNumberOfDustBins();
        const int max_dust_charge_number = dust_species_model_parameters_.GetMaxDustChargeNumber();
        const int min_dust_charge_number = dust_species_model_parameters_.GetMinDustChargeNumber();
        
        std::size_t index = species_list_.size(); // Start index from the current species list size

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
                    index++,
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

    // /**
    //  * @brief set special (specific) species index
    //  */
    // void SpeciesManager::SetSpecialSpeciesIndex() 
    // {
    //     number_of_total_species_ = species_name_list_.size();

    //     for (std::size_t ispe = 0; ispe < number_of_total_species_; ++ispe) {

    //         const std::string& species_name = species_list_[ispe]->GetName();

    //         if (species_name == kElectron) {
    //             index_electron_ = ispe;
    //         } else if (species_name == "H") {
    //             index_H_ = ispe;
    //         } else if (species_name == "H2") {
    //             index_H2_ = ispe;
    //         } else if (species_name == "He") {
    //             index_He_ = ispe;
    //         } else if (species_name == "sH") {
    //             index_sH_ = ispe;
    //         } else if (species_name == "sH2") {
    //             index_sH2_ = ispe;
    //         } else if (species_name == "mH") {
    //             index_mH_ = ispe;
    //         } else if (species_name == "mH2") {
    //             index_mH2_ = ispe;
    //         } else if (species_name == "sH2O") {
    //             index_sH2O_ = ispe;
    //         }
    //     }
    // }

    /**
     * @brief Sets the index for special species in the species list.
     * 
     * This function iterates through the list of species and checks for specific species names (e.g., "H", "He", "sH2O", etc.).
     * For each matched species, the corresponding index is assigned to a pre-defined index variable (e.g., `index_H_`, `index_He_`, etc.).
     */
    void SpeciesManager::SetSpecialSpeciesIndex() 
    {
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
            const std::string& species_name = species_list_[ispe]->GetName();  // Get species name from the list

            // Search for the species name in the special species map
            auto it = special_species_map.find(species_name);

            // If found, assign the species index to the corresponding index variable
            if (it != special_species_map.end()) {
                *(it->second) = ispe;  // Set the index for the special species
            }
        }
    }

    /**
     * @brief Updates the diffusion barriers of dust mantle species based on a reference species.
     * 
     * This function updates the diffusion barriers for all dust mantle species, adjusting their values
     * to the diffusion barrier of the "mH2O" species if certain conditions are met. Specifically, if a dust 
     * mantle species has a non-zero binding energy and its diffusion barrier is smaller than that of "mH2O",
     * the diffusion barrier is updated. Some species are excluded from the update process based on a predefined list.
     */
    void SpeciesManager::UpdateMantleSpeciesDiffusionBarriers() 
    {
        // If no dust mantle species exist or the system is not a three-phase reaction, exit the function
        if (dust_mantle_species_list_.empty() || !is_three_phase_reaction_) return;

        // List of dust mantle species to exclude from diffusion barrier update
        const std::vector<std::string> excluded_mantle_species = {"mH", "mH2", "mC", "mN", "mO"};

        // Get the diffusion barrier value for "mH2O"
        double diffusion_barrier_H2O = 0.0;
        for (const auto& dust_mantle_species : dust_mantle_species_list_) {
            // Find the species with the name "mH2O" and get its diffusion barrier
            if (dust_mantle_species->GetName() == "mH2O") {
                diffusion_barrier_H2O = dust_mantle_species->GetDiffusionBarrierOnH2Oice();
                break;  // Exit the loop once the "mH2O" species is found
            }
        }

        // Update the diffusion barrier for dust mantle species
        for (const auto& dust_mantle_species : dust_mantle_species_list_) {
            // Check if the species has a non-zero binding energy and a smaller diffusion barrier than "mH2O"
            if (dust_mantle_species->GetBindingEnergyOnH2Oice() != 0.0 &&
                dust_mantle_species->GetDiffusionBarrierOnH2Oice() < diffusion_barrier_H2O) {
                // If the species is not in the exclusion list, update its diffusion barrier
                if (!string_utils::IsInStringVector(excluded_mantle_species, dust_mantle_species->GetName())) {
                    dust_mantle_species->SetDiffusionBarrierOnH2Oice(diffusion_barrier_H2O);
                }
            }
        }
    }
    
}