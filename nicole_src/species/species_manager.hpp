/**
 * @file species_manager.hpp
 * 
 * This header file defines the SpeciesManager class, which manages various species types, including gas species, dust surface species,
 * dust mantle species, and dust species. The class provides functions for adding, searching, and retrieving information about different species.
 * It also contains methods for reading configuration files, setting up species, and calculating properties related to species.
 * 
 * @date 2025-02-12
 */

#ifndef SPECIES_MANAGER_HPP
#define SPECIES_MANAGER_HPP

#include "../nicole_defs.hpp"
#include "species.hpp"
#include "gas_species.hpp"
#include "dust_surface_species.hpp"
#include "dust_mantle_species.hpp"
#include "dust_species.hpp"

namespace nicole
{
    class SpeciesManager
    {
    public:

        /**
         * @brief Constructor for SpeciesManager.
         * 
         * Initializes a SpeciesManager object, setting up species based on the provided element manager and configuration.
         * 
         * @param ptr_element_manager Pointer to the ElementManager for handling elements.
         * @param input Configuration input.
         */
        SpeciesManager(ElementManager *ptr_element_manager, InputConfig& input);

        /**
         * @brief Constructor for SpeciesManager with a user-provided gas species list.
         * 
         * This version of the constructor allows the user to provide a list of gas species names.
         * 
         * @param ptr_element_manager Pointer to the ElementManager for handling elements.
         * @param input Configuration input.
         * @param user_gas_species_list List of gas species names.
         */
        SpeciesManager(ElementManager *ptr_element_manager, InputConfig& input, const std::vector<std::string> user_gas_species_list);

        /**
         * @brief Check the setup of species in the SpeciesManager.
         * 
         * This function checks if all species have been correctly initialized and sets up a validation report.
         * 
         * @param filename Name of the file used for checking.
         */
        void CheckSpeciesManager(const std::string& filename) const;

        // Getters: These functions retrieve specific properties of the species by their index.

        /**
         * @brief Get a species by index.
         * @param index The index of the species to retrieve from the species list.
         * @return A shared pointer to the species object, or nullptr if the index is out of bounds.
         */
        std::shared_ptr<Species> GetSpecies(const std::size_t index) const;

        /**
         * @brief Get the name of a species by index.
         * @param index The index of the species in the species list.
         * @return The name of the species, or an empty string if the index is out of bounds.
         */
        std::string GetSpeciesName(const std::size_t index) const;

        /**
         * @brief Get the charge of a species by index.
         * @param index The index of the species in the species list.
         * @return The charge of the species, or 0 if the index is out of bounds.
         */
        int GetSpeciesCharge(const std::size_t index) const;

        /**
         * @brief Get the mass of a species by index.
         * @param index The index of the species in the species list.
         * @return The mass of the species, or 0.0 if the index is out of bounds.
         */
        double GetSpeciesMass(const std::size_t index) const;

        /**
         * @brief Get the element composition of a species by index.
         * @param index The index of the species in the species list.
         * @return A reference to the vector containing the element composition, or a default composition if the species is of type Dust.
         */
        const std::vector<std::size_t>& GetSpeciesElementComposition(const std::size_t index) const;

        /**
         * @brief Get the binding energy of a species on H2O ice.
         * @param index The index of the species in the species list.
         * @return The binding energy on H2O ice, or 0.0 if the index is invalid or the species type is unsupported.
         */
        double GetSpeciesBindingEnergyOnH2Oice(const std::size_t index) const;

        /**
         * @brief Get the binding energy of a species on bare silicate.
         * @param index The index of the species in the species list.
         * @return The binding energy on bare silicate, or 0.0 if the index is invalid or the species type is unsupported.
         */
        double GetSpeciesBindingEnergyOnSilicate(const std::size_t index) const;

        /**
         * @brief Get the enthalpy of formation of a species.
         * @param index The index of the species in the species list.
         * @return The enthalpy of formation in appropriate energy units (e.g., K or eV),
         *         or 0.0 if the index is invalid or the species type does not have enthalpy information.
         */
        double GetSpeciesEnthalpyOfFormation(const std::size_t index) const;

        /**
         * @brief Get the total number of species managed by the SpeciesManager.
         * @return The total number of species.
         */
        std::size_t GetTotalNumberOfSpecies() const { return total_number_of_species_; }

        // Find functions: These functions search for species by name or specific properties.

        /**
         * @brief Find a species by its name.
         * @param name The name of the species to search for.
         * @return A shared pointer to the species if found, otherwise nullptr.
         */
        std::shared_ptr<Species> FindSpeciesByName(const std::string& name) const;

        /**
         * @brief Find a gas-phase species by its name.
         * @param species_name The name of the gas species to search for.
         * @return A shared pointer to the GasSpecies if found, otherwise nullptr.
         */
        std::shared_ptr<GasSpecies> FindGasSpeciesByName(const std::string& species_name);

        /**
         * @brief Find a dust surface species by its name.
         * @param species_name The name of the dust surface species to search for.
         * @return A shared pointer to the DustSurfaceSpecies if found, otherwise nullptr.
         */
        std::shared_ptr<DustSurfaceSpecies> FindDustSurfaceSpeciesByName(const std::string& species_name);

        /**
         * @brief Find a dust mantle species by its name.
         * @param species_name The name of the dust mantle species to search for.
         * @return A shared pointer to the DustMantleSpecies if found, otherwise nullptr.
         */
        std::shared_ptr<DustMantleSpecies> FindDustMantleSpeciesByName(const std::string& species_name);

        /**
         * @brief Find a dust species by its bin number and charge.
         * @param bin_number The bin number of the dust species.
         * @param charge The charge of the dust species.
         * @return A shared pointer to the DustSpecies if found, otherwise nullptr.
         */
        std::shared_ptr<DustSpecies> FindDustSpeciesByBinNumberAndCharge(const std::size_t bin_number, const int charge) const;

        /**
         * @brief Find the index of a species by its name.
         * @param species_name The name of the species to search for.
         * @return The index of the species if found; otherwise, `kNotFoundSpecies`.
         */
        std::size_t FindSpeciesIndex(const std::string& species_name) const;

        // Boolean functions: These functions check various conditions for species.

        /**
         * @brief Check if a species with a given name exists in the species list.
         * @param species_name The name of the species to search for.
         * @return True if the species exists, false otherwise.
         */
        bool IsChemicalSpeciesByName(const std::string& species_name) const;

        /**
         * @brief Check if a species at a given index is a dust surface species.
         * @param index The index of the species in the list.
         * @return True if the species is a dust surface species, false otherwise.
         */
        bool IsDustSurfaceSpecies(const std::size_t index) const;

        /**
         * @brief Check if a species at a given index is a mdust antle species.
         * @param index The index of the species in the list.
         * @return True if the species is a dust mantle species, false otherwise.
         */
        bool IsDustMantleSpecies(const std::size_t index) const;

        /**
         * @brief Check if a species at a given index is either a dust surface or a dust mantle species.
         * @param index The index of the species in the list.
         * @return True if the species is a dust surface or dust mantle species, false otherwise.
         */
        bool IsDustSurfaceOrDustMantleSpecies(const std::size_t index) const;

        // Friend classes that have access to SpeciesManager's private members
        friend class ReactionManager;
        friend class ReactionSimulator;
        friend class NonIdealMHDeffect;

    protected:
        
        // Setup functions for species configuration

        /**
         * @brief Set up the dust-related species from the input configuration.
         * @param input Configuration input.
         */
        void SetUpDustRelatedSpecies(InputConfig& input);

        /**
         * @brief Reads gas species from a file and adds them to the species list.
         * @param filename The name of the file containing gas species data.
         */
        void SetUpGasSpecies(const std::string& filename);

        /**
         * @brief Reads gas species from a file and adds only those specified by the user.
         * @param filename The name of the file containing gas species data.
         * @param user_gas_species_list A list of species names specified by the user.
         */
        void SetUpGasSpecies(const std::string& filename, const std::vector<std::string>& user_gas_species_list);

        /**
         * @brief Reads gas-phase species from a file and stores them in a list.
         * @param filename The file containing gas species data.
         * @param gas_species_list The list where the extracted gas species will be stored.
         */
        void ReadGasSpeciesFile(const std::string& filename, std::vector<std::shared_ptr<GasSpecies> >& gas_species_list);

        /**
         * @brief Reads a binding energy file and creates dust surface species.
         * @param filename The file containing binding energy data.
         */
        void ReadBindingEnergyFile(const std::string& filename);

        /**
         * @brief Creates and registers dust mantle species based on existing dust surface species.
         */
        void SetDustMantleSpecies();

        /**
         * @brief Reads the enthalpy of formation data from a file and assigns it to the corresponding species.
         * @param filename The file containing the enthalpy of formation information.
         */
        void ReadEnthalpyOfFormationFile(const std::string& filename);

        /**
         * @brief Generates dust species based on bin number and charge state.
         */
        void GenerateDustSpecies();

        /**
         * @brief Sets the index for special species in the species list.
         */
        void SetSpecialSpeciesIndex();

        /**
         * @brief Updates the diffusion barriers of dust mantle species based on a reference species.
         */
        void UpdateMantleSpeciesDiffusionBarriers();

        // Member variables to store pointers to species and related information
        ElementManager* ptr_element_manager_;  // Pointer to ElementManager for handling elements

        // Vectors to hold species lists
        std::vector<std::shared_ptr<Species> > species_list_;  // General list of all species
        std::vector<std::shared_ptr<GasSpecies> > gas_species_list_;  // List of gas species
        std::vector<std::shared_ptr<DustSurfaceSpecies> > dust_surface_species_list_;  // List of dust surface species
        std::vector<std::shared_ptr<DustMantleSpecies> > dust_mantle_species_list_;  // List of dust mantle species
        std::vector<std::shared_ptr<DustSpecies> > dust_species_list_;  // List of dust species

        // Lists of species names
        std::vector<std::string> species_name_list_;
        std::vector<std::string> gas_species_name_list_;

        // Dust species model parameters
        DustSpeciesModelParameters dust_species_model_parameters_;

        // Number of total species and specific categories of species
        std::size_t total_number_of_species_;
        std::size_t number_of_gas_species_;
        std::size_t number_of_dust_surface_species_;
        std::size_t number_of_dust_mantle_species_;
        std::size_t number_of_dust_species_;

        // Indices for special species
        std::size_t index_electron_;
        std::size_t index_H_;
        std::size_t index_H2_;
        std::size_t index_He_;
        std::size_t index_CO_;
        std::size_t index_sH_;
        std::size_t index_sH2_;
        std::size_t index_mH_;
        std::size_t index_mH2_;
        std::size_t index_sH2O_;

        // Flags for different types of species and reactions
        bool is_dust_species_;
        bool is_dust_surface_reaction_;
        bool is_three_phase_reaction_;
    };
}

#endif /* SPECIES_MANAGER_HPP */



// /**
//  * @file species_manager.hpp
//  */

// #ifndef SPECIES_MANAGER_HPP
// #define SPECIES_MANAGER_HPP

// #include "../nicole_defs.hpp"
// #include "species.hpp"
// #include "gas_species.hpp"
// #include "dust_surface_species.hpp"
// #include "dust_mantle_species.hpp"
// #include "dust_species.hpp"

// namespace nicole
// {
//     class SpeciesManager
//     {
//     public:

//         SpeciesManager(ElementManager *ptr_element_manager, InputConfig& input);
//         SpeciesManager(ElementManager *ptr_element_manager, InputConfig& input, const std::vector<std::string> user_gas_species_list);

//         void CheckSpeciesManager(const std::string& filename) const;

//         // Get function
//         std::shared_ptr<Species> GetSpecies(const std::size_t inex) const;
//         std::string GetSpeciesName(const std::size_t index) const;
//         int GetSpeciesCharge(const std::size_t index) const;
//         double GetSpeciesMass(const std::size_t index) const;
//         const std::vector<std::size_t>& GetSpeciesElementComposition(const std::size_t index) const;
//         double GetSpeciesBindingEnergyOnH2Oice(const std::size_t index) const;
//         double GetSpeciesBindingEnergyOnSilicate(const std::size_t index) const;
//         double GetSpeciesEnthalpyOfFormation(const std::size_t index) const;
//         std::size_t GetNumberOfTotalSpecies() const { return number_of_total_species_; }

//         // Find function
//         std::shared_ptr<Species> FindSpeciesByName(const std::string& name) const;
//         std::shared_ptr<GasSpecies> FindGasSpeciesByName(const std::string& species_name);
//         std::shared_ptr<DustSurfaceSpecies> FindDustSurfaceSpeciesByName(const std::string& species_name);
//         std::shared_ptr<DustMantleSpecies> FindDustMantleSpeciesByName(const std::string& species_name);
//         std::shared_ptr<DustSpecies> FindDustSpeciesByBinNumberAndCharge(const std::size_t bin_number, const int charge) const;
//         std::size_t FindSpeciesIndex(const std::string& species_name) const;

//         // Bool functions
//         bool IsChemicalSpeciesByName(const std::string& species_name) const;
//         bool IsSurfaceSpecies(const std::size_t index) const;
//         bool IsMantleSpecies(const std::size_t index) const;
//         bool IsSurfaceOrMantleSpecies(const std::size_t index) const;

//         friend class ReactionManager;
//         friend class ReactionSimulator;
//         friend class NonIdealMHDeffect;

//     protected:
        
//         void SetUpDustRelatedSpecies(InputConfig& input);

//         void SetUpGasSpecies(const std::string& filename);
//         void SetUpGasSpecies(const std::string& filename, const std::vector<std::string>& user_gas_species_list);
//         void ReadGasSpeciesFile(const std::string& filename, std::vector<std::shared_ptr<GasSpecies> >& gas_species_list);

//         void ReadBindingEnergyFile(const std::string& filename);
//         void SetDustMantleSpecies();
//         void ReadEnthalpyOfFormationFile(const std::string& filename);
//         void GenerateDustSpecies();
//         void SetSpecialSpeciesIndex();
//         void UpdateMantleSpeciesDiffusionBarriers();

//         ElementManager* ptr_element_manager_;

//         std::vector<std::shared_ptr<Species> > species_list_;
//         std::vector<std::shared_ptr<GasSpecies> > gas_species_list_;
//         std::vector<std::shared_ptr<DustSurfaceSpecies> > dust_surface_species_list_;
//         std::vector<std::shared_ptr<DustMantleSpecies> > dust_mantle_species_list_;
//         std::vector<std::shared_ptr<DustSpecies> > dust_species_list_;

//         std::vector<std::string> species_name_list_;
//         std::vector<std::string> gas_species_name_list_;

//         DustSpeciesModelParameters dust_species_model_parameters_;

//         std::size_t number_of_total_species_;
//         std::size_t number_of_gas_species_;
//         std::size_t number_of_dust_surface_species_;
//         std::size_t number_of_dust_mantle_species_;
//         std::size_t number_of_dust_species_;

//         std::size_t index_electron_;
//         std::size_t index_H_;
//         std::size_t index_H2_;
//         std::size_t index_He_;
//         std::size_t index_sH_;
//         std::size_t index_sH2_;
//         std::size_t index_mH_;
//         std::size_t index_mH2_;
//         std::size_t index_sH2O_;

//         bool is_dust_species_;
//         bool is_dust_surface_reaction_;
//         bool is_three_phase_reaction_;
//     };
// }

// #endif /* SPECIES_MANAGER_HPP */