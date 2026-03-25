#ifndef SPECIES_MANAGER_HPP_
#define SPECIES_MANAGER_HPP_

#include "nicole/config/input_config.hpp"
#include "nicole/element/element_manager.hpp"
#include "nicole/species/dust_mantle_species.hpp"
#include "nicole/species/dust_species.hpp"
#include "nicole/species/dust_surface_species.hpp"
#include "nicole/species/gas_species.hpp"
#include "nicole/species/species.hpp"

namespace nicole {
    class SpeciesManager {
    public:
        SpeciesManager(
            ElementManager *ptr_element_manager,
            InputConfig& config
        );

        SpeciesManager(
            ElementManager *ptr_element_manager,
            InputConfig& input,
            const std::vector<std::string> user_gas_species_list
        );

        void CheckSpeciesManager(const std::string& filename) const;

        // Getter
        std::shared_ptr<Species> GetSpecies(const std::size_t index) const;
        std::string GetSpeciesName(const std::size_t index) const;
        int GetSpeciesCharge(const std::size_t index) const;
        Real GetSpeciesMass(const std::size_t index) const;
        const std::vector<std::size_t>& GetSpeciesElementComposition(const std::size_t index) const;
        Real GetSpeciesBindingEnergyOnH2Oice(const std::size_t index) const;
        Real GetSpeciesBindingEnergyOnSilicate(const std::size_t index) const;
        Real GetSpeciesEnthalpyOfFormation(const std::size_t index) const;
        std::size_t GetTotalNumberOfSpecies() const { return total_number_of_species_; }

        // Find functions
        std::shared_ptr<Species> FindSpeciesByName(const std::string& name) const;
        std::shared_ptr<GasSpecies> FindGasSpeciesByName(const std::string& species_name);
        std::shared_ptr<DustSurfaceSpecies> FindDustSurfaceSpeciesByName(const std::string& species_name);
        std::shared_ptr<DustMantleSpecies> FindDustMantleSpeciesByName(const std::string& species_name);
        std::shared_ptr<DustSpecies> FindDustSpeciesByBinNumberAndCharge(const std::size_t bin_number, const int charge) const;
        SpeciesID FindSpeciesID(const std::string& species_name) const;

        // Boolean functions
        bool IsChemicalSpeciesByName(const std::string& species_name) const;
        bool IsDustSurfaceSpecies(const std::size_t index) const;
        bool IsDustMantleSpecies(const std::size_t index) const;
        bool IsDustSurfaceOrDustMantleSpecies(const std::size_t index) const;

        // Friend classes that have access to SpeciesManager's private members
        friend class ReactionManager;
        friend class ReactionSimulator;
        friend class NonIdealMHDeffect;

    protected:
        void ReadGasSpeciesFile(const std::string& file_path, std::vector<std::shared_ptr<GasSpecies> >& gas_species_list);
        void ReadBindingEnergyFile(const std::string& file_path);
        void ReadEnthalpyOfFormationFile(const std::string& file_path);

        void SetUpGasSpecies(InputConfig& config);
        void SetUpGasSpecies(InputConfig& config, const std::vector<std::string>& user_gas_species_list);
        
        void GenerateDustSpecies();
        void SetUpDustMantleSpecies();
        void UpdateMantleSpeciesDiffusionBarriers();
        void SetUpDustRelatedSpecies(InputConfig& input);

        void SetUpSpecialSpeciesIndex();
        void SetUpSpeciesNameList();

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

#endif /* SPECIES_MANAGER_HPP_ */
