#include "nicole/species/gas_species.hpp"

namespace nicole {
    GasSpecies::GasSpecies(
        SpeciesID id,
        const std::string& name,
        int charge,
        const std::vector<std::size_t>& element_composition,
        ElementManager* ptr_element_manager
    ) : Species(
            id,
            name,
            charge,
            element_composition,
            SpeciesType::Gas,
            ptr_element_manager
        ),
        corresponding_dust_surface_species_(nullptr),
        corresponding_dust_mantle_species_(nullptr)
    { }
}
