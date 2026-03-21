#ifndef DUST_SPECIES_HPP_
#define DUST_SPECIES_HPP_

#include "nicole/config/input_config.hpp"
#include "nicole/nicole_defs.hpp"
#include "nicole/species/species.hpp"

namespace nicole {
    class DustSpeciesModelParameters {
    public:
        DustSpeciesModelParameters(InputConfig& config);

        inline std::size_t GetNumberOfDustBins() const { return number_of_dust_bins_; }
        inline int GetMinDustChargeNumber() const { return min_dust_charge_number_; }
        inline int GetMaxDustChargeNumber() const { return max_dust_charge_number_; }
        inline Real GetDustTotalAbundance() const { return total_dust_abundances_; }
        inline Real GetDustTotalCrossSectionPernH() const { return total_dust_cross_section_per_nH_; }
        inline Real GetNumberOfSitesPerDust() const { return number_of_sites_per_a_dust_; }
        inline bool IsSizeDistributionModel() const { return is_dust_size_distribution_; }
        Real GetDustRadiusForBin(const std::size_t bin_index);
        Real GetDustCrossSectionForBin(const std::size_t bin_index);
        Real GetDustVolumeForBin(const std::size_t bin_index);
        Real GetDustAbundancesForBin(const std::size_t bin_index);

        friend class DustSpecies;

    private:
        Real GetMinimumDustRadiusForBin(const std::size_t bin_index);
        Real GetMaximumDustRadiusForBin(const std::size_t bin_index);
        Real GetNthMomentOfDustSizeDistribution(const Real moment_order);
        Real GetAverageNthMomentOfDustSizeDistribution(const Real moment_order);
        Real GetNthMomentOfDustSizeDistributionForBin(const Real moment_order, const std::size_t bin_index);
        Real GetAverageNthMomentOfDustSizeDistributionForBin(const Real moment_order, const std::size_t bin_index);

        std::size_t number_of_dust_bins_;

        // Minimum and maximum charge number of dust species
        int min_dust_charge_number_;
        int max_dust_charge_number_;

        // Internal density of the dust particles [g cm^-3]
        Real dust_internal_density_;

        // Minimum and maximum size of the dust particles [cm]
        Real dust_minimum_size_;
        Real dust_maximum_size_;

        // Power index for the dust size distribution (n(a) = C a^(power_index))
        Real power_index_;

        Real dust_to_gas_mass_ratio_;

        // Normalization factor for size distribution
        Real size_distribution_normalization_factor_;

        // Flag indicating if the size distribution model is applied
        bool is_dust_size_distribution_;

        // Total abundance of dust species
        Real total_dust_abundances_;

        // Total dust cross-section per hydrogen atom
        Real total_dust_cross_section_per_nH_;

        // Number of surface sites per dust particle
        Real number_of_sites_per_a_dust_;  
    };


    class DustSpecies : public Species {
    public:
        DustSpecies(
            SpeciesID id,
            const std::string& name,
            int charge,
            std::size_t bin_number,
            DustSpeciesModelParameters* ptr_dust_model_
        );

        // Getter
        inline Real GetRadius() const { return radius_; }
        inline Real GetVolume() const { return volume_; }
        inline Real GetCrossSection() const { return cross_section_; }
        inline std::size_t GetBinNumber() const { return bin_number_; }

        // The following methods are disabled (not used) for the DustSpecies class.
        inline const std::vector<std::size_t>& GetElementComposition() const = delete;
        inline ElementManager* GetPtrElementManager() const = delete;
        inline Real GetBindingEnergyOnH2Oice() const = delete;
        inline Real GetBindingEnergyOnBareSilicate() = delete;
        inline Real GetEnthalpyOfFormation() const = delete;
        inline void SetBindingEnergyOnH2Oice(const Real binding_energy_on_H2O_ice) = delete;
        inline void SetBindingEnergyOnBareSilicate(const Real binding_energy_on_silicate) = delete;
        inline void SetEnthalpyOfFormation(const Real enthalpy_of_formation) = delete;
        void CalculateMass() = delete;
        
    private:
        Real radius_;
        Real cross_section_;
        Real volume_;
        std::size_t bin_number_;
        DustSpeciesModelParameters* ptr_dust_model_;
    };
}

#endif /* DUST_SPECIES_HPP_ */
