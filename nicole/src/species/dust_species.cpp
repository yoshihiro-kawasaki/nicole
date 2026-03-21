#include <iostream>

#include "nicole/species/dust_species.hpp"

namespace nicole {
    DustSpeciesModelParameters::DustSpeciesModelParameters(InputConfig& input)
        : number_of_dust_bins_(input.GetInt("number_of_dust_bins")),
        min_dust_charge_number_(input.GetInt("min_dust_charge_number")),
        max_dust_charge_number_(input.GetInt("max_dust_charge_number")),
        dust_internal_density_(input.GetReal("dust_internal_density")),
        dust_minimum_size_(input.GetReal("dust_minimum_size")),
        dust_maximum_size_(input.GetReal("dust_maximum_size")),
        power_index_(input.GetReal("power_index")),
        dust_to_gas_mass_ratio_(input.GetReal("dust_to_gas_mass_ratio")),
        is_dust_size_distribution_((number_of_dust_bins_ > 1))  // Determines if the dust size distribution is considered
    {
        // Validate input parameters for the dust species model

        // Check if the number of dust bins is valid (greater than 0)
        if (number_of_dust_bins_ <= 0) {
            std::cerr << "Error: number_of_dust_bins_ <= 0 : " << number_of_dust_bins_ << std::endl;
            std::cerr << "       number_of_dust_bins_ must be > 0." << std::endl;
            throw std::runtime_error("number_of_dust_bins_ <= 0");
        }

        // Check if the minimum dust charge number is valid (must be less than 0)
        if (min_dust_charge_number_ >= 0) {
            std::cerr << "Error: min_dust_charge_number_ >= 0 : " << min_dust_charge_number_ << std::endl;
            std::cerr << "       min_dust_charge_number must be < 0." << std::endl;
            throw std::runtime_error("min_dust_charge_number >= 0");
        }

        // Check if the maximum dust charge number is valid (must be greater than 0)
        if (max_dust_charge_number_ <= 0) {
            std::cerr << "Error: max_dust_charge_number_ <= 0 : " << max_dust_charge_number_ << std::endl;
            std::cerr << "       max_dust_charge_number must be > 0." << std::endl;
            throw std::runtime_error("max_dust_charge_number <= 0");
        }

        // Check if the dust internal density is valid (must be greater than 0)
        if (dust_internal_density_ <= 0.0) {
            std::cerr << "Error: dust_internal_density <= 0 : " << dust_internal_density_ << std::endl;
            std::cerr << "       dust_internal_density must be > 0." << std::endl;
            throw std::runtime_error("dust_internal_density <= 0");
        }

        // Check if the minimum dust size is valid (must be greater than 0)
        if (dust_minimum_size_ <= 0.0) {
            std::cerr << "Error: dust_minimum_size <= 0 : " << dust_minimum_size_ << std::endl;
            std::cerr << "       dust_minimum_size must be > 0." << std::endl;
            throw std::runtime_error("dust_minimum_size <= 0");
        }

        // Check if the maximum dust size is valid (must be greater than 0)
        if (dust_maximum_size_ <= 0.0) {
            std::cerr << "Error: dust_maximum_size <= 0 : " << dust_maximum_size_ << std::endl;
            std::cerr << "       dust_maximum_size must be > 0." << std::endl;
            throw std::runtime_error("dust_maximum_size <= 0");
        }

        // Check if the dust-to-gas mass ratio is valid (must be greater than 0)
        if (dust_to_gas_mass_ratio_ <= 0.0) {
            std::cerr << "Error: dust_to_gas_mass_ratio_ <= 0 : " << dust_to_gas_mass_ratio_ << std::endl;
            std::cerr << "       dust_to_gas_mass_ratio_ must be > 0." << std::endl;
            throw std::runtime_error("dust_to_gas_mass_ratio_ <= 0");
        }

        // If dust size distribution is considered (i.e., number_of_dust_bins_ > 1)
        if (is_dust_size_distribution_) {
            // Check if the dust size parameters are consistent
            if (dust_minimum_size_ == dust_maximum_size_) {
                std::cerr << "Error: dust_minimum_size_ = dust_maximum_size_ : " << dust_minimum_size_ << std::endl;
                std::cerr << "       dust_maximum_size_ must be larger than dust_minimum_size_." << std::endl;
                throw std::runtime_error("dust_minimum_size_ = dust_maximum_size_");
            }

            if (dust_minimum_size_ > dust_maximum_size_) {
                std::cerr << "Error: dust_minimum_size_ > dust_maximum_size_ : " << std::endl;
                std::cerr << "       dust_minimum_size_ = " << dust_minimum_size_ << std::endl;
                std::cerr << "       dust_maximum_size_ = " << dust_maximum_size_ << std::endl;
                std::cerr << "       dust_minimum_size_ must be smaller than dust_maximum_size_." << std::endl;
                throw std::runtime_error("dust_minimum_size_ > dust_maximum_size_");
            }

            // Compute size distribution normalization factor using the power index and dust size range
            size_distribution_normalization_factor_ = 3.0 * (4.0 + power_index_) * dust_to_gas_mass_ratio_ * 1.4 * constants::kProtonMass
                / (4.0 * M_PI * dust_internal_density_ * (std::pow(dust_maximum_size_, 4.0 + power_index_)
                    - std::pow(dust_minimum_size_, 4.0 + power_index_)));

            // Initialize total dust abundances and cross-section per hydrogen atom (nH)
            total_dust_abundances_ = 0.0;
            total_dust_cross_section_per_nH_ = 0.0;

            // Loop through each dust bin to calculate the total abundances and cross-section
            for (std::size_t bin_number = 1; bin_number <= number_of_dust_bins_; ++bin_number) {
                Real dust_abundance = GetDustAbundancesForBin(bin_number);
                Real dust_cross_section = GetDustCrossSectionForBin(bin_number);
                total_dust_abundances_ += dust_abundance;
                total_dust_cross_section_per_nH_ += (dust_abundance * dust_cross_section);
            }

        } else {
            // Single size model (i.e., dust size distribution is not considered)

            // Calculate the dust radius, cross-section, and mass for the single dust size model
            const Real dust_radius = dust_maximum_size_;
            const Real dust_cross_section = M_PI * SQR(dust_radius);
            const Real dust_mass = (4.0 * M_PI / 3.0) * CUB(dust_radius) * dust_internal_density_;

            // Calculate the total dust abundances and cross-section per hydrogen atom (nH)
            total_dust_abundances_ = 1.4 * constants::kProtonMass * dust_to_gas_mass_ratio_ / dust_mass;
            total_dust_cross_section_per_nH_ = dust_cross_section * total_dust_abundances_;

            // Calculate the number of sites per dust grain based on surface density
            number_of_sites_per_a_dust_ = 4.0 * dust_cross_section * kDustSurfaceSitesDensity;
        }
    }


    Real DustSpeciesModelParameters::GetDustRadiusForBin(const std::size_t bin_index) {
        return GetAverageNthMomentOfDustSizeDistributionForBin(1.0, bin_index);
    }


    Real DustSpeciesModelParameters::GetDustCrossSectionForBin(const std::size_t bin_index) {
        return M_PI * GetAverageNthMomentOfDustSizeDistributionForBin(2.0, bin_index);
    }


    Real DustSpeciesModelParameters::GetDustVolumeForBin(const std::size_t bin_index) {
        return (4.0 * M_PI / 3.0) * GetAverageNthMomentOfDustSizeDistributionForBin(3.0, bin_index);
    }


    Real DustSpeciesModelParameters::GetDustAbundancesForBin(const std::size_t bin_index) {
        return GetNthMomentOfDustSizeDistributionForBin(0.0, bin_index);
    }


    Real DustSpeciesModelParameters::GetMinimumDustRadiusForBin(std::size_t bin_index) {
        Real index_s = static_cast<Real>(bin_index - 1) / static_cast<Real>(number_of_dust_bins_);
        Real dust_size = dust_minimum_size_ * std::pow(dust_maximum_size_ / dust_minimum_size_, index_s);
        return dust_size;
    }


    Real DustSpeciesModelParameters::GetMaximumDustRadiusForBin(std::size_t bin_index) {
        Real index_s = static_cast<Real>(bin_index) / static_cast<Real>(number_of_dust_bins_);
        Real dust_size = dust_minimum_size_ * std::pow(dust_maximum_size_ / dust_minimum_size_, index_s);
        return dust_size;
    }


    Real DustSpeciesModelParameters::GetNthMomentOfDustSizeDistribution(const Real moment_order) {
        return (size_distribution_normalization_factor_ / (moment_order + 1.0 + power_index_))
            * (std::pow(dust_maximum_size_, moment_order + 1.0 + power_index_) 
                - std::pow(dust_minimum_size_, moment_order + 1.0 + power_index_));
    }


    Real DustSpeciesModelParameters::GetAverageNthMomentOfDustSizeDistribution(const Real moment_order) {
        Real moment0 = GetNthMomentOfDustSizeDistribution(0.0);
        Real momentn = GetNthMomentOfDustSizeDistribution(moment_order);
        return momentn / moment0;
    }


    Real DustSpeciesModelParameters::GetNthMomentOfDustSizeDistributionForBin(const Real moment_order, const std::size_t bin_index) {
        Real dust_minimum_size_bin = GetMinimumDustRadiusForBin(bin_index);
        Real dust_maximum_size_bin = GetMaximumDustRadiusForBin(bin_index);
        return (size_distribution_normalization_factor_ / (moment_order + 1.0 + power_index_))
            * (std::pow(dust_maximum_size_bin, moment_order + 1.0 + power_index_) 
                - std::pow(dust_minimum_size_bin, moment_order + 1.0 + power_index_));
    }


    Real DustSpeciesModelParameters::GetAverageNthMomentOfDustSizeDistributionForBin(const Real moment_order, const std::size_t bin_index) {
        Real moment_an_bin = GetNthMomentOfDustSizeDistributionForBin(moment_order, bin_index);
        Real moment_a0_bin = GetNthMomentOfDustSizeDistributionForBin(0.0, bin_index);
        return moment_an_bin / moment_a0_bin;
    }


    DustSpecies::DustSpecies(
        SpeciesID id,
        const std::string& name,
        int charge,
        std::size_t bin_number,
        DustSpeciesModelParameters* ptr_dust_model
    ) : Species(
            id,
            name,
            charge,
            SpeciesType::Dust
        ),
        bin_number_(bin_number),
        ptr_dust_model_(ptr_dust_model)
    {
        // Validate bin_number_: it should be greater than zero
        if (bin_number_ <= 0) {
            std::cerr << "Error: bin_number_ <= 0 : " << bin_number_ << std::endl;
            std::cerr << "       bin_number_ must be > 0." << std::endl;
            throw std::runtime_error("bin_number_ <= 0");  // Throw error if bin_number_ is invalid
        }

        // Validate ptr_dust_model_: ensure it is not a null pointer
        if (ptr_dust_model_ == nullptr) {
            std::cerr << "Error: ptr_dust_model_ is nullptr" << std::endl;
            throw std::runtime_error("DustSpeciesModelParameters pointer is null");  // Throw error if model pointer is null
        }

        // Calculate physical properties based on the dust size distribution model
        if (ptr_dust_model_->is_dust_size_distribution_) {
            // If a dust size distribution model is used, calculate the properties for the specific bin
            radius_ = ptr_dust_model_->GetDustRadiusForBin(bin_number);  // Get the radius for the bin
            cross_section_ = ptr_dust_model_->GetDustCrossSectionForBin(bin_number);  // Get the cross section for the bin
            volume_ = ptr_dust_model_->GetDustVolumeForBin(bin_number);  // Get the volume for the bin
            mass_ = volume_ * ptr_dust_model_->dust_internal_density_;  // Calculate mass based on volume and dust internal density
        } else {
            // If a single size model is used, calculate properties using the maximum dust size
            radius_ = ptr_dust_model_->dust_maximum_size_;  // Set radius to the maximum dust size
            cross_section_ = SQR(radius_) * M_PI;  // Calculate cross section as π * radius^2
            volume_ = (4.0 * M_PI / 3.0) * CUB(radius_);  // Calculate volume as 4/3 * π * radius^3
            mass_   = volume_ * ptr_dust_model_->dust_internal_density_;  // Calculate mass based on volume and internal density
        }
    }
}
