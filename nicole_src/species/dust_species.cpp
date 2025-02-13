/**
 * @file dust_species.cpp
 * @brief Implements the functions for the DustSpecies and DustSpeciesModelParameters classes.
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#include "dust_species.hpp"

namespace nicole
{
    /**
     * @brief Constructor for DustSpeciesModelParameters.
     * 
     * Initializes various dust species model parameters from the input configuration.
     * Throws runtime errors if input values are invalid.
     * 
     * @param input Input configuration object containing the parameters.
     */
    DustSpeciesModelParameters::DustSpeciesModelParameters(InputConfig& input)
        : number_of_dust_bins_(input.GetInt("number_of_dust_bins")),
        min_dust_charge_number_(input.GetInt("min_dust_charge_number")),
        max_dust_charge_number_(input.GetInt("max_dust_charge_number")),
        dust_internal_density_(input.GetDouble("dust_internal_density")),
        dust_minimum_size_(input.GetDouble("dust_minimum_size")),
        dust_maximum_size_(input.GetDouble("dust_maximum_size")),
        power_index_(input.GetDouble("power_index")),
        dust_to_gas_mass_ratio_(input.GetDouble("dust_to_gas_mass_ratio")),
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
                double dust_abundance = GetDustAbundancesForBin(bin_number);
                double dust_cross_section = GetDustCrossSectionForBin(bin_number);
                total_dust_abundances_ += dust_abundance;
                total_dust_cross_section_per_nH_ += (dust_abundance * dust_cross_section);
            }

        } else {
            // Single size model (i.e., dust size distribution is not considered)

            // Calculate the dust radius, cross-section, and mass for the single dust size model
            const double dust_radius = dust_maximum_size_;
            const double dust_cross_section = M_PI * SQR(dust_radius);
            const double dust_mass = (4.0 * M_PI / 3.0) * CUB(dust_radius) * dust_internal_density_;

            // Calculate the total dust abundances and cross-section per hydrogen atom (nH)
            total_dust_abundances_ = 1.4 * constants::kProtonMass * dust_to_gas_mass_ratio_ / dust_mass;
            total_dust_cross_section_per_nH_ = dust_cross_section * total_dust_abundances_;

            // Calculate the number of sites per dust grain based on surface density
            number_of_sites_per_a_dust_ = 4.0 * dust_cross_section * kDustSurfaceSitesDensity;
        }
    }

    /**
     * @brief Calculates the dust radius for a given dust size bin.
     * Uses the first moment of the dust size distribution to estimate the average radius for the bin.
     * @param bin_index The index of the dust size bin.
     * @return The dust radius for the given bin.
     */
    double DustSpeciesModelParameters::GetDustRadiusForBin(const std::size_t bin_index)
    {
        // Returns the first moment (average radius) of the dust size distribution for the bin
        return GetAverageNthMomentOfDustSizeDistributionForBin(1.0, bin_index);
    }

    /**
     * @brief Calculates the dust cross-section for a given dust size bin.
     * Uses the second moment of the dust size distribution to estimate the average cross-section for the bin.
     * @param bin_index The index of the dust size bin.
     * @return The dust cross-section for the given bin.
     */
    double DustSpeciesModelParameters::GetDustCrossSectionForBin(const std::size_t bin_index)
    {
        // Returns the second moment (average cross-section) of the dust size distribution for the bin
        return M_PI * GetAverageNthMomentOfDustSizeDistributionForBin(2.0, bin_index);
    }

    /**
     * @brief Calculates the dust volume for a given dust size bin.
     * Uses the third moment of the dust size distribution to estimate the average volume for the bin.
     * @param bin_index The index of the dust size bin.
     * @return The dust volume for the given bin.
     */
    double DustSpeciesModelParameters::GetDustVolumeForBin(const std::size_t bin_index)
    {
        // Returns the third moment (average volume) of the dust size distribution for the bin
        return (4.0 * M_PI / 3.0) * GetAverageNthMomentOfDustSizeDistributionForBin(3.0, bin_index);
    }

    /**
     * @brief Calculates the dust abundance for a given dust size bin.
     * Uses the zeroth moment of the dust size distribution to estimate the dust abundance for the bin.
     * @param bin_index The index of the dust size bin.
     * @return The dust abundance for the given bin.
     */
    double DustSpeciesModelParameters::GetDustAbundancesForBin(const std::size_t bin_index)
    {
        // Returns the zeroth moment (abundance) of the dust size distribution for the bin
        return GetNthMomentOfDustSizeDistributionForBin(0.0, bin_index);
    }

    /**
     * @brief Calculates the minimum dust radius for a given dust size bin.
     * The minimum radius is calculated based on the bin index using a logarithmic scale between minimum and maximum dust size.
     * @param bin_index The index of the dust size bin.
     * @return The minimum dust radius for the given bin.
     */
    double DustSpeciesModelParameters::GetMinimumDustRadiusForBin(std::size_t bin_index)
    {
        // Calculate the dust radius at the lower bound of the bin's size range
        double index_s = static_cast<double>(bin_index - 1) / static_cast<double>(number_of_dust_bins_);
        double dust_size = dust_minimum_size_ * std::pow(dust_maximum_size_ / dust_minimum_size_, index_s);
        return dust_size;
    }

    /**
     * @brief Calculates the maximum dust radius for a given dust size bin.
     * The maximum radius is calculated based on the bin index using a logarithmic scale between minimum and maximum dust size.
     * @param bin_index The index of the dust size bin.
     * @return The maximum dust radius for the given bin.
     */
    double DustSpeciesModelParameters::GetMaximumDustRadiusForBin(std::size_t bin_index)
    {
        // Calculate the dust radius at the upper bound of the bin's size range
        double index_s = static_cast<double>(bin_index) / static_cast<double>(number_of_dust_bins_);
        double dust_size = dust_minimum_size_ * std::pow(dust_maximum_size_ / dust_minimum_size_, index_s);
        return dust_size;
    }

    /**
     * @brief Calculates the nth moment of the dust size distribution for the entire dust size range.
     * The nth moment is calculated based on the dust size distribution parameters.
     * @param moment_order The order of the moment (e.g., 0 for abundance, 1 for radius, 2 for cross-section, etc.).
     * @return The nth moment of the dust size distribution over the entire dust size range.
     */
    double DustSpeciesModelParameters::GetNthMomentOfDustSizeDistribution(const double moment_order)
    {
        // Calculate the nth moment <a^n> of the dust size distribution over the entire dust size range
        return (size_distribution_normalization_factor_ / (moment_order + 1.0 + power_index_))
            * (std::pow(dust_maximum_size_, moment_order + 1.0 + power_index_) 
                - std::pow(dust_minimum_size_, moment_order + 1.0 + power_index_));
    }

    /**
     * @brief Calculates the average nth moment of the dust size distribution for the entire dust size range.
     * The average nth moment is the ratio of the nth moment to the zeroth moment (abundance).
     * @param moment_order The order of the moment (e.g., 1 for radius, 2 for cross-section, etc.).
     * @return The average nth moment of the dust size distribution over the entire dust size range.
     */
    double DustSpeciesModelParameters::GetAverageNthMomentOfDustSizeDistribution(const double moment_order)
    {
        // Calculate the average nth moment by dividing the nth moment by the zeroth moment (abundance)
        double moment0 = GetNthMomentOfDustSizeDistribution(0.0);  // Zeroth moment (abundance)
        double momentn = GetNthMomentOfDustSizeDistribution(moment_order);  // nth moment
        return momentn / moment0;
    }

    /**
     * @brief Calculates the nth moment of the dust size distribution for a specific dust size bin.
     * This is used to calculate the moments for a given dust size bin, taking into account the bin's size range.
     * @param moment_order The order of the moment (e.g., 0 for abundance, 1 for radius, etc.).
     * @param bin_index The index of the dust size bin.
     * @return The nth moment of the dust size distribution for the specified bin.
     */
    double DustSpeciesModelParameters::GetNthMomentOfDustSizeDistributionForBin(const double moment_order, const std::size_t bin_index)
    {
        // Calculate the nth moment for the specific dust size bin, using the bin's size range
        double dust_minimum_size_bin = GetMinimumDustRadiusForBin(bin_index);
        double dust_maximum_size_bin = GetMaximumDustRadiusForBin(bin_index);
        return (size_distribution_normalization_factor_ / (moment_order + 1.0 + power_index_))
            * (std::pow(dust_maximum_size_bin, moment_order + 1.0 + power_index_) 
                - std::pow(dust_minimum_size_bin, moment_order + 1.0 + power_index_));
    }

    /**
     * @brief Calculates the average nth moment of the dust size distribution for a specific dust size bin.
     * The average nth moment for the bin is the ratio of the nth moment to the zeroth moment for that bin.
     * @param moment_order The order of the moment (e.g., 1 for radius, 2 for cross-section, etc.).
     * @param bin_index The index of the dust size bin.
     * @return The average nth moment of the dust size distribution for the specified bin.
     */
    double DustSpeciesModelParameters::GetAverageNthMomentOfDustSizeDistributionForBin(const double moment_order, const std::size_t bin_index)
    {
        // Calculate the average nth moment for the specific dust size bin
        double moment_an_bin = GetNthMomentOfDustSizeDistributionForBin(moment_order, bin_index);
        double moment_a0_bin = GetNthMomentOfDustSizeDistributionForBin(0.0, bin_index);  // Zeroth moment (abundance)
        return moment_an_bin / moment_a0_bin;
    }

    /**
     * Constructor for DustSpecies class.
     * This constructor initializes a dust species based on the provided parameters.
     * It sets the properties of the dust species depending on whether the model uses a dust size distribution or a single size model.
     * It checks the validity of the input parameters and calculates the physical properties (radius, cross section, volume, and mass).
     * 
     * @param index The index of the species.
     * @param name The name of the dust species.
     * @param charge The charge of the dust species.
     * @param bin_number The bin number corresponding to the dust size.
     * @param ptr_dust_model Pointer to the DustSpeciesModelParameters associated with the dust species.
        */
    DustSpecies::DustSpecies(
        std::size_t index,
        const std::string& name,
        int charge,
        std::size_t bin_number,
        DustSpeciesModelParameters* ptr_dust_species_parameters
    ) : Species(
            index,
            name,
            charge,
            SpeciesType::Dust
        ),
        bin_number_(bin_number),
        ptr_dust_model_(ptr_dust_species_parameters)
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