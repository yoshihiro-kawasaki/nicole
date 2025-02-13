/**
 * @file dust_species.hpp
 * @brief Defines the DustSpecies and DustSpeciesModelParameters classes for dust species management.
 *        The DustSpeciesModelParameters class manages dust-related model parameters, while the 
 *        DustSpecies class handles the properties of individual dust species.
 * @date 2025-02-12
 * @author Y. Kawasaki
 */

#ifndef DUST_SPECIES_HPP
#define DUST_SPECIES_HPP

#include "../nicole_defs.hpp"
#include "species.hpp"

namespace nicole
{
    /**
     * @class DustSpeciesModelParameters
     * @brief This class holds the parameters for the dust species model, including the dust size distribution, 
     *        charge numbers, and other physical properties such as dust internal density and dust-to-gas mass ratios.
     *        It also provides methods to calculate dust properties like radius, cross-section, volume, and abundance 
     *        based on the dust bins and the distribution model.
     */
    class DustSpeciesModelParameters
    {
    public:

        /**
         * @brief Constructor for initializing the dust species model parameters from input configurations.
         * @param input Input configuration object that provides model parameters.
         */
        DustSpeciesModelParameters(InputConfig& input);

        /**
         * @brief Getter for the number of dust bins in the model.
         * @return The number of dust bins.
         */
        inline int GetNumberOfDustBins() const { return number_of_dust_bins_; }

        /**
         * @brief Getter for the minimum dust charge number.
         * @return The minimum dust charge number.
         */
        inline int GetMinDustChargeNumber() const { return min_dust_charge_number_; }

        /**
         * @brief Getter for the maximum dust charge number.
         * @return The maximum dust charge number.
         */
        inline int GetMaxDustChargeNumber() const { return max_dust_charge_number_; }

        /**
         * @brief Getter for the total dust abundance in the model.
         * @return The total dust abundance.
         */
        inline double GetDustTotalAbundance() const { return total_dust_abundances_; }

        /**
         * @brief Getter for the total dust cross-section per hydrogen atom in the model.
         * @return The total dust cross-section per hydrogen atom.
         */
        inline double GetDustTotalCrossSectionPernH() const { return total_dust_cross_section_per_nH_; }

        /**
         * @brief Getter for the number of sites per dust particle.
         * @return The number of sites per dust particle.
         */
        inline double GetNumberOfSitesPerDust() const { return number_of_sites_per_a_dust_; }

        /**
         * @brief Check if the dust size distribution model is being used.
         * @return True if the size distribution model is used, otherwise false.
         */
        inline bool IsSizeDistributionModel() const { return is_dust_size_distribution_; }

        /**
         * @brief Get the dust radius for a specific dust bin.
         * @param bin_index The bin index for the dust size distribution.
         * @return The dust radius for the specified bin.
         */
        double GetDustRadiusForBin(const std::size_t bin_index);

        /**
         * @brief Get the dust cross-section for a specific dust bin.
         * @param bin_index The bin index for the dust size distribution.
         * @return The dust cross-section for the specified bin.
         */
        double GetDustCrossSectionForBin(const std::size_t bin_index);

        /**
         * @brief Get the dust volume for a specific dust bin.
         * @param bin_index The bin index for the dust size distribution.
         * @return The dust volume for the specified bin.
         */
        double GetDustVolumeForBin(const std::size_t bin_index);

        /**
         * @brief Get the dust abundance for a specific dust bin.
         * @param bin_index The bin index for the dust size distribution.
         * @return The dust abundance for the specified bin.
         */
        double GetDustAbundancesForBin(const std::size_t bin_index);

        // Friend class to allow DustSpecies to access private members of DustSpeciesModelParameters
        friend class DustSpecies;

    private:

        /**
         * @brief Get the minimum dust radius for a given dust bin.
         * @param bin_index The bin index for the dust size distribution.
         * @return The minimum dust radius for the specified bin.
         */
        double GetMinimumDustRadiusForBin(const std::size_t bin_index);

        /**
         * @brief Get the maximum dust radius for a given dust bin.
         * @param bin_index The bin index for the dust size distribution.
         * @return The maximum dust radius for the specified bin.
         */
        double GetMaximumDustRadiusForBin(const std::size_t bin_index);

        /**
         * @brief Calculate the nth moment of the dust size distribution.
         * @param moment_order The order of the moment to be calculated.
         * @return The nth moment of the dust size distribution.
         */
        double GetNthMomentOfDustSizeDistribution(const double moment_order);

        /**
         * @brief Calculate the average nth moment of the dust size distribution.
         * @param moment_order The order of the moment to be averaged.
         * @return The average nth moment of the dust size distribution.
         */
        double GetAverageNthMomentOfDustSizeDistribution(const double moment_order);

        /**
         * @brief Calculate the nth moment of the dust size distribution for a specific bin.
         * @param moment_order The order of the moment to be calculated.
         * @param bin_index The bin index for the dust size distribution.
         * @return The nth moment of the dust size distribution for the specified bin.
         */
        double GetNthMomentOfDustSizeDistributionForBin(const double moment_order, const std::size_t bin_index);

        /**
         * @brief Calculate the average nth moment of the dust size distribution for a specific bin.
         * @param moment_order The order of the moment to be averaged.
         * @param bin_index The bin index for the dust size distribution.
         * @return The average nth moment of the dust size distribution for the specified bin.
         */
        double GetAverageNthMomentOfDustSizeDistributionForBin(const double moment_order, const std::size_t bin_index);

        // Number of dust bins in the model
        std::size_t number_of_dust_bins_;

        // Minimum and maximum charge number of dust species
        int min_dust_charge_number_;
        int max_dust_charge_number_;

        // Internal density of the dust particles [g cm^-3]
        double dust_internal_density_;

        // Minimum and maximum size of the dust particles [cm]
        double dust_minimum_size_;
        double dust_maximum_size_;

        // Power index for the dust size distribution (n(a) = C a^(power_index))
        double power_index_;

        double dust_to_gas_mass_ratio_;

        // Normalization factor for size distribution
        double size_distribution_normalization_factor_;

        // Flag indicating if the size distribution model is applied
        bool is_dust_size_distribution_;

        // Total abundance of dust species
        double total_dust_abundances_;

        // Total dust cross-section per hydrogen atom
        double total_dust_cross_section_per_nH_;

        // Number of surface sites per dust particle
        double number_of_sites_per_a_dust_;
    };

    /**
     * @class DustSpecies
     * @brief This class represents an individual dust species, inheriting from the Species class.
     *        It includes properties such as dust radius, cross-section, and volume, and links to the
     *        corresponding DustSpeciesModelParameters object that provides the dust properties.
     */
    class DustSpecies
        : public Species
    {
    private:

        // Radius of the dust species
        double radius_;

        // Cross-sectional area of the dust species
        double cross_section_;

        // Volume of the dust species
        double volume_;

        std::size_t bin_number_;

        // Pointer to the associated DustSpeciesModelParameters
        DustSpeciesModelParameters* ptr_dust_model_;

    public:

        /**
         * @brief Constructor for DustSpecies class.
         * @param index The index of the species.
         * @param name The name of the dust species.
         * @param charge The charge of the dust species.
         * @param bin_number The bin number corresponding to the dust size.
         * @param ptr_dust_model Pointer to the DustSpeciesModelParameters associated with the dust species.
         */
        DustSpecies(
            std::size_t index,
            const std::string& name,
            int charge,
            std::size_t bin_number,
            DustSpeciesModelParameters* ptr_dust_model_
        );

        /**
         * @brief Getter for the dust radius.
         * @return The radius of the dust species.
         */
        inline double GetRadius() const { return radius_; }

        /**
         * @brief Getter for the dust volume.
         * @return The volume of the dust species.
         */
        inline double GetVolume() const { return volume_; }

        /**
         * @brief Getter for the dust cross-section.
         * @return The cross-sectional area of the dust species.
         */
        inline double GetCrossSection() const { return cross_section_; }

        /**
         * @brief Getter for the bin number associated with the dust species.
         * @return The bin number of the dust species.
         */
        inline std::size_t GetBinNumber() const { return bin_number_; }

        // The following methods are disabled (not used) for the DustSpecies class.
        inline const std::vector<std::size_t>& GetElementComposition() const = delete;
        inline ElementManager* GetPtrElementManager() const = delete;
        inline double GetBindingEnergyOnH2Oice() const = delete;
        inline double GetBindingEnergyOnBareSilicate() = delete;
        inline double GetEnthalpyOfFormation() const = delete;
        inline void SetBindingEnergyOnH2Oice(const double binding_energy_on_H2O_ice) = delete;
        inline void SetBindingEnergyOnBareSilicate(const double binding_energy_on_silicate) = delete;
        inline void SetEnthalpyOfFormation(const double enthalpy_of_formation) = delete;

        // Disable mass calculation for dust species
        void CalculateMass() = delete;

    };
}

#endif /* DUST_SPECIES_HPP */


// /**
//  * @file dust_species.hpp
//  * @brief Defines the DustSpecies and DustSpeciesModelParameters classes.
//  * @date 2025-02-12
//  * @author Y. Kawasaki
//  */

// #ifndef DUST_SPECIES_HPP
// #define DUST_SPECIES_HPP

// #include "../nicole_defs.hpp"
// #include "species.hpp"

// namespace nicole
// {
//     /**
//      * @class DustSpeciesModelParameters
//      */
//     class DustSpeciesModelParameters
//     {
//     public:

//         DustSpeciesModelParameters(InputConfig& input);

//         inline int GetNumberOfDustBins() const { return number_of_dust_bins_; }
//         inline int GetMinDustChargeNumber() const { return min_dust_charge_number_; }
//         inline int GetMaxDustChargeNumber() const { return max_dust_charge_number_; }
//         inline double GetDustTotalAbundance() const { return total_dust_abundances_; }
//         inline double GetDustTotalCrossSectionPernH() const { return total_dust_cross_section_per_nH_; }
//         inline double GetNumberOfSitesPerDust() const { return number_of_sites_per_a_dust_; }
//         inline bool IsSizeDistributionModel() const { return is_dust_size_distribution_; }

//         double GetDustRadiusForBin(const std::size_t bin_index);
//         double GetDustCrossSectionForBin(const std::size_t bin_index);
//         double GetDustVolumeForBin(const std::size_t bin_index);
//         double GetDustAbundancesForBin(const std::size_t bin_index);

//         friend class DustSpecies;

//     private:

//         double GetMinimumDustRadiusForBin(const std::size_t bin_index);
//         double GetMaximumDustRadiusForBin(const std::size_t bin_index);
//         double GetNthMomentOfDustSizeDistribution(const double moment_order);
//         double GetAverageNthMomentOfDustSizeDistribution(const double moment_order);
//         double GetNthMomentOfDustSizeDistributionForBin(const double moment_order, const std::size_t bin_index);
//         double GetAverageNthMomentOfDustSizeDistributionForBin(const double moment_order, const std::size_t bin_index);

//         std::size_t number_of_dust_bins_;
//         int min_dust_charge_number_;
//         int max_dust_charge_number_;
//         double dust_internal_density_;
//         double dust_minimum_size_;
//         double dust_maximum_size_;
//         double power_index_;
//         double dust_to_gas_mass_ratio_;
//         double size_distribution_normalization_factor_;
//         bool is_dust_size_distribution_;
//         double total_dust_abundances_;
//         double total_dust_cross_section_per_nH_;
//         double number_of_sites_per_a_dust_;
//     };

//     /**
//      * @class DustSpecies
//      */
//     class DustSpecies
//         : public Species
//     {
//     private:

//         double radius_;
//         double cross_section_;
//         double volume_;
//         std::size_t bin_number_;
//         DustSpeciesModelParameters* ptr_dust_model_;

//     public:

//         DustSpecies(
//             std::size_t index,
//             const std::string& name,
//             int charge,
//             std::size_t bin_number,
//             DustSpeciesModelParameters* ptr_dust_model_
//         );

//         // Get関数
//         inline double GetRadius() const { return radius_; }
//         inline double GetVolume() const { return volume_; }
//         inline double GetCrossSection() const { return cross_section_; }
//         inline std::size_t GetBinNumber() const { return bin_number_; }

//         // 以下の基底クラスの関数を使えなくする
//         inline const std::vector<std::size_t>& GetElementComposition() const = delete;
//         inline ElementManager* GetPtrElementManager() const = delete;
//         inline double GetBindingEnergyOnH2Oice() const = delete;
//         inline double GetBindingEnergyOnBareSilicate() = delete;
//         inline double GetEnthalpyOfFormation() const = delete;
//         inline void SetBindingEnergyOnH2Oice(const double binding_energy_on_H2O_ice) = delete;
//         inline void SetBindingEnergyOnBareSilicate(const double binding_energy_on_silicate) = delete;
//         inline void SetEnthalpyOfFormation(const double enthalpy_of_formation) = delete;

//         void CalculateMass() = delete;

//     };
// }

// #endif /* DUST_SPECIES_HPP */