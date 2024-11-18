#ifndef DUST_HPP_
#define DUST_HPP_

#include "typedefs.hpp"
#include "physical_constants.hpp"

class Dust
{
public:
    
    Dust();
    ~Dust();

    void SetDust();

    void SetSingleSizeParameter(double d_size, double dust_to_gas_mass_ratio, double dust_internal_density);
    void SetSingleSize(double d_size);
    void SetSizeDistributionParameter(double power_index_of_distribution, double min_size, double max_size,
                  double dust_to_gas_mass_ratio, double dust_internal_density);
    void SetMinDustSize(double min_size);
    void SetMaxDustSize(double min_size);

    void   SingleSizeModel();
    
    void   SizeDistributionModel();
    void   NormalizationConstantOfDistribution();

    void AggregateModel(double nH, double fdg, double N);

protected:

    double MaxDustRadiusSthBin(UnsignedInt sth_bin);
    double MinDustRadiusSthBin(UnsignedInt sth_bin);
    double MomentOfDustSizeDistribution(double n);
    double AverageMomentOfDustSizeDistribution(double n);
    double MomentOfDustSizeDistributionSthBin(double n, UnsignedInt sth_bin);
    double AverageMomentOfDustSizeDistributionSthBin(double n, UnsignedInt sth_bin);
    double DustRadiusSthBin(UnsignedInt sth_bin);
    double DustCrossSectinSthBin(UnsignedInt sth_bin);
    double DustVolumeSthBin(UnsignedInt sth_bin);

    void testDust();

    struct DustParams {
        UnsignedInt nbin_;
        double fdg_;
        double power_index_;
        double d_size_;
        double min_size_;
        double max_size_;
        double nc_;
        DoubleVector radius_;
        DoubleVector cross_sections_;
        DoubleVector mass_;
        DoubleVector abundances_;
        // DoubleVector reduced_temp_;

        double total_abundance_;
        double internal_density_;

        bool  flag_dust_evapolation_;

    } dust_p_;

private:
};

#endif // DUST_HPP_