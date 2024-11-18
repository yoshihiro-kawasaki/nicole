#include "dust.hpp"

Dust::Dust()
{
}

Dust::~Dust()
{
    
}

void Dust::SetSingleSizeParameter(double d_size, double dust_to_gas_mass_ratio, double dust_internal_density)
{
    dust_p_.d_size_ = d_size;
    dust_p_.fdg_ = dust_to_gas_mass_ratio;
    dust_p_.internal_density_ = dust_internal_density;
    return;
}

void Dust::SetSizeDistributionParameter(double power_index, double min_size, double max_size,
                  double dust_to_gas_mass_ratio, double dust_internal_density)
{
    dust_p_.power_index_ = power_index;
    dust_p_.min_size_ = min_size;
    dust_p_.max_size_ = max_size;
    dust_p_.fdg_ = dust_to_gas_mass_ratio;
    dust_p_.internal_density_ = dust_internal_density;

    dust_p_.nc_ = 3.0*(4.0 - dust_p_.power_index_)* dust_p_.fdg_ * 1.4 * HYDROGEN_MASS
        / (4*M_PI*dust_p_.internal_density_*(std::pow(dust_p_.max_size_, 4.0 - dust_p_.power_index_) 
            - std::pow(dust_p_.min_size_, 4.0 - dust_p_.power_index_)));

    return;
}

void Dust::SetSingleSize(double d_size)
{
    dust_p_.d_size_ = d_size;
    return;
}

void Dust::SetMinDustSize(double min_size)
{
    dust_p_.min_size_ = min_size;
    return;
}

void Dust::SetMaxDustSize(double max_size)
{
    dust_p_.max_size_ = max_size;
    return;
}

void Dust::SingleSizeModel()
{
    int s = 0;
    if (dust_p_.nbin_ > 1) {
        std::cerr << "error : Dust::SingleSizeModel, number of bins = " << dust_p_.nbin_ << std::endl;
        std::exit(1);
    } 
    dust_p_.radius_[s] = dust_p_.d_size_;
    dust_p_.cross_sections_[s] = M_PI * dust_p_.d_size_*dust_p_.d_size_;
    dust_p_.mass_[s] = (4.0*M_PI/3.0)*std::pow(dust_p_.d_size_, 3.0) * dust_p_.internal_density_;
    dust_p_.abundances_[s] = dust_p_.fdg_ * 1.4 * HYDROGEN_MASS / dust_p_.mass_[s];
    // reduced_temp_[s] = d_size_ * BOLTZMANN_CONSTANT * ptr_cp_->temperature_ / (CHARGE_UNIT*CHARGE_UNIT);
    dust_p_.total_abundance_ = dust_p_.abundances_[s];

    return;
}

// dust size distribution
void Dust::SizeDistributionModel()
{
    UnsignedInt s, sth_bin;
    dust_p_.total_abundance_ = 0;
    for (s = 0; s < dust_p_.nbin_; ++s) {
        sth_bin = s + 1;
        dust_p_.radius_[s] = DustRadiusSthBin(sth_bin);
        dust_p_.cross_sections_[s] = DustCrossSectinSthBin(sth_bin);
        dust_p_.mass_[s] = DustVolumeSthBin(sth_bin) * dust_p_.internal_density_;
        // reduced_temp_[s] = DustRadiusSthBin(sth_bin) * BOLTZMANN_CONSTANT * ptr_cp_->temperature_ / (CHARGE_UNIT*CHARGE_UNIT);
        dust_p_.abundances_[s] = MomentOfDustSizeDistributionSthBin(0.0, sth_bin);
        // if (flag_dust_evapolation_) {
        //     dust_grain_number_density_s_[s] = MomentOfDustSizeDistributionSthBin(0.0, sth_bin) * DustEvapolation(rv.temperature_);
        // } else {
        //     dust_grain_number_density_s_[s] = MomentOfDustSizeDistributionSthBin(0.0, sth_bin);
        // }
        dust_p_.total_abundance_ += dust_p_.abundances_[s];
    }

    return;
}

void Dust::NormalizationConstantOfDistribution()
{
    dust_p_.nc_ = 3.0*(4.0 - dust_p_.power_index_)* dust_p_.fdg_ * 1.4 * HYDROGEN_MASS
        / (4*M_PI*dust_p_.internal_density_*(std::pow(dust_p_.max_size_, 4.0 - dust_p_.power_index_) 
            - std::pow(dust_p_.min_size_, 4.0 - dust_p_.power_index_)));
}

double Dust::MaxDustRadiusSthBin(UnsignedInt sth_bin)
{
    double index_s, dust_size;
    index_s = ((double)sth_bin) / (double)dust_p_.nbin_;
    dust_size = dust_p_.min_size_ * std::pow(dust_p_.max_size_/dust_p_.min_size_, index_s);
    return dust_size;
}

double Dust::MinDustRadiusSthBin(UnsignedInt sth_bin)
{
    double index_s, dust_size;
    index_s = ((double)sth_bin - 1) / (double)dust_p_.nbin_;
    dust_size = dust_p_.min_size_ * std::pow(dust_p_.max_size_/dust_p_.min_size_, index_s);
    return dust_size;
}

double Dust::MomentOfDustSizeDistribution(double n)
{
    // ダスト半径に関するn次のモーメント <a^n>
    double moment_an;
    moment_an = (dust_p_.nc_/(n + 1.0 - dust_p_.power_index_))
                * (std::pow(dust_p_.max_size_, n + 1.0 - dust_p_.power_index_) 
                   - std::pow(dust_p_.min_size_, n +1.0 - dust_p_.power_index_));
    return moment_an;
}

double Dust::AverageMomentOfDustSizeDistribution(double n)
{
    double moment_an, moment_a0;
    moment_an = MomentOfDustSizeDistribution(n);
    moment_a0 = MomentOfDustSizeDistribution(0.0);
    return moment_an / moment_a0;
}

double Dust::MomentOfDustSizeDistributionSthBin(double n, UnsignedInt sth_bin)
{
    double momnet_an_sth;
    double dust_size_sth_min, dust_size_sth_max;

    dust_size_sth_min = MinDustRadiusSthBin(sth_bin);
    dust_size_sth_max = MaxDustRadiusSthBin(sth_bin);
    momnet_an_sth = (dust_p_.nc_/(n + 1.0 - dust_p_.power_index_))
                * (std::pow(dust_size_sth_max, n + 1.0 - dust_p_.power_index_) 
                    - std::pow(dust_size_sth_min, n + 1.0 - dust_p_.power_index_));
    
    return momnet_an_sth;
}

double Dust::AverageMomentOfDustSizeDistributionSthBin(double n, UnsignedInt sth_bin)
{
    double moment_an_sth, moment_a0_sth;
    moment_an_sth = MomentOfDustSizeDistributionSthBin(n, sth_bin);
    moment_a0_sth = MomentOfDustSizeDistributionSthBin(0.0, sth_bin);
    return moment_an_sth / moment_a0_sth;
}

double Dust::DustRadiusSthBin(UnsignedInt sth_bin)
{
    return AverageMomentOfDustSizeDistributionSthBin(1.0, sth_bin);
}

double Dust::DustCrossSectinSthBin(UnsignedInt sth_bin)
{
    return M_PI * AverageMomentOfDustSizeDistributionSthBin(2.0, sth_bin);
}

double Dust::DustVolumeSthBin(UnsignedInt sth_bin)
{
    return (4.0*M_PI/3.0) * AverageMomentOfDustSizeDistributionSthBin(3.0, sth_bin);
}


void Dust::AggregateModel(double nH, double fdg, double N)
{
    double a0, m0, sigma0, nd0, rho_d;
    double D = 2.0;
    double ad, sigma, nd;
    double rho_0 = 1.4;
    a0 = 1.0e-6;
    m0 = (4.0*M_PI/3.0) * std::pow(a0, 3.0) * rho_0;
    sigma0 = a0*a0 * M_PI;
    rho_d = fdg * 1.4*HYDROGEN_MASS*nH;
    nd0 = rho_d / m0;

    ad = a0 * std::pow(N, 1.0/D);
    sigma = sigma0 * N;
    nd = nd0 / N;

    dust_p_.radius_[0] = ad;
    dust_p_.cross_sections_[0] = sigma;
    dust_p_.mass_[0] = N * m0;
    dust_p_.abundances_[0] = nd0 / nH;
    dust_p_.total_abundance_ = dust_p_.abundances_[0];

    // std::cout << ad << std::endl;

    return;
}