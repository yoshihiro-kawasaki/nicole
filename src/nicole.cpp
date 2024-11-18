#include "nicole.hpp"

_USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>

Nicole::Nicole()
{
    chem_p_.number_density_  = 0.0;
    chem_p_.temperature_     = 0.0;
    chem_p_.ionization_rate_ = 0.0;

    relative_tolerance_ = 1.0e-5;
    absolute_tolerance_ = 1.0e-15;

    integration_time_ = 1.0e6 * SOLAR_YEAR; // 3.1557e13 [s]
    end_time_         = 0.0;

    flag_dust_evapolation_ = false;
    flag_jacobian_ = false;
    flag_error_function_ = false;
    flag_rate_coefficient_ = false;

    number_of_steps_ = 0;
    number_of_func_  = 0;
    number_of_jac_   = 0;

    flag_thermal_ionization_ = false;
    magnetic_field_ = 0.0;
}

Nicole::~Nicole()
{

}

UnsignedInt Nicole::IndexSpecies(const std::string chemical_species_name, const StringVector chemical_species_list)
{
    auto index = std::find(chemical_species_list.begin(), chemical_species_list.end(), chemical_species_name);
    if (index == chemical_species_list.end()) {
        return NOT_FOUND_SPECIES_LIST;
    } else {
        int index_num = std::distance(chemical_species_list.begin(), index);
        return index_num;
    }
}


void Nicole::ReadFile()
{
    input_.ReadFile("input/input_data.txt");
    std::string chemical_species_file_name;
    std::string chemical_reactions_file_name;
    std::string chemical_abundances_file_name;

    chemical_species_file_name = input_.Get("species_file");
    chemical_reactions_file_name = input_.Get("reaction_file");
    chemical_abundances_file_name = input_.Get("abundance_file");

    number_of_species_ = std::stoi(input_.Get("total_spcies"));
    number_of_reactions_ = std::stoi(input_.Get("total_reaction"));
    number_of_dust_grain_species_ = std::stoi(input_.Get("dust_grain_species"));

    // open file
    std::ifstream species_file(chemical_species_file_name, std::ios::in);
    FAILED_TO_OPEN(species_file, chemical_species_file_name);

    std::ifstream reaction_file(chemical_reactions_file_name, std::ios::in);
    FAILED_TO_OPEN(reaction_file, chemical_reactions_file_name);

    std::ifstream abundances_file(chemical_abundances_file_name, std::ios::in);
    FAILED_TO_OPEN(abundances_file, chemical_abundances_file_name);

    std::string str;
    StringVector line;

    // read species
    while (std::getline(species_file, str)) {
        line = Utility::Split(str, ' ');
        chem_p_.species_.push_back(line[1]);
        chem_p_.species_mass_.push_back(std::stod(line[2])*PROTON_MASS);
        chem_p_.species_charge_.push_back(std::stod(line[3]));
    }

    for (UnsignedInt i = 0; i < number_of_species_; ++i) {
        if (chem_p_.species_[i] == "e-") {
            electron_index_ = i;
        } else if (chem_p_.species_[i] == "H2") {
            molc_hydogrn_index_ = i;
            massH2_ = chem_p_.species_mass_[i];
        } else if (chem_p_.species_[i] == "He") {
            helium_index_ = i;
            massHe_ = chem_p_.species_mass_[i];
        } else if (chem_p_.species_[i] == "H") {
            hydrogen_index_ = i;
            massH_ = chem_p_.species_mass_[i];
        } else if (chem_p_.species_[i][0] == 'G') {
            chem_p_.index_dust_grains_.push_back(i);
        }
    }

    // read reactions
    while (std::getline(reaction_file, str)) {
        line = Utility::Split(str, ' ');
        chem_p_.type_of_reaction_.push_back(line[1]);
        chem_p_.index_reactant1_.push_back(IndexSpecies(line[2], chem_p_.species_));
        chem_p_.index_reactant2_.push_back(IndexSpecies(line[3], chem_p_.species_));
        chem_p_.index_product1_.push_back(IndexSpecies(line[4], chem_p_.species_));
        chem_p_.index_product2_.push_back(IndexSpecies(line[5], chem_p_.species_));
        chem_p_.index_product3_.push_back(IndexSpecies(line[6], chem_p_.species_));
        chem_p_.index_product4_.push_back(IndexSpecies(line[7], chem_p_.species_));
        chem_p_.value_for_reaction1_.push_back(std::stod(line[8]));
        chem_p_.value_for_reaction2_.push_back(std::stod(line[9]));
        chem_p_.value_for_reaction3_.push_back(std::stod(line[10]));
        chem_p_.value_for_reaction4_.push_back(std::stod(line[11]));
        chem_p_.value_for_reaction5_.push_back(std::stod(line[12]));
    }

    // read abundances
    double x;
    while (std::getline(abundances_file, str)) {
        if (str[0] == '#') continue;
        line = Utility::Split(str, ' ');
        x = std::stod(line[1]) * std::stod(line[2]);
        chem_p_.species_abundances_.Insert(line[0], x);
    }

    // close file
    abundances_file.close();
    reaction_file.close();
    species_file.close();

    // resize 
    chem_p_.reaction_rate_coefficient_.resize(number_of_reactions_);
    x_species_.resize(number_of_species_);

    // dust
    dust_p_.nbin_ = std::stoi(input_.Get("number_of_bins"));
    dust_p_.radius_.resize(dust_p_.nbin_);
    dust_p_.cross_sections_.resize(dust_p_.nbin_);
    dust_p_.mass_.resize(dust_p_.nbin_);
    dust_p_.abundances_.resize(dust_p_.nbin_);

    // hall_beta
    hall_beta_.resize(number_of_species_, 0.0);
    sigmaO_i_.resize(number_of_species_, 0.0);
    sigmaH_i_.resize(number_of_species_, 0.0);
    sigmaP_i_.resize(number_of_species_, 0.0);

    return;

} // end ReadFile()

void Nicole::ReadTest()
{
    std::string fname = "read_test.txt";
    std::ofstream file(fname, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, fname);

    // int n = chem_p_.species_.size();
    // for (int i = 0; i < n; ++i) {
    //     file << chem_p_.species_[i] << std::endl;
    // }

    int n = chem_p_.type_of_reaction_.size();
    std::cout << n << std::endl;
    for (int i = 0; i < n; ++i) {
        file << chem_p_.type_of_reaction_[i] << " " << chem_p_.species_[chem_p_.index_reactant1_[i]] << " ";
        if (chem_p_.index_reactant2_[i] == NOT_FOUND_SPECIES_LIST) {
            file << "*" << " ";
        } else {
            file << chem_p_.species_[chem_p_.index_reactant2_[i]] << " ";
        }
        file << std::endl;
    }

    std::cout << std::scientific;
    for (int i = 0; i < dust_p_.nbin_; ++i) {
        std::cout << dust_p_.radius_[i] << " " << dust_p_.cross_sections_[i] << " " << dust_p_.mass_[i]
            << " " << dust_p_.abundances_[i] << std::endl;
    }

    file.close();
    return;
}

void Nicole::SetNumberDensity(double number_density)
{
    chem_p_.number_density_ = number_density;
}


void Nicole::SetTemperature(double temperature) 
{
    chem_p_.temperature_ = temperature;
}


void Nicole::SetIonnizationRate(double ionization_rate)
{
    chem_p_.ionization_rate_ = ionization_rate;
}

void Nicole::SetRelativeTolerance(double relative_tolerance)
{
    /* defalut value : 1.0e-5 */
    relative_tolerance_ = relative_tolerance;
}


void Nicole::SetAbsoluteTolerance(double absolute_tolerance)
{
    /* defalut value : 1.0e-15 */
    absolute_tolerance_ = absolute_tolerance;
}


UnsignedInt Nicole::GetNumberOfSpecies()
{
    return number_of_species_;
}


std::string Nicole::GetSpeciesName(UnsignedInt i) {

    return chem_p_.species_[i];

}


double Nicole::GetSpeciesAbundances(UnsignedInt i) {

    return x_species_[i];

}

double Nicole::GetElectronAbundance()
{
    return x_species_[electron_index_];
}

// ############################################################################################################################# //
// ############################################################################################################################# //

void Nicole::CalculateReactionRateCoefficient()
{
    UnsignedInt i;
    UnsignedInt number_of_reaction = chem_p_.reaction_rate_coefficient_.size();
    std::string reac_type;

    for (i = 0; i < number_of_reaction; ++i) {

        reac_type = chem_p_.reaction_rate_coefficient_[i];

        if (chem_p_.type_of_reaction_[i] == IONIZATION) {

            chem_p_.reaction_rate_coefficient_[i] = chem_p_.value_for_reaction1_[i] * chem_p_.ionization_rate_;

        } else if (chem_p_.type_of_reaction_[i] == COSMIC_RAY_PROTON) {

            // direct cosmic-ray proton ionization
            chem_p_.reaction_rate_coefficient_[i] = chem_p_.value_for_reaction1_[i];

        } else if (chem_p_.type_of_reaction_[i] == COSMIC_RAY_PHOTON) {

            // direct cosmic-ray photon ionization
            chem_p_.reaction_rate_coefficient_[i] = chem_p_.value_for_reaction1_[i] * std::pow(chem_p_.temperature_/300.0, chem_p_.value_for_reaction2_[i]) 
                    * std::exp(-chem_p_.value_for_reaction3_[i]/chem_p_.temperature_);

        } else if (chem_p_.type_of_reaction_[i] == ADSORPTION_NEUTRAL) {

            chem_p_.reaction_rate_coefficient_[i] = AdsorptionNeutralParticle(chem_p_.value_for_reaction1_[i]*BOLTZMANN_CONSTANT, 
                    chem_p_.value_for_reaction2_[i]*PROTON_MASS, chem_p_.temperature_);

        } else if (chem_p_.type_of_reaction_[i] == ADSOPTION_ION_COUNTER) {

            chem_p_.reaction_rate_coefficient_[i] = ChargeParticleGrainCollision(chem_p_.value_for_reaction1_[i]*PROTON_MASS, 1.0, 
                    chem_p_.value_for_reaction2_[i], (UnsignedInt)chem_p_.value_for_reaction3_[i], chem_p_.temperature_);
            if (chem_p_.value_for_reaction4_[i] != 0.0) {
                chem_p_.reaction_rate_coefficient_[i] /= chem_p_.value_for_reaction4_[i];
            }

        } else if (chem_p_.type_of_reaction_[i] == DESORPTION_NEUTRAL) {

            chem_p_.reaction_rate_coefficient_[i] = DesorptionNeutralParticle(chem_p_.value_for_reaction1_[i]*BOLTZMANN_CONSTANT, 
                    chem_p_.value_for_reaction2_[i]*PROTON_MASS, chem_p_.temperature_);

        } else if (chem_p_.type_of_reaction_[i] == ION_GRAIN_COLLISION_COUNTER) {

            chem_p_.reaction_rate_coefficient_[i] = ChargeParticleGrainCollision(chem_p_.value_for_reaction1_[i]*PROTON_MASS, 1.0, 
                    chem_p_.value_for_reaction2_[i], (UnsignedInt)chem_p_.value_for_reaction3_[i], chem_p_.temperature_);

        } else if (chem_p_.type_of_reaction_[i] == ION_GRAIN_COLLISION_NO_COUNTER) {

            chem_p_.reaction_rate_coefficient_[i] = ChargeParticleGrainCollision(chem_p_.value_for_reaction1_[i]*PROTON_MASS, 1.0, 
                    chem_p_.value_for_reaction2_[i], (UnsignedInt)chem_p_.value_for_reaction3_[i], chem_p_.temperature_);
            if (chem_p_.value_for_reaction4_[i] != 0.0) {
                chem_p_.reaction_rate_coefficient_[i] /= chem_p_.value_for_reaction4_[i];
            }

        } else if (chem_p_.type_of_reaction_[i] == ELECTRON_GRAIN_COLLISION) {

            chem_p_.reaction_rate_coefficient_[i] = ChargeParticleGrainCollision(ELECTRON_MASS, -1.0, 
                    chem_p_.value_for_reaction1_[i], (UnsignedInt)chem_p_.value_for_reaction2_[i], chem_p_.temperature_);

        } else if (chem_p_.type_of_reaction_[i] == GRAIN_GRAIN_COLLISION) {

            chem_p_.reaction_rate_coefficient_[i] = GrainGrainCollision(chem_p_.value_for_reaction1_[i], chem_p_.value_for_reaction2_[i], 
                    (UnsignedInt)chem_p_.value_for_reaction3_[i], (UnsignedInt)chem_p_.value_for_reaction4_[i], chem_p_.temperature_);

        } else if (chem_p_.type_of_reaction_[i] == PHOTOPROCESS) {

            chem_p_.reaction_rate_coefficient_[i] = chem_p_.value_for_reaction1_[i] * std::exp(-chem_p_.value_for_reaction3_[i] * VISUAL_EXTINCTION);

        } else if (chem_p_.type_of_reaction_[i] == H2_FORMATION_ON_GRAIN_SURFACE) {

            chem_p_.reaction_rate_coefficient_[i] = H2FormationOnGrainSurfaces(chem_p_.value_for_reaction1_[i]*BOLTZMANN_CONSTANT, 
                    chem_p_.value_for_reaction2_[i]*PROTON_MASS, chem_p_.temperature_);
        
        } else {

            chem_p_.reaction_rate_coefficient_[i] = GasPhaseReaction(chem_p_.value_for_reaction1_[i], chem_p_.value_for_reaction2_[i], 
                    chem_p_.value_for_reaction3_[i], chem_p_.value_for_reaction4_[i], chem_p_.value_for_reaction5_[i], chem_p_.temperature_);
            // std::cout << std::scientific << chem_p_.reaction_rate_coefficient_[i] << std::endl;
            // std::exit(1);
        }

    } // End for 

    flag_rate_coefficient_ = true;

} // end CalculateReactionRateCoefficient

void Nicole::CheckRate()
{
    CalculateReactionRateCoefficient();
    std::string filename = "checkrate.txt";
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, filename);

    UnsignedInt n = chem_p_.reaction_rate_coefficient_.size();
    file << std::scientific;
    for (UnsignedInt i = 0; i < n; ++i) {
        // file << chem_p_.type_of_reaction_[i] << " " << chem_p_.species_[chem_p_.index_reactant1_[i]] << " ";
        // if (chem_p_.index_reactant2_[i] == NOT_FOUND_SPECIES_LIST) {
        //     file << "*" << " ";
        // } else {
        //     file << chem_p_.species_[chem_p_.index_reactant2_[i]] << " ";
        // }
        // file << chem_p_.reaction_rate_coefficient_[i] << std::endl;
        file << chem_p_.reaction_rate_coefficient_[i] << std::endl;
    }

    file.close();
    return;
}

double Nicole::GasPhaseReaction(double alpha, double beta, double gamma, double lower_temperatue_limit, 
                double upper_temperature_limit, double temperature)
{
    double rate_coefficient;
    double temperatire_range_parameter = 0.50;

    if (temperature < lower_temperatue_limit) {
        if (temperature < (1.0 - temperatire_range_parameter)*lower_temperatue_limit) {
            rate_coefficient = 0.0;
            return rate_coefficient;
        } else {
            rate_coefficient = alpha * std::pow(lower_temperatue_limit/300.0, beta) * std::exp(-gamma/lower_temperatue_limit);
            return rate_coefficient;
        }
    } else if (lower_temperatue_limit <= temperature && temperature <= upper_temperature_limit) {
        rate_coefficient = alpha * std::pow(temperature/300.0, beta) * std::exp(-gamma/temperature);
        return rate_coefficient;
    } else {
        if (temperature <= (1.0 + temperatire_range_parameter)*upper_temperature_limit) {
            rate_coefficient = alpha * std::pow(upper_temperature_limit/300.0, beta) * std::exp(-gamma/upper_temperature_limit);
            return rate_coefficient;
        } else {
            rate_coefficient = 0.0;
            return rate_coefficient;
        }
    }
} // End GasPhaseReaction

double Nicole::AdsorptionNeutralParticle(double binding_energy, double species_mass, double temperarure)
{
    // Hollenback & Salpeter 1970, Ilgner & Nelson 2006, Bai & Goodman 2009
    double rate_coefficient;
    double dust_radius, thermal_velocity, dust_cross_section;
    double c, sticking_probability;
    double kd;
    UnsignedInt s;

    c = (1.0/2.0) * std::sqrt(2*binding_energy*TRANSFERRED_ENERGY / (BOLTZMANN_CONSTANT*temperarure));
    sticking_probability = 2.0 * (c*c) * (1 - c*std::sqrt(M_PI)*std::exp(c*c)*std::erfc(c));
    thermal_velocity = std::sqrt(8.0*BOLTZMANN_CONSTANT*temperarure / (M_PI*species_mass));

    kd = 0.0;
    for (s = 0; s < dust_p_.nbin_; ++s) {
        kd += dust_p_.cross_sections_[s] * dust_p_.abundances_[s];
    }

    // rate_coefficient = thermal_velocity * dust_cross_section * sticking_probability;

    rate_coefficient = thermal_velocity * sticking_probability * kd ;

    return rate_coefficient;

} // End AdsorptionNeutralParticle


double Nicole::DesorptionNeutralParticle(double binding_energy, double species_mass, double temperature)
{
    // Hasegawa et al. 1992, Ilgner & Nelson 2006, Bai & Goodman 2009
    double rate_coefficient, nu0;

    nu0 = std::sqrt(2.0*BINDING_SITES*binding_energy / (M_PI*M_PI*species_mass));

    rate_coefficient = nu0 * std::exp(-binding_energy / (BOLTZMANN_CONSTANT*temperature));
    return rate_coefficient;

} // End DesorptionNeutralParticle

double Nicole::ChargeParticleGrainCollision(double charge_particle_mass, double particle_charge, double dust_charge, 
                                                        UnsignedInt sth_bin, double temperature)
{
    // Drain & Sutin 1987
    double rate_coefficient;
    double thermal_velocity, dust_cross_section, reduced_temperature, theta_nu;
    double sticking_probability;

    thermal_velocity = std::sqrt(8.0*BOLTZMANN_CONSTANT*temperature / (M_PI*charge_particle_mass));
    dust_cross_section = dust_p_.cross_sections_[sth_bin-1];
    reduced_temperature = dust_p_.radius_[sth_bin-1]*BOLTZMANN_CONSTANT*chem_p_.temperature_ / (CHARGE_UNIT*CHARGE_UNIT);
    theta_nu = (dust_charge == 0 ? 0.0 : dust_charge / (1 + POLARIZABILITY*std::pow(std::abs(dust_charge), -0.50)));


    if (particle_charge == -1.0) {
        // electron-grain collision
        sticking_probability = 0.6;
        if (dust_charge == 0) {
            rate_coefficient = thermal_velocity * dust_cross_section * (1.0 + POLARIZABILITY*std::pow(M_PI / (2.0*reduced_temperature), 0.50)) * sticking_probability;
        } else if (dust_charge >= 1) {
            rate_coefficient = thermal_velocity * dust_cross_section * (1.0 + dust_charge/reduced_temperature) 
                     * (1.0 + POLARIZABILITY*std::pow(2.0 / (reduced_temperature + 2.0*dust_charge), 0.50)) * sticking_probability;
        } else if (dust_charge <= -1) {
            rate_coefficient = thermal_velocity * dust_cross_section * std::exp(theta_nu/reduced_temperature)
                     * std::pow((1.0 + POLARIZABILITY*std::pow(4.0*reduced_temperature - 3.0*dust_charge, -0.50)), 2.0) * sticking_probability;
        } else {
            rate_coefficient = 0.0;
        }
    } else {
        // ion-grain collision
        sticking_probability = 1.0;
        if (dust_charge == 0) {
            rate_coefficient = thermal_velocity * dust_cross_section * (1.0 + POLARIZABILITY*std::pow(M_PI / (2.0*reduced_temperature), 0.50)) * sticking_probability;
        } else if (dust_charge >= 1.0) {
            rate_coefficient = thermal_velocity * dust_cross_section * std::exp(-theta_nu/reduced_temperature)
                     * std::pow((1.0 + POLARIZABILITY*std::pow(4.0*reduced_temperature + 3.0*dust_charge, -0.50)), 2.0) * sticking_probability;
        } else if (dust_charge <= -1.0) {
            rate_coefficient = thermal_velocity * dust_cross_section * (1.0 - dust_charge/reduced_temperature)
                     * (1.0 + POLARIZABILITY*std::pow(2.0 / (reduced_temperature - 2.0*dust_charge), 0.50)) * sticking_probability;
        } else {
            rate_coefficient = 0.0;
        }
    }

    return rate_coefficient;

} // End ChargeParticleGrainCollision

double Nicole::GrainGrainCollision(double dust_charge_1, double dust_charge_2, UnsignedInt sth_bin_1, UnsignedInt sth_bin_2, double temperature)
{
    // Umebayashi & Nakano 1990 equation (3)
    double rate_coefficient;
    double dust_radius_1, dust_radius_2, dust_mass_1, dust_mass_2, reduced_mass, thermal_velocity;

    dust_radius_1 = dust_p_.radius_[sth_bin_1 - 1];
    dust_radius_2 = dust_p_.radius_[sth_bin_2 - 1];
    dust_mass_1 = dust_p_.mass_[sth_bin_1 - 1];
    dust_mass_2 = dust_p_.mass_[sth_bin_2 - 1];
    reduced_mass = dust_mass_1*dust_mass_2 / (dust_mass_1 + dust_mass_2);
    thermal_velocity = std::sqrt(8.0*BOLTZMANN_CONSTANT*temperature / (M_PI*reduced_mass));

    rate_coefficient = M_PI*std::pow(dust_radius_1 + dust_radius_2, 2) * thermal_velocity
             * (1.0 - dust_charge_1*dust_charge_2*CHARGE_UNIT*CHARGE_UNIT / ((dust_radius_1 + dust_radius_2)*BOLTZMANN_CONSTANT*temperature));
    
    return rate_coefficient;

} // End GrainGrainCollision


double Nicole::H2FormationOnGrainSurfaces(double binding_energy, double species_mass, double temperature)
{
    double rate_coefficient;
    double dust_cross_section, thermal_velocity;
    double eta, c, sticking_probability;
    double kd;
    UnsignedInt s;
    
    thermal_velocity  = std::sqrt(8.0*BOLTZMANN_CONSTANT*temperature / (M_PI*species_mass));

    c = (1.0/2.0) * std::sqrt(2*binding_energy*TRANSFERRED_ENERGY) / (BOLTZMANN_CONSTANT*temperature);
    sticking_probability = 2.0 * (c*c) * (1 - c*std::sqrt(M_PI)*std::exp(c*c)*std::erfc(c));

    eta = 0.50; // Draine 2011

    kd = 0.0;
    for (s = 0; s < dust_p_.nbin_; ++s) {
        kd += dust_p_.cross_sections_[s] * dust_p_.abundances_[s];
    }

    rate_coefficient = 0.50 * thermal_velocity * eta * sticking_probability * kd;

    // std::cout << sticking_probability << std::endl;

    return rate_coefficient;
}

// ############################################################################################################################# //
// ############################################################################################################################# //

void Nicole::SetInitCondition(DoubleVector &x)
{
    UnsignedInt i;
    UnsignedInt s = 0;
    UnsignedInt x_length = x.size();

    for (i = 0; i < x_length; ++i) {
        // std::cout << i << std::endl;
        if (chem_p_.species_abundances_.Contains(chem_p_.species_[i])) {
            x[i] = chem_p_.species_abundances_.Get_double(chem_p_.species_[i]);
        } else {
            x[i] = 0.0;
        }

        if ((chem_p_.species_[i][0] == 'G') && (chem_p_.species_charge_[i] == 0)) {
            x[i] = dust_p_.abundances_[s];
            s += 1;
        }
    }   

    return;
}

void Nicole::SetInitConditionTest()
{
    UnsignedInt i;
    UnsignedInt ns = number_of_species_;
    UnsignedInt s = 0;

    std::cout << std::scientific;
    for (i = 0; i < ns; ++i) {
        if (chem_p_.species_abundances_.Contains(chem_p_.species_[i])) {
            std::cout << "x[" << chem_p_.species_[i] << "] = " << chem_p_.species_abundances_.Get_double(chem_p_.species_[i]) << std::endl;
        }

        if ((chem_p_.species_[i][0] == 'G') && (chem_p_.species_charge_[i] == 0)) {
            std::cout << "x[" << chem_p_.species_[i] << "] = " << dust_p_.abundances_[s] << std::endl;
            s += 1;
        }
    }

    return;
}

int Nicole::Solve()
{
    UnsignedInt number_of_equation = number_of_species_;
    double t0, tout, tout_start, t_end, tstep;
    double reltol, abstol, xmin;
    int istate, itask, iopt, itol;
    int i, j, jac_type, nstep;

    DoubleVector x(number_of_equation);

    SetInitCondition(x);

    if (!flag_rate_coefficient_) {
        CalculateReactionRateCoefficient();
    }

    reltol = relative_tolerance_;
    abstol = absolute_tolerance_;
    DoubleVector rtol(number_of_equation, reltol);
    DoubleVector atol(number_of_equation, abstol);
    IntVector iwork(20, 0);
    DoubleVector rwork(20, 0.0);
    
    t0 = 0.0;
    tout_start = 1.0 * SOLAR_YEAR * 1.0e-3;
    t_end = integration_time_;
    nstep = 30;
    tstep = std::pow(10.0, std::log10(t_end/tout_start) / ((double)nstep - 1.0));
    tout = tout_start;

    xmin = 1.0e-30;

    if (flag_jacobian_) {
        jac_type = 1;
    } else {
        jac_type = 2;
    }

    istate = 1;
    itask  = 1;
    iopt   = 0;
    itol   = 2;

    Lsoda_cpp lsoda;

    lsoda.SetTolerance(atol, rtol);
    
    for (i = 0; i < nstep; ++i) {

        lsoda.Lsoda(DifferentialEquation, number_of_equation, x, &t0, tout, itol, itask, &istate, iopt, 
                iwork, rwork, Jacobian, jac_type, &chem_p_);

        for (j = 0; j < number_of_equation; ++j) {
            if (x[j] < xmin) {
                    x[j] = xmin;
            }
        }

        if (istate < 0) {
            istate = 2;
        }
        tout *= tstep;

    }

    end_time_ = tout/SOLAR_YEAR/tstep;
    for (i = 0; i < number_of_equation; ++i) {
        x_species_[i] = x[i];
    }

    number_of_steps_ = iwork[10];
    number_of_func_  = iwork[11];
    number_of_jac_   = iwork[12];

    return 0;

} // end Solve


int Nicole::SolveAndOutput(std::string output_file)
{
    UnsignedInt number_of_equation = number_of_species_;
    double t0, tout, tout_start, t_end, tstep;
    double reltol, abstol, xmin;
    int istate, itask, iopt, itol;
    int i, j, jac_type, nstep;

    DoubleVector x(number_of_equation);

    SetInitCondition(x);

    if (!flag_rate_coefficient_) {
        CalculateReactionRateCoefficient();
    }

    reltol = relative_tolerance_;
    abstol = absolute_tolerance_;
    DoubleVector rtol(number_of_equation, reltol);
    DoubleVector atol(number_of_equation, abstol);
    IntVector iwork(20, 0);
    DoubleVector rwork(20, 0.0);
    
    t0 = 0.0;
    tout_start = 1.0 * SOLAR_YEAR * 1.0e-3;
    t_end = integration_time_;
    nstep = 30;
    tstep = std::pow(10.0, std::log10(t_end/tout_start) / ((double)nstep - 1.0));
    tout = tout_start;

    xmin = 0.0; //1.0e-35;

    if (flag_jacobian_) {
        jac_type = 1;
    } else {
        jac_type = 2;
    }

    istate = 1;
    itask  = 1;
    iopt   = 0;
    itol   = 2;

    std::ofstream file(output_file, std::ios::out);
    FAILED_TO_OPEN(file, output_file);
    for (i = 0; i < number_of_equation; ++i) {
        file << chem_p_.species_[i] << " ";
    }
    file << std::endl;
    file << std::scientific;

    Lsoda_cpp lsoda;

    lsoda.SetTolerance(atol, rtol);


    for (i = 0; i < nstep; ++i) {

        lsoda.Lsoda(DifferentialEquation, number_of_equation, x, &t0, tout, itol, itask, &istate, iopt, 
                iwork, rwork, Jacobian, jac_type, &chem_p_);
        
        file << (tout / SOLAR_YEAR);
        for (j = 0; j < number_of_equation; ++j) {
            if (x[j] < xmin) {
                x[j] = xmin;
            }
            file << " " << x[j];
        } 
        file << std::endl;
        if (istate < 0) {
            istate = 2;
        }
        tout *= tstep;

    }

    file.close();

    end_time_ = tout/SOLAR_YEAR/tstep;
    for (i = 0; i < number_of_equation; ++i) {
        x_species_[i] = x[i];
    }

    number_of_steps_ = iwork[10];
    number_of_func_  = iwork[11];
    number_of_jac_   = iwork[12];

    return 0;

} // end SolveAndOutput

void Nicole::SetIntegrationTime(double integration_time)
{
    integration_time_ = integration_time;
}

double Nicole::GetIntegrationTime()
{
    return end_time_;
}

void Nicole::OutputFile(double t, DoubleVector &x, std::string output_filename)
{
    std::ofstream file(output_filename, std::ios::out | std::ios::app);
    FAILED_TO_OPEN(file, output_filename);

    UnsignedInt i;
    UnsignedInt x_length = x.size();

    file << std::scientific;
    file << (t / SOLAR_YEAR);
    for (i = 0; i < x_length; ++i) {
        file << " " << x[i];
    }
    file << std::endl;

    file.close();

} // end output

void Nicole::ShowSpeciesAbundances()
{
    int ns = number_of_species_;
    int i;

    std::cout <<"end time = " << end_time_ << " (soler year)" << std::endl;
    std::cout << std::scientific;
    for (i = 0; i < ns; ++i) {
        std::cout << "x[" << std::setw(3) << (i+1) << " : "<< std::setw(10) << chem_p_.species_[i] << "] = " << x_species_[i] << std::endl;
    }
    std::cout << std::endl;

    return;
}

void Nicole::FlagJacobian(bool flag_jacobian)
{
    flag_jacobian_ = flag_jacobian;
}

void Nicole::ShowStatisticalNumber()
{
    std::cout << "number of steps = " << number_of_steps_ << std::endl;
    std::cout << "number of func  = " << number_of_func_  << std::endl;
    std::cout << "number of jac   = " << number_of_jac_   << std::endl;

    return;
}

void Nicole::CheckChargeState()
{
    double total_dust_number_density_result = 0.0;
    double mean_dust_cahrge = 0.0;
    double error_dust = 0.0;
    double charge = 0.0;
    
    UnsignedInt N = chem_p_.index_dust_grains_.size();
    UnsignedInt n = number_of_species_;

    for (int i = 0; i < n; ++i) {
        charge += chem_p_.species_charge_[i] * x_species_[i];
    }

    for (UnsignedInt i = 0; i < N; ++i) {
        total_dust_number_density_result += x_species_[chem_p_.index_dust_grains_[i]];
        mean_dust_cahrge += x_species_[chem_p_.index_dust_grains_[i]] * chem_p_.species_charge_[chem_p_.index_dust_grains_[i]];
        // std::cout << chemical_species_[index_dust_grains_[i]] << std::endl;
    }
    if (total_dust_number_density_result == 0.0) {
        mean_dust_cahrge = 0.0;
    } else {
        mean_dust_cahrge /= total_dust_number_density_result;
    }
    total_dust_number_denisty_ = dust_p_.total_abundance_;
    error_dust = std::abs(total_dust_number_denisty_ - total_dust_number_density_result) / total_dust_number_density_result;

    std::cout << std::endl;
    std::cout <<  std::scientific;
    std::cout << "initial total dust number density = " << total_dust_number_denisty_ << std::endl;
    std::cout << "finish total dust number density  = " << total_dust_number_density_result << std::endl;
    std::cout << "error dust number density         = " << error_dust << std::endl;
    std::cout << "mean dust charge                  = " << mean_dust_cahrge << std::endl;
    std::cout << "a = " << charge << std::endl;

    return;
}

// ############################################################################################################################# //
// ############################################################################################################################# //

void Nicole::DifferentialEquation(double t, DoubleVector &x, DoubleVector &xdot, void *data)
{
    ChemicalParams* chem_p = (ChemicalParams*)data;
    std::string top;
    UnsignedInt ir1, ir2, ip1, ip2, ip3, ip4, i;
    UnsignedInt number_of_reactions = chem_p->reaction_rate_coefficient_.size();
    UnsignedInt number_of_dust = chem_p->index_dust_grains_.size();
    UnsignedInt number_of_species = chem_p->species_.size();
    double number_density = chem_p->number_density_;
    double k = 0.0;
    

    for (i = 0; i < number_of_species; ++i) xdot[i] = 0.0;

    for (i = 0; i < number_of_reactions; ++i) {

        top = chem_p->type_of_reaction_[i];
        ir1 = chem_p->index_reactant1_[i];
        ir2 = chem_p->index_reactant2_[i];
        ip1 = chem_p->index_product1_[i];
        ip2 = chem_p->index_product2_[i];
        ip3 = chem_p->index_product3_[i];
        ip4 = chem_p->index_product4_[i];

        if ((top == IONIZATION) || (top == COSMIC_RAY_PROTON) || (top == COSMIC_RAY_PHOTON) || (top == PHOTOPROCESS)) {

            k = chem_p->reaction_rate_coefficient_[i] * x[ir1];
            xdot[ir1] -= k;
            xdot[ip1] += k;
            xdot[ip2] += k;
            if (ip3 != NOT_FOUND_SPECIES_LIST) xdot[ip3] += k;
            if (ip4 != NOT_FOUND_SPECIES_LIST) xdot[ip4] += k;

        } else if (top == ADSORPTION_NEUTRAL) {

            k = chem_p->reaction_rate_coefficient_[i] * x[ir1] * number_density;
            xdot[ir1] -= k;
            xdot[ip1] += k;

        } else if (top == DESORPTION_NEUTRAL) {

            k = chem_p->reaction_rate_coefficient_[i] * x[ir1];
            xdot[ir1] -= k;
            xdot[ip1] += k;

        } else if (top == H2_FORMATION_ON_GRAIN_SURFACE) {

            k = chem_p->reaction_rate_coefficient_[i] * x[ir1] * number_density;
            xdot[ir1] -= k;
            xdot[ir2] -= k;
            xdot[ip1] += k;

        } else {

            k = chem_p->reaction_rate_coefficient_[i] * x[ir1] * x[ir2] * number_density;
            // if (k < 1.0e-35) k = 0.0;
            xdot[ir1] -= k;
            xdot[ir2] -= k;
            xdot[ip1] += k;
            if (ip2 != NOT_FOUND_SPECIES_LIST) xdot[ip2] += k;
            if (ip3 != NOT_FOUND_SPECIES_LIST) xdot[ip3] += k;
            if (ip4 != NOT_FOUND_SPECIES_LIST) xdot[ip4] += k;

        }

    }

    return; 

} // end DifferentialEquation

void Nicole::Jacobian(double t, DoubleVector &x ,DoubleMatrix &J, void *data)
{
    ChemicalParams* chem_p = (ChemicalParams *)data;
    std::string top;
    UnsignedInt ir1, ir2, ip1, ip2, ip3, ip4, i, j;
    UnsignedInt number_of_reactions = chem_p->reaction_rate_coefficient_.size();
    UnsignedInt number_of_dust = chem_p->index_dust_grains_.size();
    UnsignedInt number_of_species = chem_p->species_.size();
    double number_density = chem_p->number_density_;
    double k = 0.0;

    // for (i = 0; i < number_of_species; ++i) {
    //     for (j = 0; j < number_of_species; ++j) {
    //         J[i][j] = 0.0;
    //     }
    // }

    for (i = 0; i < number_of_reactions; ++i) {

        top = chem_p->type_of_reaction_[i];
        ir1 = chem_p->index_reactant1_[i];
        ir2 = chem_p->index_reactant2_[i];
        ip1 = chem_p->index_product1_[i];
        ip2 = chem_p->index_product2_[i];
        ip3 = chem_p->index_product3_[i];
        ip4 = chem_p->index_product4_[i];

        if ((top == IONIZATION) || (top == COSMIC_RAY_PROTON) || (top == COSMIC_RAY_PHOTON) || (top == PHOTOPROCESS)) {

            // reactant1 derivative
            k = chem_p->reaction_rate_coefficient_[i];
            J[ir1][ir1] -= k;
            J[ip1][ir1] += k;
            J[ip2][ir1] += k;
            if (ip3 != NOT_FOUND_SPECIES_LIST) J[ip3][ir1] += k;
            if (ip4 != NOT_FOUND_SPECIES_LIST) J[ip4][ir1] += k;

        } else if (top == ADSORPTION_NEUTRAL) {

            k = chem_p->reaction_rate_coefficient_[i] * number_density;
            J[ir1][ir1] -= k; 
            J[ip1][ir1] += k;

        } else if (top == DESORPTION_NEUTRAL) {

            // reactant1 derivative
            k = chem_p->reaction_rate_coefficient_[i];
            J[ir1][ir1] -=k;
            J[ip1][ir1] +=k;

        } else if (top == H2_FORMATION_ON_GRAIN_SURFACE) {

            k = chem_p->reaction_rate_coefficient_[i] * number_density;
            J[ir1][ir1] -= k;
            J[ir2][ir1] -= k;
            J[ip1][ir1] += k;

        } else {

            // reactant1 derivative
            k = chem_p->reaction_rate_coefficient_[i] * x[ir2] * number_density;
            // if (k < 1.0e-30) k = 0.0;
            J[ir1][ir1] -= k;
            J[ir2][ir1] -= k;
            J[ip1][ir1] += k;
            if (ip2 != NOT_FOUND_SPECIES_LIST) J[ip2][ir1] += k;
            if (ip3 != NOT_FOUND_SPECIES_LIST) J[ip3][ir1] += k;
            if (ip4 != NOT_FOUND_SPECIES_LIST) J[ip4][ir1] += k;

            // Reactant2 derivative
            k = chem_p->reaction_rate_coefficient_[i] * x[ir1] * number_density;
            // if (k < 1.0e-30) k = 0.0;
            J[ir1][ir2] -= k;
            J[ir2][ir2] -= k;
            J[ip1][ir2] += k;
            if (ip2 != NOT_FOUND_SPECIES_LIST) J[ip2][ir2] += k;
            if (ip3 != NOT_FOUND_SPECIES_LIST) J[ip3][ir2] += k;
            if (ip4 != NOT_FOUND_SPECIES_LIST) J[ip4][ir2] += k;

        }

    } // end for

    return;

} // end Jacobian

// ############################################################################################################################# //
// ############################################################################################################################# //

void Nicole::SetMagneticField(double magnetic_field)
{
    magnetic_field_ = magnetic_field;
};

void Nicole::FlagThermalIonization(bool flag_thermal_ionization)
{
    flag_thermal_ionization_ = flag_thermal_ionization;
}

UnsignedInt Nicole::GetDustBinNumber(std::string dust_grain)
{
    int spos = dust_grain.find("s");
    int cpos = dust_grain.find("c");

    std::string bin_number = dust_grain.substr(spos+1, cpos-spos-1);

    return std::stoi(bin_number);
}


void Nicole::Conductivity()
{
    UnsignedInt i, s;
    UnsignedInt ns = number_of_species_;
    double rhoH2, rhoHe, rhoH;
    double T, sqrtT, logT, mass_i;
    double tau, sigmavH2, sigmavHe, sigmavH, tauH2inv, tauHeinv, tauHinv;
    double sigmavH2_coe, sigmavH_coe, sigmavHe_coe, sigmavg_coe, delta_g;
    double omega_cyc, reduced_mass, ecnB, dust_mass;

    rhoH2 = massH2_ * x_species_[molc_hydogrn_index_] * chem_p_.number_density_;
    rhoHe = massHe_ * x_species_[helium_index_] * chem_p_.number_density_;
    rhoH = massH_ * x_species_[hydrogen_index_] * chem_p_.number_density_;
    // std::cout << "H2 = " << " " << rhoH2 << std::endl;
    // std::cout << "He = " << " " << rhoHe << std::endl;
    // std::cout << "H = " << " " << rhoH << std::endl;
    T = chem_p_.temperature_;
    sqrtT = std::sqrt(chem_p_.temperature_);
    logT = std::log10(chem_p_.temperature_);
    sigmavH2_coe = 2.81e-9 * std::sqrt(POLARIZABILITY_H2);
    sigmavHe_coe = 2.81e-9 * std::sqrt(POLARIZABILITY_He);
    sigmavH_coe = 2.81e-9 * std::sqrt(POLARIZABILITY_H);
    delta_g = 1.3;
    sigmavg_coe = delta_g * std::sqrt(128.0*BOLTZMANN_CONSTANT*T/(9.0*M_PI));
    sigmaO_ = sigmaH_ = sigmaP_ = 0.0;
    // hall_beta_.resize(ns, 0.0);

    // calculate hall parameter of each species
    for (i = 0; i < ns; ++i) {

        if (chem_p_.species_charge_[i] == 0.0) continue;
        
        mass_i = chem_p_.species_mass_[i];
        if (chem_p_.species_[i] == "e-") {

            // sigmav from Pinto & Galli (2008) table1
            // electron - H2
            sigmavH2 = 1.0e-9 * sqrtT * (0.535 + 0.203*logT - 0.136*logT*logT + 0.050*std::pow(logT, 3.0));
            tauH2inv = (sigmavH2 / (ELECTRON_MASS + massH2_)) * rhoH2;
            // electron - He
            sigmavHe = 1.0e-9 * 0.428 * sqrtT;
            tauHeinv = (sigmavHe / (ELECTRON_MASS + massHe_)) * rhoHe;
            // electron - H
            sigmavH = 1.0e-9 * sqrtT * (2.841 + 0.093*logT + 0.245*logT*logT - 0.089*std::pow(logT, 3.0));
            tauHinv = (sigmavH / (ELECTRON_MASS + massH_)) * rhoH;
            // tau, cyclotron freq
            tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
            omega_cyc = chem_p_.species_charge_[i] * CHARGE_UNIT * magnetic_field_ /(ELECTRON_MASS * SPEED_OF_LIGHT);
            // hall parameter
            hall_beta_[i] = tau * omega_cyc;

        } else if (chem_p_.species_[i] == "H3+") {

            // H3+ - H2 sigmav from Pinto & Galli (2008) table1
            sigmavH2 = 1.0e-9 * (2.693 - 1.238*logT + 0.664*logT*logT - 0.089*std::pow(logT, 3.0));
            tauH2inv = (sigmavH2 / (mass_i + massH2_)) * rhoH2;
            // H3+ - He sigmav from Pinto & Galli (2008) Appendix (A.5)
            reduced_mass = mass_i * massHe_ / (mass_i + massHe_);
            sigmavHe = sigmavHe_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauHeinv = (sigmavHe / (mass_i + massHe_)) * rhoHe;
            // H3+ - H He sigmav from Pinto & Galli (2008) Appendix (A.5)
            reduced_mass = mass_i * massH_ / (mass_i + massH_);
            sigmavH = sigmavH_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauHinv = (sigmavH / (mass_i + massH_)) * rhoH;
            // tau, cyclotron freq
            tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
            omega_cyc = chem_p_.species_charge_[i] * CHARGE_UNIT * magnetic_field_ /(mass_i * SPEED_OF_LIGHT);
            // hall parameter
            hall_beta_[i] = tau * omega_cyc;

        } else if (chem_p_.species_[i] == "H+") {

            // sigmav from Pinto & Galli (2008) table1
            // H+ - H2 
            sigmavH2 = 1.0e-9 * (1.003 + 0.050*logT + 0.136*logT*logT - 0.014*std::pow(logT, 3.0));
            tauH2inv = (sigmavH2 / (mass_i + massH2_)) * rhoH2;
            // H+ - He
            sigmavHe = 1.0e-9 * (1.424 + 7.438e-6*T - 6.734e-9*T*T);
            tauHeinv = (sigmavHe / (mass_i + massHe_)) * rhoHe;
            // H+ - H 
            sigmavH = 1.0e-9 * 0.649*std::pow(T, 0.375);
            tauHinv = (sigmavH / (mass_i + massH_)) * rhoH;
            // tau, cyclotron freq
            tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
            omega_cyc = chem_p_.species_charge_[i] * CHARGE_UNIT * magnetic_field_ /(mass_i * SPEED_OF_LIGHT);
            // hall parameter
            hall_beta_[i] = tau * omega_cyc;

        } else if (chem_p_.species_[i] == "HCO+") {

            // HCO+ - H2 sigmav from Pinto & Galli (2008) table1
            sigmavH2 = 1.0e-9 * sqrtT * (1.476 - 1.409*logT + 0.555*logT*logT - 0.0775*std::pow(logT, 3.0));
            tauH2inv = (sigmavH2 / (mass_i + massH2_)) * rhoH2;
            // HCO+ - He sigmav from Pinto & Galli (2008) Appendix (A.5)
            reduced_mass = mass_i * massHe_ / (mass_i + massHe_);
            sigmavHe = sigmavHe_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauHeinv = (sigmavHe / (mass_i + massHe_)) * rhoHe;
            // HCO+ - H sigmav from Pinto & Galli (2008) Appendix (A.5)
            reduced_mass = mass_i * massH_ / (mass_i + massH_);
            sigmavH = sigmavH_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauHinv = (sigmavH / (mass_i + massH_)) * rhoH;
            // tau, cyclotron freq
            tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
            omega_cyc = chem_p_.species_charge_[i] * CHARGE_UNIT * magnetic_field_ /(mass_i * SPEED_OF_LIGHT);
            // hall parameter
            hall_beta_[i] = tau * omega_cyc;

        } else if (chem_p_.species_[i][0] == 'G') {

            // Pinto & Galli (2008) eq.(25)
            s = GetDustBinNumber(chem_p_.species_[i]);
            mass_i = dust_p_.mass_[s-1];
            // collision with H2
            sigmavH2 = dust_p_.cross_sections_[s-1] * sigmavg_coe / std::sqrt(massH2_);
            tauH2inv = (sigmavH2 / (mass_i + massH2_)) * rhoH2;
            // collision with He
            sigmavHe = dust_p_.cross_sections_[s-1] * sigmavg_coe / std::sqrt(massHe_);
            tauHeinv = (sigmavH2 / (mass_i + massHe_)) * rhoHe;
            // collision with H
            sigmavH = dust_p_.cross_sections_[s-1] * sigmavg_coe / std::sqrt(massH_);
            tauHinv = (sigmavH2 / (mass_i + massH_)) * rhoH;
            // tau, cyclotron freq
            tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
            omega_cyc = chem_p_.species_charge_[i] * CHARGE_UNIT * magnetic_field_ /(mass_i * SPEED_OF_LIGHT);
            // hall parameter
            hall_beta_[i] = tau * omega_cyc;
                             
        } else {

            // other ions sigmav from Pinto & Galli (2008) Appendix (A.5)
            // collision with
            reduced_mass = mass_i * massH2_ / (mass_i + massH2_);
            sigmavH2 = sigmavH2_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauH2inv = (sigmavH2 / (mass_i + massH2_)) * rhoH2;
            // collision with He
            reduced_mass = mass_i * massHe_ / (mass_i + massHe_);
            sigmavHe = sigmavHe_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauHeinv = (sigmavHe / (mass_i + massHe_)) * rhoHe;
            // collision with H
            reduced_mass = mass_i * massH_ / (mass_i + massH_);
            sigmavH = sigmavH_coe * std::sqrt(chem_p_.species_charge_[i]) * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
            tauHinv = (sigmavH / (mass_i + massH_)) * rhoH;
            // tau, cyclotron freq
            tau = 1.0 / (tauH2inv + tauHeinv + tauHinv);
            omega_cyc = chem_p_.species_charge_[i] * CHARGE_UNIT * magnetic_field_ /(mass_i * SPEED_OF_LIGHT);
            // hall parameter
            hall_beta_[i] = tau * omega_cyc;
                                            
        }
    } // end for 

    // calculate conductivity
    ecnB = CHARGE_UNIT * SPEED_OF_LIGHT * chem_p_.number_density_ / magnetic_field_;

    if (flag_thermal_ionization_) {
        ThermalIonization();
    }

    for (i = 0; i < ns; ++i) {
        sigmaO_ += x_species_[i] * chem_p_.species_charge_[i] * hall_beta_[i];
        sigmaO_i_[i] = x_species_[i] * chem_p_.species_charge_[i] * hall_beta_[i];
        sigmaH_ -= x_species_[i] * chem_p_.species_charge_[i] * std::pow(hall_beta_[i], 2.0) / (1.0 + std::pow(hall_beta_[i], 2.0));
        sigmaH_i_[i] = -x_species_[i] * chem_p_.species_charge_[i] * std::pow(hall_beta_[i], 2.0) / (1.0 + std::pow(hall_beta_[i], 2.0));
        sigmaP_ += x_species_[i] * chem_p_.species_charge_[i] * hall_beta_[i] / (1.0 + std::pow(hall_beta_[i], 2.0));
        sigmaP_i_[i] = x_species_[i] * chem_p_.species_charge_[i] * hall_beta_[i] / (1.0 + std::pow(hall_beta_[i], 2.0));
    }

    sigmaO_ *= ecnB;
    sigmaH_ *= ecnB;
    sigmaP_ *= ecnB;

    return;
}

void Nicole::Resistivity()
{
    double sigma_perp;
    Conductivity();
    sigma_perp = sigmaH_*sigmaH_ + sigmaP_*sigmaP_;
    etaO_ = (SPEED_OF_LIGHT*SPEED_OF_LIGHT/ (4.0*M_PI)) * (1.0/sigmaO_);
    etaH_ = (SPEED_OF_LIGHT*SPEED_OF_LIGHT/ (4.0*M_PI)) * (sigmaH_ / sigma_perp);
    etaA_ = (SPEED_OF_LIGHT*SPEED_OF_LIGHT/ (4.0*M_PI)) * (sigmaP_ / sigma_perp) - etaO_;
}

void Nicole::ThermalIonization()
{
    double xe_th, rhoH2, rhoHe;
    double sigmavH2_coe, sigmavHe_coe, sigmavH2, sigmavHe;
    double tau, tauH2inv, tauHeinv, omega_cyc;
    double potassium_charge, potassium_sigmav, potassium_hall_beta, reduced_mass;
    rhoH2 = 0.5*2.0*HYDROGEN_MASS*chem_p_.number_density_;
    rhoHe = 0.095*4.0*HYDROGEN_MASS*chem_p_.number_density_;
    sigmavH2_coe = 2.81e-9 * std::sqrt(POLARIZABILITY_H2);
    sigmavHe_coe = 2.81e-9 * std::sqrt(POLARIZABILITY_He);

    xe_th = 1.0e-12 * pow(potassium_abundances_/1.0e-7, 0.50) * pow(chem_p_.number_density_/1.0e15, -0.50) 
        * pow(chem_p_.temperature_/1.0e3, 3.0/4.0) * exp(-2.52e4/chem_p_.temperature_)/(1.14e-11);
    x_species_[electron_index_] += xe_th;

    potassium_charge = 1.0;
    reduced_mass = POTASSIUM_MASS*GAS_MOLECULAR_MASS / (POTASSIUM_MASS + GAS_MOLECULAR_MASS);
    // collision with H2
    sigmavH2 = sigmavH2_coe * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
    tauH2inv = (sigmavH2 / (POTASSIUM_MASS + massH2_)) * rhoH2;
    // collision with He
    sigmavHe = sigmavHe_coe * std::pow(reduced_mass/HYDROGEN_MASS, -0.5);
    tauH2inv = (sigmavHe / (POTASSIUM_MASS + massHe_)) * rhoHe;
    // tau
    tau = 1.0 / (tauH2inv + tauHeinv);
    omega_cyc = potassium_charge* CHARGE_UNIT * magnetic_field_ /(POTASSIUM_MASS * SPEED_OF_LIGHT);
    // hall beta
    potassium_hall_beta = tau * omega_cyc;

    sigmaO_ += xe_th * potassium_charge * potassium_hall_beta;
    sigmaH_ += xe_th * potassium_charge * std::pow(potassium_hall_beta, 2.0) / (1.0 + std::pow(potassium_hall_beta, 2.0));
    sigmaP_ += xe_th * potassium_charge * potassium_hall_beta / (1.0 + std::pow(potassium_hall_beta, 2.0));

    return;
}


void Nicole::GetResistivity(double &etaO, double &etaH, double &etaA)
{
    etaO = etaO_;
    etaH = etaH_;
    etaA = etaA_;
}

void Nicole::GetConductivity(double &sigmaO, double &sigmaH, double &sigmaP)
{
    sigmaO = sigmaO_;
    sigmaH = sigmaH_;
    sigmaP = sigmaP_;
}

double Nicole::GetHallParameter(int i)
{
    return hall_beta_[i];
}
void Nicole::GetConductivitySpecies(int i, double &sigmaOi, double &sigmaHi, double &sigmaPi) 
{
    sigmaOi = sigmaO_i_[i];
    sigmaHi = sigmaH_i_[i];
    sigmaPi = sigmaP_i_[i];
}

double Nicole::GetOhmicConductivity(int i)
{
    return sigmaO_i_[i];
}

double Nicole::GetHallConductivity(int i)
{
    return sigmaH_i_[i];
}

double Nicole::GetPedersenConductivity(int i)
{
    return sigmaP_i_[i];
}
