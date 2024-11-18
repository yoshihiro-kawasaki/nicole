#ifndef NICOLE_HPP_
#define NICOLE_HPP_

#include "typedefs.hpp"
#include "physical_constants.hpp"
#include "dust.hpp"
#include "lsoda.hpp"

class Nicole
    : public Dust
{
public:
    Nicole();
    ~Nicole();

    void ReadFile();

    void SetNumberDensity(double number_density);
    void SetTemperature(double temperature);
    void SetIonnizationRate(double ionization_rate);
    
    void SetRelativeTolerance(double ralative_tolerance);
    void SetAbsoluteTolerance(double absolute_tolerance);

    UnsignedInt GetNumberOfSpecies();
    std::string GetSpeciesName(UnsignedInt i);
    double GetSpeciesAbundances(UnsignedInt i);
    double GetElectronAbundance();

    void CalculateReactionRateCoefficient();

    int Solve();
    int SolveAndOutput(std::string output_file);
    void SetIntegrationTime(double integration_time);
    double GetIntegrationTime();
    void ShowSpeciesAbundances();
    void FlagJacobian(bool flag_jacobian);
    void ShowStatisticalNumber();
    void CheckChargeState();

    // conductivity and resistivity
    void SetMagneticField(double magnetic_field);
    void FlagThermalIonization(bool flag_thermal_ionization);
    void Conductivity();
    void Resistivity();
    void GetResistivity(double &etaO, double &etaH, double &etaA);
    void GetConductivity(double &sigmaO, double &sigmaH, double &sigmaP);
    double GetHallParameter(int i);
    void GetConductivitySpecies(int i, double &sigmaOi, double &sigmaHi, double &sigmaPi);
    double GetOhmicConductivity(int i);
    double GetHallConductivity(int i);
    double GetPedersenConductivity(int i);

private:

    UnsignedInt IndexSpecies(const std::string chemical_species_name, const StringVector chemical_species_list);
    void ReadTest();

    double GasPhaseReaction(double alpha, double beta, double gamma, double lower_temperatue_limit, double upper_temperature_limit, double temperature);
    double AdsorptionNeutralParticle(double binding_energy, double species_mass, double temperarure);
    double DesorptionNeutralParticle(double binding_energy, double species_mass, double temperature);
    double ChargeParticleGrainCollision(double charge_particle_mass, double particle_charge, double dust_charge, UnsignedInt sth_bin, double temperature);
    double GrainGrainCollision(double dust_charge_1, double dust_charge_2, UnsignedInt sth_bin_1, UnsignedInt sth_bin_2, double temperature);
    double H2FormationOnGrainSurfaces(double binding_energy, double species_mass, double temperature);
    void CheckRate();

    static void DifferentialEquation(double t, DoubleVector &x, DoubleVector &xdot, void *data);
    static void Jacobian(double t, DoubleVector &y ,DoubleMatrix &J, void *data);

    void OutputFile(double t, DoubleVector &x, std::string output_filename);
    void SetInitCondition(DoubleVector &x);
    void SetInitConditionTest();

    UnsignedInt GetDustBinNumber(std::string dust_grain);
    void ThermalIonization();

    struct ChemicalParams
    {
        double            number_density_;
        double            temperature_;
        double            ionization_rate_;
        StringVector      species_;
        Dictionary        species_abundances_;
        DoubleVector      species_mass_;
        DoubleVector      species_charge_;
        StringVector      type_of_reaction_;
        UnsignedIntVector index_reactant1_;
        UnsignedIntVector index_reactant2_;
        UnsignedIntVector index_product1_;
        UnsignedIntVector index_product2_;
        UnsignedIntVector index_product3_;
        UnsignedIntVector index_product4_;
        UnsignedIntVector index_dust_grains_;
        DoubleVector      value_for_reaction1_;
        DoubleVector      value_for_reaction2_;
        DoubleVector      value_for_reaction3_;
        DoubleVector      value_for_reaction4_;
        DoubleVector      value_for_reaction5_;
        DoubleVector      reaction_rate_coefficient_;
    } chem_p_;

    DoubleVector x_species_;

    InputConfigure input_;

    UnsignedInt number_of_species_;
    UnsignedInt number_of_reactions_;
    UnsignedInt number_of_dust_grain_species_;
    UnsignedInt electron_index_;                  // e-
    UnsignedInt molc_hydogrn_index_;              // H2
    UnsignedInt helium_index_;                    // He
    UnsignedInt hydrogen_index_;                  // H

    double massH2_;
    double massHe_;
    double massH_;

    double      dust_internal_density_;
    double      dust_to_gas_mass_ratio_;
    double      total_dust_number_denisty_;

    double      relative_tolerance_;
    double      absolute_tolerance_;
    double      integration_time_;

    bool        flag_dust_evapolation_;
    bool        flag_jacobian_;
    bool        flag_error_function_;

    bool        flag_rate_coefficient_;

    int number_of_steps_;
    int number_of_func_;
    int number_of_jac_;

    double end_time_;

    DoubleVector hall_beta_;
    DoubleVector sigmaO_i_, sigmaH_i_, sigmaP_i_;
    double magnetic_field_;
    double sigmaO_, sigmaH_, sigmaP_;
    double etaO_, etaH_, etaA_;
    bool flag_thermal_ionization_;
    double potassium_abundances_;
    
};

#endif // NICOLE_HPP_