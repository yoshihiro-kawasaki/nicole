#include "nicole.hpp"

#include <iostream>
#include <fstream>

void Test1();
void Test2();

double TempBarotropicRho(double rho)
{
    // Tsukamoto et al. 2020
    const double gamma = 7.0/5.0;
    const double rho_crit = 4.0e-14;
    double T;

    T = 10 * (1 + gamma*std::pow(rho/rho_crit, gamma - 1));

    return T;
}

int main()
{
    // Test1();
    Test2();
    return 0;
}

void Test1()
{
    Nicole nicole;
    nicole.ReadFile();
    double rho = 1.0e-15;
    double nH = rho / (1.4*HYDROGEN_MASS);
    double T = TempBarotropicRho(rho); // 10.0;
    double zeta = 1.0e-17;
    double B = 30.0e-3;
    nicole.SetNumberDensity(nH);
    nicole.SetTemperature(T);
    nicole.SetIonnizationRate(zeta);

    std::cout << std::scientific;
    std::cout << "nH   = " << nH << std::endl;
    std::cout << "T    = " << T << std::endl;
    std::cout << "zeta = " << zeta << std::endl;
    std::cout << "B    = " << B << std::endl;

    double q, rhod, amin, amax, fdg, d_size;
    q = 3.5;
    rhod = 3.0;
    amin = 5.0e-7;
    amax = 2.5e-5;
    fdg = 0.01;
    d_size = 1.0e-5;
    // nicole.SetSizeDistributionParameter(q, amin, amax, fdg, rhod);
    // nicole.SizeDistributionModel();
    nicole.SetSingleSizeParameter(d_size, fdg, rhod);
    nicole.SingleSizeModel();
    // nicole.AggregateModel(nH, fdg, 100);

    bool flag_jacobian = true;
    double reltol = 1.0e-5;
    double abstol = 1.0e-15;
    nicole.FlagJacobian(flag_jacobian);
    nicole.SetRelativeTolerance(reltol);
    nicole.SetAbsoluteTolerance(abstol);

    nicole.CalculateReactionRateCoefficient();
    double tinteg = 1.0e6 * SOLAR_YEAR;
    nicole.SetIntegrationTime(tinteg);
    nicole.Solve();
    nicole.ShowSpeciesAbundances();

    double etaO, etaH, etaA;
    nicole.SetMagneticField(B);
    nicole.Resistivity();
    nicole.GetResistivity(etaO, etaH, etaA);

    nicole.CheckChargeState();
    nicole.ShowStatisticalNumber();

    std::cout << std::scientific;
    std::cout << "etaO = " << etaO << std::endl;
    std::cout << "etaH = " << etaH << std::endl;
    std::cout << "etaA = " << etaA << std::endl;

    return;
}

void Test2()
{

    double rho, rhoi,rhof, dlnrho, nH, T, zeta, B;
    int Nrho;
    Nicole nicole;
    nicole.ReadFile();

    rhoi = 1.0e-18;
    rhof = 1.0e-9;
    Nrho = 100;
    dlnrho = (std::log10(rhof) - std::log10(rhoi)) / double(Nrho-1);
    B = 30.0e-3;
    zeta = 1.0e-17;

    double q, rhod, amin, amax, fdg, d_size;
    q = 3.5;
    rhod = 2.0;
    // amin = 5.0e-7;
    amin = 1.0e-5;
    amax = 2.5e-5;
    fdg = 0.01;
    d_size = 1.0e-5;
    // d_size = 0.035e-4;
    // nicole.SetSizeDistributionParameter(q, amin, amax, fdg, rhod);
    // nicole.SizeDistributionModel();
    nicole.SetSingleSizeParameter(d_size, fdg, rhod);
    nicole.SingleSizeModel();
    // nicole.AggregateModel(nH, fdg, 100);

    std::string filename = "test01.txt";
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, filename);

    bool flag_jacobian = true;
    double reltol = 1.0e-7;
    double abstol = 1.0e-15;
    double etaO, etaH, etaA;
    double sigmaO, sigmaH, sigmaP;
    nicole.FlagJacobian(flag_jacobian);
    nicole.SetRelativeTolerance(reltol);
    nicole.SetAbsoluteTolerance(abstol);

    nicole.SetIonnizationRate(zeta);
    nicole.SetMagneticField(B);

    int ns = nicole.GetNumberOfSpecies();
    int j;
    double tinteg;
    file << std::scientific;
    for (int i = 0; i < Nrho; ++i) {
        rho = std::pow(10.0, std::log10(rhoi)+dlnrho*i);
        nH = rho / (1.4*HYDROGEN_MASS);
        T = TempBarotropicRho(rho);
        tinteg = std::sqrt(3.0*M_PI/(32.0*GRAVITATIONAL_CONSTANT*rho)) * 10.0;
        nicole.SetNumberDensity(nH);
        nicole.SetTemperature(T);
        nicole.SetIonnizationRate(zeta);
        nicole.SetIntegrationTime(tinteg);
        // nicole.AggregateModel(nH, fdg, 1000);

        nicole.CalculateReactionRateCoefficient();
        nicole.Solve();
        nicole.Resistivity();
        nicole.GetResistivity(etaO, etaH, etaA);

        file << rho << " " << nH << " " << T << " " << zeta << " " << B << " ";
        for (j = 0; j < ns; j++) {
            file << nicole.GetSpeciesAbundances(j) << " ";
        }
        for (j = 0; j < ns; j++) {
            file << nicole.GetOhmicConductivity(j) << " ";
        }
        for (j = 0; j < ns; j++) {
            file << nicole.GetHallConductivity(j) << " ";
        }
        for (j = 0; j < ns; j++) {
            file << nicole.GetPedersenConductivity(j) << " ";
        }
        nicole.GetConductivity(sigmaO, sigmaH, sigmaP);
        file << sigmaO << " " << sigmaH << " " << sigmaP << " " << etaO << " " << etaH << " " << etaA << std::endl;
    }

    file.close();

    return;
}