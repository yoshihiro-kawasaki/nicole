#include <fstream>
#include <iostream>
#include <vector>

#include "nicole/nicole.hpp"
#include "../func.hpp"


int main() {
    // ユーザー設定の化学種
    std::vector<std::string> user_gas_species_list = {
        "H", "H2", "He", "CO", "O2", "Mg", "O", "C", "HCO", "H2O", "OH", 
        "H3+", "H2+", "H+", "HCO+", "Mg+", "He+", "C+", "O+", "O2+", "H3O+", 
        "OH+","H2O+", "e-"
    };

    // configオブジェクトの生成
    nicole::InputConfig config("input.txt");
    
    // 元素に関する管理クラスのオブジェクト生成
    nicole::ElementManager element_manager(config);

    // 化学種に関する管理クラスのオブジェクト生成
    // nicole::SpeciesManager species_manager(&element_manager, input);
    nicole::SpeciesManager species_manager(&element_manager, config, user_gas_species_list);
    std::string check_species_filename = "check_species.txt";
    species_manager.CheckSpeciesManager(check_species_filename);

    // 反応に関する管理クラスのオブジェクト生成
    nicole::ReactionManager reaction_manager(&species_manager, config);
    std::string check_reaction_file = "check_reaction_file.txt";
    reaction_manager.CheckReactionManager(check_reaction_file);

    // 化学反応計算に必要な物理パラメータの設定
    nicole::EnvironmentParameters environment_parameters;
    environment_parameters.gas_number_density = 2.0e5;
    environment_parameters.gas_temperature = 10.0;
    environment_parameters.cosmic_ray_ionization_rate = 1.0e-17;
    environment_parameters.x_rays_ionization_rate = 0.0;
    environment_parameters.visual_extinction = 15.0;
    environment_parameters.scaling_factor_uv_field = 0.0;

    // 化学反応計算のクラスのオブジェクト生成
    nicole::ReactionSimulator reaction_simulator(&species_manager, &reaction_manager, &environment_parameters, config);
    std::string check_rate_file = "check_rate_file.txt";
    reaction_simulator.CheckReactionRateCoefficient(check_rate_file);

    // 非理想磁気流体の磁気抵抗率の計算のクラスのオブジェクト生成
    nicole::NonIdealMHDeffect non_ideal_mhd_effect(&species_manager, &environment_parameters);

    // output fileの設定
    const std::string output_file = "test_resistivity.txt";
    std::ofstream file(output_file, std::ios::out | std::ios::trunc);
    if (!file.is_open()) {
        return 0;
    }

    // ファイルの１行目にラベル(時間、化学種)の記述
    const std::size_t number_of_species = species_manager.GetTotalNumberOfSpecies();
    file << std::setw(14) << "ng" << " "
         << std::setw(14) << "T"  << " "
         << std::setw(14) << "B"  << " ";
    for (int i = 0; i < number_of_species; ++i) {
        file << std::setw(14) << species_manager.GetSpeciesName(i) << " ";
    }
    file << std::setw(14) << "etaO" << " "
         << std::setw(14) << "etaH"  << " "
         << std::setw(14) << "etaA"  << std::endl;

    const double gas_number_density_start = 1.0e4;
    const double gas_number_density_end = 1.0e15;
    const int nn = 200;
    const double dlogn = (std::log10(gas_number_density_end) - std::log10(gas_number_density_start)) / static_cast<double>(nn - 1);

    for (int i = 0; i < nn; ++i) {
        // set number density and temperature
        double gas_number_density = std::pow(10.0, std::log10(gas_number_density_start) + dlogn * static_cast<double>(i));
        double gas_mass_density = 1.4 * nicole::constants::kProtonMass * gas_number_density;
        double gas_temperature = BarotropicEOS(gas_mass_density);
        double ionization_rate = IonizationRate(gas_mass_density, gas_temperature);

        environment_parameters.gas_number_density = gas_number_density;
        environment_parameters.gas_temperature = gas_temperature;
        environment_parameters.cosmic_ray_ionization_rate = ionization_rate;

        // 化学種の存在量を格納する配列の初期化
        std::vector<double> species_abundances(number_of_species);
        reaction_simulator.SetInitialSpeciesAbundances(species_abundances);

        // 積分時間に関する設定
        double t = 0.0;
        double tend = FreeFalltime(gas_mass_density) * 3.0;
        double tout = 1.0e-1 * nicole::constants::kSolarYear;
        if (tout > tend) tout = 1.0e-2 * tend;
        int nstep = 100;
        const double tstep = std::pow(10.0, std::log10(tend/tout)/static_cast<double>(nstep - 1));

        reaction_simulator.CalculateRateCoefficient();
        // 積分の実行
        for (int i = 1; i <= nstep; ++i) {
            bool success = reaction_simulator.Integrate(t, tout, species_abundances);
            tout *= tstep;
            if (!success) break;
        }

        // 磁場の設定
        const double magnetic_field = MagneticField(gas_number_density);
        environment_parameters.magnetic_field = magnetic_field;

        // 非理想磁気流体の磁気抵抗率の計算
        non_ideal_mhd_effect.CalculateHallParameters(species_abundances);
        non_ideal_mhd_effect.CalculateConductivities(species_abundances);
        non_ideal_mhd_effect.CalculateResistivities();

        // 磁気抵抗率の取得
        double etaO, etaH, etaA;
        non_ideal_mhd_effect.GetResistivite(etaO, etaH, etaA);

        file << std::scientific;
        file << std::setw(14) << gas_number_density << " " 
             << std::setw(14) << gas_temperature << " "
             << std::setw(14) << magnetic_field << " ";
        for (int i = 0; i < number_of_species; ++i) {
            file << std::setw(14) << species_abundances[i] << " ";
        }
        file << std::setw(14) << etaO << " "
             << std::setw(14) << etaH << " "
             << std::setw(14) << etaA << std::endl;
    }
    file.close();

    return 0;
}