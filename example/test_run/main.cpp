#include <fstream>
#include <iostream>
#include <vector>

#include "nicole/nicole.hpp"
#include "../func.hpp"


int main() {
    // ユーザー設定の化学種
    std::vector<std::string> user_gas_species_list = {
        "H", "H2", "He", "CO", "O2", "Mg", "O", "C", "HCO", "H2O", 
        "OH", "H3+", "H2+", "H+", "HCO+", "Mg+", "He+", "C+", "O+", 
        "CO+", "CH2+", "O2+", "H3O+", "OH+", "H2O+", "e-"
    };

    // configオブジェクトの生成
    nicole::InputConfig config("input.txt");
    
    // 元素に関する管理クラスのオブジェクト生成
    nicole::ElementManager element_manager(config);

    // 化学種に関する管理クラスのオブジェクト生成
    // nicole::SpeciesManager species_manager(&element_manager, config);
    nicole::SpeciesManager species_manager(&element_manager, config, user_gas_species_list);
    std::string check_species_filename = "check_species.txt";
    species_manager.CheckSpeciesManager(check_species_filename);

    // 反応に関する管理クラスのオブジェクト生成
    nicole::ReactionManager reaction_manager(&species_manager, config);
    std::string check_reaction_file = "check_reaction_file.txt";
    reaction_manager.CheckReactionManager(check_reaction_file);

    // 化学反応計算に必要な物理パラメータの設定
    nicole::EnvironmentParameters environment_parameters;
    environment_parameters.gas_number_density = 1.0e6;
    const double rhog = 1.4 * nicole::constants::kProtonMass * environment_parameters.gas_number_density;
    environment_parameters.gas_temperature = BarotropicEOS(rhog); // 10.0;
    environment_parameters.cosmic_ray_ionization_rate = 1.3e-17;
    environment_parameters.x_rays_ionization_rate = 0.0;
    environment_parameters.visual_extinction = 10.0;
    environment_parameters.scaling_factor_uv_field = 1.0;

    std::cout << std::scientific << "T = " << environment_parameters.gas_temperature << std::endl;

    // 化学反応計算のクラスのオブジェクト生成
    nicole::ReactionSimulator reaction_simulator(&species_manager, &reaction_manager, &environment_parameters, config);
    std::string check_rate_file = "check_rate_file.txt";
    reaction_simulator.CheckReactionRateCoefficient(check_rate_file);

    // 化学種の存在量を格納する配列の初期化
    const std::size_t number_of_species = species_manager.GetTotalNumberOfSpecies();
    std::vector<double> species_abundances(number_of_species);
    reaction_simulator.SetInitialSpeciesAbundances(species_abundances);

    // 積分時間に関する設定
    double t = 0.0;
    double tout = 1.0e-1 * nicole::constants::kSolarYear;
    const double tend = FreeFalltime(rhog) * 10.0; // 1.0e6 * constants::kSolarYear;
    std::cout << "tend = " << (tend / nicole::constants::kSolarYear) << std::endl;
    if (tout > tend) tout = 1.0e-2 * tend;
    const int nstep = 100;
    const double tstep = std::pow(10.0, std::log10(tend/tout)/static_cast<double>(nstep - 1));

    // outputの設定
    const std::string output_file = "test.txt";
    std::ofstream file(output_file, std::ios::out | std::ios::trunc);
    if (!file.is_open()) {
        return 0;
    }

    // ファイルの１行目にラベル(時間、化学種)の記述
    file << std::setw(12) << "time" << " ";
    for (int i = 0; i < number_of_species; ++i) {
        file << std::setw(12) << species_manager.GetSpeciesName(i) << " ";
    }
    file << std::endl;

    // 初期存在量の記述
    file << std::scientific << std::setw(12) << t << " ";
    for (int i = 0; i < number_of_species; ++i) {
        file << std::setw(12) << species_abundances[i] << " ";
    }
    file << std::endl;

    // 積分の実行
    for (int i = 1; i <= nstep; ++i) {
        bool success = reaction_simulator.Integrate(t, tout, species_abundances, file);
        tout *= tstep;
        if (!success) break;
    }
    // 計算結果の確認
    reaction_simulator.CheckCalculationResult(species_abundances);

    file.close();

    // 磁場の設定
    environment_parameters.magnetic_field = MagneticField(environment_parameters.gas_number_density); // 1.0e-5;

    // 非理想磁気流体の磁気抵抗率の計算のクラスのオブジェクト生成
    nicole::NonIdealMHDeffect non_ideal_mhd_effect(&species_manager, &environment_parameters);
    non_ideal_mhd_effect.CalculateHallParameters(species_abundances);
    non_ideal_mhd_effect.CalculateConductivities(species_abundances);
    non_ideal_mhd_effect.CalculateResistivities();

    // 磁気抵抗率の取得
    double etaO, etaH, etaA;
    non_ideal_mhd_effect.GetResistivite(etaO, etaH, etaA);
    std::cout << "etaO = " << std::scientific << etaO << std::endl;
    std::cout << "etaH = " << std::scientific << etaH << std::endl;
    std::cout << "etaA = " << std::scientific << etaA << std::endl;

    return 0;
}
