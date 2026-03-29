#ifndef FUNC_HPP_
#define FUNC_HPP_

#include <cmath>

#include "nicole/nicole_defs.hpp"


double FreeFalltime(const double rhog) {
    return std::sqrt(3.0 * M_PI / (32.0 * nicole::constants::kGravitationalConstant * rhog));
}


double BarotropicEOS(const double rhog) {
    // Tsukamoto et al. 2020
    const double gamma = 7.0/5.0;
    const double rho_crit = 4.0e-14;
    const double T = 10 * (1 + gamma*std::pow(rhog/rho_crit, gamma - 1));
    return T;
}


double IonizationRate(const double rhog, const double T) {
    const double SigmaCR = 96.0;
    const double Sigma = std::sqrt(nicole::constants::kBoltzmannConstant * T * rhog /
         (M_PI * nicole::constants::kGravitationalConstant * nicole::constants::kGasMolecularMass));
    const double zetaCR0 = 1.3e-17;
    const double zetaCR = zetaCR0 * std::exp(-Sigma / SigmaCR);
    const double zetaRA = 7.3e-19;
    return zetaCR + zetaRA;
}


double MagneticField(const double nH) {
    return 1.43e-7 * std::sqrt(nH);
}

#endif /* FUNC_HPP_ */