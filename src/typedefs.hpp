#ifndef TYPEDEFS_HPP_
#define TYPEDEFS_HPP_

#include "physical_constants.hpp"
#include "utils.hpp"

#include <vector>
#include <string>

#define IONIZATION                     ("IO")
#define ASSOCIATIVE_DETACHMENT         ("AD")
#define COLLISIONAL_DISSOCIATION       ("CD")
#define CHARGE_EXCHANGE                ("CE")
#define COSMIC_RAY_PROTON              ("CP")
#define COSMIC_RAY_PHOTON              ("CR")
#define DISSOCIATIVE_RECOMBINATION     ("DR")
#define ION_NEUTRAL                    ("IN")
#define MUTUAL_NEUTRALIZATION          ("MN")
#define NEUTRAL_NEUTRAL                ("NN")
#define PHOTOPROCESS                   ("PH")
#define RADIATIVE_ASSOCIATION          ("RA")
#define RADIATIVE_ELECTRON_ATTACHMENT  ("REA")
#define RADIATIVE_RECOMBINATION        ("RR")

#define ADSORPTION_NEUTRAL             ("AN")
#define ADSOPTION_ION_COUNTER          ("AIC")
#define DESORPTION_NEUTRAL             ("DN")
#define ION_GRAIN_COLLISION_COUNTER    ("CIG-C")
#define ION_GRAIN_COLLISION_NO_COUNTER ("CIG-NOC")
#define ELECTRON_GRAIN_COLLISION       ("EGC")
#define GRAIN_GRAIN_COLLISION          ("GGC")
#define H2_FORMATION_ON_GRAIN_SURFACE  ("H2FGS")

#define POLARIZABILITY                 (1.0)
#define TRANSFERRED_ENERGY             (2.0e-3 * ELECTRON_VOLT) // the amount of energy transferred to the grain particle as lattice vibrations
#define BINDING_SITES                  (1.0e15)                 // the surface number density of binding sites
#define VISUAL_EXTINCTION              (1.0)                    // the extinction by interstellar dust at visible wavelengths

#define POLARIZABILITY_H               (0.667)
#define POLARIZABILITY_H2              (0.804)
#define POLARIZABILITY_He              (0.207)

#define NOT_FOUND_SPECIES_LIST         (999)
#define FAILED_TO_OPEN(file_stream, file_name) if(!file_stream) {std::cerr << "##ERROR : Failed to open " << file_name << std::endl; exit(1); }
#define DEBUG std::cout << "ok" << std::endl;
#define DEBUG_MESSAGE(message) std::cout << message << std::endl;



typedef unsigned int              UnsignedInt;
typedef std::vector<unsigned int> UnsignedIntVector;
typedef std::vector<int>          IntVector;
typedef std::vector<double>       DoubleVector;
typedef std::vector<std::vector<double > > DoubleMatrix;
typedef std::vector<std::string>  StringVector;
typedef Utility::Dictionary       Dictionary;
typedef Utility::InputConfigure   InputConfigure;

#endif /* DEFINE_HPP_ */