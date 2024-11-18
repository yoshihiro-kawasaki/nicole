#ifndef _PHYSICALCONSTANT_HPP_
#define _PHYSICALCONSTANT_HPP_

#define _USE_MATH_DEFINES
#include <cmath>

//***************************************************************//
//************************ 物理定数 *****************************//
//***************************************************************//

// cgs Gauss unit
#define SPEED_OF_LIGHT            (2.99792458e10)                             // 光速
#define GRAVITATIONAL_CONSTANT    (6.6743e-8)                                 // 万有引力定数
#define PLANCK_CONSTANT           (6.6260755e-27)                             // プランク定数
#define DIRAC_CONSTANT            (PLANCK_CONSTANT/(2*M_PI))                  // ディラック定数
#define BOLTZMANN_CONSTANT        (1.380658e-16)                              // ボルツマン定数

#define ATOMIC_MASS_UNIT          (1.6605e-24)                                // 原子質量単位
#define PROTON_MASS               (1.6726e-24)                                // 陽子質量
#define NEUTRON_MASS              (1.6749e-14)                                // 中性子質量
#define ELECTRON_MASS             (9.1093897e-28)                             // 電子質量
#define HYDROGEN_MASS             (1.67353e-24)                               // 水素質量
#define CHARGE_UNIT               (4.803204673e-10)                           // 電荷単位
#define ELECTRON_VOLT             (1.60218e-12)                               // 1eV
#define BOHR_RADIUS               (5.2918e-9)                                 // ボーア半径
#define CLASSICAL_ELECTRON_RADIUS (2.8179e-13)                                // 古典

#define AVOGADRO_CONSTANT         (6.0221e23)                                 // アボガドロ定数
#define GAS_CONSTANT_MOL          (8.3145e7)                                  // 1molの気体定数
#define STEFAN_BOLTZMANN_CONSTANT (5.6705e-5)                                 // ステファン-ボルツマン定数
#define RADIATION_CONSTANT        (7.5646e-15)                                // 輻射定数 

#define ASTRONOMICAL_UNIT         (1.496e+13)                                 // 天文単位(cm)
#define LIGHT_YEAR                (9.4605e17)                                 // 光年(cm)
#define PARSEC                    (3.0857e18)                                 // パーセク(cm)
#define SOLAR_YEAR                (3.1557e7)                                  // 太陽年(s)

#define SOLAR_MASS                (1.989e33)                                  // 太陽質量
#define SOLAR_RADIUS              (6.960e10)                                  // 太陽（赤道）半径

#define GAS_SPECIFIC_HEAT         (1.4)                                       // ガス比熱
#define GAS_MOLECULAR_WEIGHT      (2.34)                                      // 平均分子量(ガス)
#define METAL_MOLECULAR_WEIGHT    (24.0)                                      // 平均分子量(金属)
#define GAS_MOLECULAR_MASS        (GAS_MOLECULAR_WEIGHT*PROTON_MASS)          // 平均分子質量(ガス)
#define METAL_MOLECULAR_MASS      (METAL_MOLECULAR_WEIGHT*PROTON_MASS)        // 平均分子質量(金属)
#define POTASSIUM_MASS            (39*PROTON_MASS)                            // カリウム質量
#define K_IONIZATION_POTENTIAL    (4.34*ELECTRON_VOLT)                        // カリウムのイオン化エネルギー

#endif /* _PHYSICALCONSTANT_HPP_ */