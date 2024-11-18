from mkchemicalfile import chemical_file as chf

## Ilgner & Nelson 2006
SPECIES_LIST = ["H+", "H", "H2", "H2+", "H3+", "He", "He+", "C", "C+", "CH", "CH+", "CH2", "CH2+", "N+", "N", "CH3+", "NH+", "CH3", "NH", "CH4+", "NH2+", "O", "NH2", "O+", "CH4", 
                "OH", "NH3", "NH3+", "OH+", "CH5+", "H2O", "NH4+", "H2O+", "H3O+", "Mg", "C2", "C2+", "Mg+", "C2H+", "C2H", "CN", "C2H2", "CN+", "C2H2+", "HCN","C2H3+", "HNC", 
                "C2H3", "HCN+", "H2NC+", "Si+", "HCNH+", "CO+", "CO", "Si", "N2+", "N2", "C2H4+", "SiH+", "HCO+", "HCO", "SiH", "N2H+", "NO+", "H2CO+", "SiH2", "NO", "SiH2+","H2CO", 
                "H3CO+", "SiH3", "SiH3+","S+", "CH3OH+", "SiH4", "CH3OH", "O2+", "S", "SiH4+", "O2", "HS+", "CH3OH2+", "HS", "SiH5+", "H2S", "H2S+", "H3S+", "C3", "C3+", "C3H", 
                "C3H+", "C3H2", "C3H2+", "C2N+", "CNC+", "CH2CCH", "C3H3+", "CH3CCH", "SiC+", "C3H4+", "SiC", "C2O+", "HCSi", "HC2O+", "HCSi+", "CH3CN+", "C3H5+", "CH3CN", "CH2CO", 
                "SiN+", "SiN", "CH2CO+", "CH3CNH+", "CH3CO+", "HNSi", "HNSi+", "CO2+", "SiO+", "CS", "SiO", "CO2", "CS+", "HCS", "HCO2+", "SiOH+", "HCS+", "NS+", "NS", "H2CS", "H2CS+", 
                "H3CS+", "HNS+", "C4+", "SO+", "SO", "C4", "HSO+", "C4H+", "C4H", "C3N", "C4H2+", "C3N+", "HC4H", "HC3N+", "HC3N", "C3O+", "C3O", "HC3NH+","HC3O+", "C3H2O+", 
                "H3C3O+", "C2S+", "C2S", "Fe+", "Fe", "HC2S+", "SiS", "OCS+", "SiO2", "SiS+", "OCS", "HSiS+", "HOCS+", "SO2+", "SO2", "S2", "HSO2+", "CH3C3NH+", "H2S2+", "C3S",
                "C3S+", "HC3S+", "C7+", "e-"]
                # Ilgner & Nelson (2006), Bai & Goodman (2009)

## P.Marchaned et al. 2016, koga et al. 2019
# SPECIES_LIST = ["H", "H2", "He", "C", "O", "O2", "CO", "HCO", "Mg", "H+", "H2+", "H3+", "He+", "C+", "O+", "O2+", "OH+", "O2H+", "H2O+", "CO+", "HCO+", "CH2+", "Mg+", "e-"]

## Tsukamoto et al. 2020
# SPECIES_LIST = ["H", "H2", "He", "CO", "O2", "Mg", "O", "C", "HCO", "H2O", "OH", "H3+", "H2+", "H+", "HCO+", "Mg+", "He+", "C+", "O+", "O2+", "H3O+", "OH+","H2O+", "e-"]

## P.Marchaned et al. 2016, koga et al. 2019 + Tsukamoto et al. 2020
SPECIES_LIST = ["H", "H2", "He", "CO", "O2", "Mg", "O", "C", "HCO", "H2O", "OH", "H3+", "H2+", "H+", "HCO+", "Mg+", "He+", "C+", "O+", "CO+", "CH2+", "O2+", "H3O+", "OH+", "O2H+", "H2O+", "e-"]

# SPECIES_LIST = ["H2", "H", "He", "CO", "O2", "Mg", "O", "C", "HCO", "H2O", "OH", "Si", "S", "N", "N2", "NH3", "CH4",
#                 "H3+", "H2+", "H+", "HCO+", "Mg+", "He+", "C+", "O+", "CO+", "CH2+", "O2+", "H3O+", "OH+", "O2H+", "H2O+", "Si+", "S+", "SiO+", "CH+", "N+", "N2+", "NH4+", "N2H+",
#                 "e-"]

# SPECIES_LIST = ["H2", "H", "He", "C", "CH", "CH2", "CH3", "CH4", "N", "N2", "NH", "NH2", "NH3", "O", "O2", "OH", "H2O", "CO", "CO2", "Mg", "Fe",
#                 "H+", "H2+", "H3+", "He+", "HeH+", "C+", "CH+", "CH2+", "CH3+", "CH4+", "CH5+", "N+", "N2+", "N2H+", "NH+", "NH2+", "NH3+", "NH4+", "O+", "O2+", "O2H+", "OH+", "H2O+",
#                 "H3O+", "CO+", "HCO+", "CO2+", "HCO2+", "NO+", "Mg+", "Fe+", "e-"]


chf.ASSOCIATIVE_DETACHMENT = True 
chf.COLLISIONAL_DISSOCIATION = True 
chf.CHARGE_EXCHANGE = True
chf.COSMIC_RAY_PROTON = False
chf.COSMIC_RAY_PHOTON = False
chf.DISSOCIATIVE_RECOMBINATION = True 
chf.ION_NEUTRAL = True 
chf.MUTUAL_NEUTRALIZATION = True 
chf.NEUTRAL_NEUTRAL = True 
chf.PHOTOPROCESS = False
chf.RADIATIVE_ASSOCIATION = True 
chf.RADIATIVE_ELECTRON_ATTACHMENT = True 
chf.RADIATIVE_RECOMBINATION = True 

chf.MAX_DUST_CHARGE = 1
chf.NUMBER_OF_BINS = 1
chf.GRAIN_REACTION = True
chf.FLAG_ADSORPTION = False
chf.FLAG_DESORPTION = False
chf.CHARGE_ADSOPTION = False
chf.FLAG_H2_FORMATION_ON_GRAIN_SURFACE = True

chf.LOWER_T_GPR = 10
chf.UPPER_T = True
chf.UPPER_T_GPR = 50000


## OutPut File ##
Species_File = "input/species.txt"
Rate_File = "input/reactions.txt" 

file = chf.ChemicalReaction()
file.input_species(SPECIES_LIST)
file.write_file(Species_File, Rate_File)
file.print_total()
