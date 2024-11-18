# coding: UTF-8
import numpy as np 
import sys 
import copy 
from . import const

## Umist Database
Species_File_From_UMIST = "mkchemicalfile/UMIST_database/SPECIES12.txt"
Rate_File_From_UMIST = "mkchemicalfile/UMIST_database/RATE12.dist.txt"
Binding_File_From_UMIST = "mkchemicalfile/UMIST_database/RATE12_binding_energies.dist.txt"

## Global variables
ASSOCIATIVE_DETACHMENT        = True
COLLISIONAL_DISSOCIATION      = True 
CHARGE_EXCHANGE               = True
COSMIC_RAY_PROTON             = False
COSMIC_RAY_PHOTON             = True 
DISSOCIATIVE_RECOMBINATION    = True 
ION_NEUTRAL                   = True 
MUTUAL_NEUTRALIZATION         = True 
NEUTRAL_NEUTRAL               = True 
PHOTOPROCESS                  = False 
RADIATIVE_ASSOCIATION         = True 
RADIATIVE_ELECTRON_ATTACHMENT = True 
RADIATIVE_RECOMBINATION       = True 

MAX_DUST_CHARGE = 1        ## Max of dust grain charge
NUMBER_OF_BINS  = 1

DUST_REACTION     = True
FLAG_ADSORPTION    = False
FLAG_DESORPTION    = False
CHARGE_ADSOPTION   = False # If ADSORPTION = False and  DESORPTION = False, ion colliding with dust grains is adsorbed , or 
                         # experience charge exchange and not adsorped. In the case of CHARGE_ABSOPTION == True, colliding ion is
                         # adsorbed onto dust grains,  on the other hand if CHARGE_ABSOPTION == False, instead of absorption,  
                         # ion excahnges charge with dust grains and then become neutral species .
FLAG_H2_FORMATION_ON_GRAIN_SURFACE     = False  

LOWER_T_GPR = 10
UPPER_T = True   # upper_T = Trueの場合、気相反応の上限温度はUMIST database記載のものとする。upper_T = falseの場合、全ての気相反応の上限温度を
UPPER_T_GPR = 50000 # upper_T_GPRとする

## constant global 
const.MAX_INDEX_BLANKS = 5
const.MAX_BLANKS       = 12
const.DUMMY_MASS = "0"
const.ASTERISK = "*"
const.ADSORPTION = "ADSORP" ## 吸着
const.DESORPTION = "DESORP" ## 脱着
const.ZERO = "0.00"
const.ELECTRON = "e-"

## Type of Reaction 
const.IONIZATION = "IO"                                           ## Ionization by Commisc-Ray
const.GAS_PHASE_REACTION = "GPR"                                  ## Gas phase reaction from UMIST database
const.ADSORPTION_NEUTRAL = "AN"                                   ## Adsorption of neutral particles to grains
const.ADSORPTION_ION_COUNTER = "AIC"                              ## Adsorption of positive ions to grains (positive ion-positive charged grain collision which ion has counterpart)
const.DESORPTION_NEUTRAL = "DN"                                   ## Desorption of neutral particles from grains
const.ION_NEGATIVE_GRAIN_COLLISION_COUNTER = "CIG-C"              ## Ion-negative charged Grain Collision (A Ion has a counterpart(neutral particle) in species list)
const.ION_NEGATIVE_GRAIN_COLLISION_NO_COUNTER = "CIG-NOC"         ## Ion-negative charged Grain Collision (A ion have no counterpart in species list) 
const.ELECTRON_GRAIN_COLLISION = "EGC"                            ## Electron-Grain Collision
const.GRAIN_GRAIN_COLLISION = "GGC"                               ## Grain-Grain Collison 
const.H2_FORMATION_ON_GRAIN_SURFACE = "H2FGS"

## Ionization Source List ##
const.I1 = ["H2", "CRP", "H2+", "e-", "*",  "*",  "9.70e-01"]
const.I2 = ["H2", "CRP", "H+",  "H",  "e-", "*",  "3.00e-02"]
const.I3 = ["H",  "CRP", "H+",  "e-", "*",  "*",  "5.00e-01"]
const.I4 = ["He", "CRP", "He+", "e-", "*",  "*",  "8.40e-01"]
const.IONIZATION_SOURCE = [const.I1, const.I2, const.I3, const.I4]

class ChemicalReactionBase:

    def __init__(self):
        self.input_species_list = []
        self.number_of_species = 0
        self.ions_neutral_counterpart = []
        self.ions_no_neutral_counterpart = []
        self.total_ions = []
        self.neutral_species = []

        self.species_list = []
        self.species_mass = {}
        self.mantle_species = []
        self.hydrogen_mantle = []

        self.ion_electron_reaction_no_counterpart = []
        self.reaction_type = []

        self.associative_detachment = []
        self.collisional_dissociation = []
        self.charge_excahnge = []
        self.cosmic_ray_proton = []
        self.cosmic_ray_photon = []
        self.dissociative_recombination = []
        self.ion_neutral = []
        self.mutual_neutralization = []
        self.neutral_neutral = []
        self.photoprocess = []
        self.radiative_association = []
        self.radiative_electron_attachment = []
        self.radiative_recombination = []

        self.total_species = 0
        self.total_reactions = 0

        self.species_file = None 
        self.reaction_file = None
        self.abundance_file = "input/abundances.txt"

        self.dust_array = None

    def print_total(self):
        print("species file       : " + self.species_file)
        print("reaction file      : " + self.reaction_file)
        print("abundances_file    : " + self.abundance_file)

        if DUST_REACTION == True:
            if NUMBER_OF_BINS == 1:
                print("dust model         : single_size_model")
                print("number_of_bins     : " + str(NUMBER_OF_BINS))
                print("max charge of dust : " + str(MAX_DUST_CHARGE))
                print("grain species      : " + str(2*MAX_DUST_CHARGE + 1))
            else:
                print("dust model         : size distribution_model")
                print("number_of_bins     : " + str(NUMBER_OF_BINS))
                print("max charge of dust : " + str(MAX_DUST_CHARGE))
                print("grain species      : " + str(NUMBER_OF_BINS * (2*MAX_DUST_CHARGE + 1)))
        else:
            print("dust model         : no_dust_model")
            print("number_of_bins     : " + str(0))
            print("max charge of dust : " + str(0))
            print("grain species      : " + str(0))

        print("total_spcies       : " + str(self.total_species))
        print("total reaction     : " + str(self.total_reactions))

        input_data_file = "input/input_data.txt"
        file = open(input_data_file, mode="w")
        file.write("species_file       : " + self.species_file + "\n")
        file.write("reaction_file      : " + self.reaction_file + "\n")
        file.write("abundance_file     : " + self.abundance_file + "\n")
        
        if DUST_REACTION == True:
            if NUMBER_OF_BINS == 1:
                file.write("dust_model         : single_size_model" + "\n")
                file.write("number_of_bins     : " + str(NUMBER_OF_BINS) + "\n")
                file.write("max_charge_of_dust : " + str(MAX_DUST_CHARGE) + "\n")
                file.write("dust_grain_species : " + str((2*MAX_DUST_CHARGE + 1)) + "\n")
            else:
                file.write("dust_model         : size_distribution_model" + "\n")
                file.write("number_of_bins     : " + str(NUMBER_OF_BINS) + "\n")
                file.write("max_charge_of_dust : " + str(MAX_DUST_CHARGE) + "\n")
                file.write("dust_grain_species : " + str(NUMBER_OF_BINS * (2*MAX_DUST_CHARGE + 1)) + "\n")
        else:
            file.write("dust_model         : no_dust_model" + "\n")
            file.write("number_of_bins     : " + str(0) + "\n")
            file.write("max_charge_of_dust : " + str(0) + "\n")
            file.write("dust_grain_species : " + str(0) + "\n")

        file.write("total_spcies       : " + str(self.total_species) + "\n")
        file.write("total_reaction     : " + str(self.total_reactions) + "\n")
        file.close()

class ChemicalSpecies(ChemicalReactionBase):
        
    def read_from_file(self, file_name):
        self.input_species_list = np.loadtxt(file_name, str).tolist()
        self.number_of_species = len(self.input_species_list)
        self.species_classification()
        
    def read_from_list(self, species_list):
        self.input_species_list = species_list
        self.number_of_species = len(self.input_species_list)
        self.species_classification()
        
    def species_classification(self):
        
        for species in self.input_species_list:
            if species == "e-":
                pass
            else:
                if "+" in species:
                    species_neutral = species.replace("+", '')
                    if species_neutral in self.input_species_list:
                        self.ions_neutral_counterpart.append(species)
                    else:
                        self.ions_no_neutral_counterpart.append(species)
                elif "-" in species:
                    species_neutral = species.replace("-", '')
                    if species_neutral in self.input_species_list:
                        self.ions_neutral_counterpart.append(species)
                    else:
                        self.ions_no_neutral_counterpart.append(species)
        
        self.total_ions = self.ions_neutral_counterpart + self.ions_no_neutral_counterpart
                    
    def show_species(self):
        for i in range(self.number_of_species):
            print(self.input_species_list[i])

class ReadUmistDatabase(ChemicalReactionBase):
        
    def read_species(self, species_file):
        species_data = np.loadtxt(species_file, str)
        
        for species in species_data:
            if species[1] in self.input_species_list:
                self.species_list.append(species[1])
                self.species_mass[species[1]] = species[2]
        
        diff_species = set(self.input_species_list) - set(self.species_list)
        diff_species = list(diff_species)
        if not (len(diff_species) == 0):
            print("Does not exist in UMIST database : ")
            for species in diff_species:
                print(species)
                
    def read_binding_energy(self, binding_energy_file):
        if GRAIN_REACTION == True:
            file = open(binding_energy_file, mode='rt')
            binding_energy = file.readlines()
            length_list = len(binding_energy)
            for _ in range(length_list):
                self.mantle_species.append([])

            for line in binding_energy:
                if line[0] == "#":
                    pass 
                else:
                    be_list = [i for i in line.split(' ') if i != '']
                    if be_list[0] == "C3H3":
                        be_list[0] = "CH2CCH"
                    if be_list[0] == "C4H2":
                        be_list[0] = "HC4H"
                    if be_list[0] in self.species_list:
                        del be_list[2:]
                        index_number = self.species_list.index(be_list[0])
                        be_list = [be_list[0], 'g'+be_list[0], be_list[1]]
                        self.mantle_species[index_number] = be_list

            self.mantle_species = [i for i in self.mantle_species if i != []]
            
            if FLAG_H2_FORMATION_ON_GRAIN_SURFACE == True:
                for mantle in self.mantle_species:
                    if mantle[0] == 'H':
                        self.hydrogen_mantle = mantle
                        break

            file.close()

    def read_reactions(self, reaction_rate_file):
        
        if ASSOCIATIVE_DETACHMENT == True:
            self.reaction_type.append("AD")
        if COLLISIONAL_DISSOCIATION == True:
            self.reaction_type.append("CD")
        if CHARGE_EXCHANGE == True:
            self.reaction_type.append("CE")
        if COSMIC_RAY_PROTON == True:
            self.reaction_type.append("CP")
        if COSMIC_RAY_PHOTON == True:
            self.reaction_type.append("CR")
        if DISSOCIATIVE_RECOMBINATION == True:
            self.reaction_type.append("DR")
        if ION_NEUTRAL == True:
            self.reaction_type.append("IN")
        if MUTUAL_NEUTRALIZATION == True:
            self.reaction_type.append("ME")
        if NEUTRAL_NEUTRAL == True:
            self.reaction_type.append("NN")
        if PHOTOPROCESS == True:
            self.reaction_type.append("PH")
        if RADIATIVE_ASSOCIATION == True:
            self.reaction_type.append("RA")
        if RADIATIVE_ELECTRON_ATTACHMENT == True:
            self.reaction_type.append("REA")
        if RADIATIVE_RECOMBINATION == True:
            self.reaction_type.append("RR")

        file = open(reaction_rate_file, mode='rt')
        rate_data = file.readlines()

        for reaction in rate_data:
            reaction = reaction.split(":")

            if reaction[1] in self.reaction_type:
                pass 
            else:
                continue

            if reaction[1] == "AD":
                self.append_reaction_list(reaction, self.associative_detachment)
            elif reaction[1] == "CD":
                self.append_reaction_list(reaction, self.collisional_dissociation)
            elif reaction[1] == "CE":
                self.append_reaction_list(reaction, self.charge_excahnge)
            elif reaction[1] == "CP":
                self.append_reaction_list(reaction, self.cosmic_ray_proton)
            elif reaction[1] == "CR":
                self.append_reaction_list(reaction, self.cosmic_ray_photon)
            elif reaction[1] == "DR":
                self.append_reaction_list(reaction, self.dissociative_recombination)
            elif reaction[1] == "IN":
                self.append_reaction_list(reaction, self.ion_neutral)
            elif reaction[1] == "MN":
                self.append_reaction_list(reaction, self.mutual_neutralization)
            elif reaction[1] == "NN":
                self.append_reaction_list(reaction, self.neutral_neutral)
            elif reaction[1] == "PH":
                self.append_reaction_list(reaction, self.photoprocess)
            elif reaction[1] == "RA":
                self.append_reaction_list(reaction, self.radiative_association)
            elif reaction[1] == "REA":
                self.append_reaction_list(reaction, self.radiative_electron_attachment)
            elif reaction[1] == "RR":
                self.append_reaction_list(reaction, self.radiative_recombination)

        del rate_data
        file.close()
            
    def check_species(self, species_name):
        species_list_plus = self.species_list + ["PHOTON", "CRP", "CRPHOT"]
        if (species_name == '') or (species_name == ' '):
            return True
        else:
            if species_name in species_list_plus:
                return True
            else:
                return False
    
    def list_brank(self, reaction_list):
        for i in range(len(reaction_list)):
            if (reaction_list[i] == '') or (reaction_list[i] == ' '):
                reaction_list[i] = '*'
            else:
                pass 
    
    def append_reaction_list(self, reaction, reaction_list):
        del reaction[1]
        del reaction[7]
        del reaction[12: ]
        del reaction[0]
        if ((self.check_species(reaction[0])) and (self.check_species(reaction[1])) and (self.check_species(reaction[2])) 
            and (self.check_species(reaction[3])) and (self.check_species(reaction[4])) and (self.check_species(reaction[5]))):
            self.list_brank(reaction)
            if UPPER_T == False:
                reaction[-1] = str(UPPER_T_GPR)
            reaction_list.append(reaction)

            ## イオンー電子再結合反応のうち、イオンに対応する中性粒子がSpecies Listに存在しないものをIon_Electron_Reaction_No_IonCounterに格納
            if const.ELECTRON in reaction:
                if (reaction[1] == const.ELECTRON) and (reaction[0] in self.ions_no_neutral_counterpart):
                    self.ion_electron_reaction_no_counterpart.append(reaction)
                    # print(reaction)
                elif reaction[0] == const.ELECTRON:
                    print("errer : append_reaction_list")
                    print(reaction)

class WriteFile(ChemicalReactionBase):

    def set_species_write(self, index_number, species_list):
        number_of_blanks1 = 4 - len(str(index_number))
        number_of_blanks2 = 2 * const.MAX_BLANKS - len(species_list[0]) - len(species_list[1])
        number_of_blanks3 = 8 - len(species_list[2])
        s = (" "*number_of_blanks1 + str(index_number) + "  " + species_list[0] + " "*number_of_blanks2 + species_list[1] + 
             " "*number_of_blanks3 + species_list[2] + "\n")
        return s 

    def species_charge(self, species):
        if "+" in species:
            plus_index = species.find("+")
            len_species = len(species)
            if plus_index == (len_species - 1):
                return "1"
            else:
                return species[plus_index + 1]
        elif "-" in species:
            minus_index = species.find("-")
            len_species = len(species)
            if minus_index == (len_species - 1):
                return "-1"
            else:
                return ("-" + species[minus_index + 1])
        else:
            return "0"

    def write_gas_phase_species_file(self, species_file):
        file = open(species_file, mode='a')
        index_number = 1

        ## Gas Phase Species
        for species in self.species_list:
            mass = self.species_mass[species]
            spe_ch = self.species_charge(species)
            species_list = [species, mass, spe_ch]
            s = self.set_species_write(index_number, species_list)
            file.write(s)
            index_number += 1

        return index_number

    def set_reaction_write(self, index_number, type_of_reaction, reaction_list):

        length_list = len(reaction_list)

        if length_list < 11:
            i = 11 - length_list
            for _ in range(i):
                reaction_list.append(const.ZERO)
        
        index_blanks = const.MAX_INDEX_BLANKS - len(str(index_number))
        number_of_blanks1 = const.MAX_BLANKS - len(type_of_reaction)
        number_of_blanks2 = const.MAX_BLANKS - len(reaction_list[0])
        number_of_blanks3 = const.MAX_BLANKS - len(reaction_list[1])
        number_of_blanks4 = const.MAX_BLANKS - len(reaction_list[2])
        number_of_blanks5 = const.MAX_BLANKS - len(reaction_list[3])
        number_of_blanks6 = const.MAX_BLANKS - len(reaction_list[4])
        number_of_blanks7 = 2*const.MAX_BLANKS - len(reaction_list[5]) - len(reaction_list[6])
        number_of_blanks8 = const.MAX_BLANKS - len(reaction_list[7])
        number_of_blanks9 = const.MAX_BLANKS - len(reaction_list[8])
        number_of_blanks10 = const.MAX_BLANKS - len(reaction_list[9])
        number_of_blanks11 = const.MAX_BLANKS - len(reaction_list[10])

        s = (" "*index_blanks + str(index_number) + "  " + type_of_reaction + " "*number_of_blanks1
            + reaction_list[0] + " "*number_of_blanks2 + reaction_list[1] + " "*number_of_blanks3 + reaction_list[2] + " "*number_of_blanks4
            + reaction_list[3] + " "*number_of_blanks5 + reaction_list[4] + " "*number_of_blanks6 + reaction_list[5]
            + " "*number_of_blanks7 + reaction_list[6] + " "*number_of_blanks8 + reaction_list[7] + " "*number_of_blanks9 + reaction_list[8]
            + " "*number_of_blanks10 + reaction_list[9] + " "*number_of_blanks11 + reaction_list[10] + "\n")

        return s 

    def write_gas_phase_reaction_file(self, reaction_file):
        file = open(reaction_file, mode='a')
        index_number = 1

        ## GAS PHASE REACTION ##
        if COSMIC_RAY_PROTON == True:
            for reaction in self.cosmic_ray_proton:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "CP", reaction)
                    file.write(s)
                    index_number += 1
        else:
            for reaction in const.IONIZATION_SOURCE:
                s = self.set_reaction_write(index_number, const.IONIZATION, reaction)
                file.write(s)
                index_number += 1
        
        if COSMIC_RAY_PHOTON == True:
            for reaction in self.cosmic_ray_photon:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "CR", reaction)
                    file.write(s)
                    index_number += 1
        
        if ASSOCIATIVE_DETACHMENT == True:
            for reaction in self.associative_detachment:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "AD", reaction)
                    file.write(s)
                    index_number += 1
        
        if COLLISIONAL_DISSOCIATION == True:
            for reaction in self.collisional_dissociation:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "CD", reaction)
                    file.write(s)
                    index_number += 1

        if CHARGE_EXCHANGE == True:
            for reaction in self.charge_excahnge:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "CE", reaction)
                    file.write(s)
                    index_number += 1

        if DISSOCIATIVE_RECOMBINATION == True:
            for reaction in self.dissociative_recombination:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "DR", reaction)
                    file.write(s)
                    index_number += 1

        if ION_NEUTRAL == True:
            for reaction in self.ion_neutral:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "IN", reaction)
                    file.write(s)
                    index_number += 1

        if MUTUAL_NEUTRALIZATION == True:
            for reaction in self.mutual_neutralization:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "MN", reaction)
                    file.write(s)
                    index_number += 1

        if NEUTRAL_NEUTRAL == True:
            for reaction in self.neutral_neutral:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "NN", reaction)
                    file.write(s)
                    index_number += 1

        if PHOTOPROCESS == True:
            for reaction in self.photoprocess:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "PH", reaction)
                    file.write(s)
                    index_number += 1

        if RADIATIVE_ASSOCIATION == True:
            for reaction in self.radiative_association:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "RA", reaction)
                    file.write(s)
                    index_number += 1

        if RADIATIVE_ELECTRON_ATTACHMENT == True:
            for reaction in self.radiative_electron_attachment:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "REA", reaction)
                    file.write(s)
                    index_number += 1

        if RADIATIVE_RECOMBINATION == True:
            for reaction in self.radiative_recombination:
                lower_temp = float(reaction[-2])
                if lower_temp > LOWER_T_GPR:
                    continue
                else:
                    s = self.set_reaction_write(index_number, "RR", reaction)
                    file.write(s)
                    index_number += 1

        return index_number

class Dust(WriteFile):

    def set_dust_array(self):
        self.dust_array = np.zeros([NUMBER_OF_BINS, 2*MAX_DUST_CHARGE+1], dtype=object)
        for s in range(NUMBER_OF_BINS):
            for dust_charge in range(-MAX_DUST_CHARGE, MAX_DUST_CHARGE+1):
                if dust_charge <= 0:
                    G = "Gs" + str(s+1) + "c" + str(dust_charge)
                else:
                    G = "Gs" + str(s+1) + "c+" + str(dust_charge)
                ich = MAX_DUST_CHARGE + dust_charge
                self.dust_array[s][ich] = G 

    def write_species_file(self, species_file):
        index_number = self.write_gas_phase_species_file(species_file)
        file = open(species_file, mode='a')
        ## Mantle Species ##
        if (FLAG_ADSORPTION == True) and (GRAIN_REACTION == True):
            for mantle in self.mantle_species:
                mass = self.species_mass[mantle[0]]
                spe_ch = self.species_charge(mantle[1])
                species_list = [mantle[1], mass, spe_ch]
                s = self.set_species_write(index_number, species_list)
                file.write(s)
                index_number += 1

        if GRAIN_REACTION == True:
            for s in range(NUMBER_OF_BINS):
                for ch in range(-MAX_DUST_CHARGE, MAX_DUST_CHARGE+1):
                    ich = MAX_DUST_CHARGE + ch
                    species_list = [self.dust_array[s][ich], const.DUMMY_MASS, str(ch)]
                    ss = self.set_species_write(index_number, species_list)
                    file.write(ss)
                    index_number += 1

        file.close()
        self.total_species = index_number - 1

    def write_reaction_file(self, reaction_file):
        index_number = self.write_gas_phase_reaction_file(reaction_file)

        file = open(reaction_file, mode='a')

        if DUST_REACTION == True:

            if FLAG_ADSORPTION == True:
                """
                Adsorption of neutral particles to grains 
                ガス中の中性粒子がダスト表面に吸着してマントル種形成。
                Ilgner & Nelson (2006) Table3 1-5
                ex. X + G(+) → X(G(+))
                反応係数がダストサイズに依存からbinで分けたほうがいいかも。
                それか、このままで計算する時にうまくする。
                """
                for mantle in self.mantle_species:
                    mass = self.species_mass[mantle[0]]
                    reaction = [mantle[0], const.ADSORPTION, mantle[1], const.ASTERISK, const.ASTERISK, const.ASTERISK, mantle[2], mass]
                    s = self.set_reaction_write(index_number, const.ADSORPTION_NEUTRAL, reaction)
                    file.write(s)
                    index_number += 1

            if FLAG_ADSORPTION == True:
                """
                (陽)イオンが正に帯電したダスト、あるいは帯電していないダストと衝突したとき
                イオンはダスト表面上に吸着し中性化され、mantle種を形成する。
                Ilgner & Nelson (2006) Table3 7-8
                X+ + G(+) → X(G(2+))
                X+ + G(0) → X(G(+)) 
                """
                for mantle in self.mantle_species:
                    mantle_ion = mantle[0] + "+"
                    mass = self.species_mass[mantle[0]]
                    if not (mantle_ion in self.ions_neutral_counterpart):
                        continue
                    else:
                        pass 
                    for s in range(NUMBER_OF_BINS):
                        for dust_charge in range(MAX_DUST_CHARGE):
                            ich = MAX_DUST_CHARGE + dust_charge
                            reac_grain1 = self.dust_array[s][ich]
                            prod_grain1 = self.dust_array[s][ich+1]
                            ## a = reactant1's mass, b = dust cherge, c = bin number
                            reaction = [mantle_ion, reac_grain1, prod_grain1, mantle[1], const.ASTERISK, const.ASTERISK, mass, str(dust_charge), str(s+1)] 
                            ss = self.set_reaction_write(index_number, const.ADSORPTION_ION_COUNTER, reaction)
                            file.write(ss)
                            index_number += 1
            else:
                if CHARGE_ADSOPTION == True:
                    """
                    荷電粒子がダストに吸着 ex) H+ + G0 → G+1
                    P.Marchand et al. 2016,  koga et al. 2019
                    マントル種は考えない。
                    """
                    for ion in self.total_ions:
                        mass = self.species_mass[ion]

                        for s in range(NUMBER_OF_BINS):
                            for dust_charge in range(MAX_DUST_CHARGE):
                                ich = MAX_DUST_CHARGE + dust_charge 

                                reac_grain1 = self.dust_array[s][ich]
                                prod_grain1 = self.dust_array[s][ich+1]

                                reaction = [ion, reac_grain1, prod_grain1, const.ASTERISK, const.ASTERISK, const.ASTERISK, mass, str(dust_charge), str(s+1)]
                                ss = self.set_reaction_write(index_number, const.ADSORPTION_ION_COUNTER, reaction)
                                file.write(ss)
                                index_number += 1
                else:
                    """
                    荷電粒子が(正に帯電した、あるいは中性)ダストと衝突後、ダストと電荷交換して中性種となり気相へ流出
                    Kunz & Mouschovias 2009, Ilgner and Nelson 2006 Table2 (Sano et al. 2000)
                    ex. X+ + G(0) → X + G(+)
                        X+ + G(+) → X + G(2+)
                    """
                    for mantle in self.mantle_species:
                        mantle_ion = mantle[0] + "+"
                        mass = self.species_mass[mantle[0]]
                        if not (mantle_ion in self.ions_neutral_counterpart):
                            continue
                        else:
                            pass 
                        for s in range(NUMBER_OF_BINS):
                            for dust_charge in range(MAX_DUST_CHARGE):
                                ich = MAX_DUST_CHARGE + dust_charge
                                reac_grain1 = self.dust_array[s][ich]
                                prod_grain1 = self.dust_array[s][ich+1]
                                reaction = [mantle_ion, reac_grain1, prod_grain1, mantle_ion[:-1], const.ASTERISK, const.ASTERISK, mass, str(dust_charge), str(s+1)]
                                ss = self.set_reaction_write(index_number, const.ADSORPTION_ION_COUNTER, reaction)
                                file.write(ss)
                                index_number += 1

                    for ion in self.ions_no_neutral_counterpart:
                        reaction_list = []
                        mass = self.species_mass[ion]
                        for reaction in self.ion_electron_reaction_no_counterpart:
                            if reaction[0] == ion:
                                reaction_list.append(reaction)
                        number_of_blanching = len(reaction_list)
                        for reaction in reaction_list:
                            for s in range(NUMBER_OF_BINS):
                                for dust_charge in range(MAX_DUST_CHARGE):
                                    ich = MAX_DUST_CHARGE + dust_charge
                                    reaction_copy = copy.copy(reaction)
                                    del reaction_copy[6:]
                                    reac_grain1 = self.dust_array[s][ich]
                                    prod_grain1 = self.dust_array[s][ich+1]
                                    reaction_copy[1] = reac_grain1
                                    if (reaction_copy[-2] == "*") and (reaction_copy[-1] == "*"):
                                        reaction_copy[-2] = prod_grain1
                                    elif (reaction_copy[-2] != "*") and (reaction_copy[-1] == "*"):
                                        reaction_copy[-1] = prod_grain1
                                    reaction_copy.append(mass)
                                    reaction_copy.append(str(dust_charge))
                                    reaction_copy.append(str(s+1))
                                    reaction_copy.append(str(number_of_blanching))
                                    ss = self.set_reaction_write(index_number, const.ADSORPTION_ION_COUNTER, reaction_copy)
                                    file.write(ss)
                                    index_number += 1
        
            if FLAG_DESORPTION == True:
                """
                Desorption of neutral particles from grains ##
                Ilgner & Nelson (2006) Table3 6
                脱着過程はダストサイズ(bin)に依存しない
                """
                for mantle in self.mantle_species:
                    mass = self.species_mass[mantle[0]]
                    reaction = [mantle[1], const.DESORPTION, mantle[0], const.ASTERISK, const.ASTERISK, const.ASTERISK, mantle[2], mass]
                    s = self.set_reaction_write(index_number, const.DESORPTION_NEUTRAL, reaction)
                    file.write(s)
                    index_number += 1

         #if DUST_REACTION == True:
            """
            Ion-negative charged Grain Collision (A Ion has a counterpart(neutral particle) 
            in species list)

            Ilgner & Nelson (2006) Table4 1-2
            """
            for ion in self.ions_neutral_counterpart:
                mass = self.species_mass[ion]
                neutral = ion.replace("+", '')
                for s in range(NUMBER_OF_BINS):
                    for dust_charge in range(-1, -MAX_DUST_CHARGE-1, -1):
                        ich = MAX_DUST_CHARGE + dust_charge

                        reac_grain1 = self.dust_array[s][ich]
                        prod_grain1 = self.dust_array[s][ich+1]

                        reaction = [ion, reac_grain1, prod_grain1, neutral, const.ASTERISK, const.ASTERISK, mass, str(dust_charge), str(s+1)]
                        ss = self.set_reaction_write(index_number, const.ION_NEGATIVE_GRAIN_COLLISION_COUNTER, reaction)
                        file.write(ss)
                        index_number += 1

         #if DUST_REACTION == True:
            """
            Ion-negative charged Grain Collision (A ion have no counterpart in species list) 
            Ilgner & Nelson (2006) Table4 1-2
            """
            for ion in self.ions_no_neutral_counterpart:
                reaction_list = []
                mass = self.species_mass[ion]
                for reaction in self.ion_electron_reaction_no_counterpart:
                    if reaction[0] == ion:
                        reaction_list.append(reaction)
                number_of_blanching = len(reaction_list)
                for reaction in reaction_list:
                    for s in range(NUMBER_OF_BINS):
                        for dust_charge in range(-1, -MAX_DUST_CHARGE-1, -1):
                            reaction_copy = copy.copy(reaction)
                            del reaction_copy[6:]

                            reac_grain1 = self.dust_array[s][ich]
                            prod_grain1 = self.dust_array[s][ich+1]

                            reaction_copy[1] = reac_grain1
                            if (reaction_copy[-2] == "*") and (reaction_copy[-1] == "*"):
                                reaction_copy[-2] = prod_grain1
                            elif (reaction_copy[-2] != "*") and (reaction_copy[-1] == "*"):
                                reaction_copy[-1] = prod_grain1
                            
                            reaction_copy.append(mass)
                            reaction_copy.append(str(dust_charge))
                            reaction_copy.append(str(s+1))
                            reaction_copy.append(str(number_of_blanching))
                            ss = self.set_reaction_write(index_number, const.ION_NEGATIVE_GRAIN_COLLISION_NO_COUNTER, reaction_copy)
                            file.write(ss)
                            index_number += 1

         #if DUST_REACTION == True:
            """
            Electron-Grain Collision ##
            Ilgner & Nelson (2006) Table4 3-6
            """
            for s in range(NUMBER_OF_BINS):
                for dust_charge in range(-MAX_DUST_CHARGE+1, MAX_DUST_CHARGE+1):
                    ich = MAX_DUST_CHARGE + dust_charge
                    reaction = [const.ELECTRON, self.dust_array[s][ich], self.dust_array[s][ich-1], const.ASTERISK, const.ASTERISK, const.ASTERISK, str(dust_charge), str(s+1)]
                    ss = self.set_reaction_write(index_number, const.ELECTRON_GRAIN_COLLISION, reaction)
                    file.write(ss)
                    index_number += 1

         #if DUST_REACTION == True:
            """
            Grain-Grain Collison ##
            Ilgner & Nelson (2006) Table4 7-10
            """
            for s1 in range(NUMBER_OF_BINS):
                for positive_dust_charge in range(1, MAX_DUST_CHARGE+1):
                    ipch = MAX_DUST_CHARGE + positive_dust_charge
                    reac_grain1 = self.dust_array[s1][ipch]
                    for s2 in range(NUMBER_OF_BINS):
                        for negative_dust_charge in range(-1, -MAX_DUST_CHARGE-1, -1):
                            inch = MAX_DUST_CHARGE + negative_dust_charge
                            reac_grain2 = self.dust_array[s2][inch]

                            product_charge = positive_dust_charge + negative_dust_charge

                            if product_charge < 0:
                                prod_grain1 = self.dust_array[s1][MAX_DUST_CHARGE+0]
                                prod_grain2 = self.dust_array[s2][MAX_DUST_CHARGE+product_charge] 
                            elif product_charge == 0:
                                prod_grain1 = self.dust_array[s1][MAX_DUST_CHARGE+0]
                                prod_grain2 = self.dust_array[s2][MAX_DUST_CHARGE+0] 
                            else:
                                prod_grain1 = self.dust_array[s1][MAX_DUST_CHARGE+product_charge]
                                prod_grain2 = self.dust_array[s2][MAX_DUST_CHARGE+0] 

                            reaction = [reac_grain1, reac_grain2, prod_grain1, prod_grain2, const.ASTERISK, const.ASTERISK, str(positive_dust_charge), str(negative_dust_charge), str(s1+1), str(s2+1)]
                            ss = self.set_reaction_write(index_number, const.GRAIN_GRAIN_COLLISION, reaction)
                            file.write(ss)
                            index_number += 1

            if FLAG_H2_FORMATION_ON_GRAIN_SURFACE == True:
                """
                H2 formation on grain surface ##
                """
                mass = self.species_mass[self.hydrogen_mantle[0]]
                reaction = ["H", "H", "H2", const.ASTERISK, const.ASTERISK, const.ASTERISK, self.hydrogen_mantle[2], mass]
                s = self.set_reaction_write(index_number, const.H2_FORMATION_ON_GRAIN_SURFACE, reaction)
                file.write(s)
                index_number += 1   

        ## end if DUST_REACTION == True:            

        file.close()
        self.total_reactions = index_number - 1

    

class ChemicalReaction(ChemicalSpecies, ReadUmistDatabase, Dust):

    def __init__(self): 
        super(ChemicalReaction, self).__init__()
        self.set_dust_array()

    def input_species(self, input_species_file=None, input_species_list=None):

        if (input_species_file == None) and (input_species_list == None):
            print("no input file or list")
            sys.exit()
        
        if (not (input_species_file == None)) and (input_species_list == None):
            self.read_from_file(input_species_file)
        elif (input_species_file == None) and (not (input_species_list == None)):
            self.read_from_list(input_species_list)
        else:
            print("input file and list")
            sys.exit()

    def read_species_from_umist_databae(self):
        self.read_species(Species_File_From_UMIST )
        self.read_binding_energy(Binding_File_From_UMIST)
        self.read_reactions(Rate_File_From_UMIST)

    def write_file(self, species_file, reaction_file):
        self.species_file = species_file
        self.reaction_file = reaction_file
        self.read_species_from_umist_databae()

        ## OutPut fileの既存内容　削除##
        file = open(species_file, mode='w')
        file.close()
        file = open(reaction_file, mode='w')
        file.close()

        self.write_species_file(species_file)
        self.write_reaction_file(reaction_file)

if __name__ == "__main__":
    SPECIES_LIST = ["H", "H2", "He", "CO", "O2", "Mg", "O", "C", "HCO", "H2O", "OH", "H3+", "H2+", "H+", "HCO+", "Mg+", "He+", "C+", "O+", "CO+", "CH2+", "O2+", "H3O+", "OH+", "O2H+", "H2O+", "e-"]

    ## Input File ##
    ## UMIST database http://udfa.ajmarkwick.net/  (McElroy et al. 2013)
    Species_File_From_UMIST = "../UMIST_database/SPECIES12.txt"
    Rate_File_From_UMIST = "../UMIST_database/RATE12.dist.txt"
    Binding_File_From_UMIST = "../UMIST_database/RATE12_binding_energies.dist.txt"

    ## OutPut File ##
    Species_File = "species_test.txt"
    Rate_File = "reactions_test.txt"

    test = ChemicalReaction()

    test.input_species(input_species_list=SPECIES_LIST)
    test.read_species_from_umist_databae(Species_File_From_UMIST, Rate_File_From_UMIST, Binding_File_From_UMIST)
    test.write_file(Species_File, Rate_File)
    test.print_total()