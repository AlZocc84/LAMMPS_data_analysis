#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jun 13 2022

@author: alberto_z

Reads data from:
    packmol input (optional) or input file
    multi-fragment PDB file
    lammps data file corresponding to the structure files in the  multi-fragment PDB.
Outputs:
    Lammps data file corresponding to multi-fragment pdb.
"""

#%% Basic imports
from os import getcwd
from copy import deepcopy
import itertools
import networkit as nk
from pprint import pprint
from datetime import datetime
from math import acos,degrees

#%%
now = datetime.now()
gdate = now.strftime("%d-%b-%y")

#%% Auxiliary functions

def map_mass(a_mass):
    """Maps atomic masses (in amu) to the corresponding atomic symbol."""
    atomic_masses={'1.0': 'H', '4.0': 'He', '6.9': 'Li', '9.0': 'Be', '10.8': 'B', '12.0': 'C', '14.0': 'N', '16.0': 'O', '19.0': 'F', '20.2': 'Ne', '23.0': 'Na', '24.3': 'Mg', '27.0': 'Al', '28.1': 'Si', '31.0': 'P', '32.1': 'S', '35.5': 'Cl', '39.9': 'Ar', '39.1': 'K', '40.1': 'Ca', '45.0': 'Sc', '47.9': 'Ti', '50.9': 'V', '52.0': 'Cr', '54.9': 'Mn', '55.8': 'Fe', '58.9': 'Co', '58.7': 'Ni', '63.5': 'Cu', '65.4': 'Zn', '69.7': 'Ga', '72.6': 'Ge', '74.9': 'As', '79.0': 'Se', '79.9': 'Br', '83.8': 'Kr', '85.5': 'Rb', '87.6': 'Sr', '88.9': 'Y', '91.2': 'Zr', '92.9': 'Nb', '96.0': 'Mo', '98.0': 'Tc', '101.1': 'Ru', '102.9': 'Rh', '106.4': 'Pd', '107.9': 'Ag', '112.4': 'Cd', '114.8': 'In', '118.7': 'Sn', '121.8': 'Sb', '127.6': 'Te', '126.9': 'I', '131.3': 'Xe', '132.9': 'Cs', '137.3': 'Ba', '138.9': 'La', '140.1': 'Ce', '140.9': 'Pr', '144.2': 'Nd', '145.0': 'Pm', '150.4': 'Sm', '152.0': 'Eu', '157.3': 'Gd', '158.9': 'Tb', '162.5': 'Dy', '164.9': 'Ho', '167.3': 'Er', '168.9': 'Tm', '173.1': 'Yb', '175.0': 'Lu', '178.5': 'Hf', '180.9': 'Ta', '183.8': 'W', '186.2': 'Re', '190.2': 'Os', '192.2': 'Ir', '195.1': 'Pt', '197.0': 'Au', '200.6': 'Hg', '204.4': 'Tl', '207.2': 'Pb', '209.0': 'Bi', '210.0': 'At', '222.0': 'Rn', '223.0': 'Fr', '226.0': 'Ra', '227.0': 'Ac', '232.0': 'Th', '231.0': 'Pa', '238.0': 'U', '237.0': 'Np', '244.0': 'Pu', '243.0': 'Am', '247.0': 'Bk', '251.0': 'Cf', '252.0': 'Es', '257.0': 'Fm', '258.0': 'Md', '259.0': 'No', '262.0': 'Db', '261.0': 'Rf', '266.0': 'Sg', '264.0': 'Bh', '267.0': 'Hs', '268.0': 'Mt', '271.0': 'Ds ', '272.0': 'Rg ', '285.0': 'Cn ', '284.0': 'Nh', '289.0': 'Fl', '288.0': 'Mc', '292.0': 'Lv', '295.0': 'Ts', '294.0': 'Og'}
    n_t =  atomic_masses[f"{a_mass:3.1f}".strip()]
    
    return n_t

def  remove_empty_lines(a_file_content):
    """"removes empty lines from the list a_file_content"""
    o_short_file = []
    for line in a_file_content:
        if len(line.split()) != 0:
            o_short_file.append(line)
    
    #with open("short.out","a") as out_f:
    #    out_f.writelines(o_short_file)
    #    out_f.close()
    return o_short_file



def get_sec_name(a_str):
    """Input: Lammps Data file line, returns the name of the section (without comments)
    assumes the line to already have been identified as a section header.""" 
    sec_attr_name =""
    for i,word in enumerate(a_str.split()):
        if i ==0 :
            sec_attr_name+=word
        elif i=="#":
            break
        else:
            sec_attr_name +="_"+word
    return sec_attr_name

def find_map_root_id(a_fragment_map): # Getting the structure of a fragment mao from below.
    """Finds root atom's ID from a list of dictionaries representing the system connectivity graph (produced by fragment_map)"""
    root_at_id = False
    for siblings in a_fragment_map[0].values(): # Actually useless, as there is only one.
        for child_data in siblings:
            root_at_id = child_data[0]
    if root_at_id != False:
        return root_at_id

def find_map_root_type(a_fragment_map): # Getting the structure of a fragment mao from below.
    """Finds root atom's type from a list of dictionaries representing the system connectivity graph (produced by fragment_map)"""
    root_at_type = False
    for siblings in a_fragment_map[0].values(): 
        for child_data in siblings:
            root_at_type = child_data[1]
    if root_at_type != False:
        return root_at_type

def read_atoms_dict(a_lt_sys_file):
    """Reads atom type strings from a .lt file
input: a_lt_sys_file: a .lt file to be used to read the type string for each atom.
return: out_at_t_dict: a dictionary of at_id:at_type_string. 
WARNING: this is a legacy function that has not been used for a long time. Check before use."""
    out_at_t_dict = {}
    with open(a_lt_sys_file) as lt_f:
        lines = lt_f.readlines()
        index = 0
        for i, line in enumerate(lines):
             if "write" in line and "Data Atoms" in line:
                 for line2 in lines[i+1:]:
                     if line2.strip()[0] == "}":
                         break
                     else:
                         atom_type = line2.split()[2][6:] # i.e. HS14, CAro
                         out_at_t_dict[index+1] = atom_type  
                         print(f"reading atom type {atom_type} at index {index+1}")
                         index += 1
                       
    return out_at_t_dict
 
def read_types_map(a_type_map_file):
    """WARNING: this is a legacy function that has not been used for a long time. Check before use."""
    t_types=["atoms", "bonds", "angles","dihedrals","impropers"]
    out_map_dict = {}
    with open(a_type_map_file) as type_f:
        lines = type_f.readlines()
        current=""
        for i, line in enumerate(lines):
                
             if line.split()[0] in t_types:
                 current=line.split()[0]
                 out_map_dict[current]={}
             elif len(line.split())==2:
                 out_map_dict[current][int(line.split()[0])]=int(line.split()[1])
                 
    return out_map_dict   

def read_types_n(a_n_types_file):
    """WARNING: this is a legacy function that has not been used for a long time. Check before use."""
    out_types_n_dict = {}
    with open(a_n_types_file) as type_f:
        lines = type_f.readlines()
        for line in lines:
            out_types_n_dict[f"n_{line.split()[0]}_types"] = int(line.split()[1])
    return out_types_n_dict    

#%% LammpsData

class LammpsData:
    """Class to store and manipulate information pertaining to a LAMMPS .data file. Currently assuming atom type full"""
    def __init__(self,a_lammps_data_f, is_precise=False,renumber=False,reorder=False,reassign_mols=False,name=None):
        
        if name==None:
            self.name = str(a_lammps_data_f.split(".")[0])
        else:
            self.name=name
        self.file = a_lammps_data_f
        self.ff_terms =[]
        self.precise=is_precise
            
        sections_headers = ["Pair Coeffs", "PairIJ Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs","Improper Coeffs","Masses","Atoms","Velocities","Bonds","Angles","Dihedrals","Impropers"]
        # in_section={ sec_h:False for sec_h in sections_headers}
        terms = ["Bond","Angle","Dihedral","Improper"]

        with open(a_lammps_data_f) as lammps_f:
            #skip first line
            lammps_f.readline()
            while True:
                line=lammps_f.readline()
                #print("New line: ", line)
                if not line:
                    break
    
                for sec_h in sections_headers:
                    if sec_h in line:
                        current_sec=sec_h
                        #print("Defining section: ", current_sec)
                        sec_attr_name = get_sec_name(sec_h)
                        setattr(self,sec_attr_name,{})
                        if sec_h.split()[0] in terms:
                            self.ff_terms.append(sec_h.split()[0])
                        if sec_h =="Masses" :
                            self.masses_type_labels ={}
                            self.type_to_element = {}
                        if sec_h == "Atoms":
                            self.atoms_labels = {}
                            self.atoms_per_mol = {}
                        if sec_h == "Bonds":
                            self.bonds_per_at = {}
                        line=lammps_f.readline()
                        if not line:
                            break
                            
                if len(line.split()) == 0:
                    continue
                elif "atoms" in line:
                    self.n_atoms = int(line.split()[0])
                elif "bonds" in line:
                    self.n_bonds = int(line.split()[0])
                elif "angles" in line:
                    self.n_angles = int(line.split()[0])
                elif "dihedrals" in line:
                    self.n_dihedrals = int(line.split()[0])
                elif "impropers" in line:
                    self.n_impropers = int(line.split()[0])    
                elif (len(line.split()) > 2) and ( line.split()[2] == "types"):
                    if line.split()[1] == "atom":
                        self.n_atom_types = int(line.split()[0])
                    elif line.split()[1] == "bond":
                        self.n_bond_types = int(line.split()[0])
                    elif line.split()[1] == "angle":
                        self.n_angle_types = int(line.split()[0])
                    elif line.split()[1] == "dihedral":
                        self.n_dihedral_types = int(line.split()[0])
                    elif line.split()[1] == "improper":
                        self.n_improper_types = int(line.split()[0]) 
                elif "lo" in line and "hi" in line:
                    if "xlo" in line:
                        self.xbounds = [float(x) for x in line.split()[:2]]
                    elif "ylo" in line:
                        self.ybounds = [float(x) for x in line.split()[:2]]
                    elif "zlo" in line:
                        self.zbounds = [float(x) for x in line.split()[:2]]
                elif "xy" in line:
                        self.cell_angles = [float(alpha) for alpha in line.split()[:3]]
                
                elif current_sec=="Pair Coeffs":
                    getattr(self,get_sec_name("Pair Coeffs"))[int(line.split()[0])]= line.split()[1:]

                
                elif current_sec=="PairIJ Coeffs":
                    getattr(self, get_sec_name("PairIJ coeffs"))[(int(line.split()[0]),int(line.split()[1]))] = line.split()[2:]

                elif current_sec.split()[0] in terms:                    
                    if (len(line.split("#")[0])>0 and len(line.split("#")[0].split())>2):
                        this_term_dict={}
                        if "#" in line:
                            this_term_dict["coeffs"] = line.split("#")[0].split()[1:]
                            this_term_dict["label"] = line.split("#")[1]
                        else:
                            this_term_dict["coeffs"] = line.split()[1:]

                        getattr(self, get_sec_name(current_sec))[int(line.split()[0])]= this_term_dict

                elif current_sec =="Masses":
                    this_type=int(line.split()[0])
                    getattr(self, get_sec_name(current_sec))[this_type] = float(line.split()[1])
                    try:
                        self.type_to_element[this_type]=map_mass(float(line.split()[1]))
                    except:
                        self.type_to_element[this_type]="?"
                        print(f"WARNING!! Could not define element type from mass {float(line.split()[1])}, where elements labels are needed, e.g.pdb files, questions marks will be used instead")
                    if "#" in line:
                        self.masses_type_labels[this_type] = line.split("#")[-1][:-1] # Expects mass lines to be "id mass # type"
                    else:
                        self.masses_type_labels[this_type] = self.type_to_element[this_type]
                        
                elif current_sec == "Atoms":
                    this_atom_dict={}
                    this_atom_dict["mol"] = int(line.split()[1])
                    this_atom_dict["type"] = int(line.split()[2])
                    this_atom_dict["charge"] = float(line.split()[3])
                    this_atom_dict["posX"] = float(line.split()[4])
                    this_atom_dict["posY"] = float(line.split()[5])
                    this_atom_dict["posZ"] = float(line.split()[6])
                    if len(line)==10:
                        this_atom_dict["trueX"] = int(line.split()[7])
                        this_atom_dict["trueY"] = int(line.split()[8])
                        this_atom_dict["trueZ"] = int(line.split()[9])
                        
                    if "#" in line: # Assumes atom lables are after #.
                        self.atoms_labels[line.split()[0]] = line.split("#")[1].split()[0]
                        # print(the_line.split("#"))
                        pass
                    else:
                        self.atoms_labels[line.split()[0]] = ""
                    if this_atom_dict["mol"] in self.atoms_per_mol.keys():
                        self.atoms_per_mol[this_atom_dict["mol"]].append(int(line.split()[0]))
                    else:
                        self.atoms_per_mol[this_atom_dict["mol"]] = [int(line.split()[0])]

                    getattr(self,get_sec_name(current_sec))[int(line.split()[0])] = this_atom_dict
                       
                elif current_sec == "Velocities":
                    this_vel_dict={}
                    this_vel_dict["velX"] = float(line.split()[1])
                    this_vel_dict["velY"] = float(line.split()[2])
                    this_vel_dict["velZ"] = float(line.split()[3])
                    getattr(self,get_sec_name(current_sec))[int(line.split()[0])] = this_vel_dict
                    
                elif current_sec == "Bonds":
                    this_bond_dict={}
                    this_bond_dict["type"] = int(line.split()[1])
                    this_bond_dict["atom1"] = int(line.split()[2])
                    this_bond_dict["atom2"] = int(line.split()[3])               
                    getattr(self,get_sec_name(current_sec))[line.split()[0]] = this_bond_dict
                    
                    if int(line.split()[2]) in  self.bonds_per_at.keys():
                        self.bonds_per_at[int(line.split()[2])].append(int(line.split()[3]))
                    else:
                        self.bonds_per_at[int(line.split()[2])]=[int(line.split()[3])]

                    if int(line.split()[3]) in  self.bonds_per_at.keys():
                        self.bonds_per_at[int(line.split()[3])].append(int(line.split()[2]))
                    else:
                        self.bonds_per_at[int(line.split()[3])]=[int(line.split()[2])]
                            
                elif current_sec == "Angles":
                    this_angle_dict={}
                    this_angle_dict["type"] = int(line.split()[1])
                    this_angle_dict["atom1"] = int(line.split()[2])
                    this_angle_dict["atom2"] = int(line.split()[3])
                    this_angle_dict["atom3"] = int(line.split()[4])
                    getattr(self,get_sec_name(current_sec))[line.split()[0]] = this_angle_dict   
                    
                elif current_sec == "Dihedrals":
                    this_dihedral_dict={}
                    this_dihedral_dict["type"] = int(line.split()[1])
                    this_dihedral_dict["atom1"] = int(line.split()[2])
                    this_dihedral_dict["atom2"] = int(line.split()[3])
                    this_dihedral_dict["atom3"] = int(line.split()[4])
                    this_dihedral_dict["atom4"] = int(line.split()[5])
                    getattr(self,get_sec_name(current_sec))[line.split()[0]] = this_dihedral_dict
                    
                elif current_sec == "Impropers":
                    this_improper_dict={}
                    this_improper_dict["type"] = int(line.split()[1])
                    this_improper_dict["atom1"] = int(line.split()[2])
                    this_improper_dict["atom2"] = int(line.split()[3])
                    this_improper_dict["atom3"] = int(line.split()[4])
                    this_improper_dict["atom4"] = int(line.split()[5])
                    getattr(self,get_sec_name(current_sec))[line.split()[0]] = this_improper_dict
                    
                # elif line.strip()[1] != "#":
                    # print("At line: \n", line)
                    # exit("Content of this line does not match any section.")
                
                    
                    
                
        if renumber:
            self.renumber()
        if reassign_mols:
            self.reassign_molecules()
        if reorder:
            self.reorder()
        
# ------------ Writing to files ----------------------------------------------------------------------------------------#            
    def write_molecule_template(self, a_out_lammps_mol_file):
        """Writes the Fragment info to a lammps molecule template file."""
        
        #Headers
        lammps_mol_str = "This is a molecule template # From lammps_data_analysis.py\n\n"
        lammps_mol_str += "{n_at}\t atoms\n".format(n_at=self.n_atoms)
        lammps_mol_str += "{n_b}\t bonds\n".format(n_b=self.n_bonds)
        lammps_mol_str += "{n_ang}\t angles\n".format(n_ang=self.n_angles)
        lammps_mol_str += "{n_dh}\t dihedrals\n".format(n_dh=self.n_dihedrals)
        lammps_mol_str += "{n_imp}\t impropers\n".format(n_imp=self.n_impropers)
        lammps_mol_str += "\n"
   
        #Types
        lammps_mol_str += "Types\n\n"

        for at_id,at_info in self.Atoms.items():
            
            lammps_mol_str += "{}".format(at_id)
            lammps_mol_str += "\t{}\n".format(at_info["type"])
         
        lammps_mol_str += "\n"
        
        #Charges
        lammps_mol_str += "Charges\n\n"

        for at_id,at_info in self.Atoms.items():
            
            lammps_mol_str += "{}".format(at_id)
            lammps_mol_str += "\t{:.3f}\n".format(at_info["charge"])
        lammps_mol_str += "\n"
    
        #Coords
        lammps_mol_str += "Coords\n\n"
   
        for at_id,at_info in self.Atoms.items():
            
            lammps_mol_str += "{}".format(at_id)
            lammps_mol_str += "\t{:.7f}".format(at_info["posX"])
            lammps_mol_str += "\t{:.7f}".format(at_info["posY"])
            lammps_mol_str += "\t{:.7f}\n".format(at_info["posZ"])
        lammps_mol_str += "\n"
        
        #Bonds
        lammps_mol_str += "Bonds\n\n"
   
        for bd_id, bd_info in self.Bonds.items():
            lammps_mol_str += "{}".format(bd_id)+"\t"
            lammps_mol_str += "{}".format(bd_info["type"])+"\t"
            lammps_mol_str += "{}".format(bd_info["atom1"])+"\t"
            lammps_mol_str += "{}".format(bd_info["atom2"])
            if "atoms_types" in bd_info.keys():
                lammps_mol_str += f' # {bd_info["atoms_types"]}\n'
            else:
                lammps_mol_str += "\n"
   
        lammps_mol_str += "\n"
        
        #Angles
        lammps_mol_str += "Angles\n\n"
   
        for an_id,an_info in self.Angles.items():
            lammps_mol_str += "{}".format(an_id)+"\t"
            lammps_mol_str += "{}".format(an_info["type"])+"\t"
            lammps_mol_str += "{}".format(an_info["atom1"])+"\t"
            lammps_mol_str += "{}".format(an_info["atom2"])+"\t"
            lammps_mol_str += "{}".format(an_info["atom3"])
            if "atoms_types" in an_info.keys():
                lammps_mol_str += f' # {an_info["atoms_types"]}\n'
            else:
                lammps_mol_str += "\n"
        lammps_mol_str += "\n"
        
        #Dihedrals
        lammps_mol_str += "Dihedrals\n\n"
   
        for dh_id, dh_info in self.Dihedrals.items():
            lammps_mol_str += "{}".format(dh_id)+"\t"
            lammps_mol_str += "{}".format(dh_info["type"])+"\t"
            lammps_mol_str += "{}".format(dh_info["atom1"])+"\t"
            lammps_mol_str += "{}".format(dh_info["atom2"])+"\t"
            lammps_mol_str += "{}".format(dh_info["atom3"])+"\t"
            lammps_mol_str += "{}".format(dh_info["atom4"])
            if "atoms_types" in dh_info.keys():
                lammps_mol_str += f' # {dh_info["atoms_types"]}\n'
            else:
                lammps_mol_str += "\n"
            
        lammps_mol_str += "\n"
        
        #Impropers
        lammps_mol_str += "Impropers\n\n"
   
        for imp_id, imp_info in self.Impropers.items():
            lammps_mol_str += "{}".format(imp_id)+"\t"
            lammps_mol_str += "{}".format(imp_info["type"])+"\t"
            lammps_mol_str += "{}".format(imp_info["atom1"])+"\t"
            lammps_mol_str += "{}".format(imp_info["atom2"])+"\t"
            lammps_mol_str += "{}".format(imp_info["atom3"])+"\t"
            lammps_mol_str += "{}".format(imp_info["atom4"])
            if "atoms_types" in imp_info.keys():
                lammps_mol_str += f' # {imp_info["atoms_types"]}\n'
            else:
                lammps_mol_str += "\n"
            
        lammps_mol_str += "\n"
        
        
        with open(a_out_lammps_mol_file,"w") as out_f:
            out_f.write(lammps_mol_str)
            out_f.close()
            
    def write_lammps_data(self, a_out_lammps_data_file, coeffs_out=False):
        """Writes the Fragment info to a lammps data file."""
        print("Writing output: ",a_out_lammps_data_file)
        #Headers
        lammps_data_str = "LAMMPS Description # From lammps_data_analysis.py\n\n"
        lammps_data_str += "{n_at}  atoms\n".format(n_at=self.n_atoms)
        lammps_data_str += "{n_b}  bonds\n".format(n_b=self.n_bonds)
        lammps_data_str += "{n_ang}  angles\n".format(n_ang=self.n_angles)
        lammps_data_str += "{n_dh}  dihedrals\n".format(n_dh=self.n_dihedrals)
        lammps_data_str += "{n_imp}  impropers\n".format(n_imp=self.n_impropers)
        lammps_data_str += "\n"

        #Types
        lammps_data_str += "{n_at}  atom types\n".format(n_at=self.n_atom_types)
        lammps_data_str += "{n_b}  bond types\n".format(n_b=self.n_bond_types)
        lammps_data_str += "{n_ang}  angle types\n".format(n_ang=self.n_angle_types)
        lammps_data_str += "{n_dh}  dihedral types\n".format(n_dh=self.n_dihedral_types)
        lammps_data_str += "{n_imp}  improper types\n".format(n_imp=self.n_improper_types)
        lammps_data_str += "\n"

        #Box
        lammps_data_str += "{:.3f}".format(self.xbounds[0])+" "+"{:.3f}".format(self.xbounds[1])+" xlo xhi\n"
        lammps_data_str += "{:.3f}".format(self.ybounds[0])+" "+"{:.3f}".format(self.ybounds[1])+" ylo yhi\n"
        lammps_data_str += "{:.3f}".format(self.zbounds[0])+" "+"{:.3f}".format(self.zbounds[1])+" zlo zhi\n"
        try:
            lammps_data_str += "{:.3f}".format(self.cell_angles[0])+" "+"{:.3f}".format(self.cell_angles[1])+" "+"{:.3f}".format(self.cell_angles[2])+" xy xz yz\n"
        except AttributeError:
            pass
        lammps_data_str += "\n"

        #Masses
        lammps_data_str += "Masses\n\n"
        if hasattr(self,"Masses"):
            for at_type in sorted([int(i)for i in self.Masses.keys()]):
                lammps_data_str += "{}".format(at_type)+" "+"{:.4f}".format(self.Masses[at_type])
                if  len("".join(self.masses_type_labels.values())) != 0:
                    lammps_data_str += " # {}".format(self.masses_type_labels[at_type])+"\n"
                else:
                    lammps_data_str += "\n"
            lammps_data_str += "\n"
        
        if coeffs_out == True:
            if hasattr(self,"Pair_Coeffs"):
                #Pair Coeffs
                lammps_data_str += "Pair Coeffs\n\n"
                for p_type in sorted(self.Pair_Coeffs.keys()):
                    p_Coeffs=self.Pair_Coeffs[p_type]
                    lammps_data_str += "{}".format(p_type)
                    for c in p_Coeffs:
                        lammps_data_str += " {}".format(c)
                    if hasattr(self,"atoms_names_per_type"):
                        lammps_data_str += " # {}".format(self.atoms_names_per_type[p_type])

                    lammps_data_str += "\n"
                lammps_data_str += "\n"
            if hasattr(self,"PairIJ_Coeffs"):
                #Pair Coeffs
                lammps_data_str += "PairIJ Coeffs\n\n"
                # print(self.PairIJ_Coeffs.items())
                for p_key,p_Coeffs in sorted(self.PairIJ_Coeffs.items(), key= lambda x: (x[0][0],x[0][1])):
                    # print(p_key)
                    lammps_data_str += "{}".format(p_key[0])
                    lammps_data_str += " {}".format(p_key[1])

                    for c in p_Coeffs:
                        lammps_data_str += " {}".format(c)
                    if hasattr(self,"atoms_names_per_type"):
                        lammps_data_str += " # {}".format(self.atoms_names_per_type[int(p_key[0])])
                        lammps_data_str += ",{}".format(self.atoms_names_per_type[int(p_key[1])])
                    lammps_data_str += "\n"
                lammps_data_str += "\n"   
            
            #Other Coeffs
            for term in self.ff_terms:
                if term != "Pair":
                    lammps_data_str += f"{term} Coeffs\n\n"
                    for t_type in sorted([int(i) for i in getattr(self,term+"_Coeffs").keys()]):
                        t_Coeffs=getattr(self,term+"_Coeffs")[t_type]
                        lammps_data_str += "{}".format(t_type)
                        for t_c in t_Coeffs["coeffs"]:
                            lammps_data_str += " {}".format(t_c)
                        if "label" in t_Coeffs.keys():
                            lammps_data_str += f"# {t_Coeffs['label']}"
                        lammps_data_str += "\n"
                    lammps_data_str += "\n"  
                    
        if len(self.ff_terms) == 0:
            print("Could not find any coeffs. Skipping them.")
    
        #Atoms
        lammps_data_str += "Atoms\n\n"

        for at_id,at_info in self.Atoms.items():
            
            lammps_data_str += "{}".format(at_id)
            lammps_data_str += " {}".format(at_info["mol"])
            lammps_data_str += " {}".format(at_info["type"])
            lammps_data_str += " {:.3f}".format(at_info["charge"])
            lammps_data_str += " {:.7f}".format(at_info["posX"])
            lammps_data_str += " {:.7f}".format(at_info["posY"])
            lammps_data_str += " {:.7f}".format(at_info["posZ"])
            #TODO: check why following is not working 
            if "trueX" in self.Atoms.keys():
                lammps_data_str += " {}".format(at_info["trueX"])
                lammps_data_str += " {}".format(at_info["trueY"])
                lammps_data_str += " {}".format(at_info["trueZ"])
            try:
                lammps_data_str += f" # {self.atoms_labels[at_id]}\n"
            except:
                lammps_data_str += "\n"
            else:
                lammps_data_str += ""
        lammps_data_str += "\n"
        
        if hasattr(self,"Velocities"):
            #Velocities
            lammps_data_str += "Velocities\n\n"
            
            for at_id,at_info in self.Velocities.items():
                
                lammps_data_str += "{}".format(at_id)
                lammps_data_str += " {:.7f}".format(at_info["velX"])
                lammps_data_str += " {:.7f}".format(at_info["velY"])
                lammps_data_str += " {:.7f}\n".format(at_info["velZ"])
            lammps_data_str += "\n"
        
        #Bonds
        lammps_data_str += "Bonds\n\n"

        for bd_id, bd_info in self.Bonds.items():
            lammps_data_str += "{}".format(bd_id)+" "
            lammps_data_str += "{}".format(bd_info["type"])+" "
            lammps_data_str += "{}".format(bd_info["atom1"])+" "
            lammps_data_str += "{}".format(bd_info["atom2"])
            if "atoms_types" in bd_info.keys():
                lammps_data_str += f' # {bd_info["atoms_types"]}\n'
            else:
                lammps_data_str += "\n"

        lammps_data_str += "\n"
        
        #Angles
        lammps_data_str += "Angles\n\n"

        for an_id,an_info in self.Angles.items():
            lammps_data_str += "{}".format(an_id)+" "
            lammps_data_str += "{}".format(an_info["type"])+" "
            lammps_data_str += "{}".format(an_info["atom1"])+" "
            lammps_data_str += "{}".format(an_info["atom2"])+" "
            lammps_data_str += "{}".format(an_info["atom3"])
            if "atoms_types" in an_info.keys():
                lammps_data_str += f' # {an_info["atoms_types"]}\n'
            else:
                lammps_data_str += "\n"
        lammps_data_str += "\n"
        
        #Dihedrals
        lammps_data_str += "Dihedrals\n\n"

        for dh_id, dh_info in self.Dihedrals.items():
            lammps_data_str += "{}".format(dh_id)+" "
            lammps_data_str += "{}".format(dh_info["type"])+" "
            lammps_data_str += "{}".format(dh_info["atom1"])+" "
            lammps_data_str += "{}".format(dh_info["atom2"])+" "
            lammps_data_str += "{}".format(dh_info["atom3"])+" "
            lammps_data_str += "{}".format(dh_info["atom4"])
            if "atoms_types" in dh_info.keys():
                lammps_data_str += f' # {dh_info["atoms_types"]}\n'
            else:
                lammps_data_str += "\n"
            
        lammps_data_str += "\n"
        
        #Impropers
        lammps_data_str += "Impropers\n\n"

        for imp_id, imp_info in self.Impropers.items():
            lammps_data_str += "{}".format(imp_id)+" "
            lammps_data_str += "{}".format(imp_info["type"])+" "
            lammps_data_str += "{}".format(imp_info["atom1"])+" "
            lammps_data_str += "{}".format(imp_info["atom2"])+" "
            lammps_data_str += "{}".format(imp_info["atom3"])+" "
            lammps_data_str += "{}".format(imp_info["atom4"])
            if "atoms_types" in imp_info.keys():
                lammps_data_str += f' # {imp_info["atoms_types"]}\n'
            else:
                lammps_data_str += "\n"
            
        lammps_data_str += "\n"
        
        
        with open(a_out_lammps_data_file,"w") as out_f:
            out_f.write(lammps_data_str)
            out_f.close()
            
        
    
#     @classmethod        
#     def extend_fragment(cls,a_FG, a_S,a_name = None):
#         """This function uses the information in a Synonyms object to extend a LammpsData object by adding more Atom types and duplicating all Force field terms lines to describe all possible compinations..
#         \n a_S: a Synonims object with info on atom types to be added.
#         \n a_name: name of the new combined fragment object.
#         """
        
#         #First get a new empty LammpsData obj to fill later.
#         ext_frag_obj = cls.__new__(cls)
        
       
#         #Now giving a name to the new LammpsData object
#         if a_name != None:
#             name_str = a_name
#         else:
#             name_str= "extension"
    
#         file_str = name_str + ".data"
            
#         ext_frag_obj.name = name_str
#         ext_frag_obj.file = file_str
        
        
#         #Now computing the new LammpsData attributes
        
#         #start with setting all attributes
#         for frag_attr in a_FG.__dict__.keys():
#             n_q_tot = getattr(a_FG,frag_attr)
#             setattr(ext_frag_obj ,frag_attr,n_q_tot)
            
#         # Now setting all the ones that need to be changed
        
#         # now masses
        
#         masses = {}
#         masses_type_labels = {}
        
#         # Now we extend the masses dict
                
#         n_new_at_t = a_FG.n_atom_types
#         for key,m_type in a_FG.masses_type_labels.items():
#             mass = a_FG.masses[key]
#             masses[key] = mass
#             masses_type_labels[key] = m_type 
#             if len(a_S.synonims[m_type]) > 1:
#                 syns = a_S.synonims[m_type]
#                 for syn in syns[1:]:
#                     n_new_at_t +=1
#                     new_key = str(int(n_new_at_t))
#                     masses[new_key] = mass
#                     masses_type_labels[new_key] = syn 

                
#         ext_frag_obj.masses = masses
#         ext_frag_obj.masses_type_labels = masses_type_labels
#         ext_frag_obj.n_atom_types = len(masses)
        
#         # now Pair Coeffs
        
#         pair_Coeffs = {}
        
#         # Now we extend the masses dict
        
#         n_new_at_t = a_FG.n_atom_types
#         for key,m_type in a_FG.masses_type_labels.items():
#             p_Coeffs = a_FG.Pair_Coeffs[key]
#             pair_Coeffs[key]=p_Coeffs
#             if len(a_S.synonims[m_type]) > 1:
#                 syns = a_S.synonims[m_type]
#                 for syn in syns[1:]:
#                     n_new_at_t +=1
#                     new_key = str(int(n_new_at_t))
#                     pair_Coeffs[new_key] = p_Coeffs
                
#         ext_frag_obj.Pair_Coeffs = pair_Coeffs

# #         # Now other coeffs
         
#         ffs_ls=[]
#         for i in a_FG.__dict__.keys():
#             if "_Coeffs" in i and i != "Pair_Coeffs":
#                 ffs_ls.append(i)
        
#         for j in ffs_ls:
#             setattr(ext_frag_obj, j, {})
#             n_terms = len(getattr(a_FG,j))
#             for key,val in getattr(a_FG,j).items():
#                 getattr(ext_frag_obj, j)[key] = val
#                 term_ats = val[-1].split(",")
#                 syns_term_ats = [a_S.synonims[q] for q in term_ats]
#                 full_term_ats=itertools.product(*syns_term_ats)
#                 for s in full_term_ats:
#                     if list(s) !=term_ats:
#                         n_terms +=1
#                         new_key = str(int(n_terms))
#                         new_term= ','.join(s)
#                         new_val= val[:-1]+[new_term]
#                         getattr(ext_frag_obj, j)[new_key] = new_val
#                     else:
#                         continue
#                 # add extra lines from input file, add right name at the end˚
            
        
# #         # now quantity_type type of attribute
        
        
#         for frag_attr in ["bond_","angle_","dihedral_","improper_"]:
#             attr_str=f"n_{frag_attr}types"
#             coeffs_str=frag_attr.capitalize()+"coeffs"
#             setattr(ext_frag_obj,attr_str,len(getattr(ext_frag_obj, coeffs_str)))
        

        
# #         # Now atoms
        
#         # atoms = {}
#         # at_offset = 0
#         # mol_id = 1 # 
#         # at_offset_dict = {}
#         # at_offset_dict[(mol_id)] = 0
#         # for i,frag_info in enumerate(frag_infos):
#         #     n_frag = frag_numbers[i]
#         #     for j in range(n_frag):
#         #         for key,atom in frag_info.atoms.items():
#         #             new_key = str(int(key) + at_offset)
#         #             # For now leaving the rest unchanged.
#         #             new_atom=copy.deepcopy(atom)
#         #             new_atom["type"] = atom["type"] + offsets["n_atom_types"][i]
#         #             new_atom["mol"] =  (mol_id) # -1 + atom["mol"]  # molecule_id_offset = mol_id -1. Note: this does not work properly if one of the lammps input has more than one molecule.
#         #             atoms[new_key] = new_atom
#         #         #print("now atoms is \n",atoms)
#         #         mol_id += 1
#         #         at_offset_dict[(mol_id)] = at_offset_dict[(mol_id-1)] + frag_info.n_atoms
#         #         at_offset += frag_info.n_atoms
                
#         #         #print (mol_offset,at_offset)
#         # ext_frag_obj.atoms=atoms
        
#         # print("at_offset_dict is :\n", at_offset_dict)
            
# #         #Now bonds#Can this be generalised for the others too ?
                
# #         offsets = {}        
# #         for frag_attr in ["bonds","angles","dihedrals","impropers"]:
# #             attr ={}
# #             attr_offs = offsets["n_" + frag_attr[:-1] + "_types"]
# #             mol_id = 1 # 

# #             for i,frag_info in enumerate(frag_infos):
# #                 n_frag = frag_numbers[i]
# #                 frag_offs = attr_offs[i]
# #                 for j in range(n_frag):
# #                     for key,an_attr in getattr(frag_info,frag_attr).items():
# #                           new_key = str(int(key) + frag_offs + j*getattr(frag_info,"n_"+frag_attr))
# #                           new_attr = {}
# #                           for attr_key, attr_val in an_attr.items():
# #                               if attr_key == "type":
# #                                   new_attr["type"] = an_attr["type"] + offsets["n_"+frag_attr][i]
# #                               if "atom" in attr_key:
# #                                   new_attr[attr_key] = an_attr[attr_key] + at_offset_dict[(mol_id)]
                             
# #                     attr[new_key] = new_attr
# #                     mol_id += 1 
# #             setattr(ext_frag_obj,frag_attr,attr)
        
#         return ext_frag_obj
     

    def write_types(self):
        '''Writes types.txt file for polymatic. Assumes masses, bonds, etc. have a last column with the type (e.g C-OEopt)'''
        types_data_str=""
        
        types_attr=["atom types", "bond types", "angle types", "dihedral types", "improper types"]
        
        for a_type in types_attr:
            types_data_str+=a_type+"\n"
            if a_type=="atom types":
                if len(self.masses_type_labels) ==0:
                    exit(f"Atom type label for type {a_type} could not be found as comment in the mass section. The types file cannot be completed.")
                for at_t in sorted([int(i) for i in  self.masses_type_labels.keys()]):
                    at_name = self.masses_type_labels[at_t]
                    types_data_str+=f"{at_t} {at_name}\n"
            else:
                try:
                    attr_coeff = a_type.split()[0].capitalize()+"_Coeffs"
                    for attr_t in sorted([int(i) for i in getattr(self,attr_coeff).keys()]):
                        attr_info = getattr(self,attr_coeff)[attr_t]
                        try:
                            attr_ats = attr_info["label"]
                        except:
                            exit(f"No info found on how to match atom types to {a_type.split()[0]}. The types file cannot be completed.")
                        types_data_str+=f"{attr_t} {attr_ats}\n"
                except:
                    print(f"Could not find data for {attr_coeff}. Skipping to next one.")
        
        with open("types.txt","w") as out_f:
            out_f.write(types_data_str)
            out_f.close()
        return

    def write_cassandra_charges(self):
        '''Output a charges file to be used with Cassandra'''
        if self.name == None:
            charges_f = "charges.ff"
        else:
            charges_f = f"{self.name}_charges.ff"
            
        with open(charges_f,"w") as chrg_f:
            tot_charge = 0.0
            for i,at in self.Atoms.items():
                #chrg_f.write("charge\n"+f"{i} {at['charge']}"+"\n\n")
                tot_charge += at['charge']
            print(f"Total charge is {tot_charge}")
            tot_charge_2 = 0.0
            for i,at in self.Atoms.items():
                chrg = at['charge'] - (tot_charge/self.n_atoms) # adjusted for total charge
                chrg_f.write("charge\n"+f"{i} {chrg}"+"\n\n")
                tot_charge_2 += chrg
                
            print(f"Total charge is {tot_charge_2}")
            
    def write_xyz(self,a_out_f, a_which_labels=False):
        """"Writes fragment data to xyz file."""
        
        print(f"Now writing to {a_out_f}")
        
    
        out_str = f"{self.n_atoms}"+"\n\n"
     
        for i, at in self.Atoms.items():
            try:
                label= a_which_labels[at["type"]]#+str(i)
            except:
                if a_which_labels:
                    try:
                        label = self.masses_type_labels[at["type"]]
                    except:
                        exit(f"Label for atom {at['type']} not found. Try witing the xyz with LAMMPS own numeric labels instead.")
                else:
                    label = self.type_to_element[at["type"]]

            x = at["posX"]
            y = at["posY"]
            z = at["posZ"]
            
            out_str += f"{label}     {x}      {y}      {z}\n"
        
    
        with open(a_out_f,"w") as out:
            out.write(out_str)
    
    def write_pdb(self, out_f,connect=True):
        """Writes a .pdb file from the LammpsData info."""
        with open(out_f, "w") as pdb_out: 
            #write HEADER
            header="HEADER"+4*" "+"Unclassified".ljust(40, " ")+gdate+"\n"
            title="TITLE".ljust(10, " ")+f"pdb file for {self.name}\n"
            author="AUTHOR".ljust(10, " ")+"built with lamps_data_analysis by A. Zoccante\n"
            head=header+title+author
            pdb_out.write(head)
            
            #write ATOMS section
            for i,atom in self.Atoms.items():
                if self.masses_type_labels[atom["type"]]!="":
                    name = self.masses_type_labels[atom["type"]]
                    if len(name) > 4:
                        name=name[:4]
                else:
                    name=str(atom["type"])
                res= str(atom["mol"])
                x=f"{atom['posX']:8.3f}"
                y=f"{atom['posY']:8.3f}"
                z=f"{atom['posZ']:8.3f}"
                element = self.type_to_element[atom["type"]]
                line="HETATM"+str(i).rjust(5, " ")+" "+name.ljust(4, " ")+" "+res.rjust(4, " ")+" "+"   0"+4*" "+x+y+z+f"{1.0:6.2f}"+f"{0.0:6.2f}"+10*" "+element.rjust(2," ")+2*" "+"\n"
                pdb_out.write(line)
                
            #write CONECT section
            if connect:
                for atom in self.Atoms.keys():
                    line="CONECT"+str(atom).rjust(5, " ")
                    if len(self.bonds_per_at[atom])< 5:
                        for bonded in self.bonds_per_at[atom]:
                            line+=str(bonded).rjust(5, " ")
                        line+="\n"
                    elif len(self.bonds_per_at[atom])>8:
                        exit("Coordination higher than 8 is not supported at the current momen. If suitable use option without connectivity")
                    else:
                        for bonded in self.bonds_per_at[atom][:4]:
                            line+=str(bonded).rjust(5, " ")
                        line+="\n"
                        line+="CONECT"+str(atom).rjust(5, " ")
                        for bonded in self.bonds_per_at[atom][4:]:
                            line+=str(bonded).rjust(5, " ")
                        line+="\n"
                    pdb_out.write(line)
           
            #write CRYST1
            a,b,c,alpha,beta,gamma = self.get_cell_params()
                
            line="CRYST1"+f"{a:9.3f}"+f"{b:9.3f}"+f"{c:9.3f}"+f"{alpha:7.2f}"+f"{beta:7.2f}"+f"{gamma:7.2f}"+"P 1".rjust(12, " ")+"1".rjust(4, " ")
            pdb_out.write(line)
            end="END   "
            pdb_out.write(end)

    def write_ciff(self,a_out_f,a_which_labels=False):
        """"Writes fragment data to cif file."""
        
        print(f"Now writing to {a_out_f}")
        
        preamble = """#======================================================================                                                                                                                                    

# CRYSTAL DATA

#----------------------------------------------------------------------

data_LAMMPS_phase_1


"""
        
        out_str = preamble
        
        
        name = self.name
        out_str += f"_chemical_name_common                  '{name}'\n"
        
        a,b,c,alpha,beta,gamma = self.get_cell_params()
        if (alpha != 90 or beta != 90 or gamma !=90):
            exit("Current implementation of LammpsData.write_ciff is not compatible with non-orthogonal cells.")
        out_str += f"_cell_length_a                         {a}\n"
        out_str += f"_cell_length_b                         {b}\n"
        out_str += f"_cell_length_c                         {c}\n"
        out_str += f"_cell_angle_alpha                      {alpha}\n"
        out_str += f"_cell_angle_beta                      {beta}\n"
        out_str += f"_cell_angle_gamma                      {gamma}\n"
        
        space_grp_lines = """_symmetry_cell_setting                 cubic
_space_group_name_H-M_alt              'P 1'
_symmetry_space_group_name_Hall 'P 1'
_symmetry_space_group_name_H-M  'P 1'
_symmetry_Int_Tables_number     1

loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'

"""
        out_str += space_grp_lines
        
        loop_str = """loop_
   _atom_site_label
   _atom_site_type_symbol
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_charge
   _atom_site_print_to_pdb
"""
        out_str += loop_str
        
        for i, at in self.Atoms.items():
            try:
                label= a_which_labels[at["type"]]#+str(i)
            except:
                if a_which_labels:
                    try:
                        label = self.masses_type_labels[at["type"]]
                    except:
                        exit(f"Label for atom {at['type']} not found. Try witing the xyz with LAMMPS own numeric labels instead.")
                else:
                    label = self.type_to_element[at["type"]]
            #occ = "1.0"
            fract_x = (at["posX"]-self.xbounds[0])/a
            fract_y = (at["posY"]-self.ybounds[0])/b
            fract_z = (at["posZ"]-self.zbounds[0])/c
            type_s = self.type_to_element[at["type"]]
            charge = at["charge"]
            out_str += f"   {label}    {type_s}     {fract_x}      {fract_y}      {fract_z}   {charge}   yes\n"
        

        with open(a_out_f,"w") as out:
            out.write(out_str)
            
    def write_rigid_mcf(self, mcfFile, a_which_labels=None, a_which_params=None):
            """ Prints out a .mcf file suitable for a rigid system based on the original LAMMPS data.
Arguments:
    mcfFile, string = molecular connectivity file (output)
    a_which_labels, dict = numerical to str atom types dictionary, or None for using element symbol.
    a_wich_params, dict = dictionary of atom types (str) to parameters dictionaries, or None to use LAMMPS ones if available (only correct if pair-style is lj (or related). 
returns:
    none
"""

            mcf = open(mcfFile, 'w')    

            mcf.write('!***************************************' + 
                  '****************************************\n')
            mcf.write('!Molecular connectivity file for ' + self.name + '\n')
            mcf.write('!***************************************' + 
                  '****************************************\n')

            mcf.write('!Atom Format\n')
            mcf.write('!index type element mass charge vdw_type parameters\n' + 
                  '!vdw_type="LJ", parms=epsilon sigma\n' + 
                    '!vdw_type="Mie", parms=epsilon sigma repulsion_exponent dispersion_exponent\n')
            mcf.write('\n# Atom_Info\n')
            mcf.write(str(len(self.Atoms))+'\n')
            for mcfIndex,atom in self.Atoms.items():
                mcf.write('%-4d'     % (int(mcfIndex)))
                
                try:
                    label= a_which_labels[atom["type"]]#+str(i)
                except:
                    if a_which_labels:
                        try:
                            label = self.masses_type_labels[atom["type"]]
                        except:
                            exit(f"Label for atom {atom['type']} not found. Try witing the xyz with LAMMPS own numeric labels instead.")
                    else:
                        label = self.type_to_element[atom["type"]]
                
               
                mcf.write('  %-6s'   % (label))

                mcf.write('  %-2s'   % (self.type_to_element[atom['type']]))
                mcf.write('  %7.3f'  % (self.Masses[atom['type']]))
                
                if (a_which_params == False) or not  "charge" in a_which_params[a_which_labels[atom['type']]].keys():
                    mcf.write('  %s'     % (atom['charge']))
                else:
                    mcf.write('  %s'     % ( a_which_params[a_which_labels[atom['type']]]['charge']))
                
                mcf.write('  %2s'    % "LJ")

                #mcf.write('  %2s'    % atomParms[pdbIndex]['vdw'][0])
                if a_which_params == False:
                    try:
                        mcf.write('  %8.3f'     % ( float(self.Pair_Coeffs[atom["type"]][0])/1.987204259E-03)) # Assumes LAMMPS units real are used.
                        mcf.write('  %8.3f'     % ( float(self.Pair_Coeffs[atom["type"]][1]))) # Assumes LAMMPS units real are used.
                    except:
                        try:
                            mcf.write('  %8.3f'     % ( float(self.PairIJ_Coeffs[atom["type"],atom["type"]][0])/1.987204259E-03)) # Assumes LAMMPS units real are used.
                            mcf.write('  %8.3f'     % ( float(self.PairIJ_Coeffs[atom["type"],atom["type"]][1]))) # Assumes LAMMPS units real are used.
                        except:
                            exit("Could not find LJ parameters in LAMMPS data file.")
                else:
                    mcf.write('  %8.3f'     % ( a_which_params[a_which_labels[atom['type']]]['epsilon']))
                    mcf.write('  %8.3f'     % ( a_which_params[a_which_labels[atom['type']]]['sigma']))
                # for i in range(1,len(atomParms[pdbIndex]['vdw'])):
                    # mcf.write('  %8.3f' % atomParms[pdbIndex]['vdw'][i])
                # for ring in ringList:
                    # if str(pdbIndex) in ring: mcf.write('  ring')
                mcf.write('\n')
            

            mcf.write('\n# Bond_Info\n')
            mcf.write(str(0) + '\n')

            mcf.write('\n# Angle_Info\n')
            mcf.write(str(0) + '\n')
            
            mcf.write('\n# Dihedral_Info\n')
            mcf.write(str(0) + '\n')

            mcf.write('\n# Improper_Info\n')
            mcf.write(str(0) + '\n')

            mcf.write('\n!Intra Scaling\n')
            mcf.write('!vdw_scaling    1-2 1-3 1-4 1-N\n')
            mcf.write('!charge_scaling 1-2 1-3 1-4 1-N\n')
            mcf.write('\n# Intra_Scaling\n')
            mcf.write('0. 0. %.4f 1.\n' % (0))
            mcf.write('0. 0. %.4f 1.\n' % (0))
            
            mcf.write("""!Fragment Format
    !index number_of_atoms_in_fragment branch_point other_atoms

    # Fragment_Info
    0


    # Fragment_Connectivity
    0


    END""")

            mcf.close()
            
#--------------- Utility functions -------------------------------------------------------------------------------------#      

    def get_cell_params(self):
        '''Ouputs crystallographical cell parameters (a,b,c,alpha,beta,gamma) from lammps data parameters'''''
        lx=self.xbounds[1]-self.xbounds[0]
        ly=self.ybounds[1]-self.ybounds[0]
        lz=self.zbounds[1]-self.zbounds[0]
        
        if hasattr(self, "cell_angles"):         # This applies to restricted triclinic 
            xy=self.cell_angles[0]
            xz=self.cell_angles[1]
            yz=self.cell_angles[2]
            a=lx
            b=(ly**2+xy**2)**0.5
            c=(lz**2+xz**2+yz**2)**0.5
            alpha= degrees(acos((xy*xz+ly*yz)/(b*c)))
            beta=degrees(acos(xz/c))
            gamma=degrees(acos(xy/b))
        else:
            alpha=90
            beta=90
            gamma=90
            a=lx
            b=ly
            c=lz
            
        return (a,b,c,alpha,beta,gamma)
               
      
    def add_Coeffs(self,coeffs_in_file):
        """"reads in lammps FF coeffs from a lammps input formatted file. Adds them to the LammpsData"""
        
        with open(coeffs_in_file) as coeffs_f:
            term_name=""
            for line in coeffs_f.readlines():
                if "coeff" in line:
                    new_term_name=line.split()[0].split("_")[0].capitalize()
                    if new_term_name!= term_name:
                        term_name=new_term_name
                        self.ff_terms.append(term_name)
                        if term_name == "Pair":
                            setattr(self, term_name+"IJ_Coeffs",{})
                        else:
                            setattr(self, term_name+"_Coeffs",{})
                    if term_name == "Pair":
                        getattr(self, term_name+"IJ_Coeffs")[tuple(sorted([int(i) for i in line.split()[1:3]]))]= line.split()[3:]
                    
                    else:
                        key = int(line.split()[1])
                        this_term_dict = {}
                        this_term_dict["coeffs"]= line.split()[2:]
                        getattr(self, term_name+"_Coeffs")[key]=this_term_dict

                            
    def get_system_mass(self):
        '''Returns total mass of the system in amu.'''
        o_mass =0.0
        for at in self.Atoms.values():
            o_mass += self.Masses[at["type"]]
        
        return o_mass
    
    def get_average_uma(self):
        '''Returns average mass of an atom in the system in amu.'''
        return self.get_system_mass()/self.n_atoms
    
    def get_system_volume(self):  
        '''Returns volume of the system simulation box in Aangstrom.'''
        x_l =self.xbounds[1]-self.xbounds[0]
        y_l =self.ybounds[1]-self.ybounds[0]
        z_l =self.zbounds[1]-self.zbounds[0]
        volume = x_l*y_l*z_l
        return volume
        
    def get_system_density(self):
        '''Returns density of the system in amu/Aa**3.'''
        o_mass = self.get_system_mass()
        o_vol = self.get_system_volume()
        o_dens_amuAa3= o_mass/o_vol 
        return o_dens_amuAa3
        
    def get_system_charge(self):
        '''Returns total charge of the system (note, it depends on lammps units used).'''
        o_charge = 0.0
        for at in self.Atoms.values():
            o_charge += at["charge"]
            
        return o_charge
    
    def number_atoms_per_type(self):
        '''Returns number of atoms present in the system for each atom type.'''
        n_at_t = {}
        for at in self.Atoms.values():
            if at["type"] not in n_at_t.keys():
                n_at_t[at["type"]] = 1
            else:
                n_at_t[at["type"]] += 1      
        return n_at_t
                
    def atom_IDs_for_types(self, ls_types):
        '''Returns Ids of atoms present in the system for each given atom type.'''
        #print("in types func", ls_types)
        o_at_ls_4type=[]
        #print("in atom_Ids_for_types")
        for i_at,at_info in self.Atoms.items():
            if at_info["type"] in ls_types:
                o_at_ls_4type.append(i_at) 
        #print("returning: ",o_at_ls_4type)
        return o_at_ls_4type
        
    def shift_total_charge_to_zero(self):
        '''Changes all charges by a shift equal to (total charge/number of atoms) to achieve 0 total charge'''
        tot_charge = 0.0
        for i,at in self.Atoms.items():
            #chrg_f.write("charge\n"+f"{i} {at['charge']}"+"\n\n")
            tot_charge += at['charge']
        print(f"Total charge is {tot_charge}")
        for i,at in self.Atoms.items():
            at['charge'] = at['charge'] - (tot_charge/self.n_atoms)
            
    def shift_atom_type_charge(self,a_at_type,a_charge_shift):
        '''Changes all charges of atoms of a given type by a given shift'''
        for i,at in self.Atoms.items():
            #chrg_f.write("charge\n"+f"{i} {at['charge']}"+"\n\n")
            if at["type"] == a_at_type:
                at['charge']+=a_charge_shift
                    
    def mol_sizes(self):
        '''Outputs number of molecules per size of molecule, after rechecking actual connectivity (useful as sometimes lammps react does not properly change molecule labels)'''
        #print(f"Molecules : {self.n_mols}")

        #print("Molecules per size:")
        #print(self.mol_size_info)
        return self.mol_sizes_dict
    
    def renumber(self):
        '''Renumbers atoms from 1 to n_atoms in the same order as found in the file''' 
        # renumbering atoms first
        at_c = 1
        old_to_new_at_keys = {}
        new_atoms = {}
        for key,atom_info in self.Atoms.items():
            old_to_new_at_keys[key]=at_c
            new_atoms[at_c] = atom_info
            at_c += 1
        self.Atoms = new_atoms
        
        #Now renumbering velocities
        new_vels = {}
        if hasattr(self,"Velocities"):
            for key,vel_info in self.Velocities.items():
                new_vels[new_atoms[key]]=vel_info
        
        #Now renumbering bonds etc.
        for attr in ["Bonds","Angles","Dihedrals","Impropers"]:
            #print(" at attr: ", attr)
            attr_c = 1
            new_attr = {}
            for key,attr_info in getattr(self,attr).items():
                attr_ats = []
                for key,value in attr_info.items():
                    if "atom" in key:
                        attr_ats+=[value]
                
                new_attr_info = {}
                for key2,value2 in  attr_info.items():
                    #print("at line 75: ", key2, value2)
                    if "atom" in key2:
                        new_attr_info[key2]=old_to_new_at_keys[value2]
                    else:
                        new_attr_info[key2] = value2
                new_attr[attr_c]=new_attr_info
                attr_c += 1 

            setattr(self,attr,new_attr)

            setattr(self,f"n_{attr}",len( getattr(self,attr)))

    def reorder(self):
        """Reorders the atoms by atoms index"""
        self.Atoms = {k:self.Atoms[k] for k in sorted([n for n in self.Atoms.keys()])}
        
    def reassign_molecules(self):
        '''Reassign molecules labels based on actual connectivity (which can be broken by some sofwares and by lammps react)'''
        # Create fraginfo graph for networkit
        #
        G = nk.Graph(self.n_atoms, directed=False,weighted=False)
        for key,bond_info in self.Bonds.items():
            G.addEdge(int(bond_info["atom1"])-1, int(bond_info["atom2"])-1)
        cc = nk.components.ConnectedComponents(G)
        cc.run()
        self.mol_sizes_dict = {key+1:value for key,value in cc.getComponentSizes().items()}
        self.n_mols = len(self.mol_sizes_dict)
        self.mol_size_info={}
        for mol,size in self.mol_sizes_dict.items():
            if size in self.mol_size_info.keys():
                self.mol_size_info[size]+=1
            else:
                self.mol_size_info[size]=1
        
            #test with v=atoms_number-1
        for num,at_info in self.Atoms.items():
            v=int(num)-1
            #print("Atom ", num , ", molecule: " , cc.componentOfNode(v)+1)
            at_info["mol"]= cc.componentOfNode(v)+1
            
        #Now assigning all atoms to the right molecules:
        for at_id,this_atom_dict in self.Atoms.items():
            if this_atom_dict["mol"] in self.atoms_per_mol.keys():
                self.atoms_per_mol[this_atom_dict["mol"]].append(at_id)
            else:
                self.atoms_per_mol[this_atom_dict["mol"]] = [at_id]



    @classmethod
    def remove_atoms(cls,a_frag_info,a_remove_at_ls,renumber=True):
        '''This removes a ls of atoms from the system, including their connectivity.'''
        print("This is remove_atoms")
        new_frag_info_obj = cls.__new__(cls)
        
        if hasattr(a_frag_info,"name") and a_frag_info.name != None:
            new_frag_info_obj.name = "new_" + a_frag_info.name

        new_frag_info_obj.ff_terms =a_frag_info.ff_terms
        new_frag_info_obj.precise= a_frag_info.precise
        
        new_frag_info_obj.n_atom_types = a_frag_info.n_atom_types
        new_frag_info_obj.n_bond_types = a_frag_info.n_bond_types
        new_frag_info_obj.n_angle_types = a_frag_info.n_angle_types
        new_frag_info_obj.n_dihedral_types = a_frag_info.n_dihedral_types
        new_frag_info_obj.n_improper_types = a_frag_info.n_improper_types

        new_frag_info_obj.xbounds = a_frag_info.xbounds
        new_frag_info_obj.ybounds = a_frag_info.ybounds
        new_frag_info_obj.zbounds = a_frag_info.zbounds
        
        if hasattr(a_frag_info, "cell_angles"):
                new_frag_info_obj.cell_angles = a_frag_info.cell_angles
        
        if hasattr(a_frag_info,"Masses"):
            new_frag_info_obj.Masses = a_frag_info.Masses
        if hasattr(a_frag_info,"masses_type_labels"):    
            new_frag_info_obj.masses_type_labels = a_frag_info.masses_type_labels
        if hasattr(a_frag_info,"atoms_names_per_type"):
            new_frag_info_obj.atoms_names_per_type = a_frag_info.atoms_names_per_type

        if hasattr(a_frag_info, "Pair_Coeffs"):
            new_frag_info_obj.Pair_Coeffs = a_frag_info.Pair_Coeffs
        if hasattr(a_frag_info, "PairIJ_Coeffs"):
                new_frag_info_obj.PairIJ_Coeffs = a_frag_info.PairIJ_Coeffs
        
        try:
            for term in new_frag_info_obj.ff_terms:
                if term != "Pair":
                    setattr(new_frag_info_obj,term+"_Coeffs",getattr(a_frag_info,term+"_Coeffs"))
            new_frag_info_obj.Bond_Coeffs = a_frag_info.Bond_Coeffs
            new_frag_info_obj.Angle_Coeffs = a_frag_info.Angle_Coeffs
            new_frag_info_obj.Dihedral_Coeffs = a_frag_info.Dihedral_Coeffs
            new_frag_info_obj.Improper_Coeffs = a_frag_info.Improper_Coeffs
        except AttributeError:
            pass
            
        at_c = 1
        new_frag_info_obj.Atoms = {}
        new_frag_info_obj.atoms_per_mol = {}
        old_to_new_at_keys = {}
        #print(f"atoms to be removed: {a_remove_at_ls}")
        for key,atom_info in a_frag_info.Atoms.items():
            if (key not in a_remove_at_ls) and renumber:
                old_to_new_at_keys[key]=at_c
                new_frag_info_obj.Atoms[at_c] = atom_info
                at_c += 1
            elif (key not in a_remove_at_ls) and not renumber:
                old_to_new_at_keys[key]=key
                new_frag_info_obj.Atoms[key] = atom_info
            else:
                #print(f"atom {key} to be removed:")
                continue

        #print("old_to_new_at_keys: ", old_to_new_at_keys)  
        new_frag_info_obj.n_atoms=len(new_frag_info_obj.Atoms)
            
        for attr in ["Bonds","Angles","Dihedrals","Impropers"]:
            #print(" at attr: ", attr)
            setattr(new_frag_info_obj,attr,{})
            attr_c = 1
            for key,attr_info in getattr(a_frag_info,attr).items():
                attr_ats = []
                for key,value in attr_info.items():
                    if "atom" in key:
                        attr_ats+=[value]
               # print(attr_ats)
                if len(set(a_remove_at_ls).intersection(attr_ats)) == 0:
                    #print("afte intersection: ",set(a_remove_at_ls).intersection(attr_ats))
                    new_attr_info = {}
                    for key2,value2 in  attr_info.items():
                        #print("at line 75: ", key2, value2)
                        if "atom" in key2:
                            new_attr_info[key2]=old_to_new_at_keys[value2]
                        else:
                            new_attr_info[key2] = value2
                    getattr(new_frag_info_obj,attr)[attr_c]=new_attr_info
                    attr_c += 1 
                    
            setattr(new_frag_info_obj,f"n_{attr}".lower(),len( getattr(new_frag_info_obj,attr)))
        
        if (renumber):
            new_frag_info_obj.renumber()
            new_frag_info_obj.reassign_molecules()
        return new_frag_info_obj


    @classmethod
    def remove_bonds(cls,a_frag_info,a_bond_ls):
        print("This is remove_bond")
        new_frag_info_obj = cls.__new__(cls)
        
        if hasattr(a_frag_info,"name") and a_frag_info.name != None:
            new_frag_info_obj.name = "new_" + a_frag_info.name

        new_frag_info_obj.ff_terms =a_frag_info.ff_terms
        new_frag_info_obj.precise= a_frag_info.precise
        
        new_frag_info_obj.n_atom_types = a_frag_info.n_atom_types
        new_frag_info_obj.n_bond_types = a_frag_info.n_bond_types
        new_frag_info_obj.n_angle_types = a_frag_info.n_angle_types
        new_frag_info_obj.n_dihedral_types = a_frag_info.n_dihedral_types
        new_frag_info_obj.n_improper_types = a_frag_info.n_improper_types

        new_frag_info_obj.xbounds = a_frag_info.xbounds
        new_frag_info_obj.ybounds = a_frag_info.ybounds
        new_frag_info_obj.zbounds = a_frag_info.zbounds
        
        if hasattr(a_frag_info, "cell_angles"):
                new_frag_info_obj.cell_angles = a_frag_info.cell_angles
        
        if hasattr(a_frag_info,"masses"):
            new_frag_info_obj.masses = a_frag_info.masses
        if hasattr(a_frag_info,"masses_type_labels"):    
            new_frag_info_obj.masses_type_labels = a_frag_info.masses_type_labels
        if hasattr(a_frag_info,"atoms_names_per_type"):
            new_frag_info_obj.atoms_names_per_type = a_frag_info.atoms_names_per_type

        if hasattr(a_frag_info, "Pair_Coeffs"):
            # Only assuming class1 coeffs. TODO: generalize to general coeffs case.
            new_frag_info_obj.Pair_Coeffs = a_frag_info.Pair_Coeffs
        if hasattr(a_frag_info, "PairIJ_Coeffs"):
                new_frag_info_obj.PairIJ_Coeffs = a_frag_info.PairIJ_Coeffs
        
        try:
            for term in new_frag_info_obj.ff_terms:
                if term != "Pair":
                    setattr(new_frag_info_obj,term+"_Coeffs",getattr(a_frag_info,term+"_Coeffs"))
            new_frag_info_obj.Bond_Coeffs = a_frag_info.Bond_Coeffs
            new_frag_info_obj.Angle_Coeffs = a_frag_info.Angle_Coeffs
            new_frag_info_obj.Dihedral_Coeffs = a_frag_info.Dihedral_Coeffs
            new_frag_info_obj.Improper_Coeffs = a_frag_info.Improper_Coeffs
        except AttributeError:
            pass
            
        new_frag_info_obj.Atoms = a_frag_info.Atoms
        new_frag_info_obj.atoms_per_mol = a_frag_info.atoms_per_mol

          
        new_frag_info_obj.n_atoms = a_frag_info.n_atoms
            
        for attr in ["bonds","angles","dihedrals","impropers"]:
            #print(" at attr: ", attr)
            setattr(new_frag_info_obj,attr,{})
            attr_c = 1
            for key1,attr_info in getattr(a_frag_info,attr).items():
                attr_ats = []
                for key2,value in attr_info.items():
                    if "atom" in key2:
                        attr_ats+=[value]
                #print(attr_ats)
               
                found=False
                for b in a_bond_ls:
                    b_at = b.split("-")
                    #print(set(b_at).intersection(attr_ats))
                    if len(set(b_at).intersection(attr_ats)) == 2 : 
                        found = True
                        #print("after intersection: ",set(b_at).intersection(attr_ats))
                   
                if not found:
                    getattr(new_frag_info_obj,attr)[attr_c]=getattr(a_frag_info,attr)[key1]
                    attr_c += 1 
                    
            setattr(new_frag_info_obj,f"n_{attr}",len( getattr(new_frag_info_obj,attr)))
        
        new_frag_info_obj.renumber()
        new_frag_info_obj.reassign_molecules()
        return new_frag_info_obj
        
    def atom_nbs_for_mols(self, ls_mols):
        o_at_ls_4mol=[]
        for i_at,at_info in self.Atoms.items():
            if at_info["mol"] in ls_mols:
                o_at_ls_4mol.append(i_at) 
        return o_at_ls_4mol
        
    def list_unwanted_atoms(self, a_wanted_lists):
        """Expects a list of lists of atoms to keep, returns a corresponding list of lists of atoms not to keep."""
        o_unwanted_ats_lists = []
        for frag in a_wanted_lists:
            o_unwanted_ats_lists.append([])
        for at_id in self.Atoms.keys():
            for i,frag in enumerate(a_wanted_lists):
                if int(at_id) not in frag:
                    o_unwanted_ats_lists[i].append(at_id)
        return o_unwanted_ats_lists
    
    def fragment_map(self,start_at):
        """Returns a map of a molecular system, assumed fully connected
self: the LammpsData obj representing the molecular system
start_at: the atom id to be used as the root of the system connectivity (directed) graph
returns: a list of dictionaries representing the system connectivity (directed) graph"""
        self.found_list=[]
        for i in range(0,self.n_atoms):
            self.found_list.append(False)
        self.found_list[start_at-1] = True
        
        parent_id=None
    
        frag_type_map =[{parent_id:[(start_at,self.Atoms[start_at]["type"])]}] 
        
        neigh_level = 0
        while self.found_list.count(True) != self.n_atoms:
            #print(frag_type_map)
            frag_type_map = self.find_neigh(frag_type_map,neigh_level)
    
            neigh_level+=1
        frag_type_map = self.find_neigh(frag_type_map,neigh_level)
    
        return frag_type_map
        

    def find_neigh(self, a_frag_map, a_level):
        """Adds to a fragment map the nodes that are children of the nodes at a specific level.
returns: the updated fragment map."""
        #print("a_level: ", a_level)
        neighs ={}
        level_map = a_frag_map[a_level] #check
        for parent in level_map.keys():
            if level_map[parent] != "Leaf":
                for at_id ,at_types in  level_map[parent]:    
                    #print(at_id)
                    neighs_tmp=[]
                    for neigh_at in self.bonds_per_at[at_id]:
                        #print("checking(",neigh_at,")")
                        if self.found_list[neigh_at-1] != True:
                            #print("found!")
                            neighs_tmp.append((neigh_at,self.Atoms[neigh_at]["type"]))# parent_id:[(child_id, child_type),etc.]
                            self.found_list[neigh_at-1]=True
                    if len(neighs_tmp) != 0:
                        neighs[at_id]=neighs_tmp
                    else:
                        neighs[at_id]="Leaf"
       
        a_frag_map.append(neighs)
        return a_frag_map


    def find_fragments(self,a_frag_map):
        """Finds fragments in a larger systems and return a list of lists their atoms ids."""
        
        o_fragments=[]        
        #Now looking for root atoms  
        for at_id, at_info in self.Atoms.items():
            
            root_at_id = find_map_root_id(a_frag_map)
            root_type = find_map_root_type(a_frag_map)
            
            if at_info["type"] == root_type:
                print(f"Found candidate {at_id}")
                #initializing dictionary for mapping atoms
                self.this_fragment_map={}
                
                for level in a_frag_map:
                    for children in level.values():
                        if children != "Leaf":
                            for map_at_id,map_at_type in children:
                                self.this_fragment_map[map_at_id]=False
                self.this_fragment_map[int(root_at_id)]=int(at_id) # Not false anymore (for now).
                print("initial map: ",self.this_fragment_map)
                print(self.find_map_at_id(int(at_id)))
                found_all = self.topology_overlaps(at_id, None, 0, a_frag_map)
                
                if found_all==False:
                    print(found_all)
                    pass
                else:
                    print(f"Mapping complete: the found fragment map is {self.this_fragment_map}\n")
                    res = self.get_atoms_from_map(self.this_fragment_map)
                    o_fragments.append(res) 
        if len(o_fragments)==0:
            exit("Found no candidates. Exiting script")
        else:
            return o_fragments
            

    def find_map_at_id(self,a_at_id):
        """Finds the fragment atom ID associated to a_at_id in the map.
    input: a_at_id: atom id of a system atom (tentatively) associated to a fragment map atom id (node id)
    return: map_at_id: if a corresponding fragment map atom id is found, or 
            False: if no match is found."""
        
        for map_at_id,at_id in self.this_fragment_map.items():
            if a_at_id == at_id:
                return map_at_id
        print(f"Warning: Find_map_at_id could not match {a_at_id} to any map atom.")
        return False


    def topology_overlaps(self,a_cand_id,a_parent_id, a_level,a_frag_map):
        """Checks if the local topology around a_cand_id matches the fragment. Iterative.
    input: a_cand_id: atom ID of a candidate atom (belonging to the whole system)
           a_parent_id: atom id of the candidate parent atom (belonging to the whole system)
           a_level: level of the current candidate node
           a_frag_map: the map of the fragment to be found.
    return: True: if subgraph starting at a_cand_id matches map
            False: if it does not."""
        
        # if (gDebug):
        #     print(f"In topology_overlaps, a_cand_id = {a_cand_id} ,a_parent_id= {a_parent_id}, a_level= {a_level} ,a_frag_map={a_frag_map}")
    
        # print("Finding out if we look for a leaf.")
        if a_frag_map[a_level+1][self.find_map_at_id(int(a_cand_id))] == "Leaf":
            # print("candidate is a leaf!")
            return True
            
        # Assumes current atom has been matched in previous call.
        child_l = []
        
        # print("Finding if reached terminal atom in fragment")
        for child in self.bonds_per_at[int(a_cand_id)]:
            #print(child,a_parent_id)
            #print(type(child), type(a_parent_id))
    
            if child != a_parent_id:
                #print(f"Found a child atom! child at_id {child}")
                child_l.append(self.Atoms[child]["type"])
            #print(child_l)
        if len(child_l)==0: #if the atom has no child it is a "leaf"
            #print("No children.")
            #print(a_level, a_frag_map[a_level+1])
            if a_frag_map[a_level+1][self.find_map_at_id(a_cand_id)]=="Leaf":
                return True
            else:
                return False
        
        #Now finding the types for the childs of current expected atom in the map  (for filtering)
        types_l=[]
        # print("Preparing preselection")
        map_at_id = self.find_map_at_id(int(a_cand_id))
        if map_at_id != False:
            #print("current level: ",a_level)
            #print("a_cand_id:", a_cand_id, type(a_cand_id))
            #print("map_at_id:", map_at_id, type(map_at_id))
            if a_frag_map[a_level+1][map_at_id] != "Leaf":
                for map_ats in a_frag_map[a_level+1][int(map_at_id)]: #These are the current guessed map atom childrens
                    types_l.append(map_ats[1]) 
        else:
            # print("I should never be here, HELP!")
            return False #
        #print("types_l: ", types_l)
       
        
        # Now checks if the map children types are a sublist of the actual children types
        count_c=0
        count_t=0
        for c in child_l:
            for t in types_l: 
                if t in child_l:
                        #print("count_t",count_t)
                        #print("coun_c",count_c)
                        types_l.pop(types_l.index(t))
                        #print("popping elements, now types_l is : ", types_l)
                        #print("len(types_l):",len(types_l))
                        count_t+=1
                        continue
            count_c+=1
    
            if len(types_l)==0:
                continue
            else:
                #print("Still something to do!")
                pass
        if len(types_l) != 0: #This level did not have all types in this level of map.
            # print("Could not find all atoms I was looking for.")
            return False
        else: # This is a prospective match (at this level)
            # Starting to guess assignments:
                
            found_l=[]
            parent=self.find_map_at_id(int(a_cand_id))
            frag_children = a_frag_map[a_level+1][parent]
            for sys_child_at in self.bonds_per_at[int(a_cand_id)]:
                # print(f"current sys_child_at is {sys_child_at}")
                # print(f"Have we foundit already? {sys_child_at in self.this_fragment_map.values()}")
                if sys_child_at not in self.this_fragment_map.values():
                    print(sys_child_at , a_parent_id,type(sys_child_at) , type(a_parent_id))
                    if sys_child_at != a_parent_id:
                        # print("a_parent_id: ",a_parent_id)
                        # print("CHECKING CHILD: ", sys_child_at)
        
                        # for parent,frag_children in a_frag_map[a_level+1].items():
                        #     if self.this_fragment_map[parent]==a_cand_id
                        #     print(f"looking for a match for {frag_children}")
                        #     #print("level: ",a_frag_map[a_level+1].items())
                        #     print("frag_children",frag_children)
                     
                        for frag_child in frag_children:
                            # print("frag_child: ", frag_child)
                            #print("self.Atoms[sys_child_at]['type'],  frag_child[1]: ", self.Atoms[sys_child_at]["type"] , frag_child[1])
                            if self.Atoms[sys_child_at]["type"] ==  frag_child[1]:
                                # print("atom", frag_child[0], " is a good candidate")
                                if self.this_fragment_map[int(frag_child[0])]==False:
                                    self.this_fragment_map[int(frag_child[0])]=int(sys_child_at)
                                    # print("current map is: ",self.this_fragment_map)
                                    found = self.topology_overlaps(int(sys_child_at),a_cand_id,a_level+1,a_frag_map)
                                    
                                    if not found:
                                        # print("Candidate not accepted.")
                                        self.this_fragment_map[int(frag_child[0])]=False
                                        pass
                                    else:
                                        # print(f"candidate {frag_child[0]} accepted")
                                        found_l.append(found)
                                        break
                            # print("after, child: ", child)
                        # print("current map is ",self.this_fragment_map)
            if len(found_l) == len(frag_children):
                # print("confirming candidate, current map is ",self.this_fragment_map)
                return True
            else:
                return False
                            
    
            


    def get_atoms_from_map(self, a_fragment_map):
        """Returns all system atoms matched to a fragment map.
input: a_fragment_map: a map of the fragment to be found. Assumes the match has already been achieved.
return: o_atoms_l: list of the matching atoms for the current fragment in the whole system. One found fraagment at a time."""
        o_atoms_l=list(a_fragment_map.values())
        return o_atoms_l

    def read_coords_from_xyz(self,a_xyz_file):
        """Replace coordinates of atoms with the ones from a .xyz file. Assumes systems match."""
        xyz_n_atoms=0
        xyz_coords = []
        with open(a_xyz_file) as xyz_f:
            for i,line in enumerate(xyz_f.readlines()):
                if i == 0 :
                    xyz_n_atoms=int(line.strip())
                elif i==1 :
                    continue
                else:
                    xyz_coords.append([float(x) for x in line.split()[1:]])
        if (self.n_atoms != xyz_n_atoms) or (self.n_atoms != len(xyz_coords)):
            exit("Error in get_atoms_from_xyz. Mismatch between number of atoms in .data and in .xyz")
        else:
            for i_at,atom in self.Atoms.items(): # Assumes ordering is the same.
                atom["posX"]=xyz_coords[i_at-1][0]
                atom["posY"]=xyz_coords[i_at-1][1]
                atom["posZ"]=xyz_coords[i_at-1][2]
        return

#--------------Legacy Functions, check before using !      -------------------------------------------------------------#
                
    def remap_types(self, a_types_f):
        '''Legacy function. To change types according to a given map file, built as : "section\n old_label new_labels \n ... " for any or all sections in this list:"atoms, bonds, dihedrals and/or impropers" '''
        a_types_dict=read_types_map(a_types_f)
        #renumbering atom types first
        for key,atom_info in self.Atoms.items():
            atom_info["type"]=a_types_dict["atoms"][atom_info["type"]]
        
        #Now bonds,etc.
        for attr in ["Bonds","Angles","Dihedrals","Impropers"]:
            #print(" at attr: ", attr)       
            try:
                for key,attr_info in getattr(self,attr).items():
                    attr_info["type"]=a_types_dict[attr][int(attr_info["type"])]
            except:
                pass
        #TODO. Recheck with current __init__() implementation.

    def set_types_numbers(self, a_dict_f):
        '''Legacy function'''
        a_types_n_dict=read_types_n(a_dict_f)
        for attr in a_types_n_dict.keys():
            #print(" at attr: ", attr)       
            try:
                setattr(self,attr,a_types_n_dict[attr])
            except:
                pass
            
