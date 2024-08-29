#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:57:11 2022

@author: alberto_z
"""
from lammps_data_analysis import LammpsData

#%% This will later be the body of the script.
import argparse

gTest = False

if __name__ == "__main__" and not gTest:

    parser = argparse.ArgumentParser(description=' Reads in a lammps data file and a file containing a list of atom types synonims and create a new lammps data file including all needed types.')
    parser.add_argument('-lammps', metavar='lammps_data_file',help='lammps data files')
      
    parser.add_argument('--out',help='Output lammps data file',default="data_out.lmps")
    parser.add_argument('--name', metavar='system_name',help='Name to assign to the molecular system', default=None)
    parser.add_argument('--renumber',help='Renumbers atoms and molecules for consistency',action='store_true')
    parser.add_argument('--reorder',help='Reorder atoms by atom index',action='store_true')
    parser.add_argument('--reassign_mols',help='Reassign molecules dependeing on connectivity',action='store_true')
    parser.add_argument('--print_coeffs',help="If present attempts to print FF coefficients in output .data file.",action="store_true")

    parser.add_argument('--summary',help='return a summary of the system info', action='store_true')  
   
    parser.add_argument('--tot_mass',help='return total mass of input system', action='store_true')  
    parser.add_argument('--volume',help='return system density', action='store_true') 
    parser.add_argument('--density',help='return system density', action='store_true')  
    parser.add_argument('--avemass',help='return system average mass (amu)', action='store_true')  
    parser.add_argument('--cassandra_charges',help='If present prints out charges for cassandra files praparation.', action="store_true")
    parser.add_argument('--n4type',help="If present computes the number of atoms of given type or a summary of all types (write 'all')", default=None)
    parser.add_argument('--mol_sizes',help='returns a summary of sizes for the molecules in the system', action='store_true') 
    parser.add_argument('--mol_sizes_all',help='return the size for each molecule in the system', action='store_true')   
    parser.add_argument('--charge',help='return system charge', action='store_true') 
    
    args = parser.parse_args()
    
    lammps_f = args.lammps
    
    tot_mass = args.tot_mass
    
    density = args.density
    
    avemass = args.avemass
    
    volume = args.volume
    
    charge =args.charge
    
    cassandra_charges = args.cassandra_charges
    
    
    mol_sizes = args.mol_sizes
    mol_sizes_all = args.mol_sizes_all
    
    renumber=args.renumber
    reorder=args.reorder
    reassign_mols=args.reassign_mols
        
        
    out_f = args.out
  
    
    if (lammps_f == ""):
        exit("ERROR: no input files assigned.")


    frag_info = LammpsData(lammps_f,renumber=renumber,reorder=reorder,reassign_mols=(renumber or reassign_mols),name=args.name)
    
    if args.summary:
        with open("summary.txt", "w") as summary_f:
            summary_f.write(f"Summary informations for system {frag_info.name}\n\n\n")

            summary_f.write(f"Total mass:\t{frag_info.get_system_mass():.5F} amu\n")
            summary_f.write(f"Average atomic mass:\t{frag_info.get_average_uma():.5F} amu\n")
            summary_f.write(f"Density:\t{frag_info.get_system_density():.5F} uma/ang**3\n")
            summary_f.write(f"Volume:\t{frag_info.get_system_volume():.5F} ang**3\n ONLY HOLDS FOR ORTHOGONAL CELLS\n")
            summary_f.write(f"System total charge:\t{frag_info.get_system_charge():.5F} electron charges\n\n\n")
            
            summary_f.write("Number of atoms for each type.\n")
            summary_f.write("Atom type  | Nb. of atoms | Element\n")
            summary_f.write("-----------------------------------\n")            
            n4type=list(range(1,len(frag_info.Masses)+1))
            for the_type in n4type:
                summary_f.write(f"{the_type:^10} | {len(frag_info.atom_IDs_for_types([the_type])):^12} | {frag_info.type_to_element[the_type]:^9}\n")
            summary_f.write("\n\n")
            
            frag_info.reassign_molecules()

            summary_f.write("Molecule sizes\n")
            summary_f.write("Molecule size  | nb. of molecules\n")
            summary_f.write("------------------------\n")
            for mol_size in sorted(frag_info.mol_size_info.keys()):
                summary_f.write(f" {mol_size:^14}|{frag_info.mol_size_info[mol_size]:^16}\n") 
            summary_f.write("\n\n")
            
            summary_f.write("Size of each molecule\n")
            summary_f.write("Molecule ID  | nb of atoms\n")
            summary_f.write("------------------------\n")
            size_per_mol = frag_info.mol_sizes()
            for i_mol,at_nbs in size_per_mol.items():
                summary_f.write(f" {i_mol:^12}|{at_nbs:^12}\n") 
            summary_f.write("WARNING: molecule ID refers to molecules after ressignment. May not correspond to ID in original file.\n To produce the corresponding file use flag --reassign_mols")
                
            print("Finished writing summary")

            
    
    if tot_mass:
        print(f"{frag_info.get_system_mass():.5F}")
        
    if avemass:
        print(f"{frag_info.get_average_uma():.5F}")
        
    if density:
        print(f"{frag_info.get_system_density():.5F} uma/ang**3")
         
    if volume:
        print(f"{frag_info.get_system_volume():.5F}")
        
    if charge:
        print(f"{frag_info.get_system_charge():.5F}")
      
    if cassandra_charges:
        frag_info.write_cassandra_charges()
        
    
    if args.n4type != None :
        if args.n4type != 'all': # Expecting only one integer.
            n4type=[int(args.n4type)]
        else:
            n4type=list(range(len(frag_info.Masses)))
        for the_type in n4type:
            print(f"Number of atoms of type {the_type}: {len(frag_info.atom_IDs_for_types([the_type]))}")
        
    if  mol_sizes or mol_sizes_all:
        frag_info.reassign_molecules()
        size_per_mol = frag_info.mol_sizes()
        print("Molecules per size:")
        print(frag_info.mol_size_info)

        
        if mol_sizes_all:
                
            print("# Molecule ID\tnb of atoms")
            
            for i_mol,at_nbs in size_per_mol.items():
                print(f" {i_mol}\t{at_nbs}") 

    if renumber or reorder or reassign_mols:
        print(f"renumber {renumber}, reorder {reorder}, reassign_mols {reassign_mols}")
        if renumber:
            frag_info.renumber()
        if reorder:
            frag_info.reorder()
        if reassign_mols:
            frag_info.reassign_molecules()
            
        frag_info.write_lammps_data(out_f,coeffs_out=args.print_coeffs)
                
        
exit()
