#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:12:44 2022

@author: alberto_z
"""
#%% Importing FragmentInfo from other script
from lammps_data_analysis import LammpsData
    
#%% Functions to read inputs from files.
def read_types(a_types_f):
    """Reads a types file dictionary from a polymatic style types file. Or a file with "<atom_type> <label>" on each line."""
    o_types_dict ={}
    
    with open(a_types_f, "r") as f_types:
        for line in f_types.readlines():
           
            if "atom types" in line:
                continue
            elif "bond types" in line:
                break
            else:
                key=int(line.split()[0])
                val=line.split()[1]
                o_types_dict[key]=val
    # print(o_types_dict)
    return o_types_dict

def read_params(a_params_f):
    """Reads a parameters file dictionary from a parameters file."""
    o_params_dict ={}
    
    with open(a_params_f, "r") as f_params:
        current_t = None
        start = None
        for i,line in enumerate(f_params.readlines()):
           
            if "nonbonded" in line:
                start = i
            elif start!= None:
                if i == start+1:
                    current_t = line.strip()
                    o_params_dict[current_t]={}
                elif i == start+2:
                    o_params_dict[current_t]["sigma"]=float(line.split()[1])     
                elif i == start+3:
                    o_params_dict[current_t]["epsilon"]=float(line.split()[1])        
                elif i == start+4:
                    try:
                        o_params_dict[current_t]["charge"]=float(line.split()[1])   
                    except:
                        pass
            else:
                pass
    return o_params_dict
#%% This will later be the body of the script.
import argparse

gTest = False

if __name__ == "__main__" and not gTest:

    parser = argparse.ArgumentParser(description=' Reads in a lammps data file and a file containing a list of atom types synonims and create a new lammps data file including all needed types.')
    parser.add_argument('-lammps', metavar='lammps_data_file',help='lammps data files')
    parser.add_argument('--out', metavar='out_file',help='output file name, not needed for types and cassandra charges.', default="out")
    parser.add_argument('-format',metavar='out_file_format', help='Format to output. Can be data, xzy, pdb, ciff, mcf, mol, types_list, EDIP, cassandra_charges')
    parser.add_argument('--name', metavar='system_name',help='Name to assign to the molecular system', default=None)
    parser.add_argument('--no_interactive', metavar='no_interactive',help='To skip the interactive questions use this flag with etiher <element> <own> <PATH>, see docs for details. Will print connectivity for pdb files.', default=False)
    parser.add_argument('--no_inter_mcf_p', metavar='which_params_mcf',help='To skip the interactive questions use this flag with etiher <own> or <PATH>, see docs for details', default=False)
    parser.add_argument('--print_coeffs',help="If present attempts to print FF coefficients in output .data file.",action="store_true")
    parser.add_argument('--read_coeffs',help="If present attempts to read FF coefficients from a lammps input file. Expects <PATH>", default=False)

    
    args = parser.parse_args()
    lammps_f = args.lammps
    out_f =args.out
    the_format = args.format
    read_coeffs= args.read_coeffs
    
    
    formats = ["data", "xyz", "pdb", "ciff", "mcf", "mol", "types_list", "EDIP", "cassandra_charges"]
    if the_format not in formats:
        exit(f"Format {the_format} not recognized, please check spelling. Possible formats are: {formats}")
    
    if out_f == "out":
        out_f+="."+the_format
        
        
    lammps_orig_frag=LammpsData(lammps_f,name=args.name)
    
    if read_coeffs != False:
        lammps_orig_frag.add_Coeffs(read_coeffs)
    
    if args.no_interactive=="element":
        label=False
    elif args.no_interactive=="own":
        label=True
    elif not args.no_interactive:
        pass
    else:
        label=read_types(args.no_interactive)
        
    if args.no_inter_mcf_p=="own":
        params=False
    elif not args.no_inter_mcf_p:
        if the_format == "mcf" and args.no_interactive:
            exit("Mode set to not interactive, please use flag --no_inter_mcf_p to provide parameters.")
    else:
        params=read_types(args.no_inter_mcf_p)
    
    
    match the_format:
        
        case 'data':
            lammps_orig_frag.write_lammps_data(out_f, coeffs_out=args.print_coeffs)
        
        case  'xyz':
            if not args.no_interactive:
                inp = input("Do you want to use normal .xyz element labels? (If no, you can provide labels either in the data file or as an auxiliary map file) (y/n)\n")
                if inp == "y" :
                    label= False
                else:
                    inp = input("Input file path for types map, or write 'no' to attempt using labels from the data file\n")
                    if inp == "no":
                        label= True
                    else:
                        label= read_types(inp)
            lammps_orig_frag.write_xyz(out_f,a_which_labels=label)
        
        case  'pdb':
           if not args.no_interactive:
               inp = input("Do you want to print connectivity information? (y/n) \n")
               if inp == "n" :
                   con=False
               else:
                   con=True
           else:
               con=True
           lammps_orig_frag.write_pdb(out_f,connect=con)
           
        case 'types_list':
            lammps_orig_frag.write_types()
            print("File types.txt has been written")
            
        case 'cassandra_charges':
            lammps_orig_frag.write_cassandra_charges()
            
        case 'mol':
            lammps_orig_frag.write_molecule_template(out_f)
            
        case 'ciff':
            if not args.no_interactive:
                inp = input("Do you want to use element names as .ciff atom labels? (If no, you can provide labels either in the data file or as an auxiliary map file) (y/n)\n")
                if inp == "y" :
                    label= False
                else:
                    inp = input("Input file path for types map, or write 'no' to attempt using labels from the data file\n")
                    if inp == "no":
                        label= True
                    else:
                       label= read_types(inp)
            lammps_orig_frag.write_ciff(out_f,a_which_labels=label)   
        case 'mcf':
            if not args.no_interactive:

                inp = input("Do you want to use element names as .mcf atom labels? (If no, you can provide labels either in the data file or as an auxiliary map file) (y/n)\n")
                if inp == "y" :
                    label= False
                else:
                    inp = input("Input file path for types map, or write 'no' to attempt using labels from the data file\n")
                    if inp == "no":
                        label= True
                    else:
                        label= read_types(inp)
            
            if not args.no_inter_mcf_p:
      
                inp = input("Do you want to use LAMMPS data own LJ parameters if available? (If no, you can provide parameters in an auxiliary file) (y/n)\n")
                if inp == "y" :
                    print("WARNING! This will give reasonable results only if the LAMMPS data includes pair coeffs parameters for LJ related pair_styles.")
                    params= False
                else:
                    inp = input("Input file path for parameters\n WARNING! make sure the atom labels used for the atoms matches the ones in the params file!")
                    params= read_params(inp)
            lammps_orig_frag.write_rigid_mcf(out_f,a_which_labels=label, a_which_params=params)   

    
