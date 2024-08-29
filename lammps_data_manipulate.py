#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:12:44 2022

@author: alberto_z
"""
#%% Importing FragmentInfo from other script
from lammps_data_analysis import LammpsData
    
#%% This will later be the body of the script.
import argparse

gTest = False

if __name__ == "__main__" and not gTest:

    parser = argparse.ArgumentParser(description=' Reads in a lammps data file and performs various manipulations.')
    parser.add_argument('-lammps', metavar='lammps_data_file',help='lammps data files')
    parser.add_argument('--out', metavar='out_file',help='lammps molecule output files', default="out.data")
    parser.add_argument('--name', metavar='system_name',help='Name to assign to the molecular system', default=None)
    # parser.add_argument('--no_interactive', metavar='no_interactive',help='To skip the interactive questions use this flag with etiher <element> <own> <PATH>, see docs for details. Will print connectivity for pdb files.', default=False)
    parser.add_argument('--print_coeffs',help="Prints FF coefficients in output .data file.",action="store_true")
    parser.add_argument('--insert_coeffs',help="copies coefficients data from lammps input file to .data file. Expects <path/to/input/file>",default=None)
    parser.add_argument('--renumber',help='Renumbers atoms and molecules for consistency',action='store_true')
    parser.add_argument('--reorder',help='Reorder atoms by atom index',action='store_true')
    parser.add_argument('--reassign_mols',help='Reassign molecules dependeing on connectivity',action='store_true')
    parser.add_argument('--remap_types',help='Remaps types based on a types_map file, does not change the total number of types', default=None)
    parser.add_argument('--zero_tot_q',help="If present shift all charges to achieve zero total charge", action="store_true")
    parser.add_argument('--qshift', metavar ='q_shift_list',help='list of [atom_type_A, charge shift_A, atom_type_B, charge shift_B,...] to be performed, comma separated',default=None)

    parser.add_argument('--remove', help='Remove atoms, molecules or bonds from lammps data file, according to following flags',action='store_true')
    parser.add_argument('--extract',help='Create lammps data file containing only selected atoms or molecules from lammps data file, according to following flags',action='store_true')
    parser.add_argument('--extr_separate',help='Create a separate lammps data file containing only one selected molecule each, according to following flags',action='store_true')

    parser.add_argument('--number', metavar ='number_list',help='list of atom numbers to be removed, comma separated. Ranges expressed as n_start-n_end.',default=None)
    parser.add_argument('--type',metavar='type_list',help='list of atom types to be removed from system, numerical and comma separated',default=None)
    parser.add_argument('--mol',metavar='mol_list',help='list of molecules numbers to be removed from system, numerical and comma separated. Ranges expressed as n_start-n_end.',default=None)
    parser.add_argument('--molsize',metavar="mol_size",help='Size of molecules to be removed (# of atoms).',default=None)
    parser.add_argument('--bond',metavar="bond",help='Atom pairs between wich to remove bonds. Atoms in pair dash separated, pairs comma separated (e.g. a1-a2, a6-a8)',default=None)

    
    args = parser.parse_args()
    lammps_f = args.lammps
    out_f =args.out
    coeff_inp_file = args.insert_coeffs
    
    renumber=args.renumber
    reorder=args.reorder
    reassign_mols=args.reassign_mols
        
    zero_tot_q = args.zero_tot_q
    qshift = args.qshift
  
    
        
    remap = False
    if args.remap_types != None:
        map_types_f = args.remap_types
        remap=True
        
    lammps_orig_frag=LammpsData(lammps_f,name=args.name)
    
    if coeff_inp_file != None:
        lammps_orig_frag.add_Coeffs(coeff_inp_file)
        lammps_orig_frag.write_lammps_data(out_f,coeffs_out=True)
        
    if renumber or reorder or reassign_mols:
        print(f"renumber {renumber}, reorder {reorder}, reassign_mols {reassign_mols}")
        if renumber:
            lammps_orig_frag.renumber()
        if reorder:
            lammps_orig_frag.reorder()
        if reassign_mols:
            lammps_orig_frag.reassign_molecules()
            
        lammps_orig_frag.write_lammps_data(out_f,coeffs_out=args.print_coeffs)
        
    if remap:
        lammps_orig_frag.remap_types(map_types_f)
        lammps_orig_frag.write_lammps_data(out_f,coeffs_out=args.print_coeffs)   

    if zero_tot_q:
        lammps_orig_frag.shift_total_charge_to_zero()
        if lammps_orig_frag.name == None:
            lammps_orig_frag.write_lammps_data("data_zero_tot_q.lmps", coeffs_out=args.print_coeffs)
        else:
            lammps_orig_frag.write_lammps_data(f"{lammps_orig_frag.name}_zero_tot_q.lmps", coeffs_out=args.print_coeffs)   
            
    if args.qshift != None:
        
        q_t_list = qshift.split(",")
        if len(q_t_list)%2 ==1:
            exit("qshift: Number of charges and types does not match.")
        for c in range(0,len(q_t_list)//2):
            t_c = int(q_t_list[2*c])
            q_c = float(q_t_list[2*c+1])
            lammps_orig_frag.shift_atom_type_charge(t_c,q_c)
        
        lammps_orig_frag.write_lammps_data(out_f,coeffs_out=args.print_coeffs)
    
    if args.remove or args.extract or args.extr_separate:
        
        atom_num_ls = []
        if args.extr_separate and (args.number != None or args.type != None ):
            exit("extr_separate only extract separate molecules. Use --extract or lammps_extract_fragment.py depending on your needs.")
        if (args.extr_separate or args.extract) and  args.bond != None:
            exit("Cannot extract bonds, did you mean remove?")
        if args.number != None:
            # print("in numbers")
            if "," in args.number:
                for n in args.number.split(","):
                    if "-" in n:
                        atom_num_ls += [m for m in range(int(n.split("-")[0]), int(n.split("-")[1])+1)]
                    else:
                        atom_num_ls.append(int(n))
            else:
                if "-" in args.number:
                    atom_num_ls += [n for n in range(int(args.number.split("-")[0]), int(args.number.split("-")[1])+1)]
                else:
                    atom_num_ls.append(int(args.number))
            # print("removing: ",atom_num_ls)
    
        if args.type != None:
            # print("in types")
            ls_types = [int(n) for n in args.type.split(",")]
            at_ls_from_types = lammps_orig_frag.atom_IDs_for_types(ls_types)
            #print("in remove by types: ", at_ls_from_types)
            atom_num_ls += at_ls_from_types
            # print("removing: ",atom_num_ls)
            
        if args.mol != None:
            # print("in mol")
            ls_mols =  []
            if "," in args.mol:
                for n in args.mol.split(","):
                    if "-" in n:
                        ls_mols += [int(m) for m in range(int(n.split("-")[0]), int(n.split("-")[1])+1)]
                    else:
                        ls_mols.append(int(n))
            else:
                if "-" in args.mol:
                    ls_mols += [int(n) for n in range(int(args.mol.split("-")[0]), int(args.mol.split("-")[1])+1)]
                else:
                    ls_mols.append(int(args.mol))
            # print(f"ls_mols: {ls_mols}")
            #ls_mols = [n for n in args.mol.split(",")]
            at_ls_from_mols = lammps_orig_frag.atom_nbs_for_mols(ls_mols)
            atom_num_ls += at_ls_from_mols
            # print("removing: ",atom_num_ls)
             
        if args.molsize != None:
            # print(f"Selecting all molecules with {args.molsize} atoms.")
            n_mols = 0
            # print(f"Number of molecules at start: {len(lammps_orig_frag.atoms_per_mol)}")
            for mol in lammps_orig_frag.atoms_per_mol.values():
                if len(mol) == int(args.molsize):
                    atom_num_ls += mol
                    n_mols += 1
            # print("removing: ",atom_num_ls)
            # print("Number of molecules selected: ",n_mols_selected)
        
        if args.remove:
            # print("removing: ",atom_num_ls)
            if args.bond != None:
                # print("in bond")
                bond_ls = args.bond.split(",")
                tempsys = LammpsData.remove_bonds(lammps_orig_frag,bond_ls)
                
            try:
                outsys = LammpsData.remove_atoms(tempsys,atom_num_ls)
            except:
                outsys = LammpsData.remove_atoms(lammps_orig_frag,atom_num_ls)
            outsys.reassign_molecules()
            outsys.write_lammps_data(out_f,coeffs_out=args.print_coeffs)
            
        elif args.extract:
            # print("extracting: ",atom_num_ls)
            remove_ls=[]
            for i_at,atom in lammps_orig_frag.Atoms.items():
                if i_at not in atom_num_ls:
                    remove_ls.append(i_at)
            outsys = LammpsData.remove_atoms(lammps_orig_frag,remove_ls)
            outsys.reassign_molecules()
            outsys.write_lammps_data(out_f,coeffs_out=args.print_coeffs)
            
        elif args.extr_separate:
            remove_ls=[]
            for i_at,atom in lammps_orig_frag.Atoms.items():
                if i_at not in atom_num_ls:
                    remove_ls.append(i_at)
            outsys = LammpsData.remove_atoms(lammps_orig_frag,remove_ls)
            for mol in outsys.atoms_per_mol.keys():
                remove_ls_2=[]
                for i_at,atom in outsys.Atoms.items():
                    if atom["mol"] != mol:
                        remove_ls_2.append(i_at)
                outsys2 = LammpsData.remove_atoms(outsys,remove_ls_2)
                outsys2.write_lammps_data(out_f.split(".")[0]+"_"+str(mol)+".data",coeffs_out=args.print_coeffs)
                        

