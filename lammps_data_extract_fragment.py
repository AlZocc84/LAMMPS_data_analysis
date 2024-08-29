#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:12:44 2022

@author: alberto_z
"""
#%% Importing LammpsData from other script
from lammps_data_analysis import LammpsData
from os.path import exists
from os import mkdir

#%%
def write_frag_traj(a_frags_list,a_traj_f,a_out_dir,a_out_prefix):
    
    
      
    this_frame = -1
    frames_frags = {}
    with open(a_traj_f,"r") as trr_f:
        for i,line in enumerate(trr_f.readlines()):
            if "Timestep" in line:
                this_frame +=1
                frames_frags[this_frame]={}
                start_line=i
                for i_frag,frag in enumerate(a_frags_list):
                    frames_frags[this_frame][i_frag]=[]
            elif len(line.split())>1:
                #Now looking for frag atoms:  
                    for i_frag,frag in enumerate(a_frags_list):
                        if i-start_line in frag:
                           frames_frags[this_frame][i_frag].append(line)
                       
    for i in range(0,len(a_frags_list)):
        if not exists(dir_s+f"/frag{i+1}"):
            mkdir(dir_s+f"/frag{i+1}") 
        for i_frame in range(0,this_frame+1):
            frag_xyz_data=frames_frags[i_frame][i]
            with open(dir_s+f"/frag{i+1}"+"/"+out_prefix+f"_frag{i+1}_{i_frame+1}.xyz","w") as xyz_f:
                xyz_f.write(str(len(frag_xyz_data))+"\n")
                xyz_f.write(f"frame {i_frame+1}"+"\n")
                for at_line in frag_xyz_data:
                    xyz_f.write(at_line)
                       
    print("All written.")
            
#%% Body of the script.
import argparse

gTest = False

if __name__ == "__main__" and not gTest:

    parser = argparse.ArgumentParser(description='Reads in a lammps data file representing a molecular system, and a second file representig a molecular fragment to be found and extracted from the first. \nCan read in also a trajectory file to extract fragments from all frames.')
    parser.add_argument('-lammps', metavar='lammps_data_file',help='lammps data files')
    parser.add_argument('-fragment',metavar='frag_data_file',help=' Fragment(s) to be extracted from system, lammps data file. Atom numbers are expected to be contigous and start from 1.')
    parser.add_argument('--start',metavar='start_at',help=' Starting atom for map creation.',default=None)
    parser.add_argument('--traj', metavar='xyz_trj_file',help='lammps xyz trajectory files',default=None)
    parser.add_argument('--dir', metavar='out_dir',help='output directory', default="fragments")
    parser.add_argument('--out', metavar='out_file_prefix',help='output prefix', default="fragment")
    parser.add_argument('--print_coeffs',help="If present attempts to print FF coefficients in output .data file.",action="store_true")
    parser.add_argument('--name', metavar='system_name',help='Name to assign to the molecular system', default=None)


    print("Warning: The current implementation is likely to overcount symmetric fragments.")
    args = parser.parse_args()
    lammps_f = args.lammps
    out_f =args.out
    fragment_f = args.fragment
    traj_f=args.traj
    dir_s=args.dir
    
    if not exists(dir_s):
        mkdir(dir_s)
    
    out_prefix= args.out
    
    whole_sys=LammpsData(lammps_f,name=args.name)
    
    fragment = LammpsData(fragment_f)
    if args.start == None:
        fragment.renumber()
        start_at=1
    else:
        start_at= int(args.start)
        if start_at > fragment.n_atoms:
            exit("start atom must be <= nb. atoms in fragment")
    # print("fragment.bonds_per_at: ", fragment.bonds_per_at)
    #print("whole_sys.bonds_per_at: ", whole_sys.bonds_per_at)
    
    fragment_map = fragment.fragment_map(start_at)
    
    frag_list = whole_sys.find_fragments(fragment_map)
    # print(frag_list)
    # print(len(frag_list))
    #whole_sys.write_fragments(frag_list)
    
    for i,list_remove in enumerate(whole_sys.list_unwanted_atoms(frag_list)):
        print(f"Printing the files for fragment frame {i}")
    
        new_lammps_frag=LammpsData.remove_atoms(whole_sys,list_remove,renumber=False)
        new_lammps_frag.write_lammps_data(dir_s+"/"+out_prefix+f"_{i}.data",coeffs_out=args.print_coeffs)
    if traj_f != None: 
        
        write_frag_traj(frag_list,traj_f,dir_s,out_prefix)

    exit("Done")




