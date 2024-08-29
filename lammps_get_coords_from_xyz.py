#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:02:26 2024

@author: alberto
"""

#%% Importing FragmentInfo from other script
from lammps_data_analysis import LammpsData
    
#%% This will later be the body of the script.
import argparse

gTest = False

if __name__ == "__main__" and not gTest:

    parser = argparse.ArgumentParser(description=' Reads in a lammps data file and a xyz file and writes a new lammps data file with the coordinates from the xyz file.')
    parser.add_argument('-lammps', metavar='lammps_data_file',help='lammps data files')
    parser.add_argument('--out', metavar='out_file',help='lammps molecule output files', default="out.data")
    parser.add_argument('-xyz', metavar='xyz_file',help='xyz file for coordinates')
    
    args = parser.parse_args()
    lammps_f = args.lammps
    out_f =args.out
    xyz_file = args.xyz
    
    lammps_data=LammpsData(lammps_f)
    lammps_data.read_coords_from_xyz(xyz_file)
    lammps_data.write_lammps_data(out_f)
