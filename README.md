# Analyse and manipulate lammps data files

# 

# 

## LAMMPS\_DATA\_ANALYSIS

Lammps data analysis is a set of scripts designed to easily manage, manipulate and extract information from LAMMPS data files.

This set of tools allow you to examine the contents of a lammps data file in multiple ways. To remove and extract part of the data file. And to export the data file to a number of other formats.  
These were developed to ease the use of LAMMPS in connection to other applications (such as Cassandra, Polymatic, Gromacs…)  and to simplify the production of polymers with lammps bond/react and AutoMapper (for example to fix connectivity info that can sometimes be broken during the process).  
We concentrated our effort on the most common form of LAMMPS data files. However we tried to offer a variety of options, including the possibility of including personalized labels to various elements.

All the scripts can be used as stand alone or be called by other scripts; for this purpose refer to the documentation about the  `no_interactive` flag.

## Requirements and data preparation:

### Python env:

These scripts require the use of a Python 3.6 environment with `networkit` module installed, `Pandas` is recommended but not strictly needed.  
Remember that the environment must be activated before use (e.g “conda activate \<env\_name\>)

### LAMMPS data format:

The software accept LAMMPS data files of the following format:

* atom style: full  
* units: real  
* coordinates/angles: orthogonal or restricted triclinic  
* Coefficients headers in the form: \<Coeff type\> Coeffs (e.g. Bond Coeffs, Improper Coeffs,...). These sections can be omitted.


It is possible to insert personalized labels to each line in the sections: PairIJ Coeffs, Masses, Atoms, and all the coefficients sections (Pair, PairIJ, Bond, Angle, Dihedral and Improper Coeffs). Labels need to be given in the form:

`Masses`  
`1 1.001 # LABEL`  
`2 15.996 # LABEL2`  
`…`

The label can be of any size (provided it is contained in the same line) and can contain any characters (excluding \#)  
However, in order to export the data into .pdb format the labels length need to be of max 4 characters (the software will otherwise only use the first 4 characters).

## Functionalities

The functionalities are divided in 4 categories:

* Informations about the data file: handled by `lammps_data_info` returns a variety of useful values:  
* Other data formats in which the data can be exported: handled by `lammps_data_export` returns the following formats:  
* Functions to manipulate the data file in several ways: handled by `lammps_data_manipulate:`  
* Function to extract (from a data file) all fragments of molecules corresponding to a given structure: handled by `lammps_data_extract_fragment` This script can also iterate over all frames of a trajectory keeping track of fragments across the frames.

## General syntax for all scripts

`lammps_data_<script> -lammps <PATH_to_input>` 

Plus optional arguments depending on the specific script. See manual for details or call 

`lammps_data_<script> --help`  
