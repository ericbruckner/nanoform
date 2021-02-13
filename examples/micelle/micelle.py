import numpy as np
import pandas as pd

from nanoform import main

# Make Micelle

## Filename for pdb file for MARTINI simulation
pdb = 'PA.pdb'

## Filename for corresponding MARTINI topology file
topo = 'PA.top'

## Filename simulation box .gro file
output_filename = 'PA_micelle.gro'

## Set Parameters

Lx = 8 # Box length of x dimension (nm)
Ly = 8 # Box length of y dimension (nm)
Lz = 8 # Box length of z dimension (nm)

n_radial_cross_section = 16 # number of molecules in a radial cross-section 
offset = 4 # offset of molecules from fiber core (angstroms)

main.pdb2micelle(pdb, output_filename, topfilename = topo, Lx = Lx, Ly = Ly, Lz = Lz, n_radial_cross_section = n_radial_cross_section, offset = offset)