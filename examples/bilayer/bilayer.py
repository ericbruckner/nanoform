import numpy as np
import pandas as pd

from nanoform import main

# Make Bilayer

## Filename for pdb file for MARTINI simulation
pdb = 'PA.pdb'

## Filename for corresponding MARTINI topology file
topo = 'PA.top'

## Filename simulation box .gro file
output_filename = 'PA_bilayer.gro'

## Set Parameters

Lx = 8 # Box length of x dimension (nm)
Ly = 8 # Box length of y dimension (nm)
Lz = 8 # Box length of z dimension (nm)

stacking_distance = 3.5 # distance between molecules along bilayer axes (angstroms)
offset = 4 # offset of molecules from fiber core (angstroms)

main.pdb2bilayer(pdb, output_filename, topfilename = topo, Lx = Lx, Ly = Ly, Lz = Lz, stacking_distance = stacking_distance, offset = offset)