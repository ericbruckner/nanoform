import numpy as np
import pandas as pd

from nanoform import main

# Make Fiber

## Filename for pdb file for MARTINI simulation
pdb = 'PA.pdb'

## Filename for corresponding MARTINI topology file
topo = 'PA.top'

## Filename simulation box .gro file
output_filename = 'PA_fiber.gro'

## Set Parameters

Lx = 8 # Box length of x dimension (nm)
Ly = 8 # Box length of y dimension (nm)
Lz = 8 # Box length of z dimension (nm)

n_radial_cross_section = 16 # number of molecules in a radial cross-section 
stacking_distance = 3.5 # distance between molecules along fiber length (angstroms)
offset = 4 # offset of molecules from micelle core (angstroms)

main.pdb2fiber(pdb, output_filename, topfilename = topo, Lx = Lx, Ly = Ly, Lz = Lz, n_radial_cross_section = n_radial_cross_section, stacking_distance = stacking_distance, offset = offset)