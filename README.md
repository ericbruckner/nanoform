# nanoform
A package for creating initial geometries of lipopeptide nanostructures for MARTINI simulations

    Built using python3.7

## Software prerequisites:
- Python3 (Tested on version 3.7.6)
- numpy (Tested on version 1.17.0)
- pandas (Test on version 1.1.3)
- Recommended visualization in VMD (http://www.ks.uiuc.edu/Research/vmd/)

## Description
This package contains Python scripts to generate initial geometries for MARTINI simulations of lipopeptide nanostructures (e.g. micelles, fibers, bilayers) to be run in GROMACS simulations. Initial geometries are generated from the .pdb and .top files generated by martinize.py (http://cgmartini.nl/index.php/tools2/proteins-and-bilayers/204-martinize)

main.py is the main script that lets you do all the above as shown by the examples

## Examples
An example on how to generate a fiber morphology. See "example" directory for how to make micelle and bilayer geometeries

### Import Libraries
    import numpy as np
    import pandas as pd
    from nanoform import main

### Filename for pdb file for MARTINI simulation
    pdb = 'PA.pdb'

### Filename for corresponding MARTINI topology file
    topo = 'PA.top'

### Filename for simulation box .gro file
    output_filename = 'PA_fiber.gro'

### Set Parameters

    Lx = 8 # Box length of x dimension (nm)
    Ly = 8 # Box length of y dimension (nm)
    Lz = 8 # Box length of z dimension (nm)

    n_radial_cross_section = 16 # number of molecules in a radial cross-section 
    stacking_distance = 3.5 # distance between molecules along fiber length (angstroms)
    offset = 4 # offset of molecules from micelle core (angstroms)

### Generate Structure
    main.pdb2fiber(pdb, output_filename, topfilename = topo, Lx = Lx, Ly = Ly, Lz = Lz, n_radial_cross_section = n_radial_cross_section, stacking_distance = stacking_distance, offset = offset)

### Sample Output
<img src="https://github.com/ericbruckner/nanoform/blob/main/examples/fiber/sample_output/fiber_highres.png" height="200">
