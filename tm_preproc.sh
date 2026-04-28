#!/bin/bash

### Matej Kožić | mkozic@chem.pmf.hr | 2022.11.25 15:00
### Bash script for preprocessing trajectories into a .csv matrix using TM.py
#############################################################################
#INPUTS:

topology="12_md20_2_0.67_SpoIIIEinserted_2A_10ns.gro"      #your topology file (gro, pdb, prmtop...)
trajectories="12_md20_2_0.67_SpoIIIEinserted_2A_10ns.xtc"  #your trajectory file (xtc, nc)
stride="1"			     #trajectory reading stride
saved_pdb="InsertedSpoIIIEtraj.pdb"	     #name of the created pdb
saved_csv="InsertedSpoIIIEtraj.csv"            #name of the created .csv matrix
residues="2208"                     #real number of residues (excluding missing)

#############################################################################
#############################################################################

python3 TM.py << INPUTS
traj2pdb
$topology
$trajectories
$stride
$saved_pdb
pdb2csv
$saved_pdb
$saved_csv
$residues
q
INPUTS

echo If successful trajectories processed and saved as $saved_csv

