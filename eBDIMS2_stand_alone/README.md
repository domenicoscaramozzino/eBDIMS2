#### eBDIMS2 CODE ####

Here we provide the codes to simulate transition pathways in large protein structures with eBDIMS2. Both the source code and binary codes already compiled with gfortran for LINUX and MacOs are provided.

#### SIMULATION ####

To run the eBDIMS2 analysis, you need to follow these steps:

1. Download the stand-alone code (LINUX or MacOS folders) in your working directory - or download the source code and compile it with gfortran, e.g. gfortran eBDIMS2.f90 -o eBDIMS2 -fopenmp

2. Include the PDB files of your end-state conformations in your working directory, e.g. 1ss8.pdb and 1sx4.pdb

3. Run the analysis by typing "./eBDIMS2 $PDB_ref $chain_ref $PDB_tar $chain_tar $frame_writing_frequency $DIMS_convergence" from your terminal. $PDB_ref and $PDB_tar are the name of your end-state PDB files (e.g. "1ss8" and "1sx4"), $chain_ref and $chain_tar are the chain labels of your protein model that are considered for the transition (e.g. "ABCDEFG" and "ABCDEFG"), $frame_writing_frequency represents the number of successful basing steps where an eBDIMS2 intermediate conformer is saved as a PDB file (e.g. "1000"), and $DIMS_convergence is the final requested convergence of the DIMS parameter (e.g. "99.9").

For example, if you want to analyze the GroEL transition, you can investigate the opening of the:
- Single monomer with: ./eBDIMS2 1ss8 A 1sx4 A 1000 99.9
- Single 7mer ring with: ./eBDIMS2 1ss8 ABCDEFG 1sx4 ABCDEFG 1000 99.9
- Double 7mer ring (14mer) with: ./eBDIMS2 1ss8 ABCDEFGHIJKLMN 1sx4 ABCDEFGHIJKLMN 1000 99.9

*Note that the correspondence between chain labels in multi-meric proteins is assessed based on the ordering of the chain entries, e.g. for "1ss8 ABCDEFG" and "1ss8 HIJKLMN", the transition will be run with these chain correspondences: A-H, B-I, etc.

**With this code you don't need to worry about having the same number of residues in both end-states. eBDIMS2 will find correspondences based on the numberings of residues.

***By default, the code will run on 16 OpenMP threads.

#### OUTPUTS ####

The code will generate intermediate eBDIMS2 conformers labeled with "DIMS_MD***.pdb" as well as log files reporting information about convergence and computing time.
