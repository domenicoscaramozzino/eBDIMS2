#### PCA-eBDIMS2 ####

Here we provide the codes to perform Principal Component Analysis (PCA) and simulate transition pathways in large protein structures with eBDIMS2. Binary codes are available both for LINUX and MacOs from the respective folders.

#### SIMULATION ####

To run the complete PCA-eBDIMS2 analysis, you need to follow these steps:

1. Download all codes (LINUX or MacOS folders) in your working directory. These folders include executable codes and the "eBDIMS2_pipeline_par.sh" script. This one is needed to run the whole pipeline, perform structural alignment with Gromacs*, and call all executables.

2. Include all PDB files in your working directory, as well as the text files that specify the structures in the ensemble ("ensemble.txt") and eBDIMS2 transitions ("eBDIMS2_transitions.txt"). In the "example_ACLY_input_files" folder, you can find all input files required. In "ensemble.txt", the first column must contain the 4-digit PDB code (e.g. "6pof") and the second space-separated column the corresponding chain labels ("ABCD"). To perform PCA, the correspondence between chains in multi-meric proteins is assessed based on the ordering of the chain labels in the chain entries, e.g. for "6pof ABCD" and "6hxh EFGH", chains will be aligned with these correspondences: A-E, B-F, etc. The first row in the "ensemble.txt" file must correspond to the reference structure used for alignment and PCA projections! In the "eBDIMS2_transitions.txt", each row has four space-separated columns: 4-digit PDB code of end-state 1 ("6pof"), chain labels of end-state 1 ("ABCD"), PDB code of end-state 2 ("6hxh"), chain labels of end-state 2 ("ABCD"). By inserting multiple rows, you will run multiple transitions within the same ensemble, e.g. "6uia ABCD 6o0h ABCD". For each row, both the forward ("6pof>6hxh") and backward ("6hxh>6pof") transitions will be computed.

3. Insert your input parameters in the "job.sh" script (you can find it in the "example_ACLY_input_files" folder): (1) PDB code of the reference structure (used for alignment and PC projections), e.g. "6pof"; (2) chain labels of the reference PDB, e.g. "ABCD". The PBD ID and chain labels of the reference structures must correspond to the first entry in the "ensemble.txt" file - so make sure you carefully select your reference structure!; (3) Frequency for output writing of eBDIMS2 conformers - for most transitions 1000 will work. If you want more frames, reduce this number, e.g. 100; if you want less frames, increase it, e.g. 10000; (4) DIMS convergence - for most transitions 99.9% will work. If you want a closer convergence (slower), increase it, e.g. 99.99%; if you are satisfied with lower convergence (faster), decrease it, e.g. 99%; (5) RMSD convergence - for most transitions, a default value of 0.8-1.0A will work.

4. Run the analysis typing "bash job.sh" from your terminal from your working directory.

*Since we use "gmx confmrs" to perform structural alignment, you will need Gromacs installed in your machine.

#### OUTPUTS ####

The pipeline will get all the C-alpha atoms* of your PDB models, perform structural alignments with respect to the reference PDB, perform PCA, simulate eBDIMS2 transition pathways, and project them on the PCA space. 

This will generate the following outputs:

1. CA-only and CA-only-aligned PDB models for the structural ensemble;
2. Output files from PCA (PDB model of average coordinates, variances captured by the PCs, first 10 PC eigenvectors, projection of experimental structures on the first 10 PCs);
3. Multi-model PDB files of all eBDIMS2 transitions, e.g. "transition_6pof_ABCD_6hxh_ABCD.pdb"
4. Folders for each eBDIMS2 transition path, e.g. "transition_6pof_ABCD_6hxh_ABCD", which contain: CA-only eBDIMS2 frames, two log files with information about transition convergence and computing time, and projection of each intermediate frame along the first 10 experimental PCs.

*Note that the pipeline will read only ATOM records of C-alpha atoms (with altLoc empty or equal to "A"). Before running the simulation you need to ensure that the ensemble is consistent, i.e. all models with the same number of residues, proper sequence and chain correspondences.

