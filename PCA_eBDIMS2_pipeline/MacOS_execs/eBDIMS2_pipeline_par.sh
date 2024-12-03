#!/bin/bash

ref=$1
chains_ref=$2
save_frame_freq=$3
DIMS_conv_perc=$4
RMSD_conv=$5

num_chains=${#chains_ref}
pdb_suffix=".pdb"

#Take only CA atoms from the ref struct
./write_CA $ref $chains_ref

#Read structures in the ensemble and align them to CA-only reference structure 
ref_pca=$ref"_"$chains_ref"_CA.pdb"
n_CA=$(wc -l < $ref_pca)
while read p; do
	struct=${p% *}
	chains=${p#* }
	./write_CA $struct $chains
	yes 0 | gmx confrms -f1 $ref_pca -f2 $struct"_"$chains"_CA.pdb" -one -o temp_file.pdb
	grep -a "CA" temp_file.pdb > $struct"_"$chains"_CA_ali.pdb"
	rm temp_file.pdb
	n_CA_other=$(wc -l < $struct"_"$chains"_CA_ali.pdb")
	if test $n_CA_other -ne $n_CA
	then
		echo $struct"_"$chains"_CA.pdb has a different number of CAs - check it!"
		exit 0
	fi
	./add_chain_label_ali $struct"_"$chains"_CA" $struct"_"$chains"_CA_ali"
done <ensemble.txt

#Perform PCA of the experimental structures
n_exp_conf=$(wc -l < ensemble.txt)
echo "Total number of conformations in the ensemble: $n_exp_conf"
echo " "
echo "Running PCA ... "
./run_pca_symm_red $n_CA $n_exp_conf "A"
echo " "
echo "PCA completed!"
echo " "

while read row; do

	pdb_ref="${row:0:4}"
	chain_ref="${row:5:$num_chains}"
	pdb_tar="${row:6+$num_chains:4}"
	chain_tar="${row:11+$num_chains:$num_chains}"

	#Run eBDIMS (forward transition)
	./eBDIMS2_par $pdb_ref"_"$chain_ref"_CA_ali_labeled.pdb" $pdb_tar"_"$chain_tar"_CA_ali_labeled.pdb" $save_frame_freq $DIMS_conv_perc $RMSD_conv

	#Align eBDIMS frames (forward transition) to the reference structure
	while read p; do
		p_out=${p%"$pdb_suffix"}
		p_out+="_ali.pdb"
		yes 0 | gmx confrms -f1 $ref_pca -f2 $p -one -o temp_file.pdb
		rm $p
		grep -a "CA" temp_file.pdb > $p_out
		rm temp_file.pdb
		./add_chain_label_ali ${ref_pca%"$pdb_suffix"} ${p_out%"$pdb_suffix"}
	done <eBDIMS_confrms_list.txt
	n_DIMS_conf=$(wc -l < eBDIMS_confrms_list.txt)

	#Project forward eBDIMS (aligned) frames on the PC space
	./DIMS_projection_PC $n_CA $n_DIMS_conf

	#Move all the data related to the forward transition in a "ref_tar" folder
	transition_forward="transition_"$pdb_ref"_"$chain_ref"_"$pdb_tar"_"$chain_tar".pdb"
	touch $transition_forward
	eBDIMS_forward="transition_"$pdb_ref"_"$chain_ref"_"$pdb_tar"_"$chain_tar
	mkdir $eBDIMS_forward
	mv log.txt $eBDIMS_forward
	mv log_time.txt $eBDIMS_forward
	mv DIMS_ensemble_proj_PC.txt $eBDIMS_forward
	while read p; do
		p_out=${p%"$pdb_suffix"}
		p_out+="_ali_labeled.pdb"
		echo "MODEL" >> $transition_forward
		cat $p_out >> $transition_forward
		echo "TER" >> $transition_forward
		echo "ENDMDL" >> $transition_forward
		mv $p_out $eBDIMS_forward
	done <eBDIMS_confrms_list.txt
	mv eBDIMS_confrms_list.txt $eBDIMS_forward

	#Run eBDIMS (reverse transition)
	./eBDIMS2_par $pdb_tar"_"$chain_tar"_CA_ali_labeled.pdb" $pdb_ref"_"$chain_ref"_CA_ali_labeled.pdb" $save_frame_freq $DIMS_conv_perc $RMSD_conv

	#Align eBDIMS frames (reverse transition) to the reference structure
	while read p; do
	 	p_out=${p%"$pdb_suffix"}
	  	p_out+="_ali.pdb"
		yes 0 | gmx confrms -f1 $ref_pca -f2 $p -one -o temp_file.pdb
		rm $p
		grep -a "CA" temp_file.pdb > $p_out
		rm temp_file.pdb
		./add_chain_label_ali ${ref_pca%"$pdb_suffix"} ${p_out%"$pdb_suffix"}
	done <eBDIMS_confrms_list.txt
	n_DIMS_conf=$(wc -l < eBDIMS_confrms_list.txt)

	#Project reverse eBDIMS (aligned) frames on the PC space
	./DIMS_projection_PC $n_CA $n_DIMS_conf

	#Move all the data related to the reverse transition in a "tar_ref" folder
	transition_reverse="transition_"$pdb_tar"_"$chain_tar"_"$pdb_ref"_"$chain_ref".pdb"
	touch $transition_reverse
	eBDIMS_reverse="transition_"$pdb_tar"_"$chain_tar"_"$pdb_ref"_"$chain_ref
	mkdir $eBDIMS_reverse
	mv log.txt $eBDIMS_reverse
	mv log_time.txt $eBDIMS_reverse
	mv DIMS_ensemble_proj_PC.txt $eBDIMS_reverse
	while read p; do
		p_out=${p%"$pdb_suffix"}
		p_out+="_ali_labeled.pdb"
		echo "MODEL" >> $transition_reverse
		cat $p_out >> $transition_reverse
		echo "TER" >> $transition_reverse
		echo "ENDMDL" >> $transition_reverse
		mv $p_out $eBDIMS_reverse
	done <eBDIMS_confrms_list.txt
	mv eBDIMS_confrms_list.txt $eBDIMS_reverse
	
done <eBDIMS_transitions.txt
