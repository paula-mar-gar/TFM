#!/bin/bash

for f in iPau*; do
	well=$(echo "$f" | awk -F "_" '{print $3}')
	plate=$(echo "$f" | awk -F "_" '{print $4}')
	biomass_wt=$(grep -A1 "Objective function" $f/*.Rlog | awk 'FNR==2{print $2}')
	paste <(echo "$plate") <(echo "$well") <(echo "wt") <(echo "$biomass_wt") >> biomass_mutations.txt
	echo "$well" "$plate"

	for sgd in $f/singleGeneDeletion/*.txt; do
		mut_gene=$(echo "$sgd" | grep -o "G_PA14_[0-9]*")
		biomass_mut=$(cat "$sgd" | grep "PA14_Biomass" | awk '{print $2}')
		paste <(echo "$plate") <(echo "$well") <(echo "$mut_gene") <(echo "$biomass_mut") >> biomass_mutations.txt
		echo "$mut_gene"
	done
done
