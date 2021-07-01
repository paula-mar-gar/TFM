#!/bin/bash

well=$(grep -h -A4 "Model:" iPau*/*.Rlog | grep "Model" | awk -F "constrains_" '{print $NF}' | awk -F "." '{print $1}' | awk -F "_" '{print $2}')

plate=$(grep -h -A4 "Model:" iPau*/*.Rlog | grep "Model" | awk -F "constrains_" '{print $NF}' | awk -F "." '{print $1}' | awk -F "_" '{print $1}')

biomass=$(grep -h -A1 "^Objective" iPau*/*.Rlog | grep -v "Objective" | awk -F " " '{print $NF}' | sed '/^--$/d')

paste <(echo "$plate") <(echo "$well") <(echo "$biomass") --delimiters '\t' > ../biomass.txt
