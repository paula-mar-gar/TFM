#!/bin/bash

for f in *.tab; do
	./adynFBA.R -m iPau21 	-f \
			-o \
                        -v \
                        -c "$f" \
                        --goal PA14_Biomass \
                        -V 6;
done

