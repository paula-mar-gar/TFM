cat Caz5d_mutated_genes_in_model.txt \
| while read g ; do 
	echo -n "$g	" 
        if ! grep "^Gene $g has no effect" $1 ; then 
            echo "KO" ; 
        fi
        grep "$g KO" $1 1>&2
  done > $1.genes 2> $1.genes.obj

