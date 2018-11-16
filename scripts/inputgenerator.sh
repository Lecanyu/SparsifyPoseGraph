#!/bin/bash

DATASETS=( \
	datasets/duderstadt.g2o \
	datasets/flat_duderstadt.g2o \
	datasets/intel.g2o \
	datasets/manhattan.g2o \
	datasets/sphere.g2o \
	datasets/eecs.g2o \
	datasets/flat_eecs.g2o \
	datasets/killian.g2o \
	datasets/parking_fixed.g2o \
)

KLDPERIODS=( 10 10 10 50 40 10 10 50 20 )

OUTPUT=input.txt

echo "# Parameters, in order:" > $OUTPUT
echo "#   - Algorithm (glc|[^ ]+)" >> $OUTPUT
echo "#   - Filename .g2o" >> $OUTPUT
echo "#   - Sparsification scenario (online|cluster|global)" >> $OUTPUT
echo "#   - Topology (tree|subgr|clsubgr|dense|cldense)" >> $OUTPUT
echo "#   - Linearization point (local|global)" >> $OUTPUT
echo "#   - Sparsity level ([0-9]+)" >> $OUTPUT
echo "#   - Optionally, KLD computation period ([0-9]+), default 10" >> $OUTPUT
echo "#   - Optionally, cluster size ([0-9]+), default 100" >> $OUTPUT

function print_foreach_configuration {
	for ((i=0; i<${#DATASETS[@]}; i++)); do
		FNAME=${DATASETS[$i]}
		KLDPERIOD=${KLDPERIODS[$i]}
		for SPARSITY in 2 3 4 5; do
			echo $1 $FNAME $2 $3 $4 $SPARSITY $KLDPERIOD >> $OUTPUT
		done
	done
}

print_foreach_configuration glc online tree global
print_foreach_configuration glc online dense global
print_foreach_configuration glc cluster tree global
print_foreach_configuration glc cluster dense global
print_foreach_configuration glc global tree global
print_foreach_configuration glc global dense global

print_foreach_configuration sen online tree global
print_foreach_configuration sen online subgr global
print_foreach_configuration sen online clsubgr global
print_foreach_configuration sen online dense global
print_foreach_configuration sen online cldense global
print_foreach_configuration sen online tree local
print_foreach_configuration sen online subgr local
print_foreach_configuration sen online clsubgr local
print_foreach_configuration sen online dense local
print_foreach_configuration sen online cldense local

print_foreach_configuration sen cluster tree global
print_foreach_configuration sen cluster subgr global
print_foreach_configuration sen cluster clsubgr global
print_foreach_configuration sen cluster cldense global
print_foreach_configuration sen cluster tree local
print_foreach_configuration sen cluster subgr local
print_foreach_configuration sen cluster clsubgr local
print_foreach_configuration sen cluster cldense local

print_foreach_configuration sen global tree global
print_foreach_configuration sen global subgr global
print_foreach_configuration sen global clsubgr global
print_foreach_configuration sen global cldense global
print_foreach_configuration sen global tree local
print_foreach_configuration sen global subgr local
print_foreach_configuration sen global clsubgr local
print_foreach_configuration sen global cldense local
