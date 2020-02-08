#!/bin/sh

INPUT=$1
#OUTPUT=$(basename $INPUT)
#OUTPUT=${OUTPUT%%.*}

sed '/NA/d' ${INPUT} | sed '/^26/d' | sort -k1,1n -k2,2n -k3,3 > ${INPUT}.sorted
mv ${INPUT}.sorted ${INPUT}
