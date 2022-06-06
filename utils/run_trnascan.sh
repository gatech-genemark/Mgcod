#!/bin/bash

infile=$1
output_dir=$2

if [[ ! -d ${output_dir} ]]
then
    mkdir -p ${output_dir}
fi

cat $infile | while read line
do
    out_base=$(echo $line | awk -F '/' '{print $NF}' | sed 's/.fna//')
    tRNAscan-SE -B -o ${output_dir}/${out_base}.out -m ${output_dir}/${out_base}.stats -a ${output_dir}/${out_base}.fasta $line
done
