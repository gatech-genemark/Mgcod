#!/bin/bash

input=$1 # protein sequences
database=$2 # Pfam-A.hmm
output_base=$3
threads=$4

hmmsearch -o ${output_base}.out --tblout ${output_base}.tab --cpu ${threads} ${database} ${input}