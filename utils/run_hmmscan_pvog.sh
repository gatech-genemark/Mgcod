#!/bin/bash

while getopts i:o:d:t:c: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
        d) database=${OPTARG};;
        t) table_output=${OPTARG};;
        c) cpus=${OPTARG};;
    esac
done

input_centroids=$( echo ${input} | sed 's/.faa/_centroids.faa/' )

python utils/adjust_id.py ${input}

cat ${input} | tr '[:lower:]' '[:upper:]' > ${input}1

usearch11.0.667_i86linux32 --cluster_fast ${input}1 -id 0.9 -centroids ${input_centroids}

rm ${input}1

hmmscan -o ${output} --tblout ${table_output} --max --cpu ${cpus} ${database} ${input_centroids}

