#!/bin/bash

while getopts i:o:t: flag
do
    case "${flag}" in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
        t) table_output=${OPTARG};;
    esac
done

database=/storage3/w/aaron/pvogs_hmm/AllvogHMMprofiles/all_vogs.hmm
input_centroids=$( echo ${input} | sed 's/.faa/_centroids.faa/' )

python utils/adjust_id.py ${input}

cat ${input} | tr '[:lower:]' '[:upper:]' > ${input}1

usearch11.0.667_i86linux32 --cluster_fast ${input}1 -id 0.9 -centroids ${input_centroids}

rm ${input}1

hmmscan -o ${output} --tblout ${table_output} --max --cpu 8 ${database} ${input_centroids}

