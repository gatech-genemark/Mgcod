# Brief descriptions of supporting scripts for the data analysis following Mgcod runs

For each script, a brief description of what it is doing is provided. To get information on the usage of individual scripts, please, run <script.py --help>

- `adjust_id.py` This script assigns numerical IDs to sequences in a multi-fasta file. It is called from `run_hmmscan_pvog.sh`, to ensure unique sequence IDs when functionally annotating predicted protein sequences.

- `check_circularity.py` This script checks if a contig is circular based on overlap between the ends. By defaults, a contig is considered ciruclar if ends overlap with at least 25 bp and at most 200 bp.

- `check_tandem_stop_codon_usage.py` This scripts writes out the 30 bp following a stop codon of dual-coded genes for the different genetic code models. The sequences can the be uploaded to Weblogo to create a sequence logo.

- `evaluate_runt_time.py` This scripts evaluates the run time of Mgcod in the two modes on contigs of various lengths.

- `evaluate_switching_point_prediction.py` This script evaluates the accuracy with which Mgcod predicts the switch of genetic code in simulated genomes with multiple genetic codes. It can be used in the following way:
```
./evaluate_switching_point_prediction.py -a <file_with_annotated_switch_point_coordinates> -r <directory_with_mgcod_results> -o <directory_where_to_save_figure>
```

- `get_positions_sup_trna_circular_genomes.py` This script visualizes the positions of predicted Sup-tRNAs relative to switch in genetic code. It was used to generate Figure S4 in the manuscript using the following command:
```
./get_positions_sup_trna_circular_genomes.py --path_mgcod_results <direcotry_with_mgcod_results --path_trnascan_results <directory_with_trnascan_predictions> --path_to_genomes  <directory_with_input_sequences> --target_genomes <file_listing_target_sequences> --alternative_gcode <genetic_code_id> --output_dir <directory_where_to_save_output>
```

- `parse_fastani.py` This script parses fastANI results and summarizes average ANI values.

- `parse_hmmscan_results.py` This script parses hmmscan output with functional annotation of predicted protein sequences and provides a summary of the findings.

- `run_hmmscan_pvog.sh` This script functionally annotates predicted protein sequences using hmmscan.

- `run_hmmsearch_pfam.sh` This scripts functionnally annotated predicted protein sequences using hmmsearch.

- `run_trnascan.sh` This script predicts tRNA genes in target sequences using tRNAscan-SE.

- `run_virfinder.R` This script predicts whether a target sequences is of viral origin or not using VirFinder.

- `simulate_gcodes.py` This script simulates contigs with stop codon reassignments.

- `simulate_switching_point.py` This script simulates contigs with multiple genetic codes.

- `split_genome.py` This script splits input genomes into equal sized contigs.

- `split_protein_seq_by_genetic_code.py` This script splits sequences of annotated coding sequences by genetic code.

- `summarize_predictions_metagenomes.py` This script summarizes Mgcod predictions on metagenomic contigs in terms of predicted genetic codes. It can be used in the following way:
```
./summarize_predictions_metagenomes.py -r <file_with_paths_to_mgcod_predictions> -g <file_with_paths_to_target_genomes> -o <output_file> -p <directory_where_to_save_figure>
```

- `validate_reassignment_msa.py` This script validates predicted stop codon reassignment by performing MSA of protein sequences with the same genetic code. It was used to generate Figure S3 using the following command:
```
./validate_reassignment_msa.py --input_genomes <text_files_with_target_genomes> --mgcod_results <directories_with_corresponding_mgcod_predictions> --output_dir <directories_where_to_save_msa_etc> --reassigned_stop_codon> <corresponding_reassigned_stop_codons>
```
