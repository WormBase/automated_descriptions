#!/bin/bash 


echo "total 9 scripts.";
./create_sentence_tissue_expressions_elegans.pl
echo " finishes 1";
./create_GO_sentences_elegans_species_parallel_all.pl 
echo "finishes 2"; 
./create_sentence_multiple_orthologs_species_all_parallel_all.pl
echo "finishes 3"; 
./create_sentence_multiple_orthologs_species_all_parallel_all.pl
echo "finishes 4"; 
./create_gene_regulation_expression_cluster_sentences_species_parallel_all.pl 
echo "finishes 5"; 
./create_molecule_regulation_expression_cluster_sentences_species_parallel_all.pl 
echo "finishes 6"; 
./create_anatomy_expression_cluster_sentences_species_parallel_all.pl 
echo "finishes 7"; 
./concatenate_sentences_species_parallel_all.pl
echo "finishes 8"; 
./generate_OA_concise_descriptions_parallel_all.pl
echo "finishes 9"; 
