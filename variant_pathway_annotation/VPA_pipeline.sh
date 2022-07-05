#Pipeline for functional annotation: the aim of the pipeline is to read all functional annotation input files, store   #
#them and then integrate the information they contain                                                                  #

#compare IGAP and IRIS mapping procedures, outputting a file with percentages for each locus as mapping to a given gene
python3 2_compare_iris_igap.py ../example_files_variant_pathway_annotation/1_igap_iris_mapping_genes.txt
#the output is a file called "3_merged_annotation_IGAP_IRIS_OTHER.txt" and a file called "4_merged_annotation_IGAP_IRIS_OTHER_gene_list.txt"

#We decided to be strict to these functional groups = Immune system, Amyloid metabolism, Cholesterol transport, Endocytosis and Angiogenesis

##Retrieve info about most significant GO from IGAP, this file is stored in the current folder
python3 6_parse_IGAP_annotation.py ../example_files_variant_pathway_annotation/5_IGAP_pathways_with_genes.txt
#output is set of my genes with annotation in IGAP in one of the pre-defined GO groups

##Retrieve info about GO function from literature, i.e Sleegers paper
python3 9_parse_LITERATURE_annotation.py ../example_files_variant_pathway_annotation/8_Literature_pathways_with_genes.txt 3_merged_annotation_IGAP_IRIS_OTHER.txt

##Retrieve info about GO function from Gene Ontology website
#I downloaded the full list of terms descendants of: immune system process, endocytosis, cholesterol metabolic process, cholesterol transport, beta-amyloid metabolic process, protein-lipid complex subunit organization
#I downloaded the full list of all genes associated with terms and descendants of: immune system process, endocytosis, cholesterol metabolic process, cholesterol transport, beta-amyloid metabolic process, protein-lipid complex subunit organization
#GO is a sql database, the script for the download is in 11_gene_ontology.sql
#all the downloaded files are stored in ../example_files_variant_pathway_annotation/gene_ontology/ folder
#for go in ../example_files_variant_pathway_annotation/gene_ontology/*all_genes*; do python 12_parse_GENE-ONTOLOGY_annotation.py $go; done
#printf "function\tgenes\n" > 14_assigned_functions_GENE-ONTOLOGY.txt | cat 13_assigned_function_GENE-ONTOLOGY_* >> 14_assigned_functions_GENE-ONTOLOGY.txt
#alternatively, go to amigo2 which is updated and download gene-list of all genes linked to a given go term number manually, then do the same:
for go in ../example_files_variant_pathway_annotation/gene_ontology/*all_genes*; do python3 12_parse_GENE-ONTOLOGY_annotation.py $go; done
printf "function\tgenes\n" > 14_assigned_functions_GENE-ONTOLOGY.txt | cat 13_assigned_function_GENE-ONTOLOGY_* >> 14_assigned_functions_GENE-ONTOLOGY.txt


##Retrieve info about GO function from DAVID annotation clustering
python3 16_parse_DAVID_annotation.py ../example_files_variant_pathway_annotation/15_david_annotation.txt
#input file was generated with DAVID, GO_BP + UP_KEYWORDS + UP_SEQ_FEATURE, high classification stringency, 15 clusters NB: the script can be used with other DAVID clusterings

##Make some adjustments to the files
sed 's/cholesterol/Cholesterol\/lipid/g' 14_assigned_functions_GENE-ONTOLOGY.txt | sed 's/endocytosis/Endocytosis/g' | sed 's/beta-amyloid/Beta-amyloid/g' | sed 's/immune/Immune system/g' | sed 's/angiogenesis/Angiogenesis/g' > 14_assigned_functions_GENE-ONTOLOGY.mod.txt
rm 14_assigned_functions_GENE-ONTOLOGY.txt
mv 14_assigned_functions_GENE-ONTOLOGY.mod.txt 14_assigned_functions_GENE-ONTOLOGY.txt
sed 's/Vascular\/Angiogenesis/Angiogenesis/g' 17_assigned_functions_DAVID.txt | sed 's/immune system/Immune system/g' > 17_assigned_functions_DAVID.mod.txt
rm 17_assigned_functions_DAVID.txt
mv 17_assigned_functions_DAVID.mod.txt 17_assigned_functions_DAVID.txt

##Make some stats on the functions assigned
python3 18_functional_annotation_stats.py
#no input files, no output files

##Merge information from multiple annotations into one single annotation
python3 19_merge_annotation_files.py
#no input file required, generate output file that is "19_complete_annotation_sorted.txt"

##Go to R script for going from this output to proportions of functions per locus
Rscript 21_makeProportions.R
#no input required, output is "21_final_proportions.txt"

##!!!two variants have different locus id --> change it here
sed 's/14:53391680/14:53400629/g' 22_final_proportions.txt | sed 's/11:47380340/11:47557871/g' > 23_final_proportions_adj.txt
