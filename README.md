# pathway-PRS
Pipeline to perform annotation of variants to a pre-selected set of pathways and then calculate and associate pathway-specific Polygenic Risk Scores (PRS) and typical PRS.
As such, the pipeline is divided into 2 part. The first part is the calculation of the variant-pathway weights, that is, a score ranging [0-1] describing the assignment of SNP `x` to pathway `y`. For this study, we assumed 5 pathways, therefore each SNP has a weight of association to each of the 5 pathways.

## variant-pathway mapping
The variant-pathway mapping weights can be accessed throught the Supplementary material of the published manuscript available at https://www.nature.com/articles/s41398-020-01018-7#Sec18. The Supplementary Table to be used as input for the PRS calculation is Supplementary Table S10. This file is also available in `example_files_pathwayPRS` folder. The variant-pathway annotation is a 2-step procedure:
- first, SNPs are linked to their likely affected genes
- second, genes are linked to their likely affected pathways
Finally, by looking at which genes each SNP is associated with, and which pathways each gene is associated with, it is possible to derive a weight for each SNP to each of the considered pathways. For further information how this is done, please look at Methods section in our manuscript (https://www.nature.com/articles/s41398-020-01018-7#Sec2), contact us (n.tesi@amsterdamumc.nl) or check the scripts in `variant_pathway_annotation` folder.

The calculation of variant-pathway annotation was made using ad-hoc files from previous GWASes of AD, and thus could not be replicable for a different set of SNPs and/or a different trait. If you wish to replicate out variant-pathway mapping, then the scripts and the necessary input files are located in `variant_pathway_annotation` and `example_files_variant_pathway_annotation` folders. The main file is the `VPA_pipeline.sh`, which lists all commands that should be run. This pipeline produces global and per-gene proportions of assignment of each gene and SNP to each of the considered molecular pathways.

## pathway-PRS calculation
After the variant-pathway mapping is done, it is possible to calculate PRS and pathway-PRS. Please make sure to have `PLINK2` correctly installed and accessible anywhere in your system (further info at https://www.cog-genomics.org/plink/2.0/). To execute, the script in `calculate_pathwayPRS` folder, can be used. Using the example files contained in `example_files_pathwayPRS` folder, the script can be run with:

`Rscript pathwayPRS.R ../example_files_pathwayPRS/pathwayPRS_TranslationalPsychiatry.xlsx ../example_files_pathwayPRS/plink_examples/chr1_example`

where the first argument (`../example_files_pathwayPRS/pathwayPRS_TranslationalPsychiatry.xlsx`) is the supplementary table containing the variant-pathway weights, and the second argument (`../example_files_pathwayPRS/plink_examples/chr1_example`) is the path to the chromosome 1 data (in PLINK2 format) without extension. 
In principle, you can use your own file containing the variant-pathway weights, but please keep the same format as our Supplementary table. Alternatively, you may need to run line-by-line the script. Note that the script will extract SNP dosages from PLINK2 files assuming that the SNP identifier as used in the PLINK2 files is either the SNP id (rsXXXX) or the locus id (chromosome:position). If your PLINK2 file has different IDs, this script will not work. If you still want to calculate pathway-PRS, we are happy to help trouble-shooting.

Briefly, this script will:
1. read input files
2. extract SNP dosages of the SNPs of interest from PLINK2 files
3. read dosages and merge all chromosome-specific files
4. calculate pathway-specific PRS using the variant-pathway weights
5. calculate typical PRS including all variants
6. write the PRSs as output

As such, the output of the script will be 4 PRSs:
- pathway-specific PRS including APOE variants (`PRS_perPath_apoeInc.txt`)
- pathway-specific PRS excluding APOE variants (`PRS_perPath_apoeExc.txt`)
- typical PRS including APOE variants (`PRS_AllVar_apoeInc.txt`)
- typical PRS excluding APOE variants (`PRS_AllVar_apoeExc.txt`)