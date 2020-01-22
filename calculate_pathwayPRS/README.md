# pathwayPRS
`pathwayPRS` is an executable R script that calculated pathway-specific polygenic risk scores (PRS) for the 5 pre-selected pathway associated with Alzheimer's disease (AD).

# How the tool works
`pathwayPRS` only needs two input files:
  1. the Supplementary Table file from the paper `Tesi et al., 2020`
  2. plink2 files
The supplementary table is used to take the variant-pathway annotations. Briefly, this is a weight for each of the 29 AD-associated genetic variants that quantify the involvement of each variant in the 5 pre-selected pathways.
Standard PLINK2 files (files .pvar/.psam/.pgen) are required for an internal command to extract genotype dosages.

# Where to pay attention
You should be aware of the variant identifier of the PLINK2 files. For this, you can easily check the `.var` file of a random chromosome. Make sure that the column `ID` is coded either using the variant `rsid` or `chr:pos`. The script will try to run in both cases to prevent any missing variant.
Having taken care of this, the program should run smoothly. By default PRS with and withou APOE are performed.
The program should output 4 files:
  1. pathway-specific PRS including APOE variants for each of the 5 pre-selected pathways for all the individuals specified in the PLINK2 files
  2. pathway-specific PRS excluding APOE variants for each of the 5 pre-selected pathways for all the individuals specified in the PLINK2 files
  3. total PRS including APOE variants for all the individuals specified in the PLINK2 files
  4. total PRS excluding APOE variants for all the individuals specified in the PLINK2 files

# Contact us
In case you run in any issue, please do not hesitate to contact us at `n.tesi@amsterdamumc.nl`
