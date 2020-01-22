##########################################
##########################################
# SCRIPT ON HOW TO REPRODUCE PATHWAY-PRS #
# GIVEN THE VARIANT-PATHWAY MAPPING AND  #
# A SET OF DOSAGES                       #
##########################################
##########################################

##########################################
# FILES NECESSARY TO RUN THE SCRIPT ARE: #
# 1. VARIANT-PATHWAY ANNOTATIONS         #
# 2. PLINK FILES TO EXTRACT DOSAGES      #
##########################################

##########################################
# LIBRARIES NEEDED TO RUN THE SCRIPT     #
# CHECK IF THE PACKAGES ARE INSTALLED    #
# AND IN CASE THEY ARE NOT, IT INSTALLS  #
##########################################
list.of.packages <- c("data.table", "stringr", "xlsx", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(data.table)
library(stringr)
library(xlsx)
args = commandArgs(trailingOnly=TRUE)
##########################################

##########################################
# ARGUMENTS DESCRIPTION:                 #
# 1. ARGUMENT 1 IS THE SUPPLEMENTARY     #
# TABLE THAT INCLUDES EFFECT-SIZES AND   #
# THE VARIANT MAPPINGS                   #
# 2. ARGUMENT 2 IS THE PLINK2 DOSAGES    #
# PATH TO CHROMOSOME 1 FILE              #
##########################################

#################################################
# FUNCTIONS                                     #
#################################################

# FUNCTION TO READ ARGUMENTS AND PROVIDE AND ERROR IF THEY ARE NOT ALL THERE
function.checkArguments <- function(){
  RUN <- tryCatch(
    {
      #arg. 1 is the literature file
      literature.file <- args[1]
      plink2.files <- args[2]

      RUN = TRUE
    },
    error=function(cond){
      RUN = FALSE
      print("!!! Sorry, we could not execute the program because there")
      print("!!! were issues while reading the arguments. Please fix it")
      print("!!! and run the program again. If you can not solve, please")
      print("!!! contact us at n.tesi@amsterdamumc.nl.")
    }
  )
  return(RUN)
}
#################################################

# FUNCTION TO ORDER TWO DATA FRAME IN THE SAME WAY USING THE MATCH R COMMAND
function.match <- function(df1, by.1, by.2){
  df1 <- df1[match(by.1, by.2),]
  return(df1)
}
#################################################

# FUNCTION TO CHECK AND MATCH EFFECT ALLELES IN LITERATURE AND MY DATASET
function.checkAlleles <- function(alleles, lit){
  #define locus column for alleles group
  alleles$locus <- paste(alleles$chr, alleles$pos, sep=':')

  #make sure both datasets have the same variants
  lit.flt <- subset(lit, lit$LOCUS %in% alleles$locus)

  #order same way
  alleles <- function.match(alleles, lit.flt$LOCUS, alleles$locus)

  #check alleles
  lit.flt$adj.beta <- lit.flt$BETA
  lit.flt$adj.beta[which(lit.flt$A1 != alleles$alt)] <- lit.flt$adj.beta[which(lit.flt$A1 != alleles$alt)] * (-1)

  return(lit.flt)
}
#################################################

# FUNCTION TO COMPUTE GENETIC RISK SCORE PER PATHWAY -- PRS = SUM(DOSAGE*BETA*VARIANT_PATHWAY_ANNOT)
function.grs.beta.apoeInc <- function(dosages, loci.ad.effect, xxx){
  #define output
  grs.beta.apoeInc <- matrix(data = 0, nrow = ncol(dosages), ncol = 5)
  rownames(grs.beta.apoeInc) <- colnames(dosages)
  colnames(grs.beta.apoeInc) <- colnames(xxx)[2:6]

  #main loop over pathways -- before it was looping on 3:ncol(xxx)
  for (i in 2:6){
    #select pathway of interest
    path <- xxx[, c("LOCUS", colnames(xxx)[i])]

    #user update
    print(paste("Working on ", colnames(xxx)[i], sep=''))

    #remove variants with no association to that pathway
    path <- subset(path, path[, 2] != 0)

    #isolate pathway name
    path.name <- colnames(path)[2]

    #loop over variants
    for (j in 1:nrow(path)){
      #identify snp
      snp <- as.character(path$LOCUS[j])

      #extract dosages for that snp
      code <- as.character(loci.ad.effect$V1[which(loci.ad.effect$LOCUS == snp)])
      dos <- subset(dosages, rownames(dosages) == code)

      #find proportion of pathway
      prop <- path[which(path$LOCUS == snp), path.name]

      #find beta
      beta <- loci.ad.effect[which(loci.ad.effect$LOCUS == snp), "adj.beta"]

      #compute grs here
      scores <- t(as.numeric(dos) * prop * beta)
      grs.beta.apoeInc[, path.name] <- grs.beta.apoeInc[, path.name] + scores

    }
  }
  return(grs.beta.apoeInc)
}

#function to compute PRS across all genes independently from the pathways -- normal PRS, apoe included
function.grs.AllVar.apoeIncl <- function(dosages, loci.ad.effect){
  #define output
  allVars <- matrix(data = 0, nrow = ncol(dosages), ncol = 1)
  colnames(allVars) <- "PRS"
  rownames(allVars) <- colnames(dosages)

  #order dosages and literature in the same way
  if (nrow(dosages) > 1){loci.ad.effect <- function.match(loci.ad.effect, rownames(dosages), loci.ad.effect$V1)}

  #loop over samples
  for (sam in 1:ncol(dosages)){
    #identify snps
    snps <- dosages[, sam]
    betas <- loci.ad.effect$adj.beta

    #scores
    scores <- as.numeric(snps) * betas

    #compute grs here
    grs <- sum(scores)

    #add to df
    allVars[sam, ] <- grs

  }
  return(allVars)
}

# FUNCTION TO EXTRACT DOSAGES FROM PLINK FILES GIVEN THE PATH OF CHR1 FILE ASSUMING VARIANT ID IS RSID
function.ExtractDosages <- function(literature, plink.files){
	#write temporary file with snp list
	write.table(literature$SNP, "tmp.variants", row.names=F, col.names=F, quote=F)
        write.table(literature$LOCUS, "tmp.variants2", row.names=F, col.names=F, quote=F)

	#modify plink2 file link to make it suitable for loop
	tmp.link <- str_split_fixed(plink.files, "1", 2)
	corr.link <- paste(tmp.link[, 1], "${chr}", tmp.link[, 2], sep="")

	#get the list of chromosomes to loop on
	chroms <- as.character(unique(literature$CHR))
	chrom_list <- paste0(chroms, collapse=" ")

	#make command for dosages extraction
	plink2_cmd <- paste("for chr in ", chrom_list, "; do plink2 --pfile ", corr.link, " --extract tmp.variants --export A --out chr${chr}_tmp; done", sep="")
        plink2_cmd2 <- paste("for chr in ", chrom_list, "; do plink2 --pfile ", corr.link, " --extract tmp.variants2 --export A --out chr${chr}_tmp2; done", sep="")
	system(plink2_cmd, ignore.stdout = TRUE)
        system(plink2_cmd2, ignore.stdout = TRUE)
	#clean log files
	cmd_clean <- "rm *log"
	system(cmd_clean)
	#merge files together
	cmd_bind <- "paste --delimiters='\t' *raw > chrAll_dosages.txt"
	system(cmd_bind)
	cmd_clean <- "rm *raw"
	system(cmd_clean)
	cmd_clean <- "rm *tmp*"
	system(cmd_clean)

	#read file back in R
	dos <- read.table("chrAll_dosages.txt", h=T, sep="\t", check.names=F, stringsAsFactors=F)

	#take informative columns
	x <- colnames(dos)
	col.index <- grep(":", x)
	col.index2 <- grep("rs", x)
	dos.cl <- dos[, c(col.index, col.index2)]
	dos.cl$IID <- dos$IID

	return(dos.cl)
}

#################################################
# MAIN                                          #
#################################################
set.seed(1234)

#READ ARGUMENTS AND DECIDE IF TO RUN OR NOT
RUN = function.checkArguments()

if (RUN == TRUE){
  # READ LITERATURE AND VARIANT-PATHWAY MAPPING -- THE FILE IS FROM SUPP.TABLE S11
  print("## Reading and parsing literature file..")
  literature <- read.xlsx(file = args[1], sheetIndex = 18,
                          startRow = 2, stringsAsFactors = F)
  literature$LOCUS <- paste(literature$CHR, literature$POS, sep=":")
  literature <- literature[!is.na(literature$CHOLESTEROL.LIPID),]

  # EXTRACT DOSAGES AND READ THEM BACK
  print("## Extracting dosages from PLINK2 files..")
  dosages <- function.ExtractDosages(literature, args[2])
  loci <- as.data.frame(str_split_fixed(colnames(dosages), "_", 3))
  colnames(dosages) <- loci$V1
  samples <- dosages$IID
  dosages <- t(dosages)
  colnames(dosages) <- samples
  dosages <- dosages[which(rownames(dosages) != "IID"),]

  # CHECK ALLELES
  print("## Checking alleles with literature..")
  loci.ad.effect.m1 <- merge(loci, literature, by.x='V1', by.y='LOCUS')
  loci.ad.effect.m2 <- merge(loci, literature, by.x='V1', by.y='SNP')
  loci.red1 <- loci.ad.effect.m1[, c(1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22)]
  loci.red2 <- loci.ad.effect.m2[, c(1, 2, 4, 5, 1, 6, 7, 8, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21)]
  colnames(loci.red2) <- colnames(loci.red1)
  loci.ad.effect <- rbind(loci.red1, loci.red2)
  loci.ad.effect$adj.beta <- loci.ad.effect$BETA
  loci.ad.effect$A1 <- as.character(loci.ad.effect$A1)
  loci.ad.effect$V2 <- as.character(loci.ad.effect$V2)
  loci.ad.effect$adj.beta[which(loci.ad.effect$V2 != loci.ad.effect$A1)] <- loci.ad.effect$adj.beta[which(loci.ad.effect$V2 != loci.ad.effect$A1)] * (-1)
  loci.ad.effect$LOCUS <- paste(loci.ad.effect$CHR, loci.ad.effect$POSITION, sep=":")

  # EXTRACT VARIANT-PATHWAY ANNOTATION FOR EACH VARIANT
  xxx <- literature[, c(20, 14, 15, 16, 17, 18, 19)]

  print("## Calculating PRSs..")
  #grs multiplying beta, proportion of pathways and dosages -- apoe included
  grs.beta.apoeInc <- as.data.frame(function.grs.beta.apoeInc(dosages, loci.ad.effect, xxx))

  #grs multiplying beta, proportion of pathways and dosages -- apoe excluded
  dosages.noAPOE <- dosages[which(!(rownames(dosages) %in% c("rs429358", "rs7412", "19:45412079", "19:45411941"))),]
  loci.ad.effect.noAPOE <- loci.ad.effect[which(!(loci.ad.effect$SNP %in% c("rs429358", "rs7412", "19:45412079", "19:45411941"))),]
  xxx.noAPOE <- xxx[which(!(xxx$LOCUS %in% c("19:45412079", "19:45411941", "rs429358", "rs7412"))),]
  grs.beta.apoeExc <- as.data.frame(function.grs.beta.apoeInc(dosages.noAPOE, loci.ad.effect.noAPOE, xxx.noAPOE))

  #grs including all 33 variants weighted by beta independently from the pathways -- apoe included
  AllVar.apoeIncl <- as.data.frame(function.grs.AllVar.apoeIncl(dosages, loci.ad.effect))

  #grs including all 33 variants weighted by beta independently from the pathways -- apoe excluded
  AllVar.apoeExcl <- as.data.frame(function.grs.AllVar.apoeIncl(dosages.noAPOE, loci.ad.effect.noAPOE))

  print("## Scaling PRSs..")
  AllVar.apoeIncl.scaled <- AllVar.apoeIncl
  AllVar.apoeIncl.scaled$PRS <- scale(x = AllVar.apoeIncl$PRS)
  AllVar.apoeExcl.scaled <- AllVar.apoeExcl
  AllVar.apoeExcl.scaled$PRS <- scale(x = AllVar.apoeExcl$PRS)
  cols <- colnames(grs.beta.apoeInc)
  grs.beta.apoeInc.scaled <- grs.beta.apoeInc
  grs.beta.apoeExc.scaled <- grs.beta.apoeExc
  for (i in cols){
    grs.beta.apoeInc.scaled[, i] <- scale(x = grs.beta.apoeInc[, i])
    grs.beta.apoeExc.scaled[, i] <- scale(x = grs.beta.apoeExc[, i])
  }

  print("## Writing outputs..")
  write.table(grs.beta.apoeExc.scaled, 'PRS_perPath_apoeExc.txt', quote=F, row.names = T, sep='\t')
  write.table(grs.beta.apoeInc.scaled, 'PRS_perPath_apoeInc.txt', quote=F, row.names = T, sep='\t')
  write.table(AllVar.apoeExcl.scaled, 'PRS_AllVar_apoeExc.txt', quote=F, row.names = T, sep='\t')
  write.table(AllVar.apoeIncl.scaled, 'PRS_AllVar_apoeInc.txt', quote=F, row.names = T, sep='\t')

  print("## Done")
}
