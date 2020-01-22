#compute proportions of functional annotation

library(RColorBrewer)
library(ggplot2)
rm(list=ls())

#read input file
inp <- read.table("/Users/nicco/Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/20_complete_annotation_sorted.txt", h=T, sep='\t')

### first part -- compute proportions, then exclude loci with unknown class > 90%
print("Cleaning input file...")

#define new dataset with proportions
prop <- inp

#make proportions of each functional group, unknnown included
for (i in 1:nrow(inp)){
  for (j in 5:10){
    prop[i, j] <- inp[i, j] / sum(inp[i, 5:10])
  }
}

#define new dataset with weighted proportions
weighted.prop <- prop

#make weighted proportions
weighted.prop <- prop
for (i in 1:nrow(prop)){
  for (j in 5:10){
    weighted.prop[i, j] <- prop$weight[i] * prop[i, j]
  }
}

#define max.locus as the maximum locus
max.locus <- max(weighted.prop$locus)

#define new dataset for global per-locus functional proportions
global.prop <- as.data.frame(matrix(data = NA, nrow = max.locus, ncol = ncol(weighted.prop)-2))
colnames(global.prop) <- c("locus", "pos", "immune.system", "beta.amyloid", "endocytosis", "cholesterol.lipid", "vascular", "unknown")

#loop to fill in global proportions of functions per locus
for (loc in 1:max.locus){
  subs <- weighted.prop[which(weighted.prop$locus == loc),]
  global.prop[loc, "locus"] <- loc
  global.prop[loc, "pos"] <- unique(as.vector(subs$pos))
  global.prop[loc, "immune.system"] <- sum(subs$immune.system)
  global.prop[loc, "beta.amyloid"] <- sum(subs$beta.amyloid)
  global.prop[loc, "endocytosis"] <- sum(subs$endocytosis)
  global.prop[loc, "cholesterol.lipid"] <- sum(subs$cholesterol.lipid)
  global.prop[loc, "vascular"] <- sum(subs$angiogenesis)
  global.prop[loc, "unknown"] <- sum(subs$unknown)
}

#exclude loci with >90% unknown annotations
#global.prop.flt <- subset(global.prop, global.prop$unknown < 1)
global.prop.flt <- global.prop
nrow(global.prop)
nrow(global.prop.flt)

print(paste("Before QC, there were ", nrow(global.prop), " loci."))
print(paste("After removing loci with >90% unknown function, there are ", nrow(global.prop.flt), " loci."))

### second part -- re-compute proportion for good loci
print("Now calculating proportions of function per locus")

#re-compute proportions of loci without too many unknown annotations, excluding the unknown class
inp.flt <- inp[which(inp$pos %in% global.prop.flt$pos),]

#remove unknown class and all lines (genes) with unknown annotation
inp.flt <- subset(inp.flt, inp.flt$unknown == 0)

#since I removed annotations, i need to re-compute weights otherwise it won't sum up to 1
#define new dataset
inp.flt.new <- inp.flt
#define loci
loci <- unique(inp.flt$locus)
for (i in 1:length(loci)){
  loc <- loci[i]
  subs <- inp.flt[which(inp.flt$locus == loc),]
  if (nrow(subs) == 1){
    subs$weight <- 1.00
    inp.flt.new[which(inp.flt.new$locus == loc),] <- subs
  }
  else {
    new.weights <- vector(length = nrow(subs))
    for (j in 1:nrow(subs)){
      old.w <- subs[j, "weight"]
      new.w <- old.w / sum(subs[, "weight"])
      new.weights[j] <- new.w
      }
    subs$weight <- new.weights
    inp.flt.new[which(inp.flt.new$locus == loc),] <- subs
  }
}

#define new dataset with proportions
prop.flt <- inp.flt.new
prop.flt$unknown <- NULL

#make proportions of each functional group, unknnown included
for (i in 1:nrow(inp.flt)){
  for (j in 5:9){
    if (sum(inp.flt[i, 5:9]) == 0){
      prop.flt[i, j] <- 0
    }
    else {
      prop.flt[i, j] <- inp.flt[i, j] / sum(inp.flt[i, 5:9])
    }
  }
}

#define new dataset with weighted proportions
weighted.prop.flt <- prop.flt

#make weighted proportions
for (i in 1:nrow(prop.flt)){
  for (j in 5:9){
    weighted.prop.flt[i, j] <- prop.flt$weight[i] * prop.flt[i, j]
  }
}


#define new dataset for global per-locus functional proportions
max.locus.flt <- unique(prop.flt$locus)
global.prop.flt <- as.data.frame(matrix(data = NA, nrow = length(max.locus.flt), ncol = ncol(weighted.prop.flt)-2))
colnames(global.prop.flt) <- c("locus", "pos", "immune.system", "beta.amyloid", "endocytosis", "cholesterol.lipid", "vascular")

#loop to fill in global proportions of functions per locus
for (i in 1:length(max.locus.flt)){
  loc <- max.locus.flt[i]
  subs <- weighted.prop.flt[which(weighted.prop.flt$locus == loc),]
  global.prop.flt[i, "locus"] <- loc
  global.prop.flt[i, "pos"] <- unique(as.vector(subs$pos))
  global.prop.flt[i, "immune.system"] <- sum(subs$immune.system)
  global.prop.flt[i, "beta.amyloid"] <- sum(subs$beta.amyloid)
  global.prop.flt[i, "endocytosis"] <- sum(subs$endocytosis)
  global.prop.flt[i, "cholesterol.lipid"] <- sum(subs$cholesterol.lipid)
  global.prop.flt[i, "vascular"] <- sum(subs$angiogenesis)
}
global.prop.flt # <-- this is the file without unknown with all loci (it is redundant!!!)
weighted.prop.flt # <-- this is the file without unknown with all genes

#add column for chromosome and position to both global.prop.flt and weighted.prop.flt
library(stringr)
tmp.mt <- str_split_fixed(global.prop.flt$pos, ":", 2)
global.prop.flt$chr <- as.numeric(tmp.mt[, 1])
global.prop.flt$pos <- as.numeric(tmp.mt[, 2])
global.prop.flt$locus_id <- paste(global.prop.flt$chr, global.prop.flt$pos, sep=":")
tmp.mt <- str_split_fixed(weighted.prop.flt$pos, ":", 2)
weighted.prop.flt$chr <- as.numeric(tmp.mt[, 1])
weighted.prop.flt$pos <- as.numeric(tmp.mt[, 2])
weighted.prop.flt$locus_id <- paste(weighted.prop.flt$chr, weighted.prop.flt$pos, sep=":")

#sort by chromosome and position
global.prop.flt <- global.prop.flt[order(global.prop.flt$chr, global.prop.flt$pos),]
weighted.prop.flt <- weighted.prop.flt[order(weighted.prop.flt$chr, weighted.prop.flt$pos),]

#need to add the variants with missing annotation, for this need to read a file with the total list of variants
allVars <- unique(as.character(inp$pos))
#exclude those we are not interested into
toExc <- c("16:31133100", "19:46241841", "6:41129207", "18:56189459", "15:63569902", "17:47297297", "19:1056244", "3:57226150", "6:41129252", "7:145950029")
allVars <- allVars[which(!(allVars %in% toExc))]

#now check who is missing and in case add it back with unknonwn annotation
global.props <- global.prop.flt
global.props <- global.props[!duplicated(global.props$locus_id),]
global.props <- global.props[which((global.props$locus_id %in% allVars)),]
missings <- allVars[which(!(allVars %in% global.props$locus_id))]
missings.genes <- c("NME8", "KANSL1", "HS3ST1", "CASS4", "ECHDC3", "ADAMST1")
red <- global.props[, c("immune.system", "beta.amyloid", "endocytosis", "cholesterol.lipid", "vascular", "locus_id")]
for (x in 1:length(missings)){
  tmp <- data.frame(immune.system = 0, beta.amyloid = 0, endocytosis = 0, cholesterol.lipid = 0,
                    vascular = 0, locus_id = missings[x])
  red <- rbind(red, tmp)
}
red <- red[!duplicated(red$locus_id),]
dim(red)
tmp.mt <- str_split_fixed(red$locus_id, ":", 2)
red$chr <- as.numeric(tmp.mt[, 1])
red$pos <- as.numeric(tmp.mt[, 2])
red <- red[order(red$chr, red$pos),]

#now need to do the same for the variant-gene mapping
gene.level <- weighted.prop.flt
gene.level <- gene.level[which((gene.level$locus_id %in% allVars)),]
missings <- allVars[which(!(allVars %in% gene.level$locus_id))]
missings.genes <- c("NME8", "KANSL1", "HS3ST1", "CASS4", "ECHDC3", "ADAMST1")
red.genes <- gene.level[, c("gene", "weight", "immune.system", "beta.amyloid", "endocytosis", "cholesterol.lipid", "angiogenesis", "locus_id")]
for (x in 1:length(missings)){
  tmp <- data.frame(gene = missings.genes[x], weight = 1, immune.system = 0, beta.amyloid = 0, endocytosis = 0, cholesterol.lipid = 0,
                    angiogenesis = 0, locus_id = missings[x])
  red.genes <- rbind(red.genes, tmp)
}
tmp.mt <- str_split_fixed(red.genes$locus_id, ":", 2)
red.genes$chr <- as.numeric(tmp.mt[, 1])
red.genes$pos <- as.numeric(tmp.mt[, 2])
red.genes <- red.genes[order(red.genes$chr, red.genes$pos),]

#calculate number of genes of variant -- this is done from the inp object
inp2 <- inp[which(as.character(inp$pos) %in% allVars),]
inp2$pos <- as.character(inp2$pos)
gperv <- as.data.frame(matrix(data = NA, nrow = length(unique(inp2$pos)), ncol = 2))
colnames(gperv) <- c("locus", "ngenes")
c <- 1
for (i in unique(inp2$locus)){
  #get locus
  ng <- nrow(inp2[which(inp2$locus == i),])
  
  #assign
  gperv[c, ] <- c(unique(inp2$pos[which(inp2$locus == i)]), ng)
  c <- c + 1
}

#finally the plot
dat <- red.genes
global.props <- red
function.heatmapGeneLevel.var2 <- function(dat, global.props, gperv){
  #read file with snpid and position to map
  rsid <- read.table("Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/rsid_snps.txt", h=T, 
                     dec=",", sep="\t")
  rsid$LOCUS <- paste(rsid$CHR, rsid$POSITION, sep=":")
  
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(0, 11, 17, 6))
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #define special counter for labels on the left
  tmp <- 0.25
  
  #define general counter 
  counter <- 0
  
  #general counter for heatmat on the right height
  p <- 1
  
  #set difference between labels and heatmap
  interv <- 1.5
  
  #define max locus number
  max.locus <- length(unique(dat$locus_id))
  
  #define ymax
  ymax <- nrow(dat)*h
  
  #set difference between central part and other parts -- this is for the right part
  diff <- (nrow(dat) - nrow(global.props)) / 10 * 2
  
  #set different between central part and left part -- different from above
  d2 <- (nrow(dat) - nrow(global.props)) / 10 * 2
  
  #prepare window to plot
  plot(1, xlim = c(0, 27),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = -5, y = ymax - d2 + (interv*2), labels = "Locus", adj=0.5, xpd=T, cex=1.50, font=2)
  text(x = 0, y = ymax - d2 + (interv*2), labels = "N. genes", adj=0.5, xpd=T, cex=1.50, font=2)
  text(x = 8.5, y = ymax + interv + 1, labels = "Genes", adj=1, xpd=T, cex=1.50, font=2, srt=-90)
  text(x = 9.5, y = ymax + interv + 1, labels = "Weight", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 12.5, y = ymax + interv, labels = "Immune response", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 13.5, y = ymax + interv, labels = "Beta-amyloid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 14.5, y = ymax + interv, labels = "Endocytosis", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 15.5, y = ymax + interv, labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 16.5, y = ymax + interv, labels = "Vascular dys.", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 21.5, y = ymax - diff + interv*2, labels = "Immune response", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 22.5, y = ymax - diff + interv*2, labels = "Beta-amyloid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 23.5, y = ymax - diff + interv*2, labels = "Endocytosis", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 24.5, y = ymax - diff + interv*2, labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 25.5, y = ymax - diff + interv*2, labels = "Angiogenesis", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  segments(x0 = 7.5, y0 = ymax + interv + 14.5, x1 = 10.5, y1 = ymax + interv + 14.5, lwd=4, xpd=T)
  text(x = 9, adj = 0.5, y = ymax + interv + 17.5, labels = "Variant-gene\nmapping", cex=1.50, font=4, xpd=T)
  segments(x0 = 11.5, y0 = ymax + interv + 14.5, x1 = 17.5, y1 = ymax + interv + 14.5, lwd=4, xpd=T)
  text(x = 14.5, adj = 0.5, y = ymax + interv + 17.5, labels = "Gene-pathway\nmapping", cex=1.50, font=4, xpd=T)
  segments(x0 = 20.5, y0 = ymax + 7.5, x1 = 26.5, y1 = ymax + 7.5, lwd=4, xpd=T)
  text(x = 23.5, adj = 0.5, y = ymax + 10.5, labels = "Variant-pathway\nmapping", cex=1.5, font=4, xpd=T)
  
  #define colors 
  #colors <- brewer.pal(n = 6, name = "Dark2")
  #colors <- c("grey20", "dark red", "steelblue4", "springgreen4", "palevioletred4")
  colors <- c("black", rep("dark red", 5), rep("navy", 5))
  
  #main loop across loci
  for (loc in global.props$locus_id){
    
    #select locus with all gene belonging there
    genes <- subset(dat, dat$locus_id == loc)
    
    #initialize missing annotation variable
    missingAnnot <- FALSE
    #check rowsum in order to put a bar if annotation is missing
    tt <- as.numeric(genes[1, 3:7])
    if (sum(tt) == 0){missingAnnot <- TRUE}
    
    #loop over all genes
    for (gene in 1:nrow(genes)){
      #all squares are 1x1
      #plot weight
      rect(xleft = 9, ybottom = (ymax - c), xright = 10, ytop = (ymax - c + h), col=alpha(colors[1], genes$weight[gene]), border = NULL)
      
      #plot immune system
      rect(xleft = 12, ybottom = (ymax - c), xright = 13, ytop = (ymax - c + h), col=alpha(colors[2], genes$immune[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 12.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot beta amyloid
      rect(xleft = 13, ybottom = (ymax - c), xright = 14, ytop = (ymax - c + h), col=alpha(colors[3], genes$beta.amyloid[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 13.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot endocytosis
      rect(xleft = 14, ybottom = (ymax - c), xright = 15, ytop = (ymax - c + h), col=alpha(colors[4], genes$endocytosis[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 14.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot cholesterol
      rect(xleft = 15, ybottom = (ymax - c), xright = 16, ytop = (ymax - c + h), col=alpha(colors[5], genes$cholesterol.lipid[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 15.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot angiogenesis
      rect(xleft = 16, ybottom = (ymax - c), xright = 17, ytop = (ymax - c + h), col=alpha(colors[6], genes$angiogenesis[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 16.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot gene name
      text(x = 8.5, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 4, adj = 1, xpd=T, cex=0.85)
      
      #increment c
      c <- c + h
    }
    #get locus info for number of genes linked to variant (before filters for missing annotation)
    ggg <- gperv[which(gperv$locus == loc),]
    
    #put text for locus info
    locid <- paste(genes$chr[1], genes$pos[1], sep=":")
    text(x = -6.5, y = ymax - d2 - tmp, labels = locid, cex=1.25, xpd=T, adj = 0, font = 2)
    text(x = 0.25, y = ymax - d2 - tmp, labels = ggg$ngenes, cex=1.25, xpd=T, adj = 1, font = 2)
    
    #draw line underlying locus name
    segments(x0 = -7, y0 = ymax - d2 - tmp - 0.75, x1 = 2.5, y1 = ymax - d2 - tmp - 0.75, xpd=T, lwd=1.5, col="grey60")
    
    #if first locus --> draw line at the very top
    if (loc == 1){
      #segments(x0 = 6, y0 = ymax, x1 = 15, y1 = ymax, lwd=1.5, col="grey60")
      #segments(x0 = 2.5, y0 = ymax - d2 - tmp + 0.75, x1 = 6,
      #         y1 = ymax, lwd=1.5, col="grey60")
      #segments(x0 = 15, y0 = ymax, x1 = 21, y1 = ymax - diff - counter + p, lwd=1.5, col="grey60")
    }
    
    #draw line to divide loci
    segments(x0 = 6, y0 = (ymax - c + h), x1 = 17, y1 = (ymax - c + h), xpd=T, lwd=1.5, col="grey60")
    
    #draw lines from locus column to heatmap
    segments(x0 = 2.5, y0 = ymax - d2 - tmp - 0.75, x1 = 6, y1 = (ymax - c + h), lwd=1.5, col="grey60")
    
    #connection lines to the global heatmap
    segments(x0 = 17, y0 =  (ymax - c + h), x1 = 21, y1 = ymax - diff - counter - p, 
             lwd=1.5, col="grey60")
    
    #now plot the global heatmap
    glob <- global.props[which(global.props$locus_id == loc),]
    
    #plot immune system
    rect(xleft = 21, ybottom = ymax - diff - counter - p, xright = 22, 
         ytop = ymax - diff - counter + p, col=alpha(colors[7], glob$immune.system), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 21.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot beta amyloid
    rect(xleft = 22, ybottom = ymax - diff - counter - p, xright = 23, 
         ytop = ymax - diff - counter + p, col=alpha(colors[8], glob$beta.amyloid), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 22.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot endocytosis
    rect(xleft = 23, ybottom = ymax - diff - counter - p, xright = 24, 
         ytop = ymax - diff - counter + p, col=alpha(colors[9], glob$endocytosis), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 23.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot cholesterol
    rect(xleft = 24, ybottom = ymax - diff - counter - p, xright = 25, 
         ytop = ymax - diff - counter + p, col=alpha(colors[10], glob$cholesterol.lipid), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 24.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot angiogenesis
    rect(xleft = 25, ybottom = ymax - diff - counter - p, xright = 26, 
         ytop = ymax - diff - counter + p, col=alpha(colors[11], glob$vascular), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 25.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #text on the right
    text(x = 26.5, y = ymax - diff - counter, labels = paste("~", rsid$SNP[which(rsid$LOCUS == loc)]), cex=1.25, xpd=T, adj = 0, font = 4)
    
    counter <- counter + p*2 
    tmp <- tmp + 2
    #plot gene name
    #text(x = 8.5, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 4, adj = 1, xpd=T)
    
  }
}
png('Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/Figure2.png', height = 17, width = 14, res=500, units='in')
function.heatmapGeneLevel.var2(dat, global.props, gperv)
dev.off()

write.table(red.genes, "Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/Table_S10bis.txt", quote=F, row.names = F, dec=",", sep="\t")
write.table(red, "Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/Table_S11bis.txt", quote=F, row.names = F, dec=",", sep="\t")
#########TILL HERE SSHOULD BE FINE

#write outputs
#global proportions of pathways per locus
#write.table(global.prop.flt, "Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/22_final_proportions.txt", sep='\t', quote=F, row.names = F)

#list of genes
#gene.list <- matrix(as.character(prop.flt$gene), nrow = 1)
#write.table(gene.list, 'Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/list_of_genes_ADloci.txt', row.names = F, col.names = F, quote=F, sep=',')

#print("Done. Output is '22_final_proportions.txt'")
#end

#PLOTS
#THIS IS FOR THE VARIANT-PATHWAY MAPPING PLOT
#read data -- supplementary data
#dat <- read.xlsx2("Desktop/2k19_work/SurvivalEffect_ADvariants/plots_survivalEffect/new_batch/manuscript_secondProject/20190607_paperPRS_supplementary.xlsx",
#                  sheetIndex = 7, stringsAsFactors = F, startRow = 4)
dat <- read.table("Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/22_final_proportions_perGene.txt", h=T, sep="\t", dec=",")
colnames(dat) <- c("locus_number", "chr", "position", "gene", "gene_weight", "immune", "beta-amyloid",
                   "endocytosis", "cholesterol", "vascular", "unknown", "removed")
dat$locid <- paste(dat$chr, dat$position, sep=":")
dat$removed <- NULL

function.heatmapGeneLevel.var2 <- function(dat, global.props, gperv){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(0, 11, 17, 6))
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #define special counter for labels on the left
  tmp <- 0.25
  
  #define general counter 
  counter <- 0
  
  #general counter for heatmat on the right height
  p <- 1
  
  #set difference between labels and heatmap
  interv <- 1.5
  
  #define max locus number
  max.locus <- max(dat$locus_number)
  
  #define ymax
  ymax <- nrow(dat)*h
  
  #set difference between central part and other parts -- this is for the right part
  diff <- (nrow(dat) - nrow(global.props)) / 10 * 2
  
  #set different between central part and left part -- different from above
  d2 <- (nrow(dat) - nrow(global.props)) / 10 * 2
  
  #prepare window to plot
  plot(1, xlim = c(0, 27),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = -5, y = ymax - d2 + (interv*2), labels = "Locus", adj=0.5, xpd=T, cex=1.50, font=2)
  text(x = 0, y = ymax - d2 + (interv*2), labels = "N. genes", adj=0.5, xpd=T, cex=1.50, font=2)
  text(x = 8.5, y = ymax + interv + 1, labels = "Genes", adj=1, xpd=T, cex=1.50, font=2, srt=-90)
  text(x = 9.5, y = ymax + interv + 1, labels = "Weight", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 12.5, y = ymax + interv, labels = "Immune response", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 13.5, y = ymax + interv, labels = "Beta-amyloid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 14.5, y = ymax + interv, labels = "Endocytosis", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 15.5, y = ymax + interv, labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 16.5, y = ymax + interv, labels = "Vascular dys.", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 21.5, y = ymax - diff + interv*2, labels = "Immune response", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 22.5, y = ymax - diff + interv*2, labels = "Beta-amyloid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 23.5, y = ymax - diff + interv*2, labels = "Endocytosis", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 24.5, y = ymax - diff + interv*2, labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  text(x = 25.5, y = ymax - diff + interv*2, labels = "Angiogenesis.", adj = 1, xpd=T, srt=-90, cex=1.50, font=2)
  segments(x0 = 7.5, y0 = ymax + interv + 14.5, x1 = 10.5, y1 = ymax + interv + 14.5, lwd=4, xpd=T)
  text(x = 9, adj = 0.5, y = ymax + interv + 17.5, labels = "Variant-gene\nmapping", cex=1.50, font=4, xpd=T)
  segments(x0 = 11.5, y0 = ymax + interv + 14.5, x1 = 17.5, y1 = ymax + interv + 14.5, lwd=4, xpd=T)
  text(x = 14.5, adj = 0.5, y = ymax + interv + 17.5, labels = "Gene-pathway\nmapping", cex=1.50, font=4, xpd=T)
  segments(x0 = 20.5, y0 = ymax + 7.5, x1 = 26.5, y1 = ymax + 7.5, lwd=4, xpd=T)
  text(x = 23.5, adj = 0.5, y = ymax + 10.5, labels = "Variant-pathway\nmapping", cex=1.5, font=4, xpd=T)
  
  #define colors 
  #colors <- brewer.pal(n = 6, name = "Dark2")
  #colors <- c("grey20", "dark red", "steelblue4", "springgreen4", "palevioletred4")
  colors <- c("black", rep("dark red", 5), rep("navy", 5))
  
  #main loop across loci
  for (loc in 1:29){
    
    #select locus with all gene belonging there
    genes <- subset(dat, dat$locus_number == loc)
    
    #initialize missing annotation variable
    missingAnnot <- FALSE
    #check rowsum in order to put a bar if annotation is missing
    tt <- as.numeric(genes[1, 6:10])
    if (sum(tt) == 0){missingAnnot <- TRUE}
    
    #loop over all genes
    for (gene in 1:nrow(genes)){
      #assign proportions per gene
      if (missingAnnot == FALSE){
        genes$immune[gene] <- genes$immune[gene] / sum(genes[gene, 6:10])
        genes$`beta-amyloid`[gene] <- genes$`beta-amyloid`[gene] / sum(genes[gene, 6:10])
        genes$endocytosis[gene] <- genes$endocytosis[gene] / sum(genes[gene, 6:10])
        genes$cholesterol[gene] <- genes$cholesterol[gene] / sum(genes[gene, 6:10])
        genes$vascular[gene] <- genes$vascular[gene] / sum(genes[gene, 6:10])
      }  
      
      #all squares are 1x1
      #plot weight
      rect(xleft = 9, ybottom = (ymax - c), xright = 10, ytop = (ymax - c + h), col=alpha(colors[1], genes$gene_weight[gene]), border = NULL)
      
      #plot immune system
      rect(xleft = 12, ybottom = (ymax - c), xright = 13, ytop = (ymax - c + h), col=alpha(colors[2], genes$immune[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 12.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot beta amyloid
      rect(xleft = 13, ybottom = (ymax - c), xright = 14, ytop = (ymax - c + h), col=alpha(colors[3], genes$`beta-amyloid`[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 13.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot endocytosis
      rect(xleft = 14, ybottom = (ymax - c), xright = 15, ytop = (ymax - c + h), col=alpha(colors[4], genes$endocytosis[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 14.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot cholesterol
      rect(xleft = 15, ybottom = (ymax - c), xright = 16, ytop = (ymax - c + h), col=alpha(colors[5], genes$cholesterol[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 15.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot angiogenesis
      rect(xleft = 16, ybottom = (ymax - c), xright = 17, ytop = (ymax - c + h), col=alpha(colors[6], genes$vascular[gene]), border = NULL)
      if (missingAnnot == TRUE){text(x = 16.5, y = ymax - c + h/2, labels = "X", col = "red", font = 2)}
      
      #plot gene name
      text(x = 8.5, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 4, adj = 1, xpd=T, cex=0.85)
      
      #increment c
      c <- c + h
    }
    #get locus info for number of genes linked to variant (before filters for missing annotation)
    ggg <- gperv[which(gperv$locus == genes$locid),]
    
    #put text for locus info
    locid <- paste(genes$chr[1], genes$position[1], sep=":")
    text(x = -6.5, y = ymax - d2 - tmp, labels = locid, cex=1.25, xpd=T, adj = 0, font = 2)
    text(x = 0.25, y = ymax - d2 - tmp, labels = ggg$ngenes, cex=1.25, xpd=T, adj = 1, font = 2)
    
    #draw line underlying locus name
    segments(x0 = -7, y0 = ymax - d2 - tmp - 0.75, x1 = 2.5, y1 = ymax - d2 - tmp - 0.75, xpd=T, lwd=1.5, col="grey60")
    
    #if first locus --> draw line at the very top
    if (loc == 1){
      #segments(x0 = 6, y0 = ymax, x1 = 15, y1 = ymax, lwd=1.5, col="grey60")
      #segments(x0 = 2.5, y0 = ymax - d2 - tmp + 0.75, x1 = 6,
      #         y1 = ymax, lwd=1.5, col="grey60")
      #segments(x0 = 15, y0 = ymax, x1 = 21, y1 = ymax - diff - counter + p, lwd=1.5, col="grey60")
    }
    
    #draw line to divide loci
    segments(x0 = 6, y0 = (ymax - c + h), x1 = 17, y1 = (ymax - c + h), xpd=T, lwd=1.5, col="grey60")
    
    #draw lines from locus column to heatmap
    segments(x0 = 2.5, y0 = ymax - d2 - tmp - 0.75, x1 = 6, y1 = (ymax - c + h), lwd=1.5, col="grey60")
    
    #connection lines to the global heatmap
    segments(x0 = 17, y0 =  (ymax - c + h), x1 = 21, y1 = ymax - diff - counter - p, 
             lwd=1.5, col="grey60")
    
    #now plot the global heatmap
    glob <- global.props[which(global.props$locus_id == locid),]
    
    #plot immune system
    rect(xleft = 21, ybottom = ymax - diff - counter - p, xright = 22, 
         ytop = ymax - diff - counter + p, col=alpha(colors[7], glob$immune.system), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 21.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot beta amyloid
    rect(xleft = 22, ybottom = ymax - diff - counter - p, xright = 23, 
         ytop = ymax - diff - counter + p, col=alpha(colors[8], glob$beta.amyloid), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 22.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot endocytosis
    rect(xleft = 23, ybottom = ymax - diff - counter - p, xright = 24, 
         ytop = ymax - diff - counter + p, col=alpha(colors[9], glob$endocytosis), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 23.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot cholesterol
    rect(xleft = 24, ybottom = ymax - diff - counter - p, xright = 25, 
         ytop = ymax - diff - counter + p, col=alpha(colors[10], glob$cholesterol.lipid), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 24.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #plot angiogenesis
    rect(xleft = 25, ybottom = ymax - diff - counter - p, xright = 26, 
         ytop = ymax - diff - counter + p, col=alpha(colors[11], glob$vascular), 
         border = NULL)
    if (missingAnnot == TRUE){text(x = 25.5, y = ymax - diff - counter, labels = "X", col = "red", font = 2)}
    
    #text on the right
    text(x = 26.5, y = ymax - diff - counter, labels = paste("~", locid), cex=1.25, xpd=T, adj = 0, font = 4)
    
    counter <- counter + p*2 
    tmp <- tmp + 2
    #plot gene name
    #text(x = 8.5, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 4, adj = 1, xpd=T)
    
  }
}
png('Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/heatmap_geneLevel_t5.png', height = 17, width = 14, res=500, units='in')
function.heatmapGeneLevel.var2(dat, global.props, gperv)
dev.off()


###

#heatmap 1 with annotations and weights at the single gene level including unknown, use prop for this
function.heatmapGeneLevel.includingUnk <- function(prop){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(1, 12, 6, 2))
  
  #define max locus number
  max.locus <- max(prop$locus)
  
  #define ymax
  ymax <- nrow(prop)
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #prepare window to plot
  plot(1, xlim = c(0, 7),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = 0.5, y = (ymax + 0.75), labels = "Weight", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 1.5, y = (ymax + 0.75), labels = "Immune System", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 2.5, y = (ymax + 0.75), labels = "Beta-amyloid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 3.5, y = (ymax + 0.75), labels = "Endocytosis", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 4.5, y = (ymax + 0.75), labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 5.5, y = (ymax + 0.75), labels = "Vascular", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 6.5, y = (ymax + 0.75), labels = "Unknown", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  
  #main loop across loci
  for (loc in 1:max.locus){
    #check if locus is there
    if (!(loc %in% prop$locus)){
      warning("Locus is not there, skipping!")
    } else {
      #select locus with all gene belonging there
      genes <- subset(prop, prop$locus == loc)
      
      #loop over all genes
      for (gene in 1:nrow(genes)){
        #all squares are 1x1
        #plot weight
        rect(xleft = 0, ybottom = (ymax - c), xright = 1, ytop = (ymax - c + h), col=alpha('black', genes$weight[gene]))
        
        #plot immune system
        rect(xleft = 1, ybottom = (ymax - c), xright = 2, ytop = (ymax - c + h), col=alpha('red', genes$immune.system[gene]), border = 'black')
        
        #plot beta amyloid
        rect(xleft = 2, ybottom = (ymax - c), xright = 3, ytop = (ymax - c + h), col=alpha('blue', genes$beta.amyloid[gene]), border = 'black')
        
        #plot endocytosis
        rect(xleft = 3, ybottom = (ymax - c), xright = 4, ytop = (ymax - c + h), col=alpha('green', genes$endocytosis[gene]), border = "black")
        
        #plot cholesterol
        rect(xleft = 4, ybottom = (ymax - c), xright = 5, ytop = (ymax - c + h), col=alpha('purple', genes$cholesterol.lipid[gene]), border = 'black')
        
        #plot angiogenesis
        rect(xleft = 5, ybottom = (ymax - c), xright = 6, ytop = (ymax - c + h), col=alpha('gold3', genes$angiogenesis[gene]), border = 'black')
        
        #plot unknown
        rect(xleft = 6, ybottom = (ymax - c), xright = 7, ytop = (ymax - c + h), col=alpha('dark green', genes$unknown[gene]), border = 'black')
        
        #plot gene name
        text(x = -0.25, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 3, adj = 1, xpd=T)
        
        #increment c
        c <- c + h
      }
      #put text for locus info
      text(x = -2, y = (ymax - (c - nrow(genes) - h) - (h/2 * nrow(genes))), labels = genes$pos[1], cex=1.25, font=1, xpd=T, adj = 1)
      
      #draw line to divide loci
      segments(x0 = -4, y0 = (ymax - c + h), x1 = 0, y1 = (ymax - c + h), xpd=T)
    }
  }
}
png('Desktop/plots_survivalEffect/heatmap_geneLevel_includingUnk.png', height = 38, width = 7, res=500, units='in')
function.heatmapGeneLevel.includingUnk(prop)
dev.off()

#heatmap 1 with annotations and weights at the single gene level, use prop.flt for this
function.heatmapGeneLevel <- function(prop.flt){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(1, 12, 6, 2))

  #define max locus number
  max.locus <- max(prop.flt$locus)
  
  #define ymax
  ymax <- nrow(prop.flt)
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #prepare window to plot
  plot(1, xlim = c(0, 6),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = 0.5, y = (ymax + 0.75), labels = "Weight", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 1.5, y = (ymax + 0.75), labels = "Immune System", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 2.5, y = (ymax + 0.75), labels = "Beta-amyloid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 3.5, y = (ymax + 0.75), labels = "Endocytosis", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 4.5, y = (ymax + 0.75), labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 5.5, y = (ymax + 0.75), labels = "Vascular", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  
  #main loop across loci
  for (loc in 1:max.locus){
    #check if locus is there
    if (!(loc %in% prop.flt$locus)){
      warning("Locus is not there, skipping!")
    } else {
      #select locus with all gene belonging there
      genes <- subset(prop.flt, prop.flt$locus == loc)
      
      #loop over all genes
      for (gene in 1:nrow(genes)){
        #all squares are 1x1
        #plot weight
        rect(xleft = 0, ybottom = (ymax - c), xright = 1, ytop = (ymax - c + h), col=alpha('black', genes$weight[gene]))
        
        #plot immune system
        rect(xleft = 1, ybottom = (ymax - c), xright = 2, ytop = (ymax - c + h), col=alpha('red', genes$immune.system[gene]), border = 'black')
        
        #plot beta amyloid
        rect(xleft = 2, ybottom = (ymax - c), xright = 3, ytop = (ymax - c + h), col=alpha('blue', genes$beta.amyloid[gene]), border = 'black')
        
        #plot endocytosis
        rect(xleft = 3, ybottom = (ymax - c), xright = 4, ytop = (ymax - c + h), col=alpha('green', genes$endocytosis[gene]), border = "black")
        
        #plot cholesterol
        rect(xleft = 4, ybottom = (ymax - c), xright = 5, ytop = (ymax - c + h), col=alpha('purple', genes$cholesterol.lipid[gene]), border = 'black')
        
        #plot angiogenesis
        rect(xleft = 5, ybottom = (ymax - c), xright = 6, ytop = (ymax - c + h), col=alpha('gold3', genes$angiogenesis[gene]), border = 'black')
    
        #plot gene name
        text(x = -0.25, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 3, adj = 1, xpd=T)
        
        #increment c
        c <- c + h
      }
      #put text for locus info
      text(x = -2, y = (ymax - (c - nrow(genes) - h) - (h/2 * nrow(genes))), labels = genes$pos[1], cex=1.25, font=1, xpd=T, adj = 1)
      
      #draw line to divide loci
      segments(x0 = -4, y0 = (ymax - c + h), x1 = 0, y1 = (ymax - c + h), xpd=T)
    }
  }
}
png('Desktop/plots_survivalEffect/heatmap_geneLevel.png', height = 20, width = 7, res=500, units='in')
function.heatmapGeneLevel(prop.flt)
dev.off()

#heatmap 2 with annotations at the single locus level, use prop.flt for this
function.heatmapLocusLevel <- function(prop.flt, global.prop.flt){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(1, 12, 6, 2))
  
  #define max locus number
  max.locus <- max(prop.flt$locus)
  
  #define ymax
  ymax <- nrow(prop.flt)
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #prepare window to plot
  plot(1, xlim = c(0, 5),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = 0.5, y = (ymax + 0.75), labels = "Immune System", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 1.5, y = (ymax + 0.75), labels = "Beta-amyloid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 2.5, y = (ymax + 0.75), labels = "Endocytosis", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 3.5, y = (ymax + 0.75), labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 4.5, y = (ymax + 0.75), labels = "Vascular", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  
  #main loop across loci
  for (loc in 1:max.locus){
    #check if locus is there
    if (!(loc %in% prop.flt$locus)){
      warning("Locus is not there, skipping!")
    } else {
      #select locus with all gene belonging there
      genes <- subset(prop.flt, prop.flt$locus == loc)
      glob <- subset(global.prop.flt, global.prop.flt$locus == loc)
      
      #loop over all genes
      for (gene in 1:nrow(genes)){
        #all squares are 1x1
        #plot immune system
        rect(xleft = 0, ybottom = (ymax - c), xright = 1, ytop = (ymax - c + h), col=alpha('red', glob$immune.system), border = NA)
        
        #plot beta amyloid
        rect(xleft = 1, ybottom = (ymax - c), xright = 2, ytop = (ymax - c + h), col=alpha('blue', glob$beta.amyloid), border = NA)
        
        #plot endocytosis
        rect(xleft = 2, ybottom = (ymax - c), xright = 3, ytop = (ymax - c + h), col=alpha('green', glob$endocytosis), border = NA)
        
        #plot cholesterol
        rect(xleft = 3, ybottom = (ymax - c), xright = 4, ytop = (ymax - c + h), col=alpha('purple', glob$cholesterol.lipid), border = NA)
        
        #plot angiogenesis
        rect(xleft = 4, ybottom = (ymax - c), xright = 5, ytop = (ymax - c + h), col=alpha('gold3', glob$vascular), border = NA)
        
        #plot gene name
        text(x = -0.25, y = (ymax - c + (h/2)), labels = genes$gene[gene], font = 3, adj = 1, xpd=T)
        
        #increment c
        c <- c + h
      }
      #put text for locus info
      text(x = -1.5, y = (ymax - (c - nrow(genes) - h) - (h/2 * nrow(genes))), labels = genes$pos[1], cex=1.25, font=1, xpd=T, adj = 1)
      
      #draw line to divide loci
      segments(x0 = -4, y0 = (ymax - c + h), x1 = 5, y1 = (ymax - c + h), xpd=T)
    }
  }
  #make vertical lines between pathways and first horizontal line
  for (i in 0:5){
    segments(x0 = i, y0 = ymax, x1 = i, y1 = (ymax - c + h), xpd=T)
  }
  segments(x0 = 0, y0 = ymax, x1 = 5, y1 = ymax, xpd=T)
}
png('Desktop/plots_survivalEffect/heatmap_locusLevel.png', height = 20, width = 7, res=500, units='in')
function.heatmapLocusLevel(prop.flt, global.prop.flt)
dev.off()

#heatmap 3 with proportions at the locus level, but one cell per locus
function.heatmapLocusLevelCondensed <- function(global.prop.flt){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(1, 14, 7, 2))
  
  #set map
  map <- read.table('Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/mapping_global_prop.txt', h=F, sep='\t')
  global.prop.flt <- merge(global.prop.flt, map, by.x='pos', by.y='V1')
  global.prop.flt <- global.prop.flt[!duplicated(global.prop.flt$V2),]
  global.prop.flt <- subset(global.prop.flt, global.prop.flt$V2 != "-")
  
  #define max locus number
  max.locus <- max(global.prop.flt$locus)
  
  #define ymax
  ymax <- nrow(global.prop.flt)
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #prepare window to plot
  plot(1, xlim = c(0, 5),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = 0.5, y = (ymax + 0.75), labels = "Immune System", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 1.5, y = (ymax + 0.75), labels = "Beta-amyloid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 2.5, y = (ymax + 0.75), labels = "Endocytosis", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 3.5, y = (ymax + 0.75), labels = "Cholesterol/Lipid", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = 4.5, y = (ymax + 0.75), labels = "Vascular", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  
  #main loop across loci
  for (loc in 1:max.locus){
    #check if locus is there
    if (!(loc %in% global.prop.flt$locus)){
      warning("Locus is not there, skipping!")
    } else {
      #select locus
      glob <- subset(global.prop.flt, global.prop.flt$locus == loc)
      
      #plot immune system
      rect(xleft = 0, ybottom = (ymax - c), xright = 1, ytop = (ymax - c + h), col=alpha('red', glob$immune.system), border = 'black')
        
      #plot beta amyloid
      rect(xleft = 1, ybottom = (ymax - c), xright = 2, ytop = (ymax - c + h), col=alpha('blue', glob$beta.amyloid), border = 'black')
        
      #plot endocytosis
      rect(xleft = 2, ybottom = (ymax - c), xright = 3, ytop = (ymax - c + h), col=alpha('green', glob$endocytosis), border = 'black')
        
      #plot cholesterol
      rect(xleft = 3, ybottom = (ymax - c), xright = 4, ytop = (ymax - c + h), col=alpha('purple', glob$cholesterol.lipid), border = 'black')
        
      #plot angiogenesis
      rect(xleft = 4, ybottom = (ymax - c), xright = 5, ytop = (ymax - c + h), col=alpha('gold3', glob$vascular), border = 'black')
        
      #plot locus name
      text(x = -0.25, y = (ymax - c + (h/2)), labels = glob$pos, font = 1, adj = 1, xpd=T)
      
      #plot gene name
      text(x = -1.90, y = (ymax - c + (h/2)), labels = glob$V2, font = 1, adj = 1, xpd=T)
      
      #increment c
      c <- c + h
    
      #draw line to divide loci
      segments(x0 = -2, y0 = (ymax - c + h), x1 = 5, y1 = (ymax - c + h), xpd=T)
    }
  }
}
png('Desktop/plots_survivalEffect/heatmap_locusLevel_condensed.png', height = 12, width = 7, res=500, units='in')
function.heatmapLocusLevelCondensed(global.prop.flt)
dev.off()

#take the previous heatmap and make the transpose of that
function.heatmapLocusLevelCondensed.transposed <- function(global.prop.flt){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(1, 10, 7, 2))
  
  #set map
  map <- read.table('Desktop/myPapers/ROAD_TO_SECOND_PAPER/Translational Psychiatry/rebuttal/mapping_global_prop.txt', h=F, sep='\t')
  global.prop.flt <- merge(global.prop.flt, map, by.x='pos', by.y='V1')
  global.prop.flt <- global.prop.flt[!duplicated(global.prop.flt$V2),]
  global.prop.flt <- subset(global.prop.flt, global.prop.flt$V2 != "-")
  
  #define ymax
  xmax <- nrow(global.prop.flt)

  #prepare window to plot
  plot(1, xlim = c(0, xmax),  ylim=c(0, 5), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = -1, y = 0.5, labels = "Immune System", adj = 1, xpd=T, cex=1.50, font=1)
  text(x = -1, y = 1.5, labels = "Beta-amyloid", adj = 1, xpd=T, cex=1.50, font=1)
  text(x = -1, y = 2.5, labels = "Endocytosis", adj = 1, xpd=T, cex=1.50, font=1)
  text(x = -1, y = 3.5, labels = "Cholesterol/Lipid", adj = 1, xpd=T, cex=1.50, font=1)
  text(x = -1, y = 4.5, labels = "Vascular", adj = 1, xpd=T, cex=1.50, font=1)
  
  for (i in 1:nrow(global.prop.flt)){
    #select locus
    glob <- global.prop.flt[i, ]
    
    #plot immune system
    rect(xleft = i-1, ybottom = 0, xright = i, ytop = 1, col=alpha('red', glob$immune.system), border = 'black')
    
    #plot beta amyloid
    rect(xleft = i-1, ybottom = 1, xright = i, ytop = 2, col=alpha('blue', glob$beta.amyloid), border = 'black')
    
    #plot endocytosis
    rect(xleft = i-1, ybottom = 2, xright = i, ytop = 3, col=alpha('green', glob$endocytosis), border = 'black')
    
    #plot cholesterol
    rect(xleft = i-1, ybottom = 3, xright = i, ytop = 4, col=alpha('purple', glob$cholesterol.lipid), border = 'black')
    
    #plot angiogenesis
    rect(xleft = i-1, ybottom = 4, xright = i, ytop = 5, col=alpha('gold3', glob$vascular), border = 'black')
    
    #plot locus name
    text(x = i-0.5, y = 5.25, labels = paste(glob$V2, glob$pos, sep="  "), font = 3, adj = 1, xpd=T, srt=-45)
  }
}
png('Desktop/plots_survivalEffect/heatmap_locusLevel_condensed_transposed.png', height = 4, width = 12, res=500, units='in')
function.heatmapLocusLevelCondensed.transposed(global.prop.flt)
dev.off()


#EXPRESSION DATA FROM BRAIN
#read input data -- log -> norm (per-sample) --> norm (across samples) [select indvs and genes after these steps]
inp <- read.csv('Downloads/GSE73721_Human_and_mouse_table.csv', h=T, check.names = F)

#convert gene names into upper case and manage rownames
inp$Gene <- toupper(inp$Gene)
rownames(inp) <- inp$Gene
inp$Gene <- NULL

#put colnames somewhere else also
samples <- colnames(inp)

#things to include
microglia.index <- grep("myeloid", samples)
neuron.index <- grep("neuron", samples)
oligo.index <- grep("oligo", samples)
endo.index <- grep("endo", samples)

#astrocytes is a bit more complex
astro.index <- grep("astro", samples)
astro.entries <- inp[, astro.index]
astro.final.index <- grep("ctx", colnames(astro.entries))
astro.semifinal <- astro.entries[, astro.final.index]
exclude <- grep("Fetal", colnames(astro.semifinal))
exclude.id <- colnames(astro.semifinal)[exclude]
include.id <- colnames(astro.semifinal)
include.id <- subset(include.id, !(include.id %in% exclude.id))
final.astro <- astro.semifinal[, include.id]
rm(list = c("astro.index", "astro.entries", "astro.final.index", "astro.semifinal", "exclude", "exclude.id", "include.id"))

#list of columns to add
include.cols <- c(colnames(inp)[microglia.index], colnames(inp)[neuron.index], colnames(inp)[oligo.index], colnames(inp)[endo.index], colnames(final.astro))

#gene list
gene.list <- data.frame(gene = prop.flt$gene, locus = prop.flt$locus, weight = prop.flt$weight, pos = prop.flt$pos)

#subset on samples
inp.samples <- inp[, include.cols]

#take log transformation
inp.samples.log <- log(inp.samples)
inp.samples.log.genes <- merge(inp.samples.log, gene.list, by.x='row.names', by.y='gene')
unmapped <- subset(gene.list, !(gene.list$gene %in% inp.samples.log.genes$Row.names))

#plot heatmap and distribution after log-transformation
par(new=F)
heatmap(as.matrix(inp.samples.log.genes[2:24]), col = colorRampPalette(c("blue", "white", "red"))(200), Colv = NA, Rowv = NA, scale='none', labRow = inp.samples.log.genes$Row.names)
#distribution
par(new=F)
par(mar=c(4, 6, 6, 4))
for (i in 1:ncol(inp.samples.log)){
  den <- density(inp.samples.log[, i])
  plot(den, type='l', col='blue', ylim=c(0, 1), xlim=c(-5, 10))
  par(new=T)
  
}

#scale across samples (samples as rows) 
inp.samples.log.normAcross <- t(scale(x = as.matrix(t(inp.samples.log)), center = T, scale = T))
inp.samples.log.normAcross.genes <- merge(inp.samples.log.normAcross, gene.list, by.x='row.names', by.y='gene')
inp.samples.log.normAcross.genes <- inp.samples.log.normAcross.genes[order(inp.samples.log.normAcross.genes$locus),]
unmapped <- subset(gene.list, !(gene.list$gene %in% inp.samples.log.normAcross.genes$Row.names))

#plot heatmap and distribution after normalization across samples
#distribution
par(new=F)
par(mar=c(4, 6, 6, 4))
for (i in 1:ncol(inp.samples.log.normAcross)){
  den <- density(inp.samples.log[, i])
  plot(den, type='l', col='blue', ylim=c(0, 1), xlim=c(-5, 10))
  par(new=T)
  
}

par(new=F)
par(mar=c(4, 6, 6, 4))
plot(1)
png("Desktop/plots_survivalEffect/new_batch/expressionData_microglia.png", height = 15, width = 12, res=500, units='in')
a <- heatmap(as.matrix(inp.samples.log.normAcross.genes[,2:24]), col = colorRampPalette(c("blue", "white", "red"))(200), Colv = NA, Rowv = NA, scale='none', RowSideColors = as.character(inp.samples.log.normAcross.genes$locus), labRow = inp.samples.log.normAcross.genes$Row.names)
dev.off()

toplot <- inp.samples.log.normAcross.genes
toplot.noDups <- toplot[!duplicated(toplot$Row.names),]

#heatmap 1 with annotations and weights at the single gene level, use prop.flt for this
function.heatmapGeneLevel <- function(toplot){
  #set graphical parameters
  par(mfrow=c(1,1))
  par(mar=c(1, 12, 8, 2))
  
  #remove na
  toplot <- na.omit(toplot)
  labels <- toplot[, "Row.names"]
  toplot$Row.names <- NULL
  
  #define max locus number
  max.locus <- max(toplot$locus)
  
  #define ymax
  ymax <- nrow(toplot)
  
  #define square height
  h <- 1
  
  #define counter for the plot
  c <- 1
  
  #prepare window to plot
  plot(1, xlim = c(0, 24),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', yaxt='none', col='white', bty='n')
  
  #put text on columns
  text(x = 0.5, y = (ymax + 1.75), labels = "Weight", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = mean(seq(1.5, 13, 1)), y = (ymax + 1.75), labels = "Astrocytes", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = mean(seq(13.5, 18, 1)), y = (ymax + 1.75), labels = "Oligodendrocytes", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = mean(seq(18.5, 21, 1)), y = (ymax + 1.75), labels = "Microglia", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = mean(seq(21.5, 23, 1)), y = (ymax + 1.75), labels = "Endotelial cells", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  text(x = mean(seq(23.5, 24, 1)), y = (ymax + 1.75), labels = "Neuron", adj = 1, xpd=T, srt=-45, cex=1.50, font=1)
  
  #put rectangles below samples name
  palette_pastel <- brewer.pal(n = 5, name = "Pastel1")
  rect(xleft = 1, ybottom = (ymax + 0.25), xright = 13, ytop = (ymax + 1.25), col=palette_pastel[1])
  rect(xleft = 13, ybottom = (ymax + 0.25), xright = 18, ytop = (ymax + 1.25), col=palette_pastel[2])
  rect(xleft = 18, ybottom = (ymax + 0.25), xright = 21, ytop = (ymax + 1.25), col=palette_pastel[3])
  rect(xleft = 21, ybottom = (ymax + 0.25), xright = 23, ytop = (ymax + 1.25), col=palette_pastel[4])
  rect(xleft = 23, ybottom = (ymax + 0.25), xright = 24, ytop = (ymax + 1.25), col=palette_pastel[5])
  
  #re-order toplot file according to my order
  positions <- toplot$pos
  order.my <- c("8yo ctx astro", "13yo ctx astro", "16yo ctx astro", "21yo ctx astro", "22yo ctx astro", "35yo ctx astro", "47yo ctx astro", "51yo ctx astro", "53yo ctx astro", "60yo ctx astro", "63yo ctx 1 astro", "63yo ctx 2 astro", "22yo ctx oligo", "47yo ctx oligo", "63yo ctx A oligo", "63yo ctx B oligo", "63yo ctx 3 oligo", "45yo ctx myeloid", "51yo ctx myeloid", "63 yo ctx myeloid", "13yo ctx endo", "47 yo ctx endo", "25yo ctx neuron", "locus", "weight")
  toplot.t <- t(toplot)
  toplot.t <- function.match(toplot.t, order.my, rownames(toplot.t))
  toplot <- as.data.frame(t(toplot.t))

  #convert to numbers
  for (col in 1:ncol(toplot)){
    toplot[, col] <- as.numeric(as.character(toplot[, col]))
  }
  
  #define general counter
  counter <- 0

  #main loop across loci
  for (loc in 1:max.locus){
    print(loc)
    #check if locus is there
    if (!(loc %in% toplot$locus)){
      warning("Locus is not there, skipping!")
    } else {
      #select locus with all gene belonging there
      genes <- subset(toplot, toplot$locus == loc)
    
      #nomrmalize between 0 and 1 depending on positive (red) and negative (blue) with white in the middle
      heat <- as.matrix(genes[, 1:23])
        
      #loop over all genes
      for (gene in 1:nrow(genes)){
        #increment counter
        counter <- counter + 1
        
        #all squares are 1x1
        #plot weight
        rect(xleft = 0, ybottom = (ymax - c), xright = 1, ytop = (ymax - c + h), col=alpha('black', genes$weight[gene]))
        
        for (sam in 2:24){
          #normalize expr
          expr <- heat[gene, sam-1]
          if (expr > 0){
            norm.expr <- (expr - 0) / (max(toplot[,1:23]) - 0)
            col.expr <- "dark red"
          } else {
            norm.expr <- 1 - (expr - min(toplot)) / (0 - min(toplot))
            col.expr <- "dark blue"
          }

          rect(xleft = sam - 1, ybottom = (ymax - c), xright = sam, ytop = (ymax - c + h), col=alpha(col.expr, norm.expr), border = 'black')
          }
      
        #plot gene name
        text(x = -0.25, y = (ymax - c + (h/2)), labels = labels[counter], font = 3, adj = 1, xpd=T)
      
        #increment c
        c <- c + h
      }
      #put text for locus info
      text(x = -2.5, y = (ymax - (c - nrow(genes) - h) - (h/2 * nrow(genes))), labels = positions[counter], cex=1.25, font=3, xpd=T, adj = 1)
      
      #draw line to divide loci
      segments(x0 = -5, y0 = (ymax - c + h), x1 = 24, y1 = (ymax - c + h), xpd=T, lwd=2)
    }
    segments(x0 = 1, y0 = (ymax - c + h), x1 = 1, y1 = ymax, lty = 1, lwd=2)
    segments(x0 = 13, y0 = (ymax - c + h), x1 = 13, y1 = ymax, lty = 1, lwd=2)
    segments(x0 = 18, y0 = (ymax - c + h), x1 = 18, y1 = ymax, lty = 1, lwd=2)
    segments(x0 = 21, y0 = (ymax - c + h), x1 = 21, y1 = ymax, lty = 1, lwd=2)
    segments(x0 = 23, y0 = (ymax - c + h), x1 = 23, y1 = ymax, lty = 1, lwd=2)
    segments(x0 = 24, y0 = (ymax - c + h), x1 = 24, y1 = ymax, lty = 1, lwd=2)
  }
}
png('Desktop/plots_survivalEffect/new_batch/heatmap_geneLevel_expressionData.png', height = 25, width = 12, res=500, units='in')
function.heatmapGeneLevel(toplot.noDups)
dev.off()




