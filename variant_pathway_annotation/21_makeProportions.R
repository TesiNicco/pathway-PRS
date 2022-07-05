#compute proportions of functional annotation

library(RColorBrewer)
library(ggplot2)
rm(list=ls())

#read input file
inp <- read.table("20_complete_annotation_sorted.txt", h=T, sep='\t')

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
  rsid <- read.table("rsid_snps.txt", h=T, dec=",", sep="\t")
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
png('Figure2.png', height = 17, width = 14, res=500, units='in')
function.heatmapGeneLevel.var2(dat, global.props, gperv)
dev.off()

write.table(red.genes, "Table_S10bis.txt", quote=F, row.names = F, dec=",", sep="\t")
write.table(red, "Table_S11bis.txt", quote=F, row.names = F, dec=",", sep="\t")
#########TILL HERE SSHOULD BE FINE

#write outputs
#global proportions of pathways per locus
write.table(global.prop.flt, "22_final_proportions.txt", sep='\t', quote=F, row.names = F)
write.table(weighted.prop, '22_final_proportions_perGene.txt', row.names = F, quote=F, sep='\t')

#list of genes
gene.list <- matrix(as.character(prop.flt$gene), nrow = 1)
write.table(gene.list, 'list_of_genes_ADloci.txt', row.names = F, col.names = F, quote=F, sep=',')

#print("Done. Output is '22_final_proportions.txt'")
#end

# EXPRESSION DATA FROM BRAIN
#read input data -- log -> norm (per-sample) --> norm (across samples) [select indvs and genes after these steps]
inp <- read.csv('GSE73721_Human_and_mouse_table.csv', h=T, check.names = F)

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



