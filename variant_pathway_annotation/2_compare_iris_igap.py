#!/usr/bin/python2.7

# compare IGAP and IRIS snps to gene mapping procedures

import os
import sys


# read mapping file, store genes mapped in IGAP, UK and other annotation, separately
def read_mapping_file():
    igap = {}
    iris = {}
    other = {}
    inpf = open("1_igap_iris_mapping_genes.txt").readlines()
    for line in inpf:
        genes_igap_clean = []
        genes_iris_clean = []
        genes_other_clean = []
        if line.startswith("LOCUS"):
            pass
        else:
            line = line.rstrip().split()
            locus, igap_genes, iris_genes, other_genes = line[1], line[2], line[3], line[4]
            igap_genes = igap_genes.split(",")
            iris_genes = iris_genes.split(",")
            other_genes = other_genes.split(",")
            for x in igap_genes:
                if len(x) != 0:
                    genes_igap_clean.append(x)
            for x in iris_genes:
                if len(x) != 0:
                    genes_iris_clean.append(x)
            for x in other_genes:
                if len(x) != 0:
                    genes_other_clean.append(x)
            igap[locus] = genes_igap_clean
            iris[locus] = genes_iris_clean
            other[locus] = genes_other_clean

    return igap, iris, other


# compare annotation in IGAP, UK and other annotation, to make a consensus
def compare_annotation(igap, iris, other):
    annotation = {}
    for locus in igap.keys():
        igap_genes = igap[locus]
        iris_genes = iris[locus]
        other_genes = other[locus]
        locus_genes = {}
        for gene in igap_genes:
            if gene != "-":
                locus_genes[gene] = 1
        for gene in iris_genes:
            if gene != "-":
                if gene in locus_genes.keys():
                    locus_genes[gene] = locus_genes[gene] + 1
                else:
                    locus_genes[gene] = 1
        for gene in other_genes:
            if gene != "-":
                if gene in locus_genes.keys():
                    locus_genes[gene] = locus_genes[gene] + 1
                else:
                    locus_genes[gene] = 1
        annotation[locus] = locus_genes

    sum_per_gene = []
    for locus in annotation.keys():
        print locus
        genes_n = 0
        for gene in annotation[locus].keys():
            print "\t %s : %s" % (gene, annotation[locus][gene])
            genes_n = genes_n + annotation[locus][gene]
        sum_per_gene.append(genes_n)

    #now move to proportions (percentages) instead of counts
    index = 0
    for locus in annotation.keys():
        for gene in annotation[locus].keys():
            annotation[locus][gene] = annotation[locus][gene] / float(sum_per_gene[index])
        index += 1
    print annotation


    return annotation


# write output
def write_consensus(annotation):
    out = open("3_merged_annotation_IGAP_IRIS_OTHER.txt", "w")
    out.write("locus\tpos\tannotation\n")
    locus_number = 0
    for locus in annotation.keys():
        out.write(str(locus_number))
        out.write("\t")
        locus_number += 1
        out.write(locus)
        out.write("\t")
        if len(annotation[locus]) == 1:
            gene = "".join(annotation[locus].keys())
            out.write(str(gene))
            out.write(";")
        else:
            for gene in annotation[locus].keys():
                counts = annotation[locus][gene]
                out.write(gene)
                out.write(":")
                out.write(str(counts))
                out.write(";")
        out.write("\n")
    out.close()
    return "Output is produced. Name is '3_merged_annotation_IGAP_IRIS_OTHER.txt'"

def write_geneList(annotation):
    out = open("4_merged_annotation_IGAP_IRIS_OTHER_gene_list.txt", "w")
    for locus in annotation.keys():
        for gene in annotation[locus].keys():
            out.write(gene)
            out.write("\n")
    out.close()
    return "Second output is produces. Name is '4_merged_annotation_IGAP_IRIS_OTHER.txt'"

##########MAIN

igap, iris, other = read_mapping_file()
annotation = compare_annotation(igap, iris, other)
print write_consensus(annotation)
print write_geneList(annotation)
