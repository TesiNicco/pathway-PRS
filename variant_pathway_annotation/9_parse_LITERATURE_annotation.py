#!/usr/bin/python2.7

import operator

#read input file and save annotations into a dictionary
def read_input():
    annotations = {}
    finp = open("8_Literature_pathways_with_genes.txt").readlines()
    pre_defined_groups = ["Immune system", "Beta-amyloid", "Cholesterol/lipid", "Endocytosis", "Vascular/Angiogenesis", "Unknown"]
    i = 0
    gene_list = []
    while i < len(finp):
        line = finp[i]
        line = line.rstrip()
        if i == 0:
            go = line
            gene_list = []
        elif i == len(finp) - 1:
            annotations[go] = gene_list
        elif line in pre_defined_groups:
            annotations[go] = gene_list
            go = line
            gene_list = []
        else:
            gene = line
            gene_list.append(gene)
        i += 1
    for go in annotations.keys():
        print go, annotations[go]
    return annotations

#read list of my genes
def read_my_genes():
    finp = open("3_merged_annotation_IGAP_IRIS_OTHER.txt").readlines()
    my_genes = []
    locus_genes = {}
    for line in finp:
        if line.startswith("locus"):
            pass
        else:
            line = line.rstrip().split()
            locus, pos, genes = line[0], line[1], line[2]
            genes = genes.split(";")
            for x in genes:
                x = x.split(":")
                if isinstance(x, list):
                    gene = x[0]
                    if gene in my_genes:
                        pass
                    else:
                        if gene == "":
                            pass
                        else:
                            my_genes.append(gene)
                else:
                    gene = x
                    if gene in my_genes:
                        pass
                    else:
                        if gene == "":
                            pass
                        else:
                            my_genes.append(gene)
            locus_genes[(locus, pos)] = genes
    print "\nIn my list there are %s non-duplicated genes: %s" %(len(my_genes), my_genes)
    return my_genes, locus_genes

#merge my genes and annotated genes and make summary about the merging procedure
def merge(my_genes, annotations):
    my_genes_annotated = {}
    mapped = 0
    unmapped = 0
    unknown = {}
    for gene in my_genes:
        for go in annotations.keys():
            gene_list = annotations[go]
            if gene in gene_list:
                mapped += 1
                if go in my_genes_annotated.keys():
                    if gene in my_genes_annotated[go]:
                        pass
                    else:
                        my_genes_annotated[go].append(gene)
                else:
                    my_genes_annotated[go] = []
                    my_genes_annotated[go].append(gene)
            else:
                pass
    print "\nMapped genes were %s" %(mapped)
    for go in my_genes_annotated.keys():
        print go, my_genes_annotated[go]

    #manage unknown annotations, make a dictionary of the frequency of each word and have a look at the most recurrent
    for go in annotations.keys():
        if go == "Unknown":
            gene_list = annotations[go]
            for gene in gene_list:
                gene = gene.split("(")
                gene_name, go = gene[0], gene[1].replace(")", "")
                unknown[gene_name] = go
    occurrence_terms_unknown = {}
    for gene in unknown.keys():
        function = unknown[gene]
        function = function.split()
        for word in function:
            occurrence_terms_unknown[word] = occurrence_terms_unknown.get(word, 0) + 1
    print "\nMost recurrent words in unknown category --> ", sorted(occurrence_terms_unknown.items(), key=operator.itemgetter(1))
    print "!!Might consider to add a new category!!"
    return my_genes_annotated

#write output
def write_output(my_genes_annotated):
    out = open("10_assigned_functions_LITERATURE.txt", "w")
    out.write("function\tgene\n")
    for go in my_genes_annotated.keys():
        out.write(go)
        out.write("\t")
        for gene in my_genes_annotated[go]:
            out.write(gene)
            out.write(",")
        out.write("\n")
    out.close()
    return "\nOutput is produced. Name is '10_assigned_functions_LITERATURE.txt'"

##########MAIN

annotations = read_input()
my_genes, locus_genes = read_my_genes()
my_genes_annotated = merge(my_genes, annotations)
print write_output(my_genes_annotated)