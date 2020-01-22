#/usr/bin/python2.7

import sys

#read input file and store info into a list
def read_input_file():
    finp = open(sys.argv[1]).readlines()
    line1 = finp[1]
    line1 = line1.rstrip().split()
    go = line1[0]
    gene_list = []
    print "Working on GO: %s" %(go)
    for line in finp:
        if line.startswith("superterm_name"):
            pass
        else:
            line = line.rstrip().split("\t")
            gene = line[10].upper()
            if gene in gene_list:
                pass
            else:
                gene_list.append(gene)
    print "There are %s unique genes associated with this term" %(len(gene_list))
    return gene_list, go

#modified version of function to read input file with the new gene-ontology assessment
def read_go_genes():
  finp = open(sys.argv[1]).readlines()
  print "Working on file: %s" %(sys.argv[1])
  gene_list = []
  for line in finp:
    line = line.rstrip().split("\t")
    gene = line[-1].upper()
    if gene in gene_list:
      pass
    else:
      gene_list.append(gene)
  print "There are %s unique genes associated with this term" %(len(gene_list))
  return gene_list
  
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

#merge information from annotation and my genes
def merge(my_genes, gene_list, go):
    annotated = []
    for g in my_genes:
        if g in gene_list:
            annotated.append(g)
        else:
            pass
    print "\nOut of %s genes in your list, %s were annotated to %s: %s" %(len(my_genes), len(annotated), go, annotated)
    return annotated

#write output
def write_output(annotated, go):
    init = 13
    out_name = str(init) + "_assigned_function_GENE-ONTOLOGY_" + go + ".txt"
    out = open(out_name, "a")
    out.write(go)
    out.write("\t")
    for gene in annotated:
        out.write(gene)
        out.write(",")
    out.write("\n")
    out.close()
    return "\nOutput is produced. Name is '%s'" %(out_name)

############MAIN

#gene_list, go = read_input_file()
gene_list = read_go_genes()
go = sys.argv[2]
my_genes, locus_genes = read_my_genes()
annotated = merge(my_genes, gene_list, go)
print write_output(annotated, go)
