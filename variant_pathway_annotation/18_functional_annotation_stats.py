#stats about the functional annotation procedure

import sys

def read_input(annotation):
    finp = open(annotation).readlines()
    annot = {}
    all_genes = []
    for line in finp[1:]:
        line = line.rstrip().split("\t")
        go, genes = line[0], line[1]
        genes = genes.split(",")
        for gene in genes:
            if (gene == ''):
                pass
            else:
                if gene not in all_genes:
                    all_genes.append(gene)
                if go in annot.keys():
                    annot[go].append(gene)
                else:
                    annot[go] = []
                    annot[go].append(gene)
    print ("**** Annotation -> %s" %(annotation))
    for go in annot.keys():
        print ("** GO -> %s reported %s genes" %(go, len(annot[go])))
    print ("** Total number of genes is %s" %(len(all_genes)))
    return "\n\n"

file_list = ["17_assigned_functions_DAVID.txt", "14_assigned_functions_GENE-ONTOLOGY.txt", "7_assigned_functions_IGAP.txt", "10_assigned_functions_LITERATURE.txt"]
for f in file_list:
    print (read_input(f))