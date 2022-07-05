# Libraries
import os
import sys

# read IGAP annotation file and save it into a dictionary of GO:genes
def read_annotation_file():
    finp = open(sys.argv[1]).readlines()
    i = 0
    annotation = {}
    while i < len(finp):
        line = finp[i]
        if i == 0:
            line = line.rstrip()
            go = line
            gene_list = []
        elif i == len(finp) - 1:
            annotation[(go_code, go)] = gene_list
        elif line.startswith("GO:"):
            line = line.rstrip().split()
            go_code, gene = line[0], line[1]
            if gene == "NA":
                pass
            else:
                gene_list.append(gene)
        else:
            annotation[(go_code, go)] = gene_list
            line = line.rstrip()
            go = line
            gene_list = []
        i += 1
    for go in annotation.keys():
        print ("There are %s genes associated with GO %s" % (len(annotation[go]), go))
    return annotation

# group found go terms into pre-defined set of functional groups
def group_assigned_go(annotation):
    groups = {}
    pre_defined_groups = ["Immune system", "Beta-amyloid", "Cholesterol/lipid", "Endocytosis", "Vascular/Angiogenesis",
                          "Unknown"]
    print ("\n")
    for go in annotation.keys():
        go_code, go_term = go[0], go[1]
        if ("immune" in go_term) or ("Immune" in go_term):
            print (go, " --> Immune")
            if "Immune system" in groups.keys():
                for gene in annotation[go]:
                    groups["Immune system"].append(gene)
            else:
                groups["Immune system"] = []
                for gene in annotation[go]:
                    groups["Immune system"].append(gene)

        elif ("tau" in go_term) or ("beta" in go_term) or ("amyloid" in go_term):
            print (go, " --> Amyloid")
            if "Beta-amyloid" in groups.keys():
                for gene in annotation[go]:
                    groups["Beta-amyloid"].append(gene)
            else:
                groups["Beta-amyloid"] = []
                for gene in annotation[go]:
                    groups["Beta-amyloid"].append(gene)

        elif ("lipo" in go_term) or ("lipid" in go_term) or ("cholesterol" in go_term):
            print (go, " --> Lipid")
            if "Cholesterol/lipid" in groups.keys():
                for gene in annotation[go]:
                    groups["Cholesterol/lipid"].append(gene)
            else:
                groups["Cholesterol/lipid"] = []
                for gene in annotation[go]:
                    groups["Cholesterol/lipid"].append(gene)

        else:
            print ("No groups found for GO term %s" % (go_term))
    print ("\nFunctional groups before cleaning")
    for group in groups.keys():
        print (group, groups[group])
    return groups

# exclude duplicates from each gene set
def parse_gene_list(groups):
    for group in groups.keys():
        cleaned_list = list(set(groups[group]))
        groups[group] = cleaned_list
    print ("\nFunctional groups after cleaning")
    for group in groups.keys():
        print (group, groups[group])
    return groups

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
    print ("\nIn my list there are %s non-duplicated genes: %s" %(len(my_genes), my_genes))
    return my_genes, locus_genes

#merge my genes and annotated genes and make summary about the merging procedure
def merge(my_genes, cleaned_groups):
    my_genes_annotated = {}
    mapped = 0
    unmapped = 0
    for gene in my_genes:
        for go in cleaned_groups.keys():
            gene_list = cleaned_groups[go]
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
    print ("\nMapped genes were %s" %(mapped))
    print ("Unmapped genes were %s" %(unmapped))
    for go in my_genes_annotated.keys():
        print (go, my_genes_annotated[go])
    return my_genes_annotated

#write output
def write_output(my_genes_annotated):
    out = open("7_assigned_functions_IGAP.txt", "w")
    out.write("function\tgene\n")
    for go in my_genes_annotated.keys():
        out.write(go)
        out.write("\t")
        for gene in my_genes_annotated[go]:
            out.write(gene)
            out.write(",")
        out.write("\n")
    out.close()
    return "Output is produced. Name is '6_assigned_functions_IGAP.txt'"

##########MAIN

annotation = read_annotation_file()
groups = group_assigned_go(annotation)
cleaned_groups = parse_gene_list(groups)
my_genes, locus_genes = read_my_genes()
my_genes_annotated = merge(my_genes, cleaned_groups)
print (write_output(my_genes_annotated))