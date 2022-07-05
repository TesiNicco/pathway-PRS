import os
import sys

#read mapping file, it is important for weighting more the more consistently annotated genes
def read_gene_mapping():
    my_gene = {}
    finp = open("3_merged_annotation_IGAP_IRIS_OTHER.txt")
    for line in finp:
        if line.startswith("locus"):
            pass
        else:
            line = line.rstrip().split()
            locus, pos, genes = line[0], line[1], line[2]
            genes = genes.replace(";", " ")
            genes = genes.split()
            my_gene[(locus, pos)] = genes
    return my_gene

#read generic annotation file (they all should be in the same format)
def read_anno(finp, annotation, function_list):
    finp = open(finp).readlines()
    for line in finp:
        if line.startswith("function"):
            pass
        else:
            line = line.rstrip().split("\t")
            function, genes = line[0], line[1]
            function = function.lower()
            if function not in function_list:
                print ("!!! Function was unknown !!!", function)
            else:
                counts = annotation[function]
                genes = genes.replace(",", " ")
                genes = genes.split()
                for gene in genes:
                    if gene in counts.keys():
                        counts[gene] = counts[gene] + 1
                    else:
                        counts[gene] = 1
            annotation[function] = counts

    return annotation

#write output per per-locus, per-gene
def write_out(annot_complete, my_genes):
    out = open("20_complete_annotation.txt", "w")
    out.write("locus\tpos\tgene\tweight\timmune system\tbeta-amyloid\tendocytosis\tcholesterol/lipid\tangiogenesis\tunknown\n")
    for x in my_genes.keys():
        locus, pos = int(x[0]), x[1]
        locus += 1
        for gene in my_genes[x]:
            out.write(str(locus))
            out.write("\t")
            out.write(str(pos))
            out.write("\t")
            immune_counts = 0
            beta_counts = 0
            endo_counts = 0
            chol_counts = 0
            angio_counts = 0
            flag = 0
            flag_multiMap = 0
            if ":" in gene:
                flag_multiMap = 1
                gene_spt = gene.split(":")
                gene = gene_spt[0]
                weight = gene_spt[1]
            if gene in annot_complete["immune system"].keys():
                immune_counts = annot_complete["immune system"][gene]
                flag = 1
            if gene in annot_complete["beta-amyloid"].keys():
                beta_counts = annot_complete["beta-amyloid"][gene]
                flag = 1
            if gene in annot_complete["endocytosis"].keys():
                endo_counts = annot_complete["endocytosis"][gene]
                flag = 1
            if gene in annot_complete["cholesterol/lipid"].keys():
                chol_counts = annot_complete["cholesterol/lipid"][gene]
                flag = 1
            if gene in annot_complete["angiogenesis"].keys():
                angio_counts = annot_complete["angiogenesis"][gene]
                flag = 1
            if flag == 0:
                unknown_counts = 1
            else:
                unknown_counts = 0

            if flag_multiMap == 1:
                out.write(gene)
                out.write("\t")
                out.write(str(weight))
            else:
                weight = 1.0
                out.write(gene)
                out.write("\t")
                out.write(str(weight))
            out.write("\t")
            out.write(str(immune_counts))
            out.write("\t")
            out.write(str(beta_counts))
            out.write("\t")
            out.write(str(endo_counts))
            out.write("\t")
            out.write(str(chol_counts))
            out.write("\t")
            out.write(str(angio_counts))
            out.write("\t")
            out.write(str(unknown_counts))
            out.write("\n")
    out.close()
    return "\nOutput is produced. Name is '20_complete_annotation.txt'"

################MAIN

my_genes = read_gene_mapping()

function_list = ["beta-amyloid", "immune system", "endocytosis", "cholesterol/lipid", "angiogenesis"]
all_annotations = ["10_assigned_functions_LITERATURE.txt", "7_assigned_functions_IGAP.txt", "14_assigned_functions_GENE-ONTOLOGY.txt", "17_assigned_functions_DAVID.txt"]
annotation = {"beta-amyloid" : {}, "immune system" : {}, "endocytosis" : {}, "cholesterol/lipid" : {}, "angiogenesis" : {}}

annotation_lit = read_anno("10_assigned_functions_LITERATURE.txt", annotation, function_list)
annotation_lit_igap = read_anno("7_assigned_functions_IGAP.txt", annotation_lit, function_list)
annot_lit_igap_go = read_anno("14_assigned_functions_GENE-ONTOLOGY.txt", annotation_lit_igap, function_list)
annot_complete = read_anno("17_assigned_functions_DAVID.txt", annot_lit_igap_go, function_list)
for annot in annot_complete:
    print (annot, annot_complete[annot])

print (write_out(annot_complete, my_genes))

cmd_sort = "sort -k1,1n 20_complete_annotation.txt > 20_complete_annotation_sorted.txt"
os.system(cmd_sort)
cmd_clean = "rm 20_complete_annotation.txt"
os.system(cmd_clean)

