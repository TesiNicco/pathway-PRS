import sys

#read DAVID annotation file and save content into a dictionary, without parsing the file
def read_david_file():
    david = {}
    finp = open(sys.argv[1]).readlines()
    i = 0
    while i < len(finp):
        line = finp[i]
        if i == 0:
            cluster = {}
            line = line.rstrip().split("\t")
            line_name = line[0]
            cluster_name = line_name.replace(" ", "_")
        elif i == len(finp)-1:
            david[cluster_name] = cluster
        elif line.startswith("Annotation"):
            david[cluster_name] = cluster
            cluster = {}
            line = line.rstrip().split("\t")
            line_name = line[0]
            cluster_name = line_name.replace(" ", "_")
        elif line.startswith("Category"):
            pass
        else:
            line = line.rstrip().split("\t")
            if len(line) < 2:
                pass
            else:
                category, term, genes = line[0], line[1], line[5]
                cluster[term] = [category, genes]
        i += 1
    return david

#assign a name to each cluster according to the most repetitive elements within the term names
def assign_name(david):
    annotation = {}
    groups_keywords = {"immune system" : ["immune", "complement", "antigen", "t cell"], "Beta-amyloid" : ["beta-amyloid", "alzheimer", "amyloidosis"], "Cholesterol/lipid" : ["cholesterol", "lipid", "protein-lipid"], "Endocytosis" : ["endocytosis", "proteasome", "proteasomal"], "Vascular/Angiogenesis" : ["angiogenesis", "vascular"]}
    for cluster in david.keys():
        assign_cluster = {"immune system" : 0, "Beta-amyloid" : 0, "Cholesterol/lipid" : 0, "Endocytosis" : 0, "Vascular/Angiogenesis" : 0, "Unknown" : 0}
        keywords_mapped = []
        for term in david[cluster].keys():
            term = term.replace("~", " ")
            term = term.lower()
            print (term, ";",)
            flag = 0
            for group, value in groups_keywords.items():
                for x in value:
                    if x in term:
                        assign_cluster[group] = assign_cluster[group] + 1
                        keywords_mapped.append(x)
                        flag = 1
        if flag == 0:
            assign_cluster["Unknown"] = 1

        most_rep = 0
        most_rep_lab = ""
        ex_equo = []
        for group in assign_cluster.keys():
            if assign_cluster[group] > most_rep:
                most_rep = assign_cluster[group]
                most_rep_lab = group
                annotation[cluster] = [most_rep_lab]
        for group in assign_cluster.keys():
            if assign_cluster[group] == most_rep:
                if ((group != most_rep_lab) and (group != "Unknown")):
                    ex_equo.append(group)
                    annotation[cluster].append(group)
        if len(ex_equo) > 1:
            print ("\n%s is probably associated with %s and %s" %(cluster, most_rep_lab, ex_equo))
        else:
            print ("\n%s is probably associated with %s" %(cluster, most_rep_lab))
        print ("\tEntire dictionary: %s" %(assign_cluster))
        print ("\tKeywords mapped: %s\n" %(keywords_mapped))
    return annotation

#now that I have cluster --> groups, I can translate this into genes (per cluster) --> groups
def groups_to_genes(david, annotation):
    gene_annot = {}
    for cluster in david.keys():
        for term in david[cluster].keys():
            genes = david[cluster][term][1]
            genes = genes.rstrip().replace(",", "").split()
            group = annotation[cluster]
            for x in group:
                if x in gene_annot.keys():
                    for gene in genes:
                        if gene in gene_annot[x]:
                            pass
                        else:
                            gene_annot[x].append(gene)
                else:
                    gene_annot[x] = []
                    for gene in genes:
                        gene_annot[x].append(gene)

    for group in gene_annot.keys():
        print (group, ":", ",".join(gene_annot[group]))
    return gene_annot

#write output
def write_output(gene_annot):
    out = open("17_assigned_functions_DAVID.txt", "w")
    out.write("function\tgene\n")
    for go in gene_annot.keys():
        if go != "Unknown":
            out.write(go)
            out.write("\t")
            for gene in gene_annot[go]:
                out.write(gene)
                out.write(",")
            out.write("\n")
    out.close()
    return "\nOutput is produced. Name is '17_assigned_functions_DAVID.txt'"

###########MAIN

david = read_david_file()
annotation = assign_name(david)
gene_annot = groups_to_genes(david, annotation)
print (write_output(gene_annot))