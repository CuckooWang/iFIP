def getPPInum():
    f0 = open(r"UP000002311_559292.fasta", "r")
    dict_up = {}
    for line in f0.readlines():
        if line.startswith(">"):
            sp1 = line.split("|")
            sp2 = line.split()
            upid = sp1[1]
            for i in range(len(sp2)):
                if sp2[i].startswith("GN="):
                    sp3 = sp2[i].split("=")
                    gene = sp3[1]
                    dict_up[upid] = gene

    list_at = []
    f1 = open(r"thanatos\yeast.txt","r")
    for line in f1.readlines():
        sp = line.split("\t")
        upid = sp[0]
        if upid in dict_up:
            gene = dict_up[upid]
            if gene not in list_at:
                list_at.append(gene)

    list_atg9 = []
    f2 = open(r"PPI_ATG9.txt","r")
    for line in f2.readlines():
        sp = line.strip("\n").split("\t")
        gene1 = sp[0]
        gene2 = sp[1]
        if gene1 != "ATG9" and gene1 not in list_atg9:
            list_atg9.append(gene1)
        if gene2 != "ATG9" and gene2 not in list_atg9:
            list_atg9.append(gene2)

    dict_atnum = {}
    f3 = open(r"PPI_ALL.txt", "r")
    for line in f3.readlines():
        sp = line.strip("\n").split("\t")
        gene1 = sp[0]
        gene2 = sp[1]
        if gene1 in list_atg9 and gene2 in list_at:
            if gene1 not in dict_atnum:
                list_ppi = []
                list_ppi.append(gene2)
                dict_atnum[gene1] = list_ppi
            else:
                list_ppi = dict_atnum[gene1]
                if gene2 not in list_ppi:
                    list_ppi.append(gene2)
                    dict_atnum[gene1] = list_ppi
        if gene2 in list_atg9 and gene1 in list_at:
            if gene2 not in dict_atnum:
                list_ppi = []
                list_ppi.append(gene1)
                dict_atnum[gene2] = list_ppi
            else:
                list_ppi = dict_atnum[gene2]
                if gene1 not in list_ppi:
                    list_ppi.append(gene1)
                    dict_atnum[gene2] = list_ppi
    #remove Microautophagy gene
    list_at.remove("ATG20")
    list_at.remove("CLG1")
    list_at.remove("PHO80")
    list_at.remove("PEX15")

    fw = open(r"InteratedAutophagyGenes.txt","w")
    for i in range(len(list_atg9)):
        gene = list_atg9[i]
        yes = 0
        num = 0
        if gene in list_at:
            yes = 1
        if gene in dict_atnum:
            num = len(dict_atnum[gene])
        fw.write(gene + "\t" + str(yes) + "\t" + str(num) + "\n")

    fw.flush()
    fw.close()

getPPInum()
