def Score():
    f0 = open(r"UP000002311_559292.fasta","r")
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
                    if gene not in dict_up:
                        list = []
                        list.append(upid)
                        dict_up[gene] = list
                    else:
                        list = dict_up[gene]
                        list.append(upid)
                        dict_up[gene] = list

    f1 = open(r"PPI_ATG9.txt","r")
    allatg9 = []
    for line in f1.readlines():
        sp = line.strip("\n").split("\t")
        gene1 = sp[0]
        gene2 = sp[1]
        if gene1 != "ATG9" and gene1 not in allatg9:
            allatg9.append(gene1)
        if gene2 != "ATG9" and gene2 not in allatg9:
            allatg9.append(gene2)

    trans_score = getTransScore()
    pro_score = getProScore()
    ppi_score = getPPIScore()
    fw = open("rawscore.txt","w")
    for i in range(len(allatg9)):
        gene = allatg9[i]
        all = 0.0
        trans = 1.0
        pro = 1.0
        ppi = 0.0
        if gene in trans_score:
            trans_list = trans_score[gene]
            if len(trans_list) != 3:
                print(gene + ":wrong")
            trans = (trans_list[0] + trans_list[1] + trans_list[2])/3.0
        if gene in pro_score:
            pro = pro_score[gene]
        if gene in ppi_score:
            ppi = ppi_score[gene]

        list_upid = dict_up[gene]
        fw.write(list_upid[0] + "\t" + gene + "\t" + str(all) + "\t" + str(trans) + "\t" + str(pro) + "\t" + str(ppi) + "\n")

    fw.flush()
    fw.close()

def getTransScore():
    cutoff = 0.05;
    f1 = open("gene_exp.diff","r")
    dict = {}
    for line in f1.readlines():
        if line.startswith("test_"):
            continue
        else:
            sp = line.split("\t")
            value1 = float(sp[7])
            value2 = float(sp[8])
            pval = float(sp[11])
            sp2 = sp[2].split(",")
            for i in range(len(sp2)):
                gene = sp2[i]
                tt1 = sp[4].split("_")
                tt2 = sp[5].split("_")
                if tt1[1] == tt2[1]:
                    fc = 1.0
                    if value1 != 0.0 and value2 != 0.0:
                        fc = value2/value1
                    if gene not in dict:
                        list_fcs = []
                        list_fcs.append(fc)
                        dict[gene] = list_fcs
                    else:
                        list_fcs = dict[gene]
                        list_fcs.append(fc)
                        dict[gene] = list_fcs
    return dict

def getProScore():
    f0 = open("UP000002311_559292.fasta", "r")
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
    f1 = open("Protein_PANDAout.txt")

    dict = {}
    for line in f1.readlines():
        if line.startswith("Protein"):
            continue
        else:
            sp = line.split("\t")
            sp0 = sp[0].split(";")
            for i in range(len(sp0)):
                if sp0[i].startswith("CON_") or sp0[i].startswith("REV_") or sp0[i] == "":
                    continue
                else:
                    upid = sp0[i]
                    gene = dict_up[upid]
                    atg9_0 = float(sp[1])
                    atg9_1 = float(sp[2])
                    atg9_4 = float(sp[3])
                    wt_0 = float(sp[4])
                    wt_1 = float(sp[5])
                    wt_4 = float(sp[6])

                    fc0 = atg9_0/wt_0
                    fc1 = atg9_1/wt_1
                    fc4 = atg9_4/wt_4
                    #print (str(fc0) + "\t" + str(fc1) + "\t" + str(fc2))
                    fc = (fc0+fc1+fc4)/3.0
                    #print(gene + "\t" + str(fc))
                    dict[gene] = fc
    return dict

def getPPIScore():
    f1 = open("InteratedAutophagyGenes.txt","r")
    dict = {}
    list_ppis = []
    for line in f1.readlines():
        sp = line.strip().split("\t")
        gene = sp[0]
        score = float(sp[2])
        list_ppis.append(score)
        dict[gene] = score

    max_num = max(list_ppis)
    for key in dict:
        dict[key] = dict[key]/max_num
    return dict
Score()
