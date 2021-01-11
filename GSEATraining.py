import copy
import random

import gseapy as gp
from LGTraining import predict

def perm():
    f0 = open("rawscore.txt", "r")
    dict_exp = {}
    for line in f0.readlines():
        sp = line.strip().split("\t")
        gene = sp[1]
        list_exp = []
        x1 = float(sp[3])
        x2 = float(sp[4])
        x3 = float(sp[5])
        list_exp.append(x1)
        list_exp.append(x2)
        list_exp.append(x3)
        dict_exp[gene] = list_exp
    list_coe = predict()
    max_es = ESscore(list_coe,dict_exp)
    print("Raw ES score:%s"%(max_es))

    step = 0.1
    num = 0
    max_coe = copy.deepcopy(list_coe)
    for i in range(1,10000):
        print(str(i))
        num = num + 1
        if num == 100:
            break
        tem_coe = mutation(max_coe,step)
        es = ESscore(tem_coe, dict_exp)
        if es > max_es:
            max_es = es
            print("New better ES:%s"%(max_es))
            max_coe = copy.deepcopy(tem_coe)
            f3 = open("BestES.txt","w")
            f3.write("Best ES:%s"%(max_es)+"\n")
            f3.write("Weights:" + "\t" + str(tem_coe[0]) + "\t" + str(tem_coe[1]) + "\t" + str(tem_coe[2]))
            num = 0
            f3.flush()
            f3.close()
        else:
            for j in range(1,20):
                tem_coe = mutation(tem_coe,step)
                es = ESscore(tem_coe, dict_exp)
                if es > max_es:
                    max_es = es
                    print("New better ES:%s" % (max_es))
                    max_coe = copy.deepcopy(tem_coe)
                    f3 = open("BestES.txt", "w")
                    f3.write("Best ES:%s" % (max_es) + "\n")
                    f3.write("Weights:" + "\t" + str(tem_coe[0]) + "\t" + str(tem_coe[1]) + "\t" + str(tem_coe[2]))
                    num = 0
                    f3.flush()
                    f3.close()
                    break
                else:
                    continue



def ESscore(list_coe,dict_exp):
    generateTem(list_coe,dict_exp)
    ss = gp.ssgsea(data="expression_tem.gct",
                   gene_sets="gene_sets.gmt",
                   outdir="GSEA_output",
                   sample_norm_method='rank',  # choose 'custom' for your own rank list
                   permutation_num=1000,  # skip permutation procedure, because you do need it
                   no_plot=True,  # skip plotting to speed up
                   processes=4, format='png')
    #每次循环都需要删除临时文件，因为GSEApy是个SB，会直接在文件的后面写入！！！！！

    es = readES()
    if es == -100.0:
        print("Wrong")
    return es

def generateTem(list_coe,dict_exp):
    #专门用来生成GSEA要得表达量数据的临时文件
    # 得分一直在变因此要建立临时文件
    fw = open("expression_tem.gct", "w")
    fw.write("#1.2" + "\n")
    fw.write("#201" + "\t" + "1" + "\n")
    fw.write("Gene_Symbol\tDescription\tScore\n")

    for gene in dict_exp:
        list_exp = dict_exp[gene]
        score = list_exp[0]*list_coe[0] + list_exp[1]*list_coe[1] + list_exp[2]*list_coe[2]
        fw.write(gene + "\t" + "na" + "\t" + str(score) + "\n")
    fw.flush()
    fw.close()

def readES():
    f0 = open("GSEA_output\Score\gseapy.ssgsea.gene_sets.report.txt","r")
    es = -100.0
    lines = f0.readlines()
    last_line = lines[-1]
    if last_line.startswith("ATG_Association"):
        sp = last_line.strip().split("\t")
        es = float(sp[1])
    else:
        print("Wrong")
    return es

def mutation(list_coe,step):
    pntem = random.randint(-1, 1)
    nstep = step
    if pntem < 0.0:
        nstep = 0 - step
    if pntem >= 0.0:
        nstep = step
    indtem = random.randint(0, 3)
    index = -1
    if indtem <= 1:
        index = 0
    elif indtem <= 2:
        index = 1
    elif indtem <= 3:
        index = 2
    else:
        print("Wrong")

    tem_coe = copy.deepcopy(list_coe)
    tem_coe[index] = tem_coe[index] + nstep
    return tem_coe
if __name__ == "__main__":
    perm()
