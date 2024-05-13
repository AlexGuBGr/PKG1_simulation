# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import mygene
import subprocess
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
mg = mygene.MyGeneInfo()


rng = np.random.default_rng(seed=1234)
coldict = {'E2F1':"red", 'JUN':"blue", 'NFKB1':"orange", 'RELA':"purple", 'SP3':"pink", 'STAT1':"teal", 'STAT3':"grey"}
scale = [0.03, 0.05, 0.1, 0.3, 0.5, 1]


def get_phosphoplus_data(path = "PSP-SubstrateSearch_PRKG1.xlsx"):

    # retrieving data from PhosphositePlus
    tab = pd.read_excel(path, skiprows=19)
    tab = tab[(tab['ORGANISM'] == "human")]
    prkg1_targets = tab[["GENE", "PROTEIN", "ACC#"]]
    return prkg1_targets


def get_regphos_data(path = "RegPhos_kinase_PPI_human.txt"):
    # retrieving data from RegPhos
    tab1 = pd.read_csv("RegPhos_kinase_PPI_human.txt", sep=None, engine="python")
    prkg1_targets1 = tab1[['﻿GENE_a', 'GENE_b','AC_a', 'AC_b']][tab1['﻿GENE_a'] == "PRKG1"]
    return prkg1_targets1


def ac_symbol_matcher(prkg1_targets, prkg1_targets1):
    # getting alle genes and 'core' protein accession ids from PhosphositePlus
    gene_to_prot = {}
    for i in prkg1_targets[["GENE", "PROTEIN", "ACC#"]].values:
        if i[2].split("-")[0] in gene_to_prot:
            gene_to_prot[i[2].split("-")[0]].add(i[0])
        else:
            gene_to_prot[i[2].split("-")[0]] = {i[0]}
        
    # getting all genes and 'core' protein accesion ids from RegPhos
    newones = {}
    for i in prkg1_targets1[['﻿GENE_a', 'GENE_b','AC_a', 'AC_b']].values:
        if i[2].split("-")[0] in gene_to_prot:
            gene_to_prot[i[2].split("-")[0]].add(i[0])
        else:
            if i[2].split("-")[0] in newones:
                newones[i[2].split("-")[0]].add(i[0])
            else:
                newones[i[2].split("-")[0]] = {i[0]}
    
        if i[3].split("-")[0] in gene_to_prot:
            gene_to_prot[i[3].split("-")[0]].add(i[1])
        else:
            newones[i[3].split("-")[0]] = i[1]
            if i[2].split("-")[0] in newones:
                newones[i[3].split("-")[0]].add(i[1])
            else:
                newones[i[3].split("-")[0]] = {i[1]}
                
    # checking if any ids can be mapped to multiple symbols
    for k,v in newones.items():
        if len(v) > 1:
            print(k,v)
    
    for k,v in gene_to_prot.items():
        if len(v) > 1:
            print(k,v)
    
    # checking if the additional symbols from RegPhos can be mapped to other ac_ids
    vals = set()
    for i in gene_to_prot.values():
        vals.update(i)
    
    for k,v in newones.items():
        for i in v:
            if i in vals:
                print(k,v)

    # if nothing gets printed.. all is good.



def get_unique_pkg1_targets(path1 = "PSP-SubstrateSearch_PRKG1.xlsx", path2 = "RegPhos_kinase_PPI_human.txt"):

    # retrieving unique PRKG1 substrates PhosphositePlus and RegPhos    

    prkg1_targets = get_phosphoplus_data(path1)
    prkg1_targets1 = get_regphos_data(path2)

    ac_symbol_matcher(prkg1_targets, prkg1_targets1)

    unique_pkg1_targets = set(prkg1_targets["GENE"])
    unique_pkg1_targets.update(set(prkg1_targets1['GENE_b']))

    return unique_pkg1_targets


def make_grns(grn1 = 'trrust_rawdata.human.tsv', grn2 = "RegulatoryDirections/new_kegg.human.reg.direction.txt"):
    
    ## extracting data from the TRUST gene regulatory network
    trust = pd.read_csv(grn1, sep=None, engine="python", header=None)
    
    regnet = pd.read_csv(grn2, sep=None, engine="python")
    regnet["directions"] = list(regnet["Up"])
    regnet = regnet[['#TF', 'ID', 'Target', 'ID.1','directions']]
    
    
    # combining the known edge data, removing edges with ambigous regulatory effects
    use_edges = {}
    for i in trust.values:    
        if i[2] == 'Activation':
            if (i[0], i[1]) in use_edges:
                use_edges[(i[0], i[1])] += 1
            else:
                use_edges[(i[0], i[1])] = 1
        elif i[2] == 'Repression':
            if (i[0], i[1]) in use_edges:
                use_edges[(i[0], i[1])] -= 1
            else:
                use_edges[(i[0], i[1])] = -1
    
    
    for i in regnet.values:
        if i[4] == '-->':
            if (i[0], i[2]) in use_edges:
                use_edges[(i[0], i[2])] += 1
            else:
                use_edges[(i[0], i[2])] = 1
        elif i[4] == '--|':
            if (i[0], i[2]) in use_edges:
                use_edges[(i[0], i[2])] -= 1
            else:
                use_edges[(i[0], i[2])] = -1
    
    
    #max_weight = max(np.abs(list(use_edges.values())))
    G = nx.DiGraph()
    
    for k,v in use_edges.items():
        if v > 0:
            G.add_edge(k[0], k[1], color = "green")#, weight = use_edges[k]/max_weight)
        if v < 0:
            G.add_edge(k[0], k[1], color = "red")#, weight = use_edges[k]/max_weight)

    print("Size of initial gene regulatory network:", len(G))
    return G



def get_nox_tfs(graph):
    # getting regulatory edges connected to NOX4 and NOX5
    out = {}
    for i in ["NOX4", "NOX5"]:
        out[i] = [ii[0] for ii in graph.in_edges(i)]
        
    return out
    

def draw_network(g, name):
    # network drawing function - relies on 'graphviz version 9.0.0'
    nx.drawing.nx_pydot.write_dot(g, 'my_dot.dot')
    suc = subprocess.run(["dot", "-Tsvg", 'my_dot.dot', "-o", "results/{}.svg".format(name)])
    print(suc)

def generate_pkg_to_nox_img(G, g, pkgTarg):

    # hardcoded path from PRKG1 to some NOX-TFs and NOX4/5
    # relies on 'graphviz version 9.0.0'    

    paths_for_img = [['MAPK14', 'IRF4', 'BCL6', 'NFKB1', 'AR', 'RELA'],
                    ["JUN", 'IRF1', 'E2F1']]

    to_subg = set(sum(paths_for_img, []) + ["PRKG1", "NOX4", "NOX5"])

    Gtmp = G.copy()
    for i in g.edges():
        col = g.get_edge_data(i[0], i[1])["color"]
        Gtmp.add_edge(i[0], i[1], color = col)

    for i in pkgTarg:
        Gtmp.add_edge("PRKG1", i, color = "purple")

    gtmp = nx.subgraph(Gtmp, to_subg).copy()
    
    ps = [['PRKG1'], ['MAPK14'], ['IRF4'], ['BCL6'], ['NFKB1'], ['NOX4', 'IRF1', 'AR', 'NOX5'], ['E2F1','RELA','JUN']]

    x = {}
    y = {}
    ln = 4 / 2 + 1
    for j,i in enumerate(ps):
        ln1 = len(i) / 2
        for jj,ii in enumerate(i):
            y[ii] = len(ps) - j
            x[ii] = ln + jj - ln1
    
    scle = 1
    gtmp1 = nx.DiGraph() 
    for k,v in x.items():
        if k == "NFKB1":
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5+1) * scle, y[k]*scle))
        elif k == "BCL6":
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5-1)*scle, y[k]*scle))
        elif k == "AR":
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5+0.35)*scle, y[k]*scle))
        elif k == "NOX5":
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5+0.5)*scle, y[k]))
        elif k == "IRF1":
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5-0.25)*scle, y[k]*scle))
        elif k == "JUN":
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5+1.25)*scle, y[k]*scle))
        else:
            gtmp1.add_node(k, pos = '{},{}!'.format((v*1.5)*scle, y[k]*scle))
    
    for i in gtmp.edges():
        dct = gtmp.get_edge_data(i[0], i[1])
        for k,v in dct.items():
            gtmp1.add_edge(i[0], i[1], color = v)
    nx.drawing.nx_pydot.write_dot(gtmp1, 'my_dot.dot')
    suc = subprocess.run(["dot", "-K", "neato", "-Tsvg", 'my_dot.dot', "-o", "results/{}.svg".format("PRKG1_to_NOX")])

    #return x,y, gtmp1
    
#x,y, bg = generate_pkg_to_nox_img(G, g, pkgTarg)



def prints_missing_tfs(G, symb_to_entrez):
    # prints TFs that could not be identified in the MyGene search
    for i in symb_to_entrez["missing"]:
        od = G.out_edges(i)
        if len(od) > 0:
            print(i)
            



def remove_old_symbols(g1, new_names):
    
    # removing deprecated symbols
    # All old edges of all old symbols are added to the new ones
    # the old symbols are deleted
    
    for k,v in new_names.items():
    # checking IF the current symbol already is present in the GRN..
    # if so, the new symbol inherets everything and the old is removed
        if v in g1:
            old = g1.out_edges(k)
            new = g1.out_edges(v)
            for i in old:
                if i not in new:
                    col = g1.get_edge_data(i[0],i[1])["color"]
                    g1.add_edge(v,i[1], color = col)
            assert len(set([s[1] for s in g1.out_edges(k)]).difference(set([s[1] for s in g1.out_edges(v)]))) == 0
        else:
            old = g1.out_edges(k)
            for i in old:
                col = g1.get_edge_data(i[0],i[1])["color"]
                g1.add_edge(v,i[1], color = col)
    
    g1.remove_nodes_from(list(new_names.keys()))
    
    return g1


def test_connectivity(G, pkgTarg, noxtfs_set, notinG):
    
    # Finds paths from all PRKG1 substrates to NOX-TFs
    # Outputs a list with PRKG1 substrates found in the GRN
    
    ispath = set()
    notpath = set()
    Gu = nx.to_undirected(G)
    for i in pkgTarg:
        for ii in noxtfs_set:
            try:
                tmp = nx.shortest_path(Gu, i, ii)
                ispath.add((i,ii))
            except:
                notpath.add((i,ii))
    
    dct1 = {"PKG1 target":[], "NOX4_5 TF" : []}
    dct2 = {"PKG1 target":[], "NOX4_5 TF" : []}
    dct3 = {"PKG1 target":[]}
    
    for i in sorted(list(ispath)):
        dct1["PKG1 target"].append(i[0])
        dct1["NOX4_5 TF"].append(i[1])

    for i in sorted(list(notpath)):
        dct2["PKG1 target"].append(i[0])
        dct2["NOX4_5 TF"].append(i[1])
        
    for i in notinG:
        dct3["PKG1 target"].append(i)

    with pd.ExcelWriter("results/connections_between_PKG1_substrates_and_NOX45_TFs.xlsx") as writer:
        pd.DataFrame(dct1).to_excel(writer, sheet_name="PKG1_sub_conntected", index=False)
        pd.DataFrame(dct2).to_excel(writer, sheet_name="PKG1_sub_NOT_conntected", index=False)
        pd.DataFrame(dct3).to_excel(writer, sheet_name="PKG1_sub_NOT_in_graph", index=False)


    return set([i[0] for i in ispath])



def finding_all_shortest_paths(G, pkgTarg, noxtfs):
    
    # finding all shortest paths from the PKG targets to the NOX-TFs
    no_paths = []
    pths = []
    for i in pkgTarg:
        for ii in set(sum(list(noxtfs.values()), [])):
            try:
                pths += list(nx.all_shortest_paths(G, i, ii)) 
            except:
                no_paths.append([i,ii])
    return pths, no_paths


def extracting_relevant_cc(G, pths):
    
    # Only connected components containing NOX-Tfs and PRKG1-targets are saved
    core_genes = set(sum(pths, []))
    to_remove = set()
    cc_index = []
    connect_cc = list(nx.connected_components(nx.to_undirected(G)))
    number_rem = 0
    for j,i in enumerate(connect_cc):
        if len(core_genes.intersection(i)) == 0:
            to_remove.update(i)
            number_rem += 1
        else:
            cc_index.append([j, len(i)])
    G.remove_nodes_from(to_remove)
    
    print("\nremoved {} out of {} connected components".format(number_rem, len(connect_cc)))
    print("\nSize of GRN without 'un'connected components:", len(G), "\n")
    
    return G



def remake_entrez_to_symb(ets, set_G):
    # re-compiles the entrezID-to-symbol dictionary
    out = {}
    for k,v in ets.items():
        if v in set_G:
            out[k] = v
    return out



def extracting_gene_values(G, entrez_to_symb, path = "annotated_ctrl_exprs_stacked.tsv/annotated_ctrl_exprs_stacked.tsv"):
    ## extracting expression values to create an estimate of prot abundance
    cntr_dat = pd.read_csv(path, sep="\t")
    
    # find nodes that need to be imputed
    to_impute = {}
    for i in set(list(entrez_to_symb.keys())).difference(cntr_dat["EntrezID"][cntr_dat["EntrezID"].isin(list(entrez_to_symb.keys()))]):
        if len(G.out_edges(entrez_to_symb[i])) > 0:
            to_impute[i] = entrez_to_symb[i]
    
    imputed_dct = {"EntrezID":[], "Symbol":[]}
    for k,v in to_impute.items():
        imputed_dct["EntrezID"].append(k)
        imputed_dct["Symbol"].append(v)
    pd.DataFrame(imputed_dct).to_excel("results/Imputed genes.xlsx", index=False)
    
    
    cntr_dat = cntr_dat[cntr_dat["EntrezID"].isin(list(entrez_to_symb.keys()))]
    coldat = cntr_dat.columns
    
    # mean-centering and max-scaling the data
    rep_size = int(cntr_dat.shape[0]/5)
    cntr_dat = cntr_dat.values
    entrez_in_dat = cntr_dat[:rep_size,0]
    cntr_dat= cntr_dat[:,1:].reshape((5,rep_size, 10))
    cntr_dat = cntr_dat/np.max(cntr_dat)
    cntr_dat -= np.mean(cntr_dat)
    
    mean = np.mean(cntr_dat)
    std = np.std(cntr_dat)
    
    # Imputing values of missing nodes
    imputed, entrez = impute_missing_tfs(to_impute, (5,1,10), mean, std, rng)
    cntr_dat = np.concatenate([cntr_dat, np.concatenate(imputed, axis=1)], axis=1)
    entrez_in_dat = np.concatenate([entrez_in_dat, entrez])
    
    return cntr_dat, entrez_in_dat, coldat



def impute_missing_tfs(genes, shape_for_one, mean, std, rng):

    # Imputation function
    to_add = []
    to_add_entrez = []
    for k,v in genes.items():
        dist = rng.normal(mean, std, size=shape_for_one)
        dist[dist < 0] = std
        to_add.append(dist)
        to_add_entrez.append(k)

    return to_add, to_add_entrez




def scaled_tanh(x,y):
    # A scaled tanh function - NOT USED
    out = (np.exp(x+y)-np.exp(-x+y))/(np.exp(x)+np.exp(-x))
    if np.isnan(out):
        return 1
    else:
        return out

def scaled_sigmoid(x,y):
    # the scaled sigmoid function
    return 1/(1+np.exp(-y*x))

def calc_regulation(G, start_ab, symb_to_index, pkgT, noxtfs_set, max_path_ln, rate = 0, factor = 0.3):

    """
    starting from PRKG1 substrates (pkgT), this function updates the abundance-values
    of nodes in the GRN a number of times equal to the length of  longest shortest path - 1.
    
    the function outputs:
        true_abundance:     abundance-values for NOX-TFs after ended simulation run
        change_per_step:    abundance-value-changes for all nodes per time-step

    """
    start = 1

    change_per_step = {}
    for i in noxtfs_set:
        tindx = symb_to_index[i]
        change_per_step[tindx] = np.zeros((max_path_ln)) + start_ab[tindx]

    abundance = start_ab.copy()
    tmp_targ = pkgT.copy()

    self_loops = set()
    while start < max_path_ln:
        step = {}
        new_targ = set()
        for i in tmp_targ:
            tmp_nb = list(G.out_edges(i))  
            for nb in tmp_nb:
                if nb[0] == nb[1]:
                    
                    if (nb[0], nb[1]) in self_loops:
                        #self_loops.add((nb[0], nb[1]))
                        #self_loops[(nb[0], nb[1])] = 0
                        continue
                    else:
                        self_loops.add( (nb[0], nb[1]) )
                    #    #if self_loops[(nb[0], nb[1])] == 0:
                    #    self_loops[(nb[0], nb[1])] += 1
                        #else:
                        #    self_loops[(nb[0], nb[1])] = 1#self_loops[(nb[0], nb[1])]**(start-1) 
                            #continue
                #else:
                new_targ.add(nb[1])
                tmp_val = G.get_edge_data(nb[0], nb[1])["color"]
                
                if tmp_val == "green":
                    if symb_to_index[nb[1]] in step:
                        step[symb_to_index[nb[1]]].append((scaled_sigmoid(abundance[symb_to_index[i]],factor) ))
                    else:
                        step[symb_to_index[nb[1]]] = [(scaled_sigmoid(abundance[symb_to_index[i]],factor))]
                    
                elif tmp_val == "red":
                    if symb_to_index[nb[1]] in step:
                        step[symb_to_index[nb[1]]].append(-(scaled_sigmoid(abundance[symb_to_index[i]], factor) ))
                    else:
                        step[symb_to_index[nb[1]]] = [-(scaled_sigmoid(abundance[symb_to_index[i]],factor) ) ]

        for k,v in step.items():
            abundance[k] += np.sum(v)
            if k in change_per_step:
                change_per_step[k][start] = abundance[k]
        for k in set(change_per_step).difference(step):
            change_per_step[k][start] = change_per_step[k][start - 1]

        start += 1
        tmp_targ = list(new_targ)

    #tf_reg = []
    true_abundance = []
    for i in noxtfs_set:
        #tf_reg.append([i, regu])
        true_abundance.append([i, abundance[ symb_to_index[i] ], start_ab[ symb_to_index[i] ]])

    return true_abundance, change_per_step 


def compile_means(noxtfs_set, scale, shape_, dct, val = 1):
    
    # Producing smaller dicts with averaged abundance-values
    av = {}
    for i in noxtfs_set:
        av[i] = 0


    av_per_scale = {}
    for i in scale:
        av_per_scale[i] = av.copy()


    for k,v in dct.items():
        for kk,vv in v.items():
            for i in vv:
                av[i[0]] += i[val]
                av_per_scale[k][i[0]] += i[val]

    for k,v in av.items():
        av[k] = v / (shape_[0] * shape_[2] * len(scale))


    for k,v in av_per_scale.items():
        for kk, vv in v.items():
            av_per_scale[k][kk] = vv / (shape_[0] * shape_[2])


    return av, av_per_scale



def compute_simulation(G, cntr_dat, symb_to_index, pkgTarg, noxtfs_set, max_path_ln, scale):

    """
    Function for running all simulations.
    Runs a simulation for at varying scales for ever column in all the 
    expression datasets (shape: 5, 565, 10)
    """    

    shape_ = cntr_dat.shape

    out_regu = {}
    out_true_ab = {}
    out_chst = {}
    for i in scale:
        print("scale ", i)
        out_regu[i] = {}
        out_true_ab[i] = {}
        out_chst[i] = {}
        for ii in range(shape_[0]):
            #cr = np.corrcoef(cntr_dat[ii])
            for iii in range(shape_[2]):
                tab, chst = calc_regulation(G, cntr_dat[ii][:,iii], symb_to_index, pkgTarg, noxtfs_set, max_path_ln, rate = 0, factor = i)
                #out_regu[i][(ii,iii)] = reg
                out_true_ab[i][(ii,iii)] = tab
                out_chst[i][(ii,iii)] = chst

    av_pred, av_pred_per_scale = compile_means(noxtfs_set, scale, shape_, out_true_ab, val = 1)
    av_tr, av_tr_per_scale = compile_means(noxtfs_set, scale, shape_, out_true_ab, val = 2)


    return out_true_ab, av_tr, av_tr_per_scale, av_pred, av_pred_per_scale, out_chst




##########################################



def plot_mean_expression_change(noxtfs_set, symb_to_index, max_path_ln, out_chst, entrez_to_symb, coldict, name):
    
    # plots mean expression of NOX-TFs per timestep
    
    ze = {}
    for i in sorted(list(noxtfs_set)):
        tindx = symb_to_index[i]
        ze[tindx] = np.zeros((max_path_ln))
    
    for k,v in out_chst.items():
        for kk,vv in v.items():
            for kkk,vvv in vv.items():
                ze[kkk] += vvv
    
    amount = len(out_chst) * len(out_chst[1])
    
    plt.figure(figsize=(17,10))
    ax = plt.subplot(111)
    for k,v in ze.items():
        n = entrez_to_symb[entrez_in_dat[k]]
        ax.plot(v/amount, label=n, color = coldict[n], lw=3)
    
    
    plt.legend(fontsize=20)
    plt.title("Mean Expression Across Time Steps",fontsize=24)
    plt.ylabel("Mean Expression", fontsize=20)
    plt.xlabel("Simulation Time Steps", fontsize=20)
    my_ticks = list(range( max_path_ln ))
    plt.xticks(my_ticks, [i+1 for i in my_ticks])
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.1, 1.), shadow=True, ncol=1, fontsize=15)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.savefig("results/exp_change_per_timestep_{}.svg".format(name))


##########################################
    

def plot_scaled_sigmoid(s = [0.03, 0.05, 0.1, 0.3, 0.5, 1]):
    
    # plots the curved of the scaled sigmoid function given differen scale values
    
    x = np.arange(-30,30,0.1)
    plt.figure(figsize=(17,10))
    for i in range(len(s)):
        plt.plot(x,scaled_sigmoid(x,s[i]), label="z = {}".format(np.round(s[i],2)), lw=3)
        plt.legend(fontsize=20)

    plt.title("Regulation Magnitudes",fontsize=24)
    plt.xlabel("Transcription Factor Abundance",fontsize=20)
    plt.ylabel("Regulation Size",fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.savefig("results/scaled_regulation.svg") 



def expression_change(true_exp, pred_exp, noxtfs, name, scale):
    
    # Computes and writes expression of NOX4/5 before and after simulations
    
    n4n5 = list(set(noxtfs["NOX4"]).intersection(noxtfs["NOX5"]))
    nx4only = list(set(noxtfs["NOX4"]).difference(n4n5 + noxtfs["NOX5"]))
    nx5only = list(set(noxtfs["NOX5"]).difference(n4n5 + noxtfs["NOX4"]))
    
    ntfs = dict(both = n4n5, NOX4 = nx4only, NOX5 = nx5only)

    nx4nx5 = {}
    for j,i in enumerate(scale):
        nx4 = []
        nx5 = []
        bth = []
        for ii in sorted(ntfs["NOX4"]):
            t = scaled_sigmoid(true_exp[ii][j],i)
            p = scaled_sigmoid(pred_exp[ii][j],i)
            nx4 += [t,p,p-t]
        
        for ii in sorted(ntfs["both"]):
            t = scaled_sigmoid(true_exp[ii][j],i)
            p = scaled_sigmoid(pred_exp[ii][j],i)
            bth += [t,p,p-t]
        
        for ii in sorted(ntfs["NOX5"]):
            t = scaled_sigmoid(true_exp[ii][j],i)
            p = scaled_sigmoid(pred_exp[ii][j],i)
            nx5 += [t,p,p-t]
        nx4nx5[i] = nx4 + bth + nx5
    
    df_multindex = pd.DataFrame(nx4nx5,
                            index=[["NOX4"] * len(ntfs["NOX4"]) * 3 + 
                                    ["NOX4_and_NOX5"] * len(ntfs["both"]) * 3 + 
                                    ["NOX5"] * len(ntfs["NOX5"]) * 3,
                                    sorted(ntfs["NOX4"]*3) + sorted(ntfs["both"] * 3) + 
                                    sorted(ntfs["NOX5"]*3), 
                                    ["start", "end", "difference"]*len(ntfs["NOX4"] + 
                                    ntfs["NOX5"] + 
                                    ntfs["both"])])        
    dfmT = df_multindex.T
    dfmT.index.name = "scale"
    
    n4c = len(set([i[0] for i in dfmT["NOX4"].columns]))
    n5c = len(set([i[0] for i in dfmT["NOX5"].columns]))
    bc = len(set([i[0] for i in dfmT["NOX4_and_NOX5"].columns]))
    nx4s = np.sum(dfmT["NOX4"].values.reshape((-1,n4c,3)),axis=1) + np.sum(dfmT["NOX4_and_NOX5"].values.reshape((-1,bc,3)),axis=1)
    nx5s = np.sum(dfmT["NOX5"].values.reshape((-1,n5c,3)),axis=1) + np.sum(dfmT["NOX4_and_NOX5"].values.reshape((-1,bc,3)),axis=1)

    summed_reg = pd.DataFrame(np.concatenate([nx4s, nx5s], axis=1), columns=list(zip(["NOX4"]*3 + ["NOX5"]*3 ,["start", "end", "difference"]*2)))
    summed_reg.columns = pd.MultiIndex.from_tuples(summed_reg.columns)
    summed_reg.index = scale
    summed_reg.index.name = "scale"
    
    dfmT.loc['mean'] = np.mean(dfmT.values,axis=0)    
    summed_reg.loc['mean'] = np.mean(summed_reg.values,axis=0)

    with pd.ExcelWriter("results/simulated_NOX45_expression_data_{}.xlsx".format(name)) as writer:
        dfmT.to_excel(writer, sheet_name="TF_effect_per_scale")
        summed_reg.to_excel(writer, sheet_name="summed_TF_effect_per_scale")


    #return dfmT, summed_reg
#expression_change(true_exp, pred_exp, noxtfs, "orig")



def write_data_to_excel(noxtfs_set, av_pred_per_scale, av_tr_per_scale, name, scale):
    
    # computes and writes the abundance values of the NOX-TFs upon ended simulation
    
    overall_exp_change = {"scale" : []}
    true_exp = {"scale" : []}
    pred_exp = {"scale" : []}
    
    for i in sorted(list(noxtfs_set)):
        overall_exp_change[i] = []
        true_exp[i] = []
        pred_exp[i] = []
        
    for k,v in av_pred_per_scale.items():
        overall_exp_change["scale"].append(k)
        true_exp["scale"].append(k)
        pred_exp["scale"].append(k)
        for kk,vv in v.items():
            
            overall_exp_change[kk].append(vv - av_tr_per_scale[k][kk])
            true_exp[kk].append(av_tr_per_scale[k][kk])
            pred_exp[kk].append(vv)
    
    mean_overall_exp_change = overall_exp_change.copy()
    mean_overall_exp_change.pop("scale")
    for k,v in mean_overall_exp_change.items():
        mean_overall_exp_change[k] = [np.mean(v)]
        
    
    with pd.ExcelWriter("results/simulated_expression_data_{}.xlsx".format(name)) as writer:
        pd.DataFrame(mean_overall_exp_change).to_excel(writer, sheet_name="total_mean_absolute_exp_diff", index=False)
        pd.DataFrame(overall_exp_change).to_excel(writer, sheet_name="mean_absolute_exp_diff", index=False)
        pd.DataFrame(true_exp).to_excel(writer, sheet_name="mean_true_exp_per_scale", index=False)
        pd.DataFrame(pred_exp).to_excel(writer, sheet_name="mean_pred_exp_per_scale", index=False)


    expression_change(true_exp, pred_exp, noxtfs, name, scale)


def load_data_and_make_graph(path1 = "PSP-SubstrateSearch_PRKG1.xlsx", path2 = "RegPhos_kinase_PPI_human.txt", grn1 = 'trrust_rawdata.human.tsv', grn2 = "RegulatoryDirections/new_kegg.human.reg.direction.txt"):
    
    # a 'master function' that loads all data and creates the GRN
    
    unique_pkg1_targets = get_unique_pkg1_targets(path1, path2)

    G = make_grns(grn1, grn2)
    
    noxtfs = get_nox_tfs(G)
    noxtfs_set = set(sum(list(noxtfs.values()), []))
    g = nx.subgraph(G, list(set(sum(list(noxtfs.values()),[]))) + list(noxtfs.keys())).copy()
    draw_network(g, "nox_and_tfs")
    
    symb_to_entrez = mg.querymany(list(G), scopes='symbol', species=9606, fields="entrezgene,uniprot,symbol", returnall=True)


    # Some genes (symbols) could not be identified.
    # The ones that are transcription factors are shown below
    #prints_missing_tfs(G, symb_to_entrez)

    #AES
    #ARNTL
    #ARNTL2
    #MKL1
    #RFWD2
    #SALL4A
    #STAT5
    #ZNRD1
    #U2AF1L5

    # the current symbols for these genes where ideintified by manual curation
    new_names = {"AES": "TLE5", "ARNTL":"BMAL1", "ARNTL2":"BMAL2", 
                 "MKL1":"MRTFA", "RFWD2":"COP1", "SALL4A":"SALL4", 
                 "STAT5":"STAT5A", "ZNRD1":"POLR1H", 
                 "U2AF1L5":"LOC102724594"}
    
    print("\nThe following symbols could not be found via MyGene, but have been replaced by manual curation:")
    print("old:\t\t\tnew:")
    for k,v in new_names.items():
        print(k,"\t\t\t",v)
    
    
    #SALL4A is an isoform of SALL4
    nn_dct = {"old names":[], "new names":[]}
    for k,v in new_names.items():
        nn_dct["old names"].append(k)
        nn_dct["new names"].append(v)    
    pd.DataFrame(nn_dct).to_excel("results/Gene symbols replaced.xlsx", index = False)

    G = remove_old_symbols(G, new_names)

    # getting entrezID for the new symbols
    symb_to_entrez_addon = mg.querymany(list(new_names.values()), scopes='symbol', species=9606, fields="entrezgene,uniprot,symbol", returnall=True)

    # adding it together
    for k,v in symb_to_entrez_addon.items():
        symb_to_entrez[k] += v

    # making a entrezID to symbol dictionary
    entrez_to_symb = {}
    for i in symb_to_entrez["out"]:
        if "entrezgene" in i:
            entrez_to_symb[int(i["entrezgene"])] = i["query"]
    
    
    # here I remove non-tf nodes
    to_remove = set()
    for i in list(G.out_degree(list(G))):
        #if i[0] == "NOX4" or i[0] == "NOX5":
        #    continue
        if i[1] == 0:
            to_remove.add(i[0])
    
    G.remove_nodes_from(to_remove)
    print("\nSize of gene regulatory network only with tfs:", len(G))
    
    # getting the PKG1 targets available in the GRN
    pkgTarg = unique_pkg1_targets.intersection(set(G))
    
    # finding PKG1 substrates connected to any NOX4/5 TFs
    pkgTarg = test_connectivity(G, pkgTarg, noxtfs_set, unique_pkg1_targets.difference(pkgTarg))


    pths, no_paths = finding_all_shortest_paths(G, pkgTarg, noxtfs)
    max_path_ln = len(max(pths,key=len))
    
    # Only connected components containing NOX-Tfs and PKG1-targets are saved
    G = extracting_relevant_cc(G, pths)


    entrez_to_symb = remake_entrez_to_symb(entrez_to_symb, set(G))

    assert len(entrez_to_symb) == len(G)

    return G, pkgTarg, entrez_to_symb, pths, max_path_ln, noxtfs, g, unique_pkg1_targets, noxtfs_set


def shufler(a, sq = False):
    # Simple shuffle-function
    aa = a.flatten()
    
    if sq == True:
        aa = np.square(aa)
    
    rng.shuffle(aa)
    aa = aa.reshape(a.shape)
    return aa


def random_sampler(size, rng):
    # simple function for random sampling
    return rng.normal(size=size)


def get_tf_w_pths_to_noxtfs(G, noxtfs_set, pkgTarg, rng):
    
    
    lg = list(set(G).difference(noxtfs_set.union(pkgTarg)))
    pths1, no_paths1 = finding_all_shortest_paths(G, lg, noxtfs)
    
    nxtfs = {}
    tf = ""
    collected = set()
    for j,i in enumerate(pths1):
        if j == 0:
            tf = i[0]
            nxtfs[tf] = []

        if i[0] == tf:
            nxtfs[tf].append(i)
            collected.add(i[-1])

        else:
            if collected != noxtfs_set or len(set(sum(nxtfs[tf],[])).intersection(pkgTarg)) > 0:
                nxtfs.pop(tf)
            
            tf = i[0]
            nxtfs[tf] = [i]
            collected = set()
            collected.add(i[-1])

    lnxtfs = sorted(list(nxtfs))
    rng.shuffle(lnxtfs)

    return nxtfs, lnxtfs

def append_double_dict(inst, add):
    for k in inst:
        for kk in inst[k]:
            inst[k][kk].append( add[k][kk] )
    return inst
    
def append_triple_dict(inst, add):
    for k in inst:
        for kk in inst[k]:
            for kkk in inst[k][kk]:
                inst[k][kk][kkk].append( add[k][kk][kkk] )
    return inst

def make_like3(use):
    out = {}
    for k in use:
        out[k] = {}
        for kk in use[k]:
            out[k][kk] = {}
            for kkk in use[k][kk]:
                out[k][kk][kkk] = []

    return out

def make_like2(use):
    out = {}
    for k in use:
        out[k] = {}
        for kk in use[k]:
            out[k][kk] = []

    return out



def combine_runs(oc, times):
    b = {}
    for k in oc:
        b[k] = {}
        for kk in oc[k]:
            b[k][kk] = {}
            for kkk in oc[k][kk]:
                zz = np.zeros((len(max(oc[k][kk][kkk],key=len))))
                zzc = zz.copy()
                for i in oc[k][kk][kkk]:
                    zz[:len(i)] += i
                    zzc[:len(i)] += 1
                b[k][kk][kkk] = zz/times #* (zzc / np.max(zzc))
    return b


def run_multiple_times(G, use_dat, symb_to_index, pkgTarg, noxtfs_set, scale, nxtfs, lnxtfs, max_path_ln):
    
    lnpkgt = len(pkgTarg)
    #times = len(lnxtfs) - lnpkgt + 1
    times = len(lnxtfs) // lnpkgt
    all_stfs = []
    for i in range(times):
        stf = lnxtfs[lnpkgt*i:lnpkgt*(i+1)]
        print(i+1, stf)
        all_stfs.append(stf)
        start = sum([nxtfs[k] for k in stf],[])
        #max_path_ln = len(max(start, key=len))
        gg = nx.subgraph(G, set(sum(start,[])))
        out_true_ab, av_tr, av_tr_per_scale, av_pred, av_pred_per_scale, out_chst = compute_simulation(gg, use_dat, symb_to_index, stf, noxtfs_set, max_path_ln, scale)

        if i == 0:
            oc = make_like3(out_chst)
            atps = make_like2(av_tr_per_scale)
            apps = make_like2(av_pred_per_scale)
            
        oc = append_triple_dict(oc, out_chst)
        atps = append_double_dict(atps, av_tr_per_scale)
        apps = append_double_dict(apps, av_pred_per_scale)
        
    return atps, apps, oc, combine_runs(oc, times)


def diff_chst(noxtfs_set, symb_to_index, max_path_ln, out_chst, oc_comb, entrez_to_symb, coldict, name):
    ze = {}
    ze1 = {}
    for i in sorted(list(noxtfs_set)):
        tindx = symb_to_index[i]
        ze[tindx] = np.zeros((max_path_ln))
        ze1[tindx] = np.zeros((max_path_ln))
    
    for k,v in out_chst.items():
        for kk,vv in v.items():
            for kkk,vvv in vv.items():
                ze[kkk] += vvv
                ze1[kkk] += oc_comb[k][kk][kkk]
    
    amount = len(out_chst) * len(out_chst[1])
    
    plt.figure(figsize=(17,10))
    ax = plt.subplot(111)
    for k,v in ze.items():
        n = entrez_to_symb[entrez_in_dat[k]]
        ax.plot(v/amount - ze1[k]/amount, label=n, color = coldict[n], lw=3)
        
    plt.legend(fontsize=20)
    plt.title("Mean Expression Difference Across Time Steps",fontsize=24)
    plt.ylabel("Mean Expression", fontsize=20)
    plt.xlabel("Simulation Time Steps", fontsize=20)
    my_ticks = list(range( max_path_ln ))
    plt.xticks(my_ticks, [i+1 for i in my_ticks])
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.1, 1.), shadow=True, ncol=1, fontsize=15)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.savefig("results/diff_exp_change_per_timestep_{}.svg".format(name))


#diff_chst(noxtfs_set, symb_to_index, max_path_ln, out_chst, b, entrez_to_symb, coldict, "test")

def diff_chst_subplots(noxtfs_set, symb_to_index, max_path_ln, out_chst, oc_comb, entrez_to_symb, coldict, name):
    ze = {}
    ze1 = {}
    ze2 = {}
    amount = len(out_chst) * len(out_chst[1])
    for i in sorted(list(noxtfs_set)):
        tindx = symb_to_index[i]
        ze[tindx] = np.zeros((max_path_ln))
        ze1[tindx] = np.zeros((max_path_ln))
        ze2[tindx] = np.zeros((max_path_ln))

    for k,v in out_chst.items():
        for kk,vv in v.items():
            for kkk,vvv in vv.items():
                ze[kkk] += vvv/amount
                ze1[kkk] += oc_comb[k][kk][kkk]/amount
                ze2[kkk] += vvv/amount - oc_comb[k][kk][kkk]/amount
    

    fig, axes = plt.subplots(1, 3, figsize=(12*3, 12), sharey=True)
    for j,i in enumerate([ze, ze1, ze2]):
        for k,v in i.items():
            n = entrez_to_symb[entrez_in_dat[k]]
            axes[j].plot(v, label=n, color = coldict[n], lw=3)
    
        if j == 0:
            axes[j].set_ylabel("Mean Expression", fontsize=34)
        axes[j].set_xlabel("Simulation Time Steps", fontsize=34)
        my_ticks = list(range( max_path_ln ))
        axes[j].set_xticks(my_ticks, [i+1 for i in my_ticks], fontsize=30)
        if j == 2:
            axes[j].legend(fontsize=30)
        axes[j].tick_params(axis='y', labelsize=30)
        
    fig.tight_layout()
    plt.savefig("results/diff_exp_change_per_timestep_{}.svg".format(name))


def make_small_dct(noxtfs_set, nms, scales, lst):
    out_dct = {"scale":[]}
    for i in noxtfs_set:
        out_dct[i] = []
    
    for i in range(len(nms)):
        out_dct[nms[i]].append(lst[i])
        if scales[i] not in out_dct["scale"]:
            out_dct["scale"].append(scales[i])
            
    return out_dct


def calc_pvals(out_chst, oc, noxtfs_set):
    real = {}
    bb = {}
    real_molc = {}
    bb_molc = {}
    for k in out_chst:
        real[k] = {}
        bb[k] = {}
        for kk in out_chst[1][(0,0)]:
            real[k][kk] = []
            bb[k][kk] = []
            real_molc[kk] = []
            bb_molc[kk] = []
        
    for k in oc:
        for kk in oc[k]:
            for kkk in out_chst[k][kk]:
                real[k][kkk] += list(out_chst[k][kk][kkk])
                real_molc[kkk] += list(out_chst[k][kk][kkk])
            
    for k in oc:
        for kk in oc[k]:
            for kkk in oc[k][kk]:
                bb[k][kkk] += sum([list(i) for i in oc[k][kk][kkk]],[])
                bb_molc[kkk] += sum([list(i) for i in oc[k][kk][kkk]],[])
                
    ps = []
    scales = []
    nms = []
    for k,v in real.items():
        for kk,vv in v.items():
            _,p = mannwhitneyu(vv, bb[k][kk])
            n = entrez_to_symb[entrez_in_dat[kk]]
            nms.append(n)
            ps.append(p)
            scales.append(k)
            
    for k,v in real_molc.items():
        _,p = mannwhitneyu(v, bb_molc[k])
        n = entrez_to_symb[entrez_in_dat[k]]
        nms.append(n)
        ps.append(p)
        scales.append("all_combined")
    
    cp = multipletests(ps, alpha=0.01, method="bonferroni")
    
    dp = make_small_dct(noxtfs_set, nms, scales, ps)
    dft = make_small_dct(noxtfs_set, nms, scales, cp[0])
    dq = make_small_dct(noxtfs_set, nms, scales, cp[1])
    

    with pd.ExcelWriter("results/simulated_expression_data_diff.xlsx") as writer:
        pd.DataFrame(dp).to_excel(writer, sheet_name="p_values", index=False)
        pd.DataFrame(dft).to_excel(writer, sheet_name="bonf_true_false", index=False)
        pd.DataFrame(dq).to_excel(writer, sheet_name="bonf_corrected", index=False)
    

    #return bb
#calc_pvals(out_chst, oc, noxtfs_set)



def av_of_double_dct(dct):
    out = {}
    for k in dct:
        out[k] = {}
        for kk in dct[k]:
            out[k][kk] = np.mean(dct[k][kk])
    return out


def make_fig_data(g, tfs={"NOX4", "NOX5"}):
    
    org = pd.read_excel("results/simulated_expression_data_orig.xlsx", sheet_name="mean_absolute_exp_diff")
    org = pd.concat([org, pd.DataFrame(np.mean(org.values,axis=0, keepdims=True), columns=org.columns)], axis=0)
    sc = list(org["scale"])
    sc[-1] = "Mean:"
    org["scale"] = sc


    ra = pd.read_excel("results/simulated_expression_data_random_start.xlsx", sheet_name="mean_absolute_exp_diff")
    ra = pd.concat([ra, pd.DataFrame(np.mean(ra.values,axis=0, keepdims=True), columns=ra.columns)], axis=0)
    sc = list(ra["scale"])
    sc[-1] = "Mean:"
    ra["scale"] = sc
    
    
    cl = list(org.columns[1:][org.values[-1][1:] < 0])
    
    t1 = org.copy()
    #t1[cl] -= ra[cl]
    t1[cl] = org[cl] - ra[cl]#(ra[cl] - org[cl]) / (ra[cl] + org[cl])
    t1 = t1[["scale"] + cl]
    
    deg = []
    for i in cl:
        deg.append(np.sum([1 for i in g.out_edges(i) if i[1] in tfs]) / len(tfs))
    
    t2 = pd.DataFrame({ "TF" : cl, "score" : t1[cl].values[-1] * np.array(deg)})
    t2 = t2.sort_values(by="score").reset_index(drop=True)
    t2["ranking"] = np.arange(len(t2)) + 1


    with pd.ExcelWriter("results/figure_data.xlsx") as writer:
        org.to_excel(writer, sheet_name="A", index=False)
        t1.to_excel(writer, sheet_name="B", index=False)
        t2.to_excel(writer, sheet_name="C", index=False)



if __name__ == "__main__":
    
    G, pkgTarg, entrez_to_symb, pths, max_path_ln, noxtfs, g, unique_pkg1_targets, noxtfs_set = load_data_and_make_graph()
    
    nxtfs, lnxtfs = get_tf_w_pths_to_noxtfs(G, noxtfs_set, pkgTarg, rng)
    
    generate_pkg_to_nox_img(G, g, pkgTarg)
    
    cntr_dat, entrez_in_dat, coldat = extracting_gene_values(G, entrez_to_symb, path = "annotated_ctrl_exprs_stacked.tsv/annotated_ctrl_exprs_stacked.tsv")
    
    symb_to_index = {}
    for j,i in enumerate(entrez_in_dat):
        symb_to_index[entrez_to_symb[i]] = j


    names = ["orig", "random_start", "c1_shuf", "c2_zeroed", "c3_random"]
    
    for i in range(len(names)):

        print("running {}".format(names[i]))
            
        if i != 1:
            if i == 0:
                use_dat = cntr_dat.copy()
                plot_scaled_sigmoid(s = scale)
                gg = nx.subgraph(G, set(sum(pths, [])))
            elif i == 2:
                use_dat = shufler(cntr_dat.copy(), sq = False)
            elif i == 3:
                use_dat = np.zeros(cntr_dat.shape)
            elif i == 4:
                use_dat = random_sampler(cntr_dat.shape, rng)
            
            out_true_ab, av_tr, av_tr_per_scale, av_pred, av_pred_per_scale, out_chst = compute_simulation(gg, use_dat, symb_to_index, pkgTarg, noxtfs_set, max_path_ln, scale)

            plot_mean_expression_change(noxtfs_set, symb_to_index, max_path_ln, out_chst, entrez_to_symb, coldict, names[i])
    
            write_data_to_excel(noxtfs_set, av_pred_per_scale, av_tr_per_scale, names[i], scale)
        
        else:
            use_dat = cntr_dat.copy()
            atps, apps, oc, oc_comb = run_multiple_times(G, use_dat, symb_to_index, pkgTarg, noxtfs_set, scale, nxtfs, lnxtfs, max_path_ln)
            calc_pvals(out_chst, oc, noxtfs_set)
            plot_mean_expression_change(noxtfs_set, symb_to_index, max_path_ln, oc_comb, entrez_to_symb, coldict, names[i])
            diff_chst_subplots(noxtfs_set, symb_to_index, max_path_ln, out_chst, oc_comb, entrez_to_symb, coldict, "orig_randomstart_diff")
            write_data_to_excel(noxtfs_set, av_of_double_dct(apps), av_of_double_dct(atps), names[i], scale)
        
        print("End of simulation\n")
    
    make_fig_data(g)
        
        