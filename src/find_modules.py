import networkx as nx
import pandas as pd
import numpy as np
import random
import sknetwork as skn
import infomap
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'arial'

# community colors
colors = ['#00a679', '#edab00', '#7b68b5', '#e74e00', '#48af17','#ac7500', '#666666', 
        '#f60000', '#277cbb', '#00b94a', '#a533a3', '#ff7100', '#f9ff00', '#b04e17', '#ff68bd',
        '#469990', '#e07099', '#3cb44b', '#808000', '#800000', '#ffe119', '#4363d8', '#f58231', 
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', 
        '#fffac8', '#aaffc3', '#ffd8b1', '#000075', '#808080', '#000000']

def skmp_to_nxp(skg, nxg, labels):
    '''scikit network partition format to networkx partition format'''
    partition={}
    for n in range(len(labels)):
        partition[list(nxg.nodes())[n]] = labels[list(skg.names).index(list(nxg.nodes())[n])]
    return partition

def nxtup_to_nxp(nxg, comp):
    '''networkx Girvan Newman partition to networkx partition format'''
    tup = tuple(sorted(c) for c in next(comp))
    partition={}
    for n in range(len(list(nxg.nodes()))):
        for c in range(len(tup)):
            if list(nxg.nodes())[n] in tup[c]:   
                partition[list(nxg.nodes())[n]] = c
    return partition

def info_nxp(nxg, directed=False):
    '''perform informap module detection
    infomap partition format to networkx partition format'''
    if directed:
        im = infomap.Infomap('--two-level --directed')
    else:
        im = infomap.Infomap()
    
    mapps = dict(list(zip(list(nxg.nodes()), range(len(nxg.nodes())))))
    imnxg = nx.relabel_nodes(nxg, mapps)
    im.add_links(list(imnxg.edges.data('weight')))
    im.run()
    partition = im.get_modules()
    return {k: partition[nid] for k, nid in mapps.items()}

def plot_degree(G, n_bins = 10):
    '''compute graph degree distribution binned by log10'''
    degree = np.array(list(dict(nx.degree(G)).values()))
    minx = np.log10(np.min(degree)) if np.min(degree) >=1 else 0.0
    maxx = np.log10(np.max(degree))
    bins = np.logspace(minx, maxx, n_bins, base=10)
    y,_ = np.histogram(degree, bins = bins, density=True)
    x = bins[1:] - np.diff(bins)/2.0
    return x, y

def louvain()

def jaccard_idx(lsta, lstb):
    '''calculate jaccard index'''
    a = set(lsta)
    b = set(lstb)
    c = a.intersection(b)
    idx = float(len(c)/(len(a) + len(b) - len(c)))
    return idx

def get_cmap(partition):
    '''assign colors to communities'''
    # use 1111 to represent absent nodes, set color to white
    if 1111 not in list(partition.values()):
        ncolor = len(np.unique(np.array(list(partition.values()))))
        colormap = dict(zip(list(np.unique(np.array(list(partition.values())))), 
                            (colors*(int(ncolor/22)+1))[:ncolor]))
    else:
        ncolor = len(np.unique(np.array(list(partition.values()))))-1
        colormap = dict(zip(list(filter(lambda x:x !=1111, list(np.unique(np.array(list(partition.values())))))), 
                            (colors*(int(ncolor/22)+1))[:ncolor]))
        colormap[1111]='#ffffff'
    return colormap

def get_greymap(partition, module):
    '''highlight the color of selected modules'''
    ncolor = len(np.unique(np.array(list(partition.values()))))
    colormap = dict(zip(list(np.unique(np.array(list(partition.values())))), 
                        (colors*(int(ncolor/22)+1))[:ncolor]))
    for n in colormap:
        if n not in module:
            colormap[n] = '#999999'
    return colormap

def goplot(file):
    '''generate gene ontology bar plot'''
    dfannots = pd.read_csv(file)
    dfannots = dfannots[dfannots['p.adjust']<0.01]
    dfannots['refgs'] = [int(n.split('/')[0]) for n in dfannots['BgRatio']]
    dfannots = dfannots[(dfannots['refgs']<1500)&(dfannots['refgs']>10)]
    dfannots = dfannots.drop_duplicates(subset='geneID', keep='first')
    dfannots = dfannots[:10]
    dfannots = dfannots.sort_values(by=['p.adjust'], ascending=False)
    if len(dfannots)>=1:
        if len(dfannots)>1:
            fig = plt.figure(figsize=(2, 1.5*np.log(len(dfannots))))
        else:
            fig = plt.figure(figsize=(2, 0.4))
        ax = fig.add_subplot(1, 1, 1)
        plt.barh(list(dfannots['Description']), -np.log10(dfannots['p.adjust']), 
                 color=colors[int(file[-5])], alpha=0.2, edgecolor=colors[int(file[-5])], linewidth=2)
        plt.barh(list(dfannots['Description']), -np.log10(dfannots['p.adjust']), 
                 color='none', alpha=0.7, edgecolor=colors[int(file[-5])], linewidth=2)
        plt.yticks(fontsize=13)
        plt.title(file, fontsize=13)
        plt.xlabel('-log'+r'$_{10}$'+'(FDR)', fontsize=13)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
    else:
        print(file)
    return


