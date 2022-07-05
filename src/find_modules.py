import random
import networkx as nx
import pandas as pd
import numpy as np
import sknetwork as skn
import infomap
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'arial'

# community colors
colors = ['#00a679', '#edab00', '#7b68b5', '#e74e00', '#48af17','#ac7500', 
        '#277cbb', '#00b94a', '#f60000', '#a533a3', '#ff7100', '#f9ff00', '#b04e17', '#ff68bd',
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

def louvain_cluster(G, directed = True, random_state=False, random_s=10):
    '''perform louvain clustering and return a dict of module memberships'''
    skg = skn.data.convert_edge_list(list(G.edges.data('weight')), directed=directed)
    if random_state:
        louvain = skn.clustering.Louvain(random_state=random_s)
    else:
        louvain = skn.clustering.Louvain()
    labels = louvain.fit_transform(skg.adjacency)
    partition = skmp_to_nxp(skg, G, labels)
    return partition

def jaccard_idx(lsta, lstb):
    '''calculate jaccard index'''
    a = set(lsta)
    b = set(lstb)
    c = a.intersection(b)
    idx = float(len(c)/(len(a) + len(b) - len(c)))
    return idx

def _position_nodes(g, partition, scale, **kwargs):
    # positions of nodes within communities
    communities = dict()
    for node, community in partition.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]
            
    pos = dict()
    for ci, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, k=scale)
        pos.update(pos_subgraph)
    return pos

def _find_between_community_edges(g, partition):
    edges = dict()
    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]
        
        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]
    return edges

def _position_communities(g, partition, **kwargs):
    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(g, partition)
    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))
        
    # find layout for communities
    pos_communities = nx.circular_layout(hypergraph, **kwargs)
    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]
    return pos

def community_layout(g, partition, scales_c=5., scales_n=1.):
    # modified from https://stackoverflow.com/questions/43541376/how-to-draw-communities-with-networkx/43541777
    """
    Compute the layout for a modular graph.
    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance
        graph to plot

    partition -- dict mapping int node -> int community
        graph partitions
    Returns:
    --------
    pos -- dict mapping int node -> (float x, float y)
        node positions
    """
    pos_communities = _position_communities(g, partition, scale=scales_c)
    pos_nodes = _position_nodes(g, partition, scale=scales_n)
    # combine positions
    pos = dict()
    for node in g.nodes():
        pos[node] = pos_communities[node] + pos_nodes[node]
    return pos

def get_cmap(partition, highlight=False, modulelist=[], alpha=0.7):
    '''assign colors to communities'''
    ncolor = len(np.unique(np.array(list(partition.values()))))
    colormap = dict(zip(list(np.unique(np.array(list(partition.values())))), 
                        (colors*2)[:ncolor]))
    if highlight:
        for c in colormap:
            if c not in modulelist:
                c_arry = np.array(mpl.colors.hex2color(colormap[c]))
                colormap[c]=tuple(alpha+c_arry*(1-alpha))
            else:
                colormap[c]=mpl.colors.hex2color(colormap[c])
    return colormap

def highlight_cmap(list1, list2, highlight = '#fc0303', backgroud='#696969'):
    '''return background colors of list1 with highlighted colors if the element appears in list2'''
    colormap=[]
    for n in list1:
        if n in list2:
            colormap.append(highlight)
        else:
            colormap.append(backgroud)
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


