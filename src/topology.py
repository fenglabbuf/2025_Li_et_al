import networkx as nx
import pandas as pd
import numpy as np
import powerlaw
from multiprocessing import Process
import os.path

def run_with_limited_time(func, args, kwargs, time):
    '''runs a function with time limit'''
    p = Process(target=func, args=args, kwargs=kwargs)
    p.start()
    p.join(time)
    if p.is_alive():
        p.terminate()
        return False
    return True

def read_directed_grn(infile):
    '''read file as networkx directed graph'''
    net = pd.read_csv(infile, sep='\t', header=None)
    DG = nx.DiGraph()
    DG.add_weighted_edges_from(list(zip(net[0], net[1], net[2])))
    return DG

def remove_edge(nxg, threshold, component='all'):
    '''remove edges in graph weighted smaller than edge weight threshold (absolute value) and 
    keep the largest connected component in the remaining graph. Component has to be one of 
    all, positive or negative to keep all, positively weighted edges or negatively weighted 
    edges in the graph'''
    if component=='all':
        edgelst = [edge for edge, weight in dict(nxg.edges()).items() \
                    if abs(weight['weight']) < threshold]
    elif component=='positive':
        edgelst = [edge for edge, weight in dict(nxg.edges()).items() \
                    if weight['weight'] < threshold]
    elif component=='negative':
        edgelst = [edge for edge, weight in dict(nxg.edges()).items() \
                    if -weight['weight'] < threshold]
    nxgc = nxg.copy()
    nxgc.remove_edges_from(edgelst)
    
    nxug = nx.Graph()
    nxug.add_edges_from(nxgc.edges())
    if not nx.is_connected(nxug):
        nodes = max(nx.connected_components(nxug), key=len)
        removen = [n for n in list(nxgc.nodes()) if n not in nodes]
        nxgc.remove_nodes_from(removen)
    nxgc.remove_nodes_from(list(nx.isolates(nxgc)))
    return nxgc

def nxedge_to_positive(DGo):
    '''transform negative to positive edge weights in a negatively weighted graph'''
    DG = nx.DiGraph()
    DG.add_weighted_edges_from([(a, b, abs(c)) for (a, b, c) in list(DGo.edges.data('weight'))])
    return DG

def compute_topology_parameters(DG, outfile, treatment, run, threshold, component):
    '''return number of nodes, number of edges, small world coefficient sigma, 
    number of nodes in graph dominating set, powerlaw parameters of a graph to text file'''
    nodes = len(DG.nodes())
    edges = len(DG.edges())
    # transform directed to undirected graph as 
    # sigma has to be computed in an undirected graph
    G = nx.Graph()
    G.add_edges_from(DG.edges())
    sigma = nx.sigma(G, niter=1, nrand=1)
    dmn = len(nx.dominating_set(DG))

    # compute powerlaw parameters
    data = np.array(list(dict(nx.degree(DG)).values()))
    results = powerlaw.Fit(data, discrete=True)
    (exp_r, exp_p) = results.distribution_compare('power_law', 'exponential')
    (logn_r, logn_p) = results.distribution_compare('power_law', 'lognormal')
    alpha = results.power_law.alpha
    xmin = results.power_law.xmin

    d_params = {'treatment':treatment, 'run':run, 'threshold':threshold, 'component':component}
    d_params.update({'nodes':nodes, 'edges':edges, 'sigma':sigma, 'dominating_size':dmn})
    d_params.update({'alpha':alpha, 'xmin':xmin, 'exp_r': exp_r, 'exp_p':exp_p})
    d_params.update({'logn_r':logn_r, 'logn_p':logn_p})

    if not os.path.exists(outfile):
        with open(outfile, 'w') as fw:
            fw.write('\t'.join(list(d_params.keys())))
            fw.write('\n')

    with open(outfile, 'a') as f:
        for value in list(d_params.values()):
            f.write(f'{value}\t')
        f.write('\n')
    return 