from find_modules import *
from topology import *

random_s = 10
random.seed(random_s)

def compute_grn_centrality(DG, modules, random_s=10, count=False):
    '''return number of nodes, gene names, global and module centrality to pandas dataframe 
    from input graph and dict of community membership'''
    df = pd.DataFrame.from_dict(modules, orient='index', columns=['community'])
    d_btw = nx.betweenness_centrality(DG, weight='weight', seed=random_s)
    d_odeg = nx.out_degree_centrality(DG)
    df.insert(0, 'betweenness_centrality', pd.Series(d_btw.values(), index=d_btw.keys()), allow_duplicates=True)
    df.insert(0, 'outdegree_centrality', pd.Series(d_odeg.values(), index=d_odeg.keys()), allow_duplicates=True)
    disconnect = [(a,b) for (a,b) in list(DG.edges()) if modules[a]!=modules[b]]
    DG.remove_edges_from(disconnect)
    d_btw = nx.betweenness_centrality(DG, weight='weight', seed=random_s)
    d_odeg = nx.out_degree_centrality(DG)
    df.insert(0, 'local_betweenness_centrality', pd.Series(d_btw.values(), index=d_btw.keys()), allow_duplicates=True)
    df.insert(0, 'local_outdegree_centrality', pd.Series(d_odeg.values(), index=d_odeg.keys()), allow_duplicates=True)
    if count:
        dfcount = df.groupby('community', as_index=False).count()[['community', 'betweenness_centrality']]
        dfcount = dfcount.rename(columns={'betweenness_centrality':'community_size'})
        df = df.reset_index().merge(dfcount, left_on='community', right_on='community', how='left').set_index('index')
    return pd.concat([df], keys=[len(modules)], names=['nodes', 'geneid'])