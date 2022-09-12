import numpy as np
from scipy import stats


def hypergeometric_test(total_genes_expressed, n_genes_of_interest, n_genes_picked, n_overlap):
    '''
    Hypergeometric test: does a list of genes significantly overlap with another list of genes?

    Parameters
    ----------
    total_genes_expressed : int
    The total number of genes expressed in a dataset.
    n_genes_of_interest : int
    The number of genes in a pathway of interest.
    n_genes_picked : int 
    The number of genes differentially expressed in a cluster.
    n_overlap : int
    The number of genes that are both in the gene list of the pathway and the differential genes list.

    Returns
    -------
    The probability that there are at least `n_overlap` genes in common.
    '''
    hypergeom_model = stats.hypergeom(M=total_genes_expressed, n=n_genes_of_interest, N=n_genes_picked)
    p_val = hypergeom_model.sf(n_overlap-1)
    return p_val


def adjust_p_value_fdr(p_val_list):
    '''
    Benjamini and Hochberg method: adjusting p-value for multiple tests
    
    Parameters
    ----------
    p_val_list : list
    A list of p-values to be adjusted.

    Returns
    -------
    The adjusted p-values.
    '''
    m = len(p_val_list)
    k = np.argsort(p_val_list) + 1
    p_val_adj = np.array(p_val_list) * m/k
    p_val_adj[p_val_adj>=1] = 1
    return p_val_adj
