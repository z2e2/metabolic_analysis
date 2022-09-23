import numpy as np
from scipy import stats
from statsmodels.stats import multitest
from scipy import interpolate

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

def adjust_p_value_qval(pv, pi_0=1):
    '''
    Storey method: adjusting p-value for multiple tests
    
    Parameters
    ----------
    p_val_list : list
    A list of p-values to be adjusted.
    pi_0 : float
    Estimated proportion of true null hypothesis

    Returns
    -------
    The adjusted p-values.
    '''
    m = len(pv)
    if pi_0 is None:
        m0 = []
        lam = np.arange(0, 0.9, 0.01)
        for i in lam:
            m0.append(np.sum(pv > i)/(1-i))
        pi_0 = np.array(m0)/m
        tck = interpolate.splrep(lam, pi_0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        pi_0 = pi_0[-1] if pi_0[-1] < 1 else 1

    sorted_index = np.argsort(pv)
    qv_le = np.array([np.sum(np.array(pv)<=item) for item in pv])
    fdr = pv*pi_0*m/(qv_le)
    for i in range(m):
        idx = sorted_index[::-1][i]
        if i == 0:
            fdr[idx] = min(fdr[idx], 1)
        else:
            fdr[idx] = min(fdr[idx], fdr[sorted_index[::-1][i-1]])
    return fdr

def adjust_p_value_bh(p_val_list):
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
    res = multitest.multipletests(p_val_list, method='fdr_bh')
    return res[1]

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
    k = np.array([np.sum(np.array(p_val_list)<=item) for item in p_val_list])
    p_val_adj = np.array(p_val_list) * m/k
    p_val_adj[p_val_adj>=1] = 1
    return p_val_adj

def get_overlapping_met_genes(met_dict, all_genes, verb = False):
    '''
    Get the genes in the metabolism pathway that were also captured in a dataset

    Parameters
    ----------
    met_dict : dict
    For each metabolism term of interest, the value is a list of genes associated with that term.
    all_genes : list
    All of the genes captured in a dataset
    verb : bool
    Optionally print out the number of genes that are overlapping between the metabolism terms and data

    Returns
    ----------
    overlap_res : A filtered metabolism dict with non-shared genes removed.
    n_removed_terms : for each metabolism term, the number of genes removed
    '''
    overlap_res = {}
    n_removed_terms = {}
    for k,v in met_dict.items():
        overlap_res[k] = [i for i in v if i in all_genes]
        if verb:
            print(f'{k}\toriginal: {len(v)} | overlapping: {len(overlap_res[k])} | removed: {len(v) - len(overlap_res[k])}')
        n_removed_terms[k] = len(v) - len(overlap_res[k])
    return overlap_res, n_removed_terms
