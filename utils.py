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
