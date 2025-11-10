import pandas as pd
import anndata as ad


def get_repair_df(adata,
                  pathways=None,
                  groupby = ['cell_line','timepoint'],
                 agg_fun = 'mean'):

    if pathways is None:
        if 'pathways' not in adata.uns.keys():
            raise('Error! If pathways is not a key of adata.uns then it must be specified')
        pathways = adata.uns['pathways']
    
    df = adata.obs[adata.obs['transfected']].\
            groupby(groupby,observed=True)[pathways].agg('mean').\
    reset_index().melt(id_vars=groupby,
                       var_name='pathway',
                       value_name='repair')

    counts = adata.obs[adata.obs['transfected']].\
            groupby(groupby,observed=True)[pathways].count().\
    reset_index().melt(id_vars=groupby,
                       var_name='pathway',
                       value_name='n')
    merge_keys = list(groupby) + ['pathway']
    df = df.join(counts.set_index(merge_keys)['n'],on=merge_keys)

    return(df)


