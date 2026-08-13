import pandas as pd
import scanpy as sc

from hcrseq.common.quantify import HCRseqQuantifier
from hcrseq.common.ref import Reference

def quantify_repair(h5_file,bam,
                    ref_path,
                    min_mapq = 5):
    """
    Iterates over bam to quantify repair metrics, then adds them to anndata object in h5_file
    """

    if h5_file.endswith('.h5'):
        adata = sc.read_10x_h5(h5_file)
    else:
        adata = sc.read(h5_file)
    
    ref = Reference.load(ref_path)
    reporter_names = [reporter.name for reporter in ref.reporters]

    plasmid_idx = adata.var.index.isin(reporter_names)
    adata.obs[adata.var.index[plasmid_idx]] = adata[:,plasmid_idx].to_df()
    adata = adata[:,~plasmid_idx]

    # Perform the reporter counting and quantification
    q = HCRseqQuantifier(ref,("CB","UB"))
    q.count_umis(bam,min_mapq=min_mapq)

    # Store the reporter count information (one obsm matrix of cell x reporter per metric)
    cell_counts = q.get_cell_counts()
    counts_df = pd.concat(
        {reporter_name: pd.DataFrame.from_dict(counts_by_cb, orient='index')
         for reporter_name, counts_by_cb in cell_counts.items()},
        names=['reporter','CB']).fillna(0)

    for metric in counts_df.columns:
        adata.obsm[metric + '_counts'] = adata.obs[[]].join(
            counts_df[metric].unstack('reporter'),how='left').fillna(0)

    # Map per-cell repair pathway measurements into adata.obs
    repair_measurements = q.quantify_repair_per_cell()
    for pathway_name, values_by_cb in repair_measurements.items():
        pathway_df = pd.DataFrame(values_by_cb).T
        adata.obs[pathway_name] = pathway_df['value'].reindex(adata.obs.index)
        adata.obs[pathway_name + '_sd'] = pathway_df['sd'].reindex(adata.obs.index)

    # Record pathway names in adata.uns
    adata.uns['pathways'] = list(repair_measurements.keys())

    return(adata)

