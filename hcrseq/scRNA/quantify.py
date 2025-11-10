import pysam
import pandas as pd
import scanpy as sc

def count_umis_at_lesion(bam, lesion_info,min_mapq = 5):

    L = pd.read_csv(lesion_info, sep='\t', index_col='reporter')

    
    counts = dict()
    with pysam.AlignmentFile(bam) as bam_in:
        for contig,row in L.iterrows():
            
            counts[contig] = dict()
            
            for pileupcolumn in bam_in.pileup(contig,row['position']-1,row['position'],
                                              truncate=True,
                                              min_mapping_quality=min_mapq,
                                             max_depth=1e9):

                print(f"Checking contig={contig}, pos={row['position']}")
                i=0
                for pileupread in pileupcolumn.pileups:
                    i+=1
                    tags = dict(pileupread.alignment.get_tags())
                
                    if ('CB' not in tags) or ('UB' not in tags):
                        continue
                        
                    if pileupread.is_del:
                        key = 'del'
                    elif not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        
                        if base == L.loc[contig, 'repaired_base']:
                            key = 'repaired'
                        elif base == L.loc[contig, 'unrepaired_base']:
                            key = 'unrepaired'
                        else:
                            key = 'other'

                    cb = tags['CB']
                    ub = tags['UB']
                    
                    if cb not in counts[contig]:
                        counts[contig][cb] = dict()
                    
                    if ub not in counts[contig][cb]:
                        counts[contig][cb][ub] = dict({'repaired':0,'unrepaired':0,'del' : 0,'other':0})
                    counts[contig][cb][ub][key] +=1
                print(f'analyzed {i} reads')

                 
    umis = get_umi_counts(counts)
    return(umis)

def get_umi_counts(counts):
    umis = dict()
    for contig,v in counts.items():
        for cb,vv in v.items():
            for ub,vvv in vv.items():
                key = pd.Series(vvv).idxmax()
            
                if (contig,cb) not in umis:
                    umis[(contig,cb)] = dict({'repaired':0,'unrepaired':0,'del' : 0,'other':0})
                umis[(contig,cb)][key] += 1
    return(umis)

def quantify_repair(h5_file,bam,
                    lesion_info,
                    pathway_info,
                    min_mapq = 5):

    if h5_file.endswith('.h5'):
        adata = sc.read_10x_h5(h5_file)
    else:
        adata = sc.read(h5_file)
    
    plasmid_idx = adata.var.index.str.match('GFP_')
    adata.obs[adata.var.index[plasmid_idx]] = adata[:,plasmid_idx].to_df()
    adata = adata[:,~plasmid_idx]
    
    # Collect site-specific counts from bam
    repair_counts = count_umis_at_lesion(bam,lesion_info,min_mapq=min_mapq)
    
    df = pd.DataFrame(repair_counts).T

    for k in df.columns:
        adata.obsm[k + '_counts'] = adata.obs[[]].join(df[k].unstack().T,how='left')

    df['f'] = df['repaired'] / (df['repaired']+df['unrepaired'])
        
    df = df.reset_index().pivot(index='level_1',columns='level_0')
    df.columns = df.columns.swaplevel(0,1)
    df = df.sort_index(axis=1)

    # Filter to called cells
    df = df[df.index.isin(adata.obs.index)]

    P = pd.read_csv(pathway_info,sep='\t')

    ## Identify control plasmid
    idx = P['pathway'].str.match('control',case=False)
    if sum(idx)>1:
        raise("Error, multiple control plasmids specified!")
    if sum(idx)==1:
        control = P.loc[idx,'reporter'].iloc[0]
    if sum(idx)==0:
        print('Warning! No control plasmid specified:')
        control = 'dummy'

    
    # Quantify pathways besides control    
    for ind,row in P[~idx].iterrows():
        
        if row['metric']=='fraction_repaired':
            adata.obs[row['pathway']] = df[row['reporter']]['f']
        elif row['metric']=='abundance':
            adata.obs[row['pathway']] = adata.obs[row['reporter']] /  adata.obs[control]

    return(adata)

