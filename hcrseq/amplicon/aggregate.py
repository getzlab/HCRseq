import pandas as pd
import yaml
import json

def aggregate_repair_outputs(files,ids):
    """
    Given list of pathway repair metrics
    """
    xs = list()
    for id,file in zip(ids,files):
        with open(file, 'rt') as f:
            x = pd.DataFrame(yaml.safe_load(f)).T
            x.index.name  = 'pathway'
            x['sample'] = id
            xs.append(x)
    x = pd.concat(xs,axis=0).reset_index()

    return(x)


def aggregate_count_outputs(files, ids):
    """
    Merges count output files into single dataframe
    """
    xs = list()
    for sid, file in zip(ids, files):
        with open(file, 'rt') as f:
            x = yaml.safe_load(f)
            x = pd.DataFrame(x).reset_index().melt(id_vars='index')
            x.columns = ['metric', 'reporter', 'count']
            x['sample'] = sid
            xs.append(x)
    x = pd.concat(xs, axis=0).reset_index(drop=True)

    return (x)

def aggregate_cutadapt_outputs(files,ids):

    xs = list()
    for sid, file in zip(ids, files):
        with open(file, 'rt') as f:
            x = json.load(f)
            x = {'sample' : sid,'input_reads' : x['read_counts']['input'],
             'reads_post_cutadapt' :  x['read_counts']['output']}
            xs.append(x)
    df = pd.DataFrame(xs).set_index('sample')

    return(df)

def aggregate_results(cutadapt_files,count_files,repair_files,ids,outstem):
    R = aggregate_repair_outputs(repair_files,ids)
    C = aggregate_count_outputs(count_files,ids)
    Q = aggregate_cutadapt_outputs(cutadapt_files,ids)

    Q['final_umis'] = C[C['metric']=='total'].groupby('sample')['count'].sum()

    R.to_csv(f'{outstem}.repair_measurements.csv')
    C.to_csv(f'{outstem}.reporter_counts.csv')
    Q.to_csv(f'{outstem}.qc_metrics.csv')



