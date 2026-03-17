from hcrseq.common.ref import Reference, Reporter, Pathway 

if __name__ == '__main__':
    plasmid_info = 'plasmid_info.txt'
    reporter_info = 'reporter_info.txt'
    pathway_info = 'pathway_calculations.txt'
    
    # Now when you build and save, the pickle will remember the source 
    # as 'hcrseq.common.ref.Reference' instead of '__main__.Reference'
    r = Reference.build(plasmid_info, reporter_info, pathway_info)
    r.write_pickle('HCRseq_v0.4.pkl')
