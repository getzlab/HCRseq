

BP = {'A' : 'T',
      'C' : 'G',
      'G' : 'C',
      'T' : 'A',
      'N' : 'N'}

def rc(seq):
    return("".join([BP[b] for b in reversed(seq)]))

