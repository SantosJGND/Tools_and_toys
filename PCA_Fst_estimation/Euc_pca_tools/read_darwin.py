import numpy as np
import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


def read_refs(index_file):
    indxs = recursively_default_dict()
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs[int(line[0])][line[1]] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())]



def read_Darwin(darwin_file,miss= '9'):
    d= 0
    Input= open(darwin_file,'r')
    Names= []
    gen= []
    for line in Input:
        line= line.rstrip().split('\t')
        
        if d== 0:
            d+= 1
            Nsnps= max([int(x) for x in line[1:]])
            continue

        art= line[1:]
        art= [[x,'0'][int(x == ' ')] for x in art]

        Names.append(line[0])

        art= [int(x) for x in art]
        art= [[x,np.nan][int(x == miss)] for x in art]
        
        if len(art) < Nsnps:
                        art.extend([0] * (Nsnps - len(art)))
        
        
        gen.append(art)
        d += 1
    Input.close()
    
    return gen, Names

