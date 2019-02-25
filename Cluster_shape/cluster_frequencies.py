# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 13:01:08 2017

@author: jgarcia
"""

import itertools as it

import numpy as np
import pandas as pd
import scipy
import collections

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import scale


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


import StructE_tools as Ste
from sklearn.metrics import silhouette_samples, silhouette_score
from AMOVA_func import AMOVA_FM42, amova_cofactor

########## START HERE #############

import os
import argparse
parser = argparse.ArgumentParser()

### optional arguments
parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--bim",help = "snp information bim format.")

parser.add_argument("--amova",action= 'store_true',help= 'perform AMOVA')

parser.add_argument("--aims",type= str,default= '',help = "genetic window file file")

parser.add_argument("--supervised",action= 'store_true',help= 'use ref groups for supervised analysis.')

parser.add_argument("--mrg",type = int,default = 1000,help = "Margin to genic regions provided in aims file.")

parser.add_argument("--ref",help = "reference accessions indexes in genofile.")

parser.add_argument("--admx",help = "admixed accessions indexes in genofile.")

parser.add_argument("--cl_freqs",type = int,default = 35,help = "Minimum cluster size from which to calculate freqs")
###
parser.add_argument("--clsize",type = int,default = 15,help = "Minimum cluster size")
###
parser.add_argument("--out",type= str,default= '',help = "output directory")

parser.add_argument("--id",type= str,default= 'StructE',help = "Give your analysis an ID.")
### 
parser.add_argument("--dr",default = 'PCA',help = "Dimensionality reduction. options: PCA, NMF")
### 
parser.add_argument("--ncomp",type = int,default = 4,help = "Number of components kept in case of PCA reduction")
### 
parser.add_argument("--random",action='store_true',help= 'random regions')

parser.add_argument("--randN",type= int,default= 100,help= 'Number of random regions.')
###
parser.add_argument("-w",type = int,default = 150, help = "random region size")

args = parser.parse_args()



########## Some new functions.

####
#### read block and profile files
####

print('To begin reading from: ')
print(args.books)

Books= Ste.read_geno_books(args.books)

print(Books)

Home= args.out

########## Complementary files.
#### read names (Fam file in same order as genotype data)
Fam = Ste.FAMread(args.fam)

#### read SNP info (.bim file, same order and number as in geno file)
MissG, Gindex = Ste.BIMread(args.bim)


#### read reference ind files.
refs_lib, Parents, absent_refs  = Ste.read_refs(args.ref,Fam)

if absent_refs:
    print(",".join([str(x) for x in absent_refs]) + ' absent from provided refs.')


admx_lib, Crossed, absent_admx  = Ste.read_refs(args.admx,Fam)

if len(absent_admx) > (0.5 * len(Crossed)):
    print('over half the admixed missing from fam file.')
    
admx_lib.update(refs_lib)

Geneo = admx_lib

print('Population labels: {}'.format(Geneo.keys()))

print('Population Sizes: {}'.format([len(x) for x in Geneo.values()]))


####
####

#### Select regions to analyse.

Bandwidth_split = 30
Names_store = Fam

### Select regions to study
## read from provided file
##

## Generate random regions
if args.random:
    Genes = Ste.Gen_rand(MissG,Books.Chr,args.randN,args.w)
    print('random analysis selected. {} regions of size {} generated.'.format(args.randN,args.w))
else:
    if args.aims == '':
        print('--random nor --Aims provided.')
    ## /gs7k1/home/jgarcia/miRNA/MirBase_clean.txt
    Genes = Ste.read_selected(args.aims,mrg= args.mrg)

## remove regions with no SNPs
Genes = {x:Genes[x] for x in Genes.keys() if len(Genes[x]) > 0}

## procure sequences:
Geneo_order= list(Geneo.keys())
ref_order= list(refs_lib.keys())

Whose= [z for z in it.chain(*[Geneo[x] for x in Geneo_order])]
Sup_labels= np.repeat(ref_order,[len(refs_lib[x]) for x in ref_order])

SequenceStore = Ste.Extract_to_dict(Genes,MissG,Whose,Books)

refs_lib= {x:[Whose.index(z) for z in refs_lib[x]] for x in refs_lib.keys()}

###########################
### Prepare libraries
###########################
### Define parameters and libraries of analyses.
### Define parameters and libraries of analyses.

if args.amova:
    print('amova selected.')
Bandwidth_split= 30
clsize= args.clsize

control_N= 100

Results = {x:recursively_default_dict() for x in SequenceStore.keys()}

Frequencies = {x:recursively_default_dict() for x in SequenceStore.keys()}

Clover= []
Coordinates= []
Clusters_coords= []

Construct = recursively_default_dict()
PC_var= recursively_default_dict()


for CHR in SequenceStore.keys():
    print('going on CHR: '+ str(CHR))
    for c in SequenceStore[CHR].keys():
        
        ### PCA and MeanShift of information from each window copied from *FM36_Galaxy.py.
        Sequences= [SequenceStore[CHR][c][x] for x in Whose]
        Sequences= np.array(Sequences)
        
        if Sequences.shape[1] <= 3:
            Results[CHR][c] = [0,0]
            continue
        
        Sequences= np.nan_to_num(Sequences)
        
        pca = PCA(n_components=args.ncomp, whiten=False,svd_solver='randomized').fit(Sequences)
        data = pca.transform(Sequences)
        PC_var[CHR][c]= [x for x in pca.explained_variance_]
        
        params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
        grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

        ######################################
        ####### TEST global Likelihood #######
        ######################################
        Focus_labels = list(range(Sequences.shape[0]))
        
        #### Mean Shift approach
        ## from sklearn.cluster import MeanShift, estimate_bandwidth

        bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=len(Focus_labels))
        if bandwidth <= 1e-3:
            bandwidth = 0.1

        ms = MeanShift(bandwidth=bandwidth, cluster_all=False, min_bin_freq=clsize)
        ms.fit(data[Focus_labels,:])
        labels = ms.labels_
        
        Tree = {x:[Focus_labels[y] for y in range(len(labels)) if labels[y] == x] for x in [g for g in list(set(labels)) if g != -1]}
        Keep= [x for x in Tree.keys() if len(Tree[x]) > clsize]

        Tree= {x:Tree[x] for x in Keep}
        Ngps= len(Tree)
        SpaceX = {x:data[Tree[x],:] for x in Tree.keys()}

        ### Extract MScluster likelihood by sample
        
        for hill in SpaceX.keys():
            
            if len(Tree[hill]) >= args.cl_freqs:
                cl_seqs= Sequences[Tree[hill],:]
                
                
                freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]
                
                Frequencies[CHR][c][hill]= freq_vector
            
            grid.fit(data[Tree[hill],:])
            
            # use the best estimator to compute the kernel density estimate
            kde = grid.best_estimator_
            
            P_dist = kde.score_samples(data[Tree[hill],:])
            Dist = kde.score_samples(data)
            P_dist= np.nan_to_num(P_dist)
            Dist= np.nan_to_num(Dist)
            if np.std(P_dist) == 0:
                Dist= np.array([int(Dist[x] in P_dist) for x in range(len(Dist))])
            else:
                Dist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Dist)
            Dist= np.nan_to_num(Dist)
            Construct[CHR][c][hill] = Dist
        
            ######################################### 
        ############# AMOVA ################
            #########################################
        if args.supervised:
            labels= Sup_labels
            Who= [z for z in it.chain(*[refs_lib[x] for x in ref_order])]
            Ngps= len(refs_lib)
            
            for hill in refs_lib.keys():
                
                if len(refs_lib[hill]) >= args.cl_freqs:
                    cl_seqs= Sequences[refs_lib[hill],:]
                    
                    freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]
                    
                    Frequencies[CHR][c][hill]= freq_vector
        else:
            Who = [x for x in range(len(labels)) if labels[x] != -1 and labels[x] in Keep]
            labels = [labels[x] for x in Who]
            Who= [Focus_labels[x] for x in Who]
            for hill in SpaceX.keys():
                
                if len(Tree[hill]) >= args.cl_freqs:
                    cl_seqs= Sequences[Tree[hill],:]
                    
                    
                    freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]
                    
                    Frequencies[CHR][c][hill]= freq_vector
        
        if len(list(set(labels))) == 1:
            Results[CHR][c]= [0,1]
            continue
        
        if args.amova:
            AMOVA,Cig = AMOVA_FM42(data[Who,:],labels,n_boot=0,metric= 'euclidean')
            print('couting: {}, Ngps: {}'.format(AMOVA,Ngps))
            Results[CHR][c] = [AMOVA,Ngps]

print("MS PCA.")

print('print..')
if Home:
    Home += '/'

##
filename= Home + args.id + '_KDE.txt'
os.makedirs(os.path.dirname(filename), exist_ok=True)
Output= open(filename,'w')

for CHR in Frequencies.keys():
    for bl in Frequencies[CHR].keys():
        for hill in Frequencies[CHR][bl].keys():
            Output.write('\t'.join([str(x) for x in [CHR,bl,hill]]) + '\t')
            Output.write('\t'.join([str(x) for x in Construct[CHR][bl][hill]]))
            Output.write('\n')

Output.close()

###
###

filename= Home + args.id + '_freqs.txt'
os.makedirs(os.path.dirname(filename), exist_ok=True)
Output= open(filename,'w')

for CHR in Frequencies.keys():
    for bl in Frequencies[CHR].keys():
        for hill in Frequencies[CHR][bl].keys():
            Output.write('\t'.join([str(x) for x in [CHR,bl,hill]]) + '\t')
            Output.write('\t'.join([str(x) for x in Frequencies[CHR][bl][hill]]))
            Output.write('\n')

Output.close()



##
## AMOVA

if args.amova:
    filename= Home + args.id + '_AMOVA.txt'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    Output= open(filename,'w')

    Output.write('\t'.join(['CHR','start','end','id','AMOVA','n_clusters']))
    Output.write('\n')

    for CHR in Results.keys():
        for bl in Results[CHR].keys():
            line= [CHR,bl,Genes[CHR][bl][0],Genes[CHR][bl][1]]
            line.extend(Results[CHR][bl])
            print(line)
            
            Output.write('\t'.join([str(x) for x in line]))
            Output.write('\n')

    Output.close()

print("Done.")