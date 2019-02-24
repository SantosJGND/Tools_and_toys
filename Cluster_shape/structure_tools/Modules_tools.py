import numpy as np
import pandas as pd
import itertools as it

import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

import re
import matplotlib.pyplot as plt

from matplotlib.collections import BrokenBarHCollection

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


###################################################################################
####  Load data  #######################################################################
## ##############


def read_refs(index_file):
    indxs = recursively_default_dict()
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs[int(line[0])][line[1]] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())]



def read_Darwin(darwin_file):
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
        #art= [[x,np.nan][int(x == 9)] for x in art]
        
        if len(art) < Nsnps:
                        art.extend([0] * (Nsnps - len(art)))
        
        art= [x / 2 for x in art]
        gen.append(art)
        d += 1
    Input.close()
    
    
    return gen, Names


def OriginbySNMF(Geno_Q,t):
    """
    Classes individuals according to Group assignment by SNMF
    using user provided threshold (.8 advised). returns dict.
    """
    Geneo = open(Geno_Q,"r")
    Ind = 0
    Groups = recursively_default_dict()
    for line in Geneo:
        line= line.split()
        line= [float(x.strip("\n")) for x in line]
        line= [x / sum(line) for x in line]
        if Ind == 0:
            Groups = {x:[] for x in range(len(line) + 1)}
        bagged = 0
        for value in range(len(line)):
            if line[value] >= t:
                Groups[value].append(line)
                bagged += 1
        if bagged == 0:
            Groups[len(line)].append(line)
        Ind += 1
    return Groups



##########################################################################################
### Fst
### Calculate pairwise Fst based on frequency vectors selected.
### return total Fst
def return_fsts(vector_lib,pops):
    
    H= {pop: [1-(vector_lib[pop,x]**2 + (1 - vector_lib[pop,x])**2) for x in range(vector_lib.shape[1])] for pop in pops}
    Store= []
    for comb in it.combinations(pops,2):
        P= [sum([vector_lib[x,i] for x in comb]) / len(comb) for i in range(vector_lib.shape[1])]
        HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
        Fst= np.mean([(HT[x] - np.mean([H[p][x] for p in comb])) / HT[x] for x in range(len(P))])
        
        Store.append([comb,Fst])
    
    ### total fst:
    P= [sum([vector_lib[x,i] for x in pops]) / len(pops) for i in range(vector_lib.shape[1])]
    HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
    FST= np.mean([(HT[x] - np.mean([H[p][x] for p in pops])) / HT[x] for x in range(len(P))])
    
    return pd.DataFrame(Store,columns= ['pops','fst']),FST


### Dont return global Fst
def return_fsts2(freq_array):
    pops= range(freq_array.shape[0])
    H= {pop: [1-(freq_array[pop,x]**2 + (1 - freq_array[pop,x])**2) for x in range(freq_array.shape[1])] for pop in range(freq_array.shape[0])}
    Store= []

    for comb in it.combinations(H.keys(),2):
        P= [sum([freq_array[x,i] for x in comb]) / len(comb) for i in range(freq_array.shape[1])]
        HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
        per_locus_fst= [[(HT[x] - np.mean([H[p][x] for p in comb])) / HT[x],0][int(HT[x] == 0)] for x in range(len(P))]
        per_locus_fst= np.nan_to_num(per_locus_fst)
        Fst= np.mean(per_locus_fst)

        Store.append([comb,Fst])
    
    return pd.DataFrame(Store,columns= ['pops','fst'])



###################################################################""""
### Local sampling correct (Notebook 8.)


def local_sampling_correct(data_now,n_comp):
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data_now)
    feats= pca.transform(data_now)
    
    N= 50
    bandwidth = estimate_bandwidth(feats, quantile=0.15)
    params = {'bandwidth': np.linspace(np.min(feats), np.max(feats),30)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    ## perform MeanShift clustering.
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=False, min_bin_freq=5)
    ms.fit(feats)
    labels1 = ms.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}

    ## Extract the KDE of each cluster identified by MS.
    Proxy_data= []

    for lab in label_select.keys():
        if len(label_select[lab]) < 3:
            continue
            
        Quanted_set= feats[label_select[lab],:]
        grid.fit(Quanted_set)

        kde = grid.best_estimator_
        Extract= kde.sample(N)
        Return= pca.inverse_transform(Extract)
        
        Proxy_data.extend(Return)
    
    Proxy_data= np.array(Proxy_data)
    
    print([len(x) for x in label_select.values()])
    pca2 = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(Proxy_data)
    var_comp= pca2.explained_variance_ratio_
    
    New_features= pca2.transform(data_now)# * var_comp
    return New_features, var_comp


################################################
### Data generation

#### Union - calculate 3D overlap from KDE estimation
#### Function to calculate overlap given coordinates matrix and dictionary of indicies.
def extract_profiles_union(global_data,target_ind_dict,threshold,P):
    ## estimate the bandwith
    params = {'bandwidth': np.linspace(np.min(global_data), np.max(global_data),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

    ## perform MeanShift clustering.
    combine= {}
    for bull in target_ind_dict.keys():
        grid.fit(global_data[target_ind_dict[bull],:])
        combine[bull]= grid.best_estimator_    

    Stats= recursively_default_dict()

    for combo in it.combinations(target_ind_dict.keys(),2):
        pop1= combo[0]
        pop2= combo[1]

        All_coords= [x for x in it.chain(*[target_ind_dict[z] for z in combo])]

        Quanted_set= global_data[All_coords,:]

        i_coords, j_coords, z_coords = np.meshgrid(np.linspace(min(Quanted_set[:,0]),max(Quanted_set[:,0]),P),
                              np.linspace(min(Quanted_set[:,1]),max(Quanted_set[:,1]),P),
                                np.linspace(min(Quanted_set[:,2]),max(Quanted_set[:,2]),P), indexing= 'ij')


        traces= [x for x in it.product(range(P),range(P),range(P))]

        background= np.array([i_coords,j_coords,z_coords])

        background= [background[:,c[0],c[1],c[2]] for c in traces]

        background=np.array(background)

        pop1_fist= combine[pop1].score_samples(background)
        #pop1_fist= np.exp(pop1_fist)
        P_dist_pop1= combine[pop1].score_samples(global_data[target_ind_dict[pop1],:])
        pop1_fist = scipy.stats.norm(np.mean(P_dist_pop1),np.std(P_dist_pop1)).cdf(pop1_fist)
        pop1_fist= [int(x >= threshold) for x in pop1_fist]
        
        pop2_fist= combine[pop2].score_samples(background)
        #pop2_fist= np.exp(pop2_fist)
        P_dist_pop2= combine[pop2].score_samples(global_data[target_ind_dict[pop2],:])
        pop2_fist = scipy.stats.norm(np.mean(P_dist_pop2),np.std(P_dist_pop2)).cdf(pop2_fist)
        pop2_fist= [int(x >= threshold) for x in pop2_fist]

        
        pop1_and_2= len([x for x in range(background.shape[0]) if pop1_fist[x] == 1 and pop2_fist[x] == 1])
        pop1_I_pop2= pop1_and_2 / float(sum(pop1_fist))
        pop2_I_pop1= pop1_and_2 / float(sum(pop2_fist))
        
        total_overlap= pop1_and_2 / float(sum(pop1_fist) + sum(pop2_fist) - pop1_and_2)
        
        empty_space= 1 - (sum(pop1_fist) + sum(pop2_fist) - pop1_and_2) / background.shape[0]
        
        Stats[combo][pop1]= pop1_I_pop2
        Stats[combo][pop2]= pop2_I_pop1
        Stats[combo]['empty']= empty_space
        Stats[combo]['PU']= total_overlap
        
    
    return Stats


#### reference KDE

def extract_profiles(global_data,target_ind_dict):
    ## estimate the bandwith
    params = {'bandwidth': np.linspace(np.min(global_data), np.max(global_data),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    cluster_profiles= {x:[] for x in target_ind_dict.keys()}
    
    
    ## perform MeanShift clustering.
    combine= {}
    for bull in target_ind_dict.keys():
        Quanted_set= global_data[target_ind_dict[bull],:]
        grid.fit(Quanted_set)
        kde = grid.best_estimator_

        P_dist = kde.score_samples(Quanted_set)
        Fist = kde.score_samples(global_data)

        ## Normalizing log-likelihood estimates by those of the reference set and extracting their cdf.
        Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
        cluster_profiles[bull].append(Fist)

    
    return cluster_profiles


def extract_profiles_class(global_data,target_ind_dict):
    '''
    copy of the previous function. change of name to deal with local 
    function similarity.
    '''
    ## estimate the bandwith
    cluster_profiles= recursively_default_dict()
    params = {'bandwidth': np.linspace(np.min(global_data), np.max(global_data),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    combine= {}
    for bull in target_ind_dict.keys():
        
        Quanted_set= global_data[target_ind_dict[bull],:]
        grid.fit(Quanted_set)
        kde = grid.best_estimator_
        
        P_dist = kde.score_samples(Quanted_set)
        Fist = kde.score_samples(global_data)

        ## Normalizing log-likelihood estimates by those of the reference set and extracting their cdf.
        Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
        cluster_profiles[bull]=Fist

    
    return cluster_profiles
    

    
    
####################################################
### Classification of KDE collections

def Merge_class(Ref_profiles,focus_indicies,Out,Diff_threshold,X_threshold):
    Blocks_genome = recursively_default_dict()
    
    for CHR in Ref_profiles.keys():
        
        Points = sorted(Out[CHR].keys())
        Likes = Ref_profiles[CHR]
        N_pops= len(Likes[[x for x in Likes.keys()][0]])
        Pop_labels= Likes[[x for x in Likes.keys()][0]].keys()
        
        Likes = {x:[Likes[bl][x] for bl in sorted(Likes.keys())] for x in Pop_labels}
        Likes = {x:np.array([y[0] for y in Likes[x]]) for x in Likes.keys()}

        Topo = []
        
        #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
        #range_Crossed = [x for x in range(Aro.shape[0])]
        
        for acc in focus_indicies:
            Guys = np.array([Likes[x][:,acc] for x in Pop_labels])
            Guys = np.nan_to_num(Guys)
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            
            Test = [int(x <= X_threshold) for x in np.amax(np.array(Guys),axis = 0)]     
            
            #
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
            
            #
            Consex = [x for x in it.combinations(range(N_pops),2)]
            if Consex:
                for h in range(len(maxim)):
                    CL = []
                    for j in Consex:
                        Diff = Guys[h,j]
                        if maxim[h] not in j or len([x for x in Diff if x < X_threshold]) > 0:
                            continue
                        if max(Diff) <= X_threshold:
                            Diff = 0
                        else:
#                            Diff = int(len([x for x in Diff if x <= Diff_threshold]) == 1)
                            Diff = abs(max(Diff)) / abs(min(Diff))
                            Diff = int(Diff > Diff_threshold)
                        
                        if Diff == 0:
                            CL.append(j)
                    
                    if len(CL) == 2:
                        maxim[h] = 7
                    if len(CL) == 1:
                        maxim[h] = sum(CL[0]) + N_pops
            
            maxim[where_X] = N_pops
            
            if not Consex:
                for h in range(len(maxim)):
                    #maxim[h] = int(10*Guys[h,0])
                    maxim[h]= int(Guys[h,0] <= X_threshold)
            
            
            Topo.append(maxim + 1)
        
        
        Topo = np.array(Topo).T
        
        Clove = {CHR:{Points[x]:Topo[x,] for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome


def Merge_class2(Ref_profiles,target_indx,focus_indicies,Out,Diff_threshold,X_threshold):
    Blocks_genome = recursively_default_dict()
    
    for CHR in Ref_profiles.keys():
        print(CHR)
        Points = sorted(Out[CHR].keys())
        Likes = Ref_profiles[CHR]
        
        N_pops= len(Likes[[x for x in Likes.keys()][0]])
        Pop_labels= Likes[[x for x in Likes.keys()][0]].keys()
        print("number of reference populations: {0}".format(N_pops))
        #Likes = {x:[Likes[bl][x] for bl in sorted(Likes.keys())] for x in Pop_labels}
        #Likes = {x:np.array([y[0] for y in Likes[x]]) for x in Likes.keys()}

        Topo = []
        
        #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
        #range_Crossed = [x for x in range(Aro.shape[0])]
        
        for bl in Likes.keys():
            maxim= []
            
            
            for acc in focus_indicies:
                Test= [Likes[bl][x][0][acc] for x in sorted(Likes[bl].keys())]
                
                Guys= [Likes[bl][x][0][acc] for x in sorted(Likes[bl].keys())]
                
                if max(Test) <= X_threshold:
                    maxim.append(N_pops)
                else:
                    Consex = [x for x in it.combinations(range(N_pops),2)]
                    CL = []
                    for j in Consex:
                        Diff = [Guys[i] for i in j]
                        
                        if Guys.index(max(Guys)) not in j or len([x for x in Diff if x < X_threshold]) > 0:
                            continue
                        if max(Diff) <= X_threshold:
                            Diff = 0
                        else:
#                            Diff = int(len([x for x in Diff if x <= Diff_threshold]) == 1)
                            Diff = abs(max(Diff)) / abs(min(Diff))
                            Diff = int(Diff > Diff_threshold)
                        
                        if Diff == 0:
                            CL.append(j)
                    
                    if len(CL) == 2:
                        maxim.append(7)
                    if len(CL) == 1:
                        maxim.append(sum(CL[0]) + N_pops)
                    else:
                        maxim.append(Guys.index(max(Guys)))
            
            Topo.append([maxim[c] + 1 for c in range(len(maxim))]) 
        
        Topo = np.array(Topo)
        print(Topo.shape)
        Clove = {CHR:{Points[x]:Topo[x,] for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome


def Merge_class_mix(Ref_profiles,focus_indicies,Out,Diff_threshold,X_threshold):
    Blocks_genome = recursively_default_dict()
    Clove= recursively_default_dict()
    for CHR in Ref_profiles.keys():
        print(CHR)
        Points = sorted(Out[CHR].keys())
        Likes = Ref_profiles[CHR]
        N_pops= len(Likes[[x for x in Likes.keys()][0]])
        Pop_labels= Likes[[x for x in Likes.keys()][0]].keys()
        print("number of reference populations: {0}".format(N_pops))
        
        Topo = []
        
        for bl in Likes.keys():
            Guys = np.array([Likes[bl][x] for x in Pop_labels])
            Guys = np.nan_to_num(Guys)
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            
            Test = [int(x <= X_threshold) for x in np.amax(np.array(Guys),axis = 0)]     
            
            #
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
            
            #
            Consex = [x for x in it.combinations(range(N_pops),2)]
            if Consex:
                for h in range(len(maxim)):
                    CL = []
                    for j in Consex:
                        Diff = Guys[h,j]
                        if maxim[h] not in j or len([x for x in Diff if x < X_threshold]) > 0:
                            continue
                        if max(Diff) <= X_threshold:
                            Diff = 0
                        else:
                            Diff = abs(max(Diff)) / abs(min(Diff))
                            Diff = int(Diff > Diff_threshold)
                        
                        if Diff == 0:
                            CL.append(j)
                    
                    if len(CL) == 2:
                        maxim[h] = 7
                    if len(CL) == 1:
                        maxim[h] = sum(CL[0]) + N_pops
            
            maxim[where_X] = N_pops
            
            if not Consex:
                for h in range(len(maxim)):
                    maxim[h] = int(10*Guys[h,0])    
            
            Clove[CHR][bl]= maxim + 1
            Topo.append(maxim + 1)
        
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome


#########################################################
### Ideogram Processing

def compress_ideo(df,Out,chromosome_list):
    
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        
        Chr = int(re.search('Region_(.+?)_',chromosome_list[CHR]).group(1))
        sub = df[df.chrom == chromosome_list[CHR]]
        Coordinates = sorted(sub.start)
        Size = sub.shape[0]
        start = min(df.start)
        First = sub.gieStain.iloc[0]
        for index in range(len(Coordinates)):
            row = sub[sub.start == Coordinates[index]]
            if index == 0:
                continue
            if index == (Size - 1):
                if row.gieStain.iloc[0] == First:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
                else:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
            else:
                if row.gieStain.iloc[0] == First:
                    continue
                else:
                    new_set.append([chromosome_list[CHR],start,row.start.iloc[0]-1,First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
    
    new_set = pd.DataFrame(new_set,columns = ['chrom', 'start', 'end', 'gieStain'])
    return new_set



def compress_ideo_vII(df,Out,chromosome_list):
    
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        
        Chr = int(re.search('Region_(.+?)_',chromosome_list[CHR]).group(1))
        sub = df[df.chrom == chromosome_list[CHR]]
        Coordinates = sorted(sub.start)
        Size = sub.shape[0]
        start = min(df.start)
        First = sub.gieStain.iloc[0]
        for index in range(len(Coordinates)):
            row = sub[sub.start == Coordinates[index]]
            if index == 0:
                continue
            if index == (Size - 1):
                if row.gieStain.iloc[0] == First:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
                else:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
            else:
                if row.gieStain.iloc[0] == First:
                    continue
                else:
                    new_set.append([chromosome_list[CHR],start,row.start.iloc[0]-1,First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
    
    new_set = pd.DataFrame(new_set,columns = ['chrom', 'start', 'end', 'gieStain'])
    return new_set



# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def return_ideogram(ideo, chromosome_list, Comparison_threshold, Outlier_threshold, out= True):
    # Height of each ideogram
    chrom_height = 1

    # Spacing between consecutive ideograms
    chrom_spacing = 0

    # Height of the gene track. Should be smaller than `chrom_spacing` in order to
    # fit correctly
    gene_height = 0.0

    # Padding between the top of a gene track and its corresponding ideogram
    gene_padding = 0.0


    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing



    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing
    

    # Colors for different chromosome stains
    color_lookup = {
        'red': [255, 0, 0],
        'yellow': [255, 255, 0],
        'blue': [0, 0, 255],
        'orange': [255, 165, 0],
        'green': [50, 205, 50],
        'black': [0, 0, 0],
        'purple': [128, 0, 128],
        'silver': [211, 211, 211],
    }

    # Add a new column for colors
    
    ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Width, height (in inches)
    figsize = (10, 30)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the ideogram data...
    print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_xticks([x for x in range(min(ideo.start),max(ideo.end),int(10000))])
    ax.set_xticklabels([round(x / float(10000),3) for x in range(min(ideo.start),max(ideo.end),int(10000))])
    plt.xticks(fontsize = 5,rotation = 90)
    ax.tick_params(axis = 'x',pad = 10)

    ax.tick_params(axis='y', which='major', pad=30)
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list, fontsize = 5)
    ax.axis('tight')
    if out == True:
        plt.savefig('Ideo_step_' + '_OutlierTh' + str(Outlier_threshold) + '_Z' +str(Comparison_threshold)+ '.png',bbox_inches = 'tight')
    return fig


