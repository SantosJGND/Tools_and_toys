import pandas as pd
import numpy as np
import itertools as it
import collections
import re
import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

import plotly.graph_objs as go
from plotly import tools
from plotly.offline import iplot

from IPython.display import clear_output

from structure_tools.Modules_tools import return_fsts2

from structure_tools.AMOVA_func import AMOVA_FM42, amova_cofactor


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)



def read_geno_nanum(filename, row_info= 6,header_info= 9,phased= False):

    info_summ= {}
    info_save= list(range(row_info))

    header_len= header_info
    summary= []

    Miss= recursively_default_dict()
    
    Input= open(filename)

    genotype= []
    d= 0

    for line in Input:    
        line= line.strip()

        if d in info_save:

            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split('=')
            info_summ[line[0]] = ''.join(line[1:])
            d += 1
            continue

        if d == (len(info_save)):
            print(info_summ)
            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split()

            columns= line[:header_len]
            Names= line[header_len:]

            d += 1
            continue

        if d > (len(info_save)):
            line= line.split()
            seq= []

            info= line[:header_len]
            chrom= re.search(r'\d+', line[0]).group()
            info[0]= chrom

            summary.append(info)

            for ind in range(header_len,len(line)):
                locus= line[ind]
                #print(locus)
                alleles= locus.split(':')[0]
                #print(alleles)
                if '.' in alleles:
                    alleles= ''.join([[x,'0'][int(x == '.')] for x in list(alleles)])
                alleles= list(map(int, re.findall(r'\d+', alleles)))
                if len(alleles) != 2:
                    print(alleles)
                if phased:
                    seq.extend(alleles)
                else:
                    seq.append(sum(alleles))

            genotype.append(seq)
            d += 1

    Input.close()

    summary= np.array(summary)
    summary= pd.DataFrame(summary,columns= columns)
    genotype= np.array(genotype).T

    return genotype, summary, Names


def simple_read_vcf(filename,row_info= 5,header_info= 9,phased= False):

    Input= open(filename)

    info_summ= {}
    info_save= list(range(row_info))

    phased= False
    header_len= header_info
    summary= []

    Miss= recursively_default_dict()

    genotype= []
    d= 0

    for line in Input:    
        line= line.strip()

        if d in info_save:

            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split('=')
            info_summ[line[0]] = ''.join(line[1:])
            d += 1
            continue
        
        if d == (len(info_save)):
            line= ''.join(filter(lambda ch: ch not in "#", line))
            line= line.split()

            columns= line[:header_len]

            Fam= {
                line[x]: x for x in range(header_len,len(line))
            }

            d += 1
            continue

        if d > (len(info_save)):
            line= line.split()
            seq= []

            info= line[:header_len]
            chrom= re.search(r'\d+', line[0]).group()
            info[0]= chrom

            summary.append(info)
            
            for ind in range(header_len,len(line)):
                locus= line[ind]
                #print(locus)
                alleles= locus.split(':')[0]
                #print(alleles)
                if '.' in alleles:
                    alleles= ''.join([[x,'0'][int(x == '.')] for x in list(alleles)])
                alleles= list(map(int, re.findall(r'\d+', alleles)))
                if len(alleles) != 2:
                    alleles
                if phased:
                    seq.extend(alleles)
                else:
                    seq.append(sum(alleles))

            genotype.append(seq)
            d += 1

    Input.close()

    summary= np.array(summary)
    summary= pd.DataFrame(summary,columns= columns)
    genotype= np.array(genotype).T

    return genotype, summary, info_save


def check_densities(vector_lib_2,N):
    
    who= np.random.choice(list(range(vector_lib_2.shape[0])),N,replace= False)
    
    freqs= []
    
    for pop in who:
        
        freq_vector= vector_lib_2[pop,:]

        X_plot = np.linspace(0, 1, 100)

        kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(np.array(freq_vector).reshape(-1,1))

        log_dens= kde.score_samples(X_plot.reshape(-1,1))
                        
        freqs.append(np.exp(log_dens))
    
    freqs= np.array(freqs)
    
    return freqs


def geno_window_split(genotype,summary,Steps= 25,window_size=100):

    window_starts= list(np.arange(0,genotype.shape[1],Steps))
    if window_starts[-1] != genotype.shape[1]:
        window_starts.append(genotype.shape[1])

    Windows= recursively_default_dict()
    Out= recursively_default_dict()

    lengths_winds= []

    for splyt in range(len(window_starts) - 1):
        IN= window_starts[splyt]
        OUT= window_starts[splyt] + window_size

        if OUT > genotype.shape[1]:
            OUT= genotype.shape[1] - 1
        range_window= [IN,OUT]

        lengths_winds.append(OUT-IN)

        chrom= int(summary.CHROM[range_window[0]])

        start= int(summary.POS[range_window[0]])
        end= int(summary.POS[range_window[1]])

        Windows[chrom][start]= genotype[:,range_window[0]:range_window[1]]
        Out[chrom][start]= end

        if OUT-IN < window_size:
            break

    return Windows, Out


def window_fst_sup(Windows,ref_labels,labels1,Chr= 1,ncomp= 4,range_sample= [],rand_sample= 0):
    
    kde_class_labels= labels1
    kde_label_dict= {
        z:[x for x in range(len(kde_class_labels)) if kde_class_labels[x] == z] for z in list(set(kde_class_labels))
    }
    
    if rand_sample:
        sample= rand_sample
        sample_range= [0,sample]
        Freq_extract= {
            Chr:{
                bl:Windows[Chr][bl] for bl in np.random.choice(list(Windows[Chr].keys()),sample,replace= True)
            }
        }

    if range_sample:
        sample_range= range_sample
        Freq_extract= {
            Chr:{
                bl:Windows[Chr][bl] for bl in list(sorted(Windows[Chr].keys()))[sample_range[0]:sample_range[1]]
            }
        }
    
    sim_fst= []
    
    for c in Freq_extract[Chr].keys():
        Sequences= Windows[Chr][c]

        if Sequences.shape[1] <= 3:
            Results[Chr][c] = [0,0]
            print('hi')
            continue

        Sequences= np.nan_to_num(Sequences)

        pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized').fit(Sequences)
        data = pca.transform(Sequences)
        Ngps= len(ref_labels)
        these_freqs= []
        
        for hill in ref_labels:
            cl_seqs= Sequences[kde_label_dict[hill],:]

            freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]
            these_freqs.append(freq_vector)
            
        Pairwise= return_fsts2(np.array(these_freqs))
        sim_fst.append(list(Pairwise.fst))
    
    ###
        
    return sim_fst


def window_analysis(Windows,ref_labels,labels1,Chr= 1,ncomp= 4,amova= True,supervised= True,include_who= [],
                    range_sample= [130,600],rand_sample= 0,clsize= 15,cl_freqs= 5,Bandwidth_split= 20):
    

    kde_class_labels= labels1
    kde_label_dict= {
        z:[x for x in range(len(kde_class_labels)) if kde_class_labels[x] == z] for z in list(set(kde_class_labels))
    }

    if include_who:
        include= [x for x in range(len(kde_class_labels)) if kde_class_labels[x] in include_who]
        ref_labels= include_who
        kde_class_labels= [kde_class_labels[x] for x in include]

        kde_label_dict= {
            z:[x for x in range(len(kde_class_labels)) if kde_class_labels[x] == z] for z in include_who
        }
    
    
    if rand_sample:
        sample= rand_sample
        sample_range= [0,sample]
        Freq_extract= {
            Chr:{
                bl:Windows[Chr][bl] for bl in np.random.choice(list(Windows[Chr].keys()),sample,replace= True)
            }
        }

    if range_sample:
        sample_range= range_sample
        Freq_extract= {
            Chr:{
                bl:Windows[Chr][bl] for bl in list(sorted(Windows[Chr].keys()))[sample_range[0]:sample_range[1]]
            }
        }


    Results = {
        'header': ['Chr','window'],
        'info': []
    }


    Frequencies = {
        'header': ['Chr','window','cl'],
        'coords': [],
        'info': []
    }


    Construct =  {
        'header': ['Chr','window','cl'],
        'coords': [],
        'info': []
    }


    PC_var=  {
        'header': ['Chr','window'],
        'coords': [],
        'info': []
    }


    pc_density= []
    pc_coords= []

    sim_fst= []

    for c in Freq_extract[Chr].keys():
        Sequences= Windows[Chr][c]

        if Sequences.shape[1] <= 3:
            Results[Chr][c] = [0,0]
            print('hi')
            continue

        Sequences= np.nan_to_num(Sequences)

        pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized').fit(Sequences)
        data = pca.transform(Sequences)

        if include_who:
            data= data[include,:]

        ##### PC density
        PC= 0

        pc_places= data[:,PC]

        X_plot = np.linspace(-8, 8, 100)

        kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(np.array(pc_places).reshape(-1,1))

        log_dens= kde.score_samples(X_plot.reshape(-1,1))

        pc_density.append(np.exp(log_dens))
        pc_coords.append(pc_places)

        PC_var['coords'].append([Chr,c])
        PC_var['info'].append([x for x in pca.explained_variance_])
        ###
        params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
        grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

        ######################################
        ####### TEST global Likelihood #######
        ######################################
        Focus_labels = list(range(data.shape[0]))

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

        these_freqs= []
        ### Extract MScluster likelihood by sample

        for hill in SpaceX.keys():

            if len(Tree[hill]) >= cl_freqs:
                if supervised == False:
                    print('hi')
                    cl_seqs= Sequences[Tree[hill],:]

                    freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]

                    Frequencies['coords'].append([Chr,c,hill])
                    Frequencies['info'].append(freq_vector)
                    these_freqs.append(freq_vector)

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
            Construct['coords'].append([Chr,c,hill])
            Construct['info'].append(Dist)

            ######################################### 
        ############# AMOVA ################
            #########################################

        if supervised:
            labels= [x for x in kde_class_labels if x in ref_labels]
            Who= [z for z in it.chain(*[kde_label_dict[x] for x in ref_labels])]
            Ngps= len(ref_labels)
            
            
            print(ref_labels)
            for hill in ref_labels:

                if len(kde_label_dict[hill]) >= cl_freqs:
                    if include_who:
                        Seq_specific= Sequences[include,:]
                        
                    cl_seqs= Seq_specific[kde_label_dict[hill],:]

                    freq_vector= [float(x) / (cl_seqs.shape[0] * 2) for x in np.sum(cl_seqs,axis= 0)]

                    Frequencies['coords'].append([Chr,c,hill])
                    Frequencies['info'].append(freq_vector)
                    these_freqs.append(freq_vector)

        else:
            Who = [x for x in range(len(labels)) if labels[x] != -1 and labels[x] in Keep]
            labels = [labels[x] for x in Who]
            Who= [Focus_labels[x] for x in Who]

        #
        Pairwise= return_fsts2(np.array(these_freqs))
        sim_fst.extend(Pairwise.fst)
        


        if len(list(set(labels))) == 1:
            Results['coords'].append([Chr,c])
            Results['info'].append([AMOVA,Ngps])
            continue

        if amova:
            clear_output()
            AMOVA,Cig = AMOVA_FM42(data[Who,:],labels,n_boot=0,metric= 'euclidean')
            print('counting: {}, Ngps: {}'.format(AMOVA,Ngps))
            Results['info'].append([Chr,c,AMOVA,Ngps])

    Results['info']= pd.DataFrame(np.array(Results['info']),columns= ['chrom','window','AMOVA','Ngps'])
    
    X_plot = np.linspace(0, .3, 100)

    freq_kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(np.array(sim_fst).reshape(-1,1))

    log_dens = freq_kde.score_samples(X_plot.reshape(-1,1))

    fig_roost_dens= [go.Scatter(x=X_plot, y=np.exp(log_dens), 
                                mode='lines', fill='tozeroy', name= '',
                                line=dict(color='blue', width=2))]
    ##

    layout= go.Layout(
        title= 'allele frequency distribution across clusters',
        yaxis= dict(
            title= 'density'
        ),
        xaxis= dict(
            title= 'fst'
        )
    )

    fig = go.Figure(data=fig_roost_dens, layout= layout)

    return Frequencies, sim_fst, Results, Construct, pc_density, pc_coords, fig
