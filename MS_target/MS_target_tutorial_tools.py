import numpy as np 
import collections
import itertools as it

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth

def Windows_KDE_amova(SequenceStore,admx_lib,refs_lib):

    Geneo= admx_lib

    Geneo_order= list(Geneo.keys())
    ref_order= list(refs_lib.keys())

    Whose= [z for z in it.chain(*[Geneo[x] for x in Geneo_order])]
    Sup_labels= list(np.repeat(Geneo_order,[len(Geneo[x]) for x in Geneo_order]))

    ### Define parameters and libraries of analyses.

    Results = {x:recursively_default_dict() for x in SequenceStore.keys()}

    Construct = recursively_default_dict()
    PC_var= recursively_default_dict()


    for CHR in SequenceStore.keys():
        print('going on CHR: '+ str(CHR))
        for c in SequenceStore[CHR].keys():

            ### PCA and MeanShift of information from each window copied from *FM36_Galaxy.py.
            Sequences= [SequenceStore[CHR][c][x] for x in Whose]
            Sequences= np.array(Sequences) 

            Sequences= np.nan_to_num(Sequences)

            pca = PCA(n_components=KDE_comps, whiten=False,svd_solver='randomized').fit(Sequences)
            data = pca.transform(Sequences)
            PC_var[CHR][c]= [x for x in pca.explained_variance_]

            params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
            grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

            ######################################
            ####### TEST global Likelihood #######
            ######################################
            Focus_labels = [z for z in it.chain(*refs_lib.values())]

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

                grid.fit(data[Tree[hill],:])

                # use the best estimator to compute the kernel density estimate
                kde = grid.best_estimator_

                # normalize kde derived log-likelihoods, derive sample p-values
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
            if supervised:
                labels= Sup_labels
                Who= list(range(Sequences.shape[0]))
                Ngps= len(refs_lib)

            else:
                Who = [x for x in range(len(labels)) if labels[x] != -1 and labels[x] in Keep]
                labels = [labels[x] for x in Who]
                Who= [Focus_labels[x] for x in Who]

            if amova:
                clear_output()
                Bool_set= Sequences[Who,:].astype(bool)
                print('chr {}, where: {}, supervised: {}, n clusters: {}'.format(CHR,c,str(supervised),Ngps))
                Amova1,Ciggy= AMOVA_FM42(Bool_set,labels,n_boot=0,metric= 'jaccard')
                Amova2,Ciggy= AMOVA_FM42(data[Who,:],labels,n_boot=0,metric= 'euclidean')
                Amova3,Ciggy= AMOVA_FM42(Bool_set,labels,n_boot=0,metric= 'hamming')
                print('old: ; jaccard: {}; PCA euc: {}; nHam: {}'.format(Amova1,Amova2,Amova3))
                Results[CHR][c] = [Ngps,Amova1,Amova2,Amova3]
    
    return Results, Construct, PC_var


def KDE_pcaPlot(gp):
    from plotly import tools
    import plotly.graph_objs as go
    from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
    
    fig_pca_subplots = tools.make_subplots(rows=2, cols=2)
    
    for subp in range(4):

        n_plot= subp
        
        gradient= Cameo[:,gp-1]
        
        print(gradient.shape)
        
        if subp >= 2:
            subp -= 2

            trace= go.Scatter(
            x = feats[:,0],
            y = feats[:,subp + 1],
            mode= "markers",
            marker= {
                'color': gradient,
                'colorscale':'Viridis',
                'line': {'width': 0},
                'size': 6,
                'symbol': 'circle',
                "opacity": .6
            })

            fig_pca_subplots.append_trace(trace, int(n_plot/float(2)) + 1, subp + 1)

        else:
            coords= {z:[x for x in Subset if label_vector[x] == z] for z in list(set(label_vector))}
            Col_vec= ['red','yellow','blue','black','green']
            for i in coords.keys():
                if coords[i]:
                    trace= go.Scatter(
                    x = feats[coords[i],0],
                    y = feats[coords[i],subp + 1],
                    mode= "markers",
                    name= str(i),
                    marker= {
                    'color': Col_vec[i],
                    'line': {'width': 0},
                    'size': 6,
                    'symbol': 'circle',
                    "opacity": .8})

                    fig_pca_subplots.append_trace(trace, int(n_plot/float(2)) + 1, subp + 1)

        fig_pca_subplots['layout']['yaxis' + str(n_plot + 1)].update(title='PC{}'.format(subp + 2))
        fig_pca_subplots['layout']['xaxis' + str(n_plot + 1)].update(title='PC1')
    
    layout = go.Layout()
    
    fig= go.Figure(data=fig_pca_subplots, layout=layout)
    iplot(fig)