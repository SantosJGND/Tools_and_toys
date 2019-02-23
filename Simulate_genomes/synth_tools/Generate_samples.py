import numpy as np
import itertools as it

import plotly.graph_objs as go
import plotly
from plotly.offline import download_plotlyjs, plot, iplot

from sklearn.decomposition import PCA
import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


import synth_tools.StructE_tools as Ste


###

def OriginbySNMF(Geno_Q,t):
    """
    Classes individuals according to Group assignment by SNMF
    using user provided threshold (.8 advised). returns dict.
    """
    Geneo = open(Geno_Q,"r")
    Ind = 0
    Inds= {}
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
        Inds[Ind]= line
    return Groups,Inds


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



def Gen_samples(Pops,Sizes,vector_lib,prior_func,prior_kwargs,return_pca= False,n_comp= 100,prior= 'sinusoid',range_diff= [0,100],steps= 100):
    
    print('...')
    
    Npops= len(Sizes)
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(vector_lib)
    features = pca.transform(vector_lib)
    
    target= [0,1]

    Windows= recursively_default_dict()
    
    threshold= .005
    P= 30

    Fst_labels= []

    Fst_crawl= []

    Fst_windows= recursively_default_dict()
    
    
    if len(Sizes) > Npops:
        print('Size vector longer than N pops. using first {}.'.format(Npops))
        Sizes= Sizes[:Npops]

    for angle in np.arange(range_diff[0],range_diff[1]):
        coords= features[Pops,:]
        vector2= coords[target[0]] - coords[target[1]]
        
        coords= prior_func(coords,target,vector2,angle,**prior_kwargs)
        
        new_freqs= pca.inverse_transform(coords[:Npops])
        
        N_pops= len(Sizes)
        

        data= []

        for k in range(N_pops):

            probs= new_freqs[k,:]

            probs[(probs > 1)]= 1
            probs[(probs < 0)]= 0
            m= Sizes[k]
            Haps= [[np.random.choice([1,0],p= [1-probs[x],probs[x]]) for x in range(vector_lib.shape[1])] for acc in range(m)]

            data.extend(Haps)

        data= np.array(data)

        if return_pca:
                pca2 = PCA(n_components=n_comp, whiten=False,svd_solver='randomized')

                data= pca2.fit_transform(data)

        Windows[int(angle)]= data

        ##### FSTs
        Pairwise= Ste.return_fsts2(new_freqs)
        Pairwise['angle']= [angle] * Pairwise.shape[0]



        Fst_windows[int(angle)]= Pairwise

        ### store stuff.
        Windows[int(angle)]= data
    
    
    Windows= {1:Windows}
    Fst_windows= {1:Fst_windows}
    print('Done.')
    
    return Windows, Fst_windows


def Gen_samples_II(Pops,Sizes,vector_lib,label_package,Origins,prior_func,prior_kwargs,Cop_choice,window_size= 5000,Chr= 1,return_pca= False,n_comp= 100,range_diff= [0,100],steps= 100):
    
    labels= label_package['labels']
    Whose= label_package['Whose']
    ind_to_group= label_package['ind_to_group']
    
    print('...')
    color_ref= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']
    
    
    Npops= len(Sizes)
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(vector_lib)
    features = pca.transform(vector_lib)
    
    Windows= recursively_default_dict()
    Out= {Chr:{}}
    
    threshold= .005
    P= 30

    Fst_labels= []

    Fst_crawl= []

    Fst_windows= recursively_default_dict()
    
    Ideo= []
    current= recursively_default_dict()
    
    
    if len(Sizes) > Npops:
        print('Size vector longer than N pops. using first {}.'.format(Npops))
        Sizes= Sizes[:Npops]

    for angle in np.arange(range_diff[0],range_diff[1]):
        
        bl= int(angle*window_size)
        end= int(bl+ window_size - 1)
        Out[1][bl]= end
        
        Cop_local= {}
        
        for gp in Cop_choice.keys():
            cop_func= Cop_choice[gp]['cop_func']
            
            cop, ID= cop_func(angle,range_diff,**Cop_choice[gp]['cop_kwargs'])
            Cop_local[gp]= cop
            
        
        coords= features[Pops,:]
        
        coords= prior_func(coords,angle,range_diff,**prior_kwargs)
        
        new_freqs= pca.inverse_transform(coords[:Npops]).T
        np.random.shuffle(new_freqs)
        new_freqs= new_freqs.T
        
        N_pops= len(Sizes)
        
        data= []

        for acc in range(len(Whose)):
            Subject = 'sample' + str(acc)
            
            transition_p= Origins[ind_to_group[acc][0]][ind_to_group[acc][1]]
            COp= Cop_local[ind_to_group[acc][0]]
            
            if acc in current.keys():
                cross_over= np.random.choice([0,1], p=[1-COp,COp])
                if cross_over == 1:
                    k= np.random.choice(labels, p=transition_p)
                    current[acc]= k
                else:
                    k= current[acc]
            else:
                
                k= np.random.choice(labels, p=transition_p) #ind_to_group[acc][0]
                current[acc]= k
            
            probs= new_freqs[k,:]

            probs[(probs > 1)]= 1
            probs[(probs < 0)]= 0

            Haps= [np.random.choice([1,0],p= [1-probs[x],probs[x]]) for x in range(vector_lib.shape[1])]

            Stock = ['Region_chr'+str(Chr)+ '_' + Subject,bl,end,color_ref[k]]
            Ideo.append(Stock)
            data.append(Haps)

        data= np.array(data)
        
        if return_pca:
                pca2 = PCA(n_components=n_comp, whiten=False,svd_solver='randomized')

                data= pca2.fit_transform(data)

        Windows[int(angle)]= data

        ##### FSTs
        Pairwise= Ste.return_fsts2(new_freqs)
        Pairwise['angle']= [angle] * Pairwise.shape[0]
        Fst_labels.extend(Pairwise.pops)

        Fst_crawl.extend(Pairwise.fst)



        Fst_windows[int(angle)]= Pairwise

        ### store stuff.
        Windows[int(angle)]= data
    
    
    Windows= {Chr:Windows}
    Fst_windows= {1:Fst_windows}
    print('Done.')
    
    return Windows, Fst_windows, Ideo, Out


def Check_Path(Npops,vector_lib,prior_func,prior_kwargs,Pops= [], random= True,n_comp= 100,range_diff= [0,100],steps= 100):
    
    pca = PCA(n_components=100, whiten=False,svd_solver='randomized').fit(vector_lib)
    features = pca.transform(vector_lib)# * pca.explained_variance_ratio_
    
    
    if random:
        Pops= np.random.choice(vector_lib.shape[0],Npops,replace= False)
    
    
    Windows= recursively_default_dict()
    
    Fst_labels= []
    
    Fst_crawl= []
    
    Fst_windows= recursively_default_dict()
    
    
    for angle in np.arange(range_diff[0],range_diff[1]):
        
        coords= features[Pops,:]
        
        coords, prior= prior_func(coords,angle,range_diff,passport= True,**prior_kwargs)
        
        new_freqs= pca.inverse_transform(coords)
        new_freqs[new_freqs > 1] = 1
        new_freqs[new_freqs < 0] = 0
        
        Pairwise= Ste.return_fsts2(new_freqs)
        Pairwise['angle']= [angle] * Pairwise.shape[0]
        #
        Fst_labels.extend(Pairwise.pops)

        Fst_crawl.extend(Pairwise.fst)
        

        Fst_windows[int(angle)]= Pairwise

        ### store stuff
    
    Fst_windows= {1:Fst_windows}

    fig_data= [go.Scatter(
    x= [x for x in Fst_windows[1].keys()],
    y= [Fst_windows[1][x].fst[i] for x in Fst_windows[1].keys()],
    mode= 'markers',
    name= '{}'.format([x for x in it.combinations(range(Npops),2)][i])
    ) for i in range(len([x for x in it.combinations(range(Npops),2)]))
    ]

    layout = go.Layout(
        title= 'Fst across sets. prior: {}'.format(prior),
        yaxis=dict(
            title='fst',
            range= [0,.5]),
        xaxis=dict(
            title='Proxy genome position')
    )

    fig= go.Figure(data=fig_data, layout=layout)
    
    if random:
        return fig, Pops, prior
    else: return fig, prior



def plot_GenFst(Fst_lib,Npops,Chr):
    
    fig_data= [go.Scatter(
    x= [x for x in Fst_lib[1].keys()],
    y= [Fst_lib[1][x].fst[i] for x in Fst_lib[1].keys()],
    mode= 'markers',
    name= '{}'.format([x for x in it.combinations(range(Npops),2)][i])
    ) for i in range(len([x for x in it.combinations(range(Npops),2)]))
    ]

    layout = go.Layout(
        title= 'Fst vs. distance in vector feature space',
        yaxis=dict(
            title='fsts',
            range= [0,.5]),
        xaxis=dict(
            title='eucledian distance in feature space')
    )

    fig= go.Figure(data=fig_data, layout=layout)

    iplot(fig)


