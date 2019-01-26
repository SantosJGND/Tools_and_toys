

import numpy as np
import plotly.graph_objs as go


def Check_cop(Cop_func,Cop_kwargs,range_diff= [-10,10],steps= 100):
    
    Cop_windows= {}
    
    for angle in np.linspace(range_diff[0],range_diff[1],steps):
        
        Cop, ID= Cop_func(angle,range_diff,**Cop_kwargs)
        
        Cop_windows[angle]= Cop
    
    
    probs= [x for x in Cop_windows.keys()]
    winds= [Cop_windows[x] for x in Cop_windows.keys()]
    
    fig_data= [go.Scatter(
    x= probs,
    y= winds,
    mode= 'markers'
    )]
    
    return fig_data



#######
#######

import scipy
import numpy as np
import math
from sklearn.neighbors import KernelDensity


def Rec_prob_uniform(angle,range_windows,Cop):
    '''
    Uniform recombination across windows.
    - Cop: recombination rate.
    '''
    ID= 'uniform'
    return Cop, ID

def Rec_prob_sinusoid(angle,range_windows,freq= 1,c = 0,multiplier= 1):
    '''
    Sinusoid recombination rate. Takes (sin + 1) of progress along window range given.
    - Cop_range: multiply (sin(x) + 1) by a given factor. 
    '''
    ID= 'sinusoid'
    
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    #print(progress)
    
    Cop= math.sin(math.pi * freq * progress + c) * multiplier / 2 + multiplier / 2
    
    return Cop, ID

def Rec_prob_modal(angle,range_windows,modes,multiplier,N= 100,bandwidth= .2):
    '''
    Modal recombination rate. Select location of each mode in proportion to the range given ([0,1]).
    Distribution is constructed using KDE of values given.
    - modes: position of modes. vector. 
    - N= 100: number of samples for KDE estimation.
    - bandwidth: bandwidth of KDE estimation.
    - multiplier: multiply by a given factor.
    '''
    ID= 'modal. modes: {}'.format(modes)
    N_modes= len(modes)
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    sam= np.repeat(modes,N)
    sam= np.array(sam).reshape(-1,1)
    
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(sam)
    
    log_dens = kde.score_samples(progress)
    
    Cop= np.exp(log_dens)[0] * multiplier
    #print(Cop)
    return Cop, ID

def Rec_prob_region(angle,range_windows,region,prob_basal,prob):
    '''
    Range specific probability of recombination. Select Prob rec in and outside of range.
    - region: range. list.
    - prob_basal: Prob rec outside range.
    - prob: Prob rec inside range.
    '''
    ID= 'region'
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        Cop= prob
    else:
        Cop= prob_basal
    
    return Cop, ID
    

def Rec_prob_rdist(angle,range_windows,multiplier,c,loc,scale):
    '''
    Rdist distribution of recombination probability using range provided as [0,1]. scipy.stats.rdist.pdf used.
    - c, loc and scale: scipy rdist shape, mean and sd parameters.
    - multiplier: multiply pdf by given factor.
    '''
    ID= 'rdist'
    
    #Size= (range_windows[1] - range_windows[0])
    
    progress= (angle - range_windows[0]) / ((range_windows[1] - range_windows[0]) / 2) - 1 
    #print(progress)
    
    Cop= multiplier * scipy.stats.rdist.pdf(progress,c,loc,scale)
    if Cop > 1:
        Cop= 1
    if Cop < 0:
        Cop= 0
    
    return Cop, ID

