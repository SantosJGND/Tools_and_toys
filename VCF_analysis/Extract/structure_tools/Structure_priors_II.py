import math
import numpy as np
from structure_tools.Generate_samples import fst_select

### Functions to manipulate genetic structure.

def sin_prior(vector_lib,Fsts_test,angle,range_windows,target= [0,1],freq= 2,range_fst= [0,0.15],passport= False):
    '''
    Sinusoid differentiation between targets.
    - index of target pop vetors in Pops matrix. list.
    - fst_max: range of divergence pattern in fsts.
    - passport: return function ID. Boolean.
    '''
    ID= 'sinusoid'
    
    fst_middle= sum(range_fst) / float(2)
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    fst_wanted= np.sin(progress * np.pi * freq) * ((range_fst[1] - range_fst[0]) / 2) + fst_middle
    
    who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
    
    new_freqs= vector_lib[who,:]
    
    if passport:
        return new_freqs,ID
    else:
        return new_freqs


def linear_prior(vector_lib,Fsts_test,angle,range_windows,target= [0,1],region= [-5,5],slope= 1,range_fst= [0,0.15],passport= False):
    '''
    Linear differentiation between target populations.
    - region: span of differentiation pattern in prop to range provided.
    - fst_max: range of divergence pattern in fsts.
    - passport: return function ID. Boolean.
    '''
    ID= 'linear'
    
    fst_middle= sum(range_fst) / float(2)
    
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    fst_wanted= range_fst[0] + progress * (range_fst[1] - range_fst[0]) * slope
    
    who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
    
    new_freqs= vector_lib[who,:]
    
    if passport:
        return new_freqs,ID
    else:
        return new_freqs


def introgression_prior(vector_lib,Fsts_test,angle,range_windows,target= [0,1],region= [0,1],passport= False):
    '''
    Use the same vector for two populations at a given range.
    - region: span of differentiation pattern in prop to range provided.
    '''
    ID= 'introgression'
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        fst_wanted= 0
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    else:
        fst_wanted= .1
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    
    if passport:
        return new_freqs,ID
    else:
        return new_freqs
    

def alien_prior_I(vector_lib,Fsts_test,angle,range_windows,target= [0,1],fst= .2,region= [-5,5],passport= False):
    ID= 'alien I'
    
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        fst_wanted= fst
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    else:
        fst_wanted= .1
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    
    if passport:
        return new_freqs,ID
    else:
        return new_freqs


def alien_prior_II(vector_lib,Fsts_test,angle,range_windows,target= [0,1],fst= .2,region= [-5,5],passport= False):
    ID= 'alien II'
    
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        fst_wanted= fst
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    else:
        fst_wanted= 0
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    
    if passport:
        return new_freqs,ID
    else:
        return new_freqs




def alien_prior_III(vector_lib,Fsts_test,angle,range_windows,target= [0,1],fst_a= 0.2,fst_b= .2,region= [-5,5],passport= False):
    ID= 'alien III'

    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        fst_wanted= fst_a
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    else:
        fst_wanted= fst_b
        who= fst_select(Fsts_test,fst_wanted,range_allow= 0.01)
        
        new_freqs= vector_lib[who,:]
    
    if passport:
        return new_freqs,ID
    else:
        return new_freqs



