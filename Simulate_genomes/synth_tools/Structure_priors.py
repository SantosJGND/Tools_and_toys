import math

### Functions to manipulate genetic structure.

### Functions to manipulate genetic structure.

def sin_prior(coords,angle,range_windows,target= [0,1],freq= 2,fst_max= 0.2,passport= False):
    '''
    Sinusoid differentiation between targets.
    - index of target pop vetors in Pops matrix. list.
    - fst_max: range of divergence pattern in fsts.
    - passport: return function ID. Boolean.
    '''
    
    vector2= coords[target[0]] - coords[target[1]]
    
    ID= 'sinusoid'
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    coords[target[0]] = coords[target[0]] - [(fst_max * 10 - 1) * math.sin(math.pi * freq * progress) * x for x in vector2]
    
    if passport:
        return coords,ID
    else:
        return coords


def linear_prior(coords,angle,range_windows,target= [0,1],region= [-5,5],slope= 1,passport= False):
    '''
    Linear differentiation between target populations.
    - index of target pop vetors in Pops matrix. list.
    - region: span of differentiation pattern in prop to range provided.
    - fst_max: range of divergence pattern in fsts.
    - passport: return function ID. Boolean.
    '''
    
    vector2= coords[target[0]] - coords[target[1]]
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    ID= 'linear'
    
    if progress >= region[0] and progress <= region[1]:
        
        coords[target[0]] = coords[target[0]] + [progress * x * slope for x in vector2]
    
    if passport:
        return coords,ID
    else:
        return coords


def introgression_prior(coords,angle,range_windows,target= [0,1],region= [0,1],passport= False):
    '''
    Use the same vector for two populations at a given range.
    - region: span of differentiation pattern in prop to range provided.
    '''
    ID= 'introgression'
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        coords[target[0]] = coords[target[1]]
    
    if passport:
        return coords,ID
    else:
        return coords
    

def alien_prior_I(coords,angle,range_windows,target= [0,1],fst= .2,region= [-5,5],passport= False):
    ID= 'alien I'
    
    vector2= coords[target[0]] - coords[target[1]]
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        coords[target[0]] = coords[target[0]] + [(10 * fst - 1) * x for x in vector2]
    else:
        coords[target[0]] = coords[target[1]]
    
    if passport:
        return coords,ID
    else:
        return coords


def alien_prior_II(coords,angle,range_windows,target= [0,1],fst= .2,region= [-5,5],passport= False):
    ID= 'alien II'
    
    vector2= coords[target[0]] - coords[target[1]]
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        coords[target[0]] = coords[target[0]] + [(10 * fst - 1) * x for x in vector2]
    
    if passport:
        return coords,ID
    else:
        return coords


def alien_prior_III(coords,angle,range_windows,target= [0,1],fst_a= 0.2,fst_b= .2,region= [-5,5],passport= False):
    ID= 'alien III'
    vector2= coords[target[0]] - coords[target[1]]
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    if progress >= region[0] and progress <= region[1]:
        coords[target[0]] = coords[target[0]] + [(10 * fst_b-1) * x for x in vector2]
        
    else:
        coords[target[0]] = coords[target[0]] + [(10 * fst_a-1) * x for x in vector2]
    
    if passport:
        return coords,ID
    else:
        return coords

