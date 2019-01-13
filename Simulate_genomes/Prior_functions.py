from math import sin

### Functions to manipulate genetic structure.

def sin_prior(coords,target,vector2,angle,fst_max= 0.2):
    
    coords[target[0]] = coords[target[0]] - [(fst_max * 10 - 1) * sin(angle) * x for x in vector2]
    
    return coords


def linear_prior(coords,target,vector2,angle,region= [-5,5],slope= 1):
    
    if angle >= region[0] and angle <= region[1]:
        progression= abs(angle - region[0]) / (region[1] - region[0])
        coords[target[0]] = coords[target[0]] + [progression * x * slope for x in vector2]
    
    return coords


def introgression_prior(coords,target,vector2,angle, region):
    
    if angle >= region[0] and angle <= region[1]:
        coords[target[0]] = coords[1]
    
    return coords


def alien_prior_I(coords,target,vector2,angle,fst= .2,region= [-5,5]):
    
    if angle >= region[0] and angle <= region[1]:
        coords[target[0]] = coords[target[0]] + [(10 * fst - 1) for x in vector2]
        #coords[target[0]] = coords[target[0]] + [sin(angle) * x for x in vector2] remember to try this though
        #coords[target[0]] = coords[len(coords) - 1]
    else:
        coords[target[0]] = coords[target[1]]
    
    return coords


def alien_prior_II(coords,target,vector2,angle,fst= .2,region= [-5,5]):
    
    if angle >= region[0] and angle <= region[1]:
        coords[target[0]] = coords[target[0]] + [(10 * fst - 1) for x in vector2]
        #coords[target[0]] = coords[target[0]] + [sin(angle) * x for x in vector2] remember to try this though
        #coords[target[0]] = coords[len(coords) - 1]
    
    return coords



def alien_prior_III(coords,target,vector2,angle,fst_a= 0.2,fst_b= .2,region= [-5,5]):
    
    if angle >= region[0] and angle <= region[1]:
        coords[target[0]] = coords[target[0]] + [(10 * fst_b - 1) for x in vector2]
        #coords[target[0]] = coords[target[0]] + [sin(angle) * x for x in vector2] remember to try this though
        #coords[target[0]] = coords[len(coords) - 1]
    else:
        coords[target[0]] = coords[target[0]] + [(10 * fst_a - 1) for x in vector2]
    
    return coords
