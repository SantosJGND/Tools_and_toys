from random import choices
from scipy.stats import expon
from functools import reduce
from math import factorial
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
import numpy as np
import itertools as it

def nCk(n,k):
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )


def get_time(k= 4,Nt=10):
    
    data= []
    
    while len(data) < Nt:
        um= expon.rvs(size= 1,scale= 1 / nCk(k,2))[0]
        if um > 0:
            data.append(um)
    
    return data


def sim_coales(k= 4):
    
    Tcs= []
    
    for co in range(k-1):
        
        td= get_time(k= k - co,Nt=1)
        
        Tcs.extend(td)
        
        
    return Tcs

def get_time_II(k= 4,Theta= 1,Nt=10):
    
    par= k * (k - 1 + Theta) / 2
    
    data= []
    
    while len(data) < Nt:
        um= expon.rvs(size= 1,scale= 1 / par)[0]
        if um > 0:
            data.append(um)
    
    return data


def toincoss(k= 4,Theta= 1):
    
    coal= (k - 1) / (k - 1 + Theta)
    
    probs= [1- coal, coal]
    
    head= np.random.choice([0,1],1,p= probs)
    
    return head



def SimI_net(k= 4):

    times= sim_coales(k= k)
    times_acc= [sum(times[:(x+1)]) for x in range(len(times))]
    #print(times_acc)
    #times_acc= times_acc[::-1]
    

    node_ord= np.random.choice(list(range(k)), k, replace= False)

    sim_keys= {z: [] for z in range(k)}
        
    edges= []
    times_dict= {}
    
    leaves= {z: [z] for z in range(k)}
    
    surface= list(range(k))
    #random.shuffle(surface)
    
    for cl in range(len(times)):
        
        if cl == len(times) - 1:
            new_nd= -1
        
        else:
            new_nd= cl + k
        
        pair_idx= np.random.choice(list(range(len(surface))),2,replace= False)
        
        pair= tuple([surface[x] for x in pair_idx])
        
        new_edges= [tuple([new_nd,x]) for x in pair]
        edges.extend(new_edges)
        
        surface= [surface[x] for x in range(len(surface)) if x not in pair_idx]        
        surface.append(new_nd)
        
        sim_keys[new_nd]= list(pair)
        
        ti= times_acc[cl]
        times_dict[new_nd]= ti
        
        ti= round(ti,3)
        leaves[new_nd]= ['t',ti]#['t: {}'.format(ti)]
    
    return sim_keys, times_acc, leaves, edges, times_dict


def SimIII(k= 4, Theta= 1):
    sim_keys, times_acc, leaves, edges, times_dict= SimI_net(k= k)
    
    proxy_edges= list(edges)
    
    for ed in proxy_edges:

        t0= times_dict[ed[0]]
        
        if ed[1] in times_dict.keys():
            t1= times_dict[ed[1]]
        else:
            t1= 0
        
        tt= t0 - t1

        Pexp= tt * Theta / 2

        muts= np.random.poisson(Pexp,1)[0]

        if muts > 0:
            
            ## times
            times_b= np.random.uniform(size= muts) * tt
            times_b= np.array(sorted(times_b)) + t1
            times_b= np.round(times_b,3)
            times_b= times_b[::-1]
            db= 0
            ##
            curDictL= len(sim_keys)

            new_dict= {}
            new_edges= []
            new_leaves= {}

            for z in range(curDictL,curDictL + muts - 1):
                new_dict[z]= [z+1]
                #new_leaves[z]= ['mut t: {}'.format(times_b[db])]
                new_leaves[z]= ['mut',times_b[db]]
                new_edges.append((z,z+1))
                
                db += 1
            
            new_dict.update({curDictL+muts-1: [ed[1]]})
            new_leaves[curDictL+muts-1]= ['mut',times_b[db]] #['mut t: {}'.format(times_b[db])]
            
            new_edges.append((ed[0],curDictL))
            new_edges.append((curDictL+muts-1,ed[1]))

            sim_keys[ed[0]]= [x for x in sim_keys[ed[0]] if x != ed[1]]
            sim_keys[ed[0]].append(curDictL)

            edges= [x for x in edges if x != ed]

            leaves.update(new_leaves)
            sim_keys.update(new_dict)
            edges.extend(new_edges)
    
    return sim_keys, leaves, edges


def get_PA(leaves,Ne,MRCA= 2,mut_tag= 'mut'):
    '''
    get mutations occuring since the MRCA.
    '''
    
    PA_times= [g[1] for g in leaves.values() if mut_tag in g]
    PA_times= [x for x in PA_times if x <= MRCA]
    
    return PA_times
