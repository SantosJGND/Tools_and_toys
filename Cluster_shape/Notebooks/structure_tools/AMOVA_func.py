import scipy
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import pairwise_kernels
import numpy as np
import random

def amova_cofactor(L,allele_profiles,populations,total_populations,metric= 'hamming'):
    population_list = list(set(populations))
    Npops= len(population_list)
    
    
    
    coords= {
        z:[x for x in range(len(allele_profiles)) if populations[x] == z] for z in population_list
    }
    
    nu_total= np.triu_indices(allele_profiles.shape[0],1)
    SSTOT= pairwise_distances(allele_profiles, metric=metric)
    SSTOT= SSTOT[nu_total]
    sstot_length= len(SSTOT)
    SSTOT= np.sum(SSTOT)
    SSTOT = SSTOT/float(L)
    
    SSWP_each = 0
    SSWP_divisor = 0
    SSWP = 0
    
    for population in population_list:
        data1= allele_profiles[coords[population],:]
        nu_gp= np.triu_indices(data1.shape[0],1)
        
        gp_dist= pairwise_distances(data1, metric=metric)
        gp_dist= gp_dist[nu_gp]
        SSWP_each= np.sum(gp_dist)
        
        SSWP_divisor= len(gp_dist)
        
        SSWP_divisor = (2*SSWP_divisor+0.25)**0.5 - 0.5
        
        if SSWP_each != 0:
            SSWP += float(SSWP_each)/float(SSWP_divisor)
    
    SSAP = SSTOT - SSWP
    
    squared_count_sum = float(sum([list(populations).count(x)**2 for x in list(set(populations))]))
    
    total_samples = float(L)
    total_pops = float(total_populations)
    dfAP = total_populations - 1
    dfWP = L - total_populations
    MSAP = float(SSAP/dfAP)
    MSWP = float(SSWP/dfWP)
    N0 = float((total_samples - float(squared_count_sum/total_samples)) * float(1/(total_pops-1)))
    VAP = float((MSAP - MSWP)/N0)
    if VAP + MSWP == 0:
        PhiPT = 0
    else:
        PhiPT = float(VAP/(VAP + MSWP))
    return PhiPT



def AMOVA_FM42(allele_profiles,populations,n_boot,metric= 'hamming'):
    '''
    allele_profiles: list of haplotype vectors (numeric, string doesnt matter).
    populations: list of population assignment of accessions in allele_profiles.
                -> same length as allele_profiles, same order.
    will treat NA's as alleles.. either pre-remove them or modify CombDiffrecv3
    '''
    different_populations = list(set(populations))
    population_list = different_populations
    total_populations = len(different_populations)
    
    #allele_profiles = [''.join([str(x) for x in y]) for y in allele_profiles]
    
    PhiPT = amova_cofactor(len(allele_profiles),allele_profiles,populations,total_populations,metric)
    District = []
    
    for IT in range(n_boot):
        District.append(amova_cofactor(len(allele_profiles),allele_profiles,random.sample(populations,len(populations))),total_populations,metric)
    Sign = 0
    if District:
        District= [x[0] for x in District]
        Sign= scipy.stats.norm(np.mean(District),np.std(District)).cdf(PhiPT)
        
    return PhiPT,Sign
# 