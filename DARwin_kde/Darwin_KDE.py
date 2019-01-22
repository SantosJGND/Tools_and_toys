
from Kernel_tools import *
########## START HERE #############

import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--CHR",type=int,help = "Chromosome number")

parser.add_argument("--geno",help = "Chromosome file to work on")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--bim",help = "snp information bim format.")

parser.add_argument("--ref",help = "reference accessions indexes in genofile.")

parser.add_argument("--admx",help = "admixed accession indexes in geno file")

parser.add_argument("--start",type= int,default= 0,help = "where to start, in snps")

parser.add_argument("--end",type= int,default= 0,help = "where to end, in snps")

parser.add_argument("--bornes",type= int,default= 1000,help = "bornes")

parser.add_argument("--out",type= str,default= '',help = "output directory")

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter [savgol filter]")
###
parser.add_argument("--MSprint",action= "store_true",help = "if given prints cluster stats.")
###
parser.add_argument("--VARprint",action= "store_true",help = "if given prints PC explained variance per window. If PCA is not chosen just prints out 0's")
###
parser.add_argument("--id",type= str,default= 'gene_rand',help = "Give your analysis an ID. default is set to integer 2")
###
parser.add_argument("-c",action = "store_true",help = "specific accession choice file")
###
parser.add_argument("-w",type = int,default = 200, help = "Window size - markers")
### 
parser.add_argument("--threshold",type = float,default = 5,help = "Intermediate classification threshold")
###
parser.add_argument("--outlier",type = float,default = 0.001,help = "Outlier threshold")
### 
parser.add_argument("--clustmethod",default = "MeanShift",choices = ["MeanShift","DBscan","HDBscan"],help = "Clustering method to extract reference specific clusters. MS, dbscan and hdbsan available. MS preferred")
###
parser.add_argument("--het",type = float,default = 5e-2, help = "Heterozygosity filter")
### 
parser.add_argument("--dr",default = 'NMF',help = "Dimensionality reduction. options: PCA, NMF")
### 
parser.add_argument("--ncomp",type = int,default = 4,help = "Number of components kept in case of PCA reduction")
### 
parser.add_argument("--outmethod",default = "None",help = "Outlier filter of population refs method. options: DBSCAN, NMF, Perc, Z, Isolation_forest.")
###
parser.add_argument("--overlap",type= int,default = 100,help = "Overlap between windows, in snps")
###

args = parser.parse_args()



########## Complementary files.

def read_refs(index_file,Fam_lib):
    indxs = recursively_default_dict()
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs[int(line[0])][Fam_lib[line[1]]] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())]


####

Fam = FAMread(args.fam)

MissG, Gindex = BIMread(args.bim)


GenoSUF = args.geno

admx_lib, Crossed = read_refs(args.admx,Fam)
refs_lib, Parents = read_refs(args.ref,Fam)
N_pops= len(refs_lib)

admx_lib.update(refs_lib)

Geneo = admx_lib

CHR = args.CHR
BIN = args.bin

Home = args.out

if len(args.out) > 0:
    Home= args.out + '/'

print('Number of markers in this chr: {}'.format(len(MissG[CHR])))

print('Population labels: {}'.format(Geneo.keys()))

print('Population Sizes: {}'.format([len(x) for x in Geneo.values()]))
####
####


def Set_up(Chr0,Chr1,Sub,MissG):
    Settings= []
    for Chr in range(Chr0,Chr1):
        Size= len(MissG[Chr])
        Steps = [x for x in range(0,Size+1,int(Size/float(Sub)))]
        Steps[-1] = Size
        for i in range(len(Steps)-1):
            Settings.append([Chr,Steps[i],Steps[i+1]-1])
    return Settings


def local_sampling_correct(data):
    ncomp= data.shape[1]
    pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized')
    features= pca.fit_transform(data)
    
    N= 100
    bandwidth = estimate_bandwidth(features, quantile=0.2)
    if bandwidth <= 1e-3:
        bandwidth = 0.1
    
    params = {'bandwidth': np.linspace(np.min(features), np.max(features),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    ## perform MeanShift clustering.
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False, min_bin_freq=20)
    ms.fit(features)
    labels1 = ms.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
    
    ## Extract the KDE of each cluster identified by MS.
    Proxy_data= []
    
    for lab in label_select.keys():
        if len(label_select[lab]) < 3:
            continue
        Quanted_set= features[label_select[lab],:]
        
        grid.fit(Quanted_set)
        
        kde = grid.best_estimator_
        Extract= kde.sample(N)
        Return= pca.inverse_transform(Extract)
        Proxy_data.extend(Return)
    
    Proxy_data= np.array(Proxy_data)
    
    pca = PCA(n_components=ncomp, whiten=False,svd_solver='randomized').fit(Proxy_data)
    New_features= pca.transform(data)
    
    return New_features



   
#### Fam
#### Whose
#### Miss
#### Geneo
#### Parents
#### GenoFile
Miss= MissG[CHR]
GenoFile= GenoSUF 

######################
######################

### Here define some things
Window = args.w

#### Intermediate classification threshold
Diff_threshold = args.threshold
X_threshold= args.outlier
Filter_Het = args.het


## Dimensionality reduction: PCA, NMF
DIMr = args.dr
n_comp = args.ncomp
PC_var_lim = .01
##
## savgol filter parameter
BIN = args.bin

## Chose clustering method
Method = args.clustmethod

## Outlier filter method: DBSCAN, NMF, Perc, Z, Isolation_forest.
Outlier_method = args.outmethod


quantile_filter_1 = 20 ### quantile filter for Percentile based outlier.

## KDE estimation tool
## Attention: sklearn uses GridsearchCV. While this implies more parameters, it 
## nonetheless ensures the construction of a KDE when scipy KDE would break down.
## This is because scipy KDE will always estimate bandwidth itself, thus using the points given only,
## while sklearn KDE allows the bandwidth to be passed. This allows us in turn to set a minimum,
## and ensure it to never come to 0.
## sklearn, scipy

KDE_tool = 'sklearn'
Bandwidth_split = 30

### Output normalization: CDF, Max
normalize = 'CDF'

### Resampling size for CDF estimates.
KDE_samples = 1000

# Controlling for global likelihood.
Control = False

start= args.start - args.bornes

end= args.end + args.bornes

##
##
##
##
### These things define themselves


Crossed = [x for x in Geneo.keys() if x not in Parents]

Whose = list(it.chain(*Geneo.values()))
SequenceStore = {fy:[] for fy in Whose}

Likes = {x:[] for x in range(len(Parents))}

Accuracy = []

Geno = open(GenoFile,"r")
Points = []
Points_out = []

PC_var= []

Win = 0
Index = 0
Intervals = []

Construct = recursively_default_dict()

for line in Geno:
    Codes = [0,nan,2,0,0,0,0,0,0,nan]
    d = Miss[Index][0]
    if Index > end:
        break
    if len([x for x in Whose if line[x] =='1']) / float(len(Whose)) > Filter_Het:
        Index += 1
        continue
    if d >= start and d <= end:
        for judas in SequenceStore.keys():
            SequenceStore[judas].append(Codes[int(line[judas])])
        Win += 1
    Index += 1


Geno.close()


#### N polymorphic loci

Sentence= np.array([SequenceStore[x] for x in Whose])

polyM= np.nansum(Sentence,axis= 0)
var_index= [x for x in range(len(polyM)) if polyM[x] > 0]

print('{} polymorphic sites.'.format(str(len(var_index))))

###
NA_Ind= [sum(np.isnan(x)) for x in Sentence]
ind_exclude= [x for x in range(len(NA_Ind)) if NA_Ind[x] > .5 * len(var_index)]

NA_mark= [sum(np.isnan(x)) for x in Sentence.T]
snp_exclude= [x for x in range(len(NA_mark)) if NA_mark[x] > .5 * len(NA_mark)]
#####


s1 = time.time()
window_start = Index - Window + 1
Seq = SequenceStore

Aro = np.nan_to_num(np.array([Seq[x] for x in it.chain(*[Geneo[x] for x in Crossed])]))

Daddy = np.nan_to_num(np.array([Seq[x] for x in it.chain(*[Geneo[x] for x in Parents])]))

b = np.concatenate((Aro,Daddy),axis = 0)

b= b[:,var_index]

data = np.zeros((b.shape[0],b.shape[1]+1))
data[:,:-1] = b

pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data)
data = pca.transform(data)
PC_var.append([x for x in pca.explained_variance_])
data= local_sampling_correct(data)

Accurate = []
params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)


######################################### 
############# TEST #####################
#########################################

for D in range(len(Parents)):
    Where = [sum([len(Geneo[x]) for x in Parents[0:D]])+Aro.shape[0],sum([len(Geneo[x]) for x in Parents[0:(D+1)]])+Aro.shape[0]]
    Where = [int(x) for x in Where]
    Quanted_set = data[Where[0]:Where[1],:]
    Below = []
    ###################            ###################
    ## KDE estimation on filtered parental data set ##
    ###################            ###################
    
    Indexes = range(Quanted_set.shape[0])
    
    Indexes = [x for x in Indexes if x not in Below]
    
    Set_ref = np.vstack({tuple(row) for row in Quanted_set[Indexes,:]})
    
    grid.fit(Quanted_set[Indexes,:])
    kde = grid.best_estimator_
    #Quanted_set = kde.sample(KDE_samples)
    P_dist = kde.score_samples(Quanted_set)
    Fist = kde.score_samples(data)
    P_dist = np.nan_to_num(P_dist)
    Fist = np.nan_to_num(Fist)
    
    if sum(np.isnan(P_dist)) == len(P_dist):
        if len(Likes[D]) == 0:
            Likes[D].append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
            Accurate.append([int(x in range(Where[0],Where[1])) for x in range(data.shape[0])])
        else:
            Likes[D].append(Likes[D][-1])
            Accurate.append(Likes[D][-1])
        continue
    
    ### Neutrality tests of filtered reference pop KDE derived log-Likelihoods.
    #Normality.append(scipy.stats.mstats.normaltest(P_dist)[1])
    ######
    if normalize == 'CDF':
        if np.std(P_dist)== 0:
            Fist= np.array([int(Fist[x] in P_dist) for x in range(len(Fist))])
        else:
            Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
    else:
        Fist = kde(data.T)
        Fist = Fist / max(Fist)
        Fist = [round(x,5) for x in Fist]
    if Control == True:
        Fist = Fist * Mortal
    
    Accurate.append(Fist)

Accurate= np.array(Accurate).T
Guys= Accurate
Test = [int(x <= X_threshold) for x in np.amax(np.array(Guys),axis = 1)]

maxim = np.argmax(Guys,axis = 1)
where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]

Consex = [x for x in it.combinations(range(N_pops),2)]

for h in range(len(maxim)):
    CL = []
    for j in Consex:
        Diff = Guys[h,j]
        if maxim[h] not in j or len([x for x in Diff if x < X_threshold]) > 0:
            continue
        if max(Diff) <= X_threshold:
            Diff = 0
        else:
#                            Diff = 1 - int(len([x for x in Diff if x >= X_threshold]) == 2)
            Diff = abs(max(Diff)) / abs(min(Diff))
            Diff = int(Diff > Diff_threshold)
        
        if Diff == 0:
            CL.append(j)
    
    if len(CL) == 2:
        maxim[h] = 7
    if len(CL) == 1:
        maxim[h] = sum(CL[0]) + N_pops

maxim[where_X] = N_pops


Original_labels= np.repeat(Crossed + Parents,[len(Geneo[x]) for x in Crossed + Parents])

Indexes= [y for y in it.chain(*[Geneo[x] for x in Crossed + Parents])]
Names= [Fam[y] for y in Indexes]

ID= args.id
L= len(var_index)

Home= ID + "_to_Darwin/"
filename= Home + ID + "_IDs.txt"

os.makedirs(os.path.dirname(filename), exist_ok=True)
Output= open(filename,"w")

for x in range(len(Whose)):
    Output.write("\t".join([str(x) for x in [Names[x],Original_labels[x],maxim[x]]]) + "\n")

Output.close()

filename= Home + ID + "_DataMatrix.txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output= open(filename,"w")

Output.write("Vars" + "\t")
for unit in range(L):
    Output.write(str(unit + 1) + "\t")

Output.write('\n')

for x in range(len(Whose)):
    Output.write(Names[x] + '\t')
    seq= SequenceStore[Indexes[x]]
    seq= [seq[x] for x in var_index]
    
    nantes= np.isnan(seq)
    nantes= [x for x in range(len(nantes)) if nantes[x] == True]
    seq= [str(x) for x in seq]
    
    for ney in nantes:
        seq[ney] = " "
    
    Output.write('\t'.join(seq) + "\n")

Output.close()


filename= Home + ID + "_Info.txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output= open(filename,"w")

Output.write('{} polymorphic snps.'.format(str(len(var_index))) + '\n')

Output.write('\n')

Output.write('Length: {}.'.format(str(end - start)) + '\n')


Output.close()

print('Done.')

