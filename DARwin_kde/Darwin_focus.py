
from Kernel_tools import *
########## START HERE #############

import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--CHR",type=int,help = "Chromosome number")

parser.add_argument("--geno",help = "Chromosome file to work on")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--bim",help = "snp information bim format.")

parser.add_argument("--focus",help = "reference accessions indexes in genofile.")

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
parser.add_argument("--ncomp",type = int,default = 4,help = "Number of components kept in case of PCA reduction")
### 
parser.add_argument("--outmethod",default = "None",help = "Outlier filter of population refs method. options: DBSCAN, NMF, Perc, Z, Isolation_forest.")
###
parser.add_argument("--overlap",type= int,default = 100,help = "Overlap between windows, in snps")
###

args = parser.parse_args()



########## Complementary files.

def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs


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

###
Focus= read_focus(args.focus)
Whose= [Fam[x] for x in Focus]

###

CHR = args.CHR
BIN = args.bin

Home = args.out

if len(args.out) > 0:
    Home= args.out + '/'

print('Number of markers in this chr: {}'.format(len(MissG[CHR])))

####
####

   
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
n_comp = args.ncomp
##
## savgol filter parameter
BIN = args.bin

## Chose clustering method
Method = args.clustmethod


## KDE estimation tool
## Attention: sklearn uses GridsearchCV. While this implies more parameters, it 
## nonetheless ensures the construction of a KDE when scipy KDE would break down.
## This is because scipy KDE will always estimate bandwidth itself, thus using the points given only,
## while sklearn KDE allows the bandwidth to be passed. This allows us in turn to set a minimum,
## and ensure it to never come to 0.
## sklearn, scipy

Bandwidth_split = 30

start= args.start - args.bornes

end= args.end + args.bornes

##
##
##
##
### These things define themselves
###

SequenceStore = {fy:[] for fy in Whose}

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
    if d > end:
        break
    if len([x for x in Whose if line[x] =='1']) / float(len(Whose)) > Filter_Het:
        Index += 1
        continue
    if d >= start and d <= end:
        for judas in Whose:
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

Sequences= [Seq[x] for x in Whose]
Sequences= np.array(Sequences)
Sequences= np.nan_to_num(Sequences)

b = Sequences

b= b[:,var_index]

data = np.zeros((b.shape[0],b.shape[1]+1))
data[:,:-1] = b


Names= [Fam[y] for y in Whose]

ID= args.id
L= len(var_index)

Home= ID + "_to_Darwin/"
filename= Home + ID + "_IDs.txt"

os.makedirs(os.path.dirname(filename), exist_ok=True)
Output= open(filename,"w")

for x in range(len(Whose)):
    Output.write("\n".join([str(x) for x in Names]))

Output.close()

filename= Home + ID + "_DataMatrix.txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output= open(filename,"w")

Output.write("Vars" + "\t")
Output.write('\t'.join([str(unit + 1) for unit in range(L)]))
Output.write('\n')

for x in range(len(Whose)):
    Output.write(Names[x] + '\t')
    seq= Sentence[x]
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

