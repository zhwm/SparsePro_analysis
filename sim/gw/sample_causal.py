import pandas as pd
import numpy as np
import argparse
from scipy.special import softmax
import os

def sample_causal(sdir,snp,K,seed,beta,anno=None,aw=None):
    '''Return causal settings with plink bim file, number of causal variants, RNG seed, effect size, annotation and annotation enrichment weight as input'''
    np.random.seed(seed)
    rs=pd.read_csv(snp,sep='\t',header=None)
    if anno is None and aw is None:
        print("No annotation information provided.")
        prob = np.ones(rs.shape[0]) * 1/rs.shape[0]
    elif anno is not None and aw is not None:
        anno=pd.read_csv(anno,sep='\t',index_col=0)
        W=np.array(aw)
#        print("Use annotation {} with enrichment weight {} for simulation.".format(anno.columns,aw))
        ano=anno.loc[rs[0]]
        prob=softmax(np.dot(ano.values,W))
    else:
        print("Please provide both annotation and enrichment weight.")
    idxall=np.random.choice(rs[0], K, replace=False, p=prob)
    causal=pd.DataFrame({'causal':idxall})
    causal['beta']=np.random.choice(beta,K)
    causal.to_csv(os.path.join(sdir,'C{}.causal'.format(seed)),index=False,header=False,sep='\t')
#    if anno is not None:
#        print(anno.loc[causal['causal']])
    return causal

parser = argparse.ArgumentParser(description='Get Causal ID')
parser.add_argument('--sdir', type=str, default=None, help='folder for save output', required=True)
parser.add_argument('--snp', type=str, default=None, help='pruned in file', required=True)
parser.add_argument('--K', type=int, default=None, help='number of causal variants', required=True)
parser.add_argument('--seed', type=int, default=None, help='seed', required=True)
parser.add_argument('--anno', type=str, default=None, help='annotation files')
parser.add_argument('--aw', type=float, default=None, nargs='+', help='enrichment weights')
parser.add_argument('--beta', type=float, default=None,nargs='+', help='effect size')
args = parser.parse_args()

causal=sample_causal(args.sdir,args.snp,args.K,args.seed,np.array(args.beta),args.anno,np.array(args.aw))
#print(causal)
