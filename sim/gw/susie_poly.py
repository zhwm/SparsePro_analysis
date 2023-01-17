import pandas as pd
import argparse
import time
import os
import numpy as np
from scipy.special import softmax
from scipy.stats import chi2
from scipy.stats import chi2_contingency
from scipy.stats import entropy
import scipy.sparse as sparse

import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
susieR = importr('susieR')

np.set_printoptions(precision=4, linewidth=200)

#obtain from https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/readme_ld.txt
def load_ld_npz(ld_prefix):
    #load the SNPs metadata
    gz_file = '%s.gz'%(ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps['A1'] + '.' + df_ld_snps['A2']
    #load the LD matrix
    npz_file = '%s.npz'%(ld_prefix)
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s'%(npz_file))
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps

def get_XX_XtX_ytX_z(LD,Z,N):
    XX = np.ones(len(Z)) * N
    XtX = LD * N
    ytX = Z * np.sqrt(N)
    return XX,XtX,ytX

def get_HESS_h2_z(LD,Z,N,LDthres=0.2,gwp=1e-5):
    '''calculate local heritabilities'''
    zsquare= Z**2
    maxzid = np.argmax(zsquare)
    pvalmin = 1 - chi2.cdf(zsquare[maxzid],1)
    var_b = zsquare[maxzid]/N
    idx_retain = []
    idx_exclude = [i for i in range(len(Z))]
    while len(idx_exclude)>0:
        maxid = idx_exclude[np.argmax(zsquare[idx_exclude])] #find the idx with smallest p-value
        pval = 1 - chi2.cdf(zsquare[maxid],1)
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(np.abs(LD[maxid,:])>LDthres)[0]]
        if pval>gwp:
            break
    Indidx = np.sort(idx_retain) #obtain independent signals through P + T
    P = len(Indidx)
    LD_id = LD[np.ix_(Indidx,Indidx)]
    R_inv = np.linalg.inv(LD_id)
    vec_id = Z[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(),R_inv),vec_id)-P)/(N-P)
    if h2_hess<0.0001:
        h2_hess = 0.0001
    if h2_hess>0.9:
        h2_hess = 0.9
    return h2_hess, var_b, pvalmin, P

def ukb(args):
    z = pd.read_csv(args.zdir,sep="\s+",header=None,index_col=0) #no header, first column idx(chr.pos.ref.alt)
    snpvar = pd.read_csv(args.snpvar,sep='\t')
    snpvar.index=['22.'+'.'.join(i.split('.')[1:]) for i in snpvar['SNP']] #fake 22:snpvar.index = snpvar['SNP']
    minvar = min(snpvar['SNPVAR'])
    z['prob'] = minvar
    z.loc[snpvar.index,'prob'] = snpvar['SNPVAR']
    print("summary statistics loaded at {}".format(time.strftime("%Y-%m-%d %H:%M")))
    ldlists=pd.read_csv(args.ukb,sep='\s+') #header with 3 columns
    pip = []
    pip_name = []
    cs = []
    cs_pip = []
    z_vec = []
    for ite in range(len(ldlists)):
        ld,start,end = ldlists.iloc[ite,0:3]
        ldfile = ld.replace('.npz','')
        df_R, df_ld_snps = load_ld_npz(os.path.join(args.LDdir,ldfile))
        idx = df_R.index.intersection(z.index)
        LD = df_R.loc[idx,idx]
        if len(idx)<10:
            print("less than 10 variants matched, skipping")
            continue
        pos = [int(i.split('.')[1]) for i in idx]
        effidx = [i for i in range(len(idx)) if ((pos[i] >= start) & (pos[i] < end))]
        effnum = len(effidx)
        print('{} variants in the range of {} to {}'.format(effnum, start, end))
        if effnum <=10:
            print('less than 10 effective variants, skipping')
            continue
        h2_hess, var_b, p_chi, K=get_HESS_h2_z(LD.values, z.loc[idx,1].values, args.N, gwp=args.gwp)
        print("{} variants loaded from {} explaining {:2.2%} of trait heritability with smallest p-value {:.2e} and deemed at most {} effects  \n".format(LD.shape[1], ldfile, h2_hess, p_chi, K))
        if p_chi>args.gwp:
            print('Smallest p-value in this region is larger than defined threshold, continue to next block\n')
            continue
        if args.K is not None:
            K = args.K
        susierss = susieR.susie_rss(z.loc[idx,1].values.reshape(-1,1),LD.values,L=K,n=args.N,prior_weights=z.loc[idx,'prob'].values.reshape(-1,1))
        pip_vec = [i for i in susieR.susie_get_pip(susierss)]
        ssets = susieR.susie_get_cs(susierss,Xcorr=LD.values,coverage = 0.95, min_abs_corr = 0.5)
        susie_sets = ssets[0]
        if susie_sets == ro.rinterface.NULL:
            print("No effect detected")
        else:
            scoverage = [i for i in ssets[3]]
            mcs = {}
            for set_i, susie_set in enumerate(susie_sets):
                print(set_i)
                print(susie_set)
                mcs[set_i]=np.array(susie_set)-1
            print("Detected k = {}\n".format(len(mcs)))
            for e in mcs:
                if mcs[e][0] in effidx:
                    mcs_idx = [idx[j] for j in mcs[e]]
                    print('The {}-th effect contains effective variants:'.format(e))
                    print('causal variants: {}'.format(mcs_idx))
                    print('posterior coverage: {}'.format(scoverage[e]))
                    print()
                    cs.append(mcs_idx)
                    cs_pip.append(scoverage[e])

        pip.extend([pip_vec[i] for i in effidx])
        z_vec.extend([z.loc[idx,1].values[i] for i in effidx])
        pip_name.extend([idx[i] for i in effidx])
    allPIP = pd.DataFrame({"idx":pip_name,'z':z_vec, "pip":pip})
    allPIP.to_csv(os.path.join(args.save,"{}.pip".format(args.prefix)),sep='\t',header=False,index=False)
    allcs = pd.DataFrame({"cs":['/'.join(i) for i in cs],"pip":cs_pip})
    allcs.to_csv(os.path.join(args.save,"{}.cs".format(args.prefix)),sep='\t',header=True,index=False)
    ldlists.to_csv(os.path.join(args.save,"{}.h2".format(args.prefix)),sep='\t',header=True,index=False)

parser = argparse.ArgumentParser(description='SparsePro Commands:')
parser.add_argument('--ukb', type=str, default=None, help='genome-wide finemapping mode: path to LD lists')
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--LDdir', type=str, default=None, help='path to ld files', required=True)
parser.add_argument('--snpvar', type=str, default=None, help='path to polyfun snpvar',required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--K', type=int, default=None, help='largest number of effect')
parser.add_argument('--gwp', type=int, default=1e-5, help='p-value threshold for heritability estimate')

args = parser.parse_args()
if not os.path.exists(args.save):
    os.makedirs(args.save)

ukb(args)
