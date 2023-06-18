import pandas as pd
import argparse
import os
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri

ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
susieR = importr('susieR')

np.set_printoptions(precision=4, linewidth=200)


def get_XX_XtX_ytX_z(LD, Z, N):
    XX = np.ones(len(Z)) * N
    XtX = LD * N
    ytX = Z * np.sqrt(N)
    return XX, XtX, ytX


def zld(args):
    ldlists = pd.read_csv(args.zld, sep='\s+')  # z ld anno
    print("LD list with {} LD blocks loaded\n".format(len(ldlists)))
    ldlists['h2'] = 0.0
    ldlists['pval'] = 1.0
    ldlists['varb'] = 0.0
    for ite in range(len(ldlists)):
        cs = []
        cs_pip = []
        zfile, ldfile = ldlists.iloc[ite, 0:2]
        LD = pd.read_csv(os.path.join(args.zdir, ldfile), sep='\s+', header=None)  # no header square symmetric matrix
        z = pd.read_csv(os.path.join(args.zdir, zfile), sep='\t', header=None, index_col=0)  # no header matched vector
        assert len(LD) == len(z), 'LD file and Z file must match!'
        susierss = susieR.susie_rss(z[1].values.reshape(-1, 1), LD.values, L=args.K, n=args.N)
        pip_vec = [i for i in susieR.susie_get_pip(susierss)]
        ssets = susieR.susie_get_cs(susierss, Xcorr=LD.values, coverage=0.95, min_abs_corr=0.5)
        z['pip'] = pip_vec
        susie_sets = ssets[0]
        if susie_sets == ro.rinterface.NULL:
            print("No effect detected")
        else:
            scoverage = [i for i in ssets[3]]
            mcs = {}
            for set_i, susie_set in enumerate(susie_sets):
                print(set_i)
                print(susie_set)
                mcs[set_i] = np.array(susie_set) - 1
            print("Detected k = {}\n".format(len(mcs)))
            for e in mcs:
                mcs_idx = [z.index[j] for j in mcs[e]]
                print('The {}-th effect contains effective variants:'.format(e))
                print('causal variants: {}'.format(mcs_idx))
                print('posterior coverage: {}'.format(scoverage[e]))
                print()
                cs.append(mcs_idx)
                cs_pip.append(scoverage[e])
        z.to_csv(os.path.join(args.save, "{}.pip".format(zfile)), sep='\t', header=False, index=True)
        allcs = pd.DataFrame({"cs": ['/'.join(i) for i in cs], "pip": cs_pip})
        allcs.to_csv(os.path.join(args.save, "{}.cs".format(zfile)), sep='\t', header=True, index=False)
    ldlists.to_csv(os.path.join(args.save, "{}.h2".format(args.prefix)), sep='\t', header=True, index=False)


parser = argparse.ArgumentParser(description='SparsePro Commands:')
parser.add_argument('--zld', type=str, default=None,
                    help='locus finemapping mode: path to matched zscore and ld lists (or annotation file)')
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--K', type=int, default=None, help='largest number of effect', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')

args = parser.parse_args()

if not os.path.exists(args.save):
    os.makedirs(args.save)

zld(args)
