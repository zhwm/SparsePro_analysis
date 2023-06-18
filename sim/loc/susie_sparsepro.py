import pandas as pd
import argparse
import os
import numpy as np
from scipy.stats import chi2, entropy
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri

ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
susieR = importr('susieR')


def get_XX_XtX_ytX_z(LD, Z, N):
    XX = np.ones(len(Z)) * N
    XtX = LD * N
    ytX = Z * np.sqrt(N)
    return XX, XtX, ytX


def get_HESS_h2_z(LD, Z, N, ptLD=0.2, ptp=1e-5):
    """calculate local heritabilities"""
    zsquare = Z ** 2
    pvalmin = chi2.sf(max(zsquare), 1)
    idx_retain = []
    idx_exclude = [i for i in range(len(Z))]
    while len(idx_exclude) > 0:  # P + T
        maxid = idx_exclude[np.argmax(zsquare[idx_exclude])]  # find the idx with smallest p-value
        pval = chi2.sf(zsquare[maxid], 1)
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(np.abs(LD[maxid, :]) > ptLD)[0]]
        if pval > ptp:
            break
    Indidx = np.sort(idx_retain)  # obtain independent signals
    P = len(Indidx)
    LD_id = LD[np.ix_(Indidx, Indidx)]
    R_inv = np.linalg.inv(LD_id)
    vec_id = Z[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(), R_inv), vec_id) - P) / (N - P)
    var_b = np.max(zsquare[Indidx]) / N
    if h2_hess < 0.0001:
        h2_hess = 0.0001
    if h2_hess > 0.2:
        h2_hess = 0.2
        print("Heritability estimates is exceeding 0.2, which is rare for a locus, please check!!")
    return h2_hess, var_b, pvalmin, P


def get_sp_pip(gamma):
    return np.max((gamma), axis=1)


def get_sp_effect(gamma, cthres=0.95, ethres=20.0):
    vidx = np.argsort(-gamma, axis=1)
    matidx = np.argsort(-gamma, axis=0)
    pt, kt = gamma.shape
    mat_eff = np.zeros((pt, kt))  # effective gamma
    for p in range(pt):
        mat_eff[p, vidx[p, 0]] = gamma[p, vidx[p, 0]]
    csum = mat_eff.sum(axis=0).round(2)
    print("Attainable coverage for effect groups: {}".format(csum))  # statistical evidence
    eff = {}
    eff_gamma = {}
    for k in range(kt):
        if csum[k] > cthres:
            if entropy(mat_eff[:, k]) < np.log(ethres):
                for p in range(pt):
                    if np.sum(mat_eff[matidx[0:p, k], k]) > cthres * csum[k] or mat_eff[matidx[p, k], k] < 0.01:
                        eff[k] = matidx[0:p, k]
                        eff_gamma[k] = mat_eff[eff[k], k].round(4)
                        break
    return eff, eff_gamma


def zld(args):
    ldlists = pd.read_csv(args.zld, sep='\s+')  # z ld anno
    print("LD list with {} LD blocks loaded\n".format(len(ldlists)))
    for ite in range(len(ldlists)):
        cs = []
        cs_pip = []
        spcs = []
        spcs_pip = []
        zfile, ldfile = ldlists.iloc[ite, 0:2]
        LD = pd.read_csv(os.path.join(args.zdir, ldfile), sep='\s+', header=None)  # no header square symmetric matrix
        z = pd.read_csv(os.path.join(args.zdir, zfile), sep='\t', header=None, index_col=0)  # no header matched vector
        assert len(LD) == len(z), 'LD file and Z file must match!'
        h2_hess, var_b, p_chi, K = get_HESS_h2_z(LD.values, z[1].values, args.N)
        susiehess = susieR.susie_rss(z[1].values.reshape(-1, 1),
                                     LD.values,
                                     L=args.K,
                                     n=args.N,
                                     scaled_prior_variance=var_b,
                                     residual_variance=1-h2_hess,
                                     estimate_prior_variance=False,
                                     estimate_residual_variance=False)
        alpha = ro.FloatVector(susiehess.rx2('alpha'))
        alpha_np = np.array(alpha).transpose()
        sp_pip = get_sp_pip(alpha_np)
        pip_vec = [i for i in susieR.susie_get_pip(susiehess)]
        ssets = susieR.susie_get_cs(susiehess, Xcorr=LD.values, coverage=0.95, min_abs_corr=0.5)
        mcs, eff_gamma = get_sp_effect(alpha_np)
        z['pip'] = pip_vec
        z['sp_pip'] = sp_pip
        if len(mcs) == 0:
            print("No effect detected")
        else:
            print("Detected k = {}\n".format(len(mcs)))
            for e in mcs:
                mcs_idx = [z.index[j] for j in mcs[e]]
                print('The {}-th effect contains effective variants:'.format(e))
                print('causal variants: {}'.format(mcs_idx))
                print('posterior inclusion probabilities: {}'.format(eff_gamma[e]))
                print()
                spcs.append(mcs_idx)
                spcs_pip.append(eff_gamma[e])
        susie_sets = ssets[0]
        if susie_sets == ro.rinterface.NULL:
            print("No effect detected")
        else:
            scoverage = [i for i in ssets[3]]
            mcs = {}
            for set_i, susie_set in enumerate(susie_sets):
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
        print("Updated cs")
        z.to_csv(os.path.join(args.save, "{}.pip".format(zfile)), sep='\t', header=False, index=True)
        allcs = pd.DataFrame({"cs": ['/'.join(i) for i in cs], "pip": cs_pip})
        allcs.to_csv(os.path.join(args.save, "{}.cs".format(zfile)), sep='\t', header=True, index=False)
        allspcs = pd.DataFrame({"cs": ['/'.join(i) for i in spcs],
                                "pip": ['/'.join([str(j) for j in i]) for i in spcs_pip]})
        allspcs.to_csv(os.path.join(args.save, "{}.sp.cs".format(zfile)), sep='\t', header=True, index=False)


parser = argparse.ArgumentParser(description='SparsePro Commands:')
parser.add_argument('--zld', type=str, default=None,
                    help='locus finemapping mode: path to matched zscore and ld lists (or annotation file)')
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--K', type=int, default=None, help='largest number of effect', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--gwp', type=int, default=1e-5, help='p-value threshold for heritability estimate')


args = parser.parse_args()

if not os.path.exists(args.save):
    os.makedirs(args.save)

zld(args)
