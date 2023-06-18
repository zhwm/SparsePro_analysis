import pandas as pd
import numpy as np
from scipy.special import softmax
from scipy.stats import chi2_contingency
import argparse


def chi_squared(row):
    if any(row[0:4] == 0):
        return '1.00e+00'
    chi2, p, dof, expctd = chi2_contingency(np.array([[row[0], row[1]], [row[2], row[3]]]), lambda_="log-likelihood")
    return '{:.2e}'.format(p)


def get_uni_enrich(annomat, pip):
    """univariate test"""
    W_new = pd.Series(0.0, index=annomat.columns)
    df = pd.DataFrame(0.0, columns=['k0', 'k1', 'r0', 'r1'], index=annomat.columns)
    for a in annomat.columns:
        df.loc[a, 'k0'] = (1 - annomat[a]).sum()
        df.loc[a, 'k1'] = annomat[a].sum()
        df.loc[a, 'r0'] = np.dot(pip, (1 - annomat[a]).values)
        df.loc[a, 'r1'] = np.dot(pip, annomat[a].values)
        W_new[a] = np.log(df.loc[a, 'r1'] * df.loc[a, 'k0'] / df.loc[a, 'r0'] / df.loc[a, 'k1'])
    df = df.round(2)
    df['W'] = W_new.round(2)
    df['W_se'] = np.sqrt(1 / df['r1'] + 1 / df['r0'] - 1 / df['k1'] - 1 / df['k0']).round(2)
    df['p'] = df.apply(chi_squared, axis=1)
    return df


def get_all_enrich(annomat, pip, maxite=20):
    """multivariate test"""
    W_new = pd.Series(0.0, index=annomat.columns)
    for ite in range(maxite):
        W_old = W_new.copy()
        df = pd.DataFrame(0.0, columns=['k0', 'k1', 'r0', 'r1'], index=annomat.columns)
        for a in annomat.columns:
            wt = softmax(np.dot(annomat.drop(a, axis=1).values, W_new.drop(a))).squeeze()
            df.loc[a, 'k0'] = np.dot(wt, (1 - annomat[a]).values) * annomat.shape[0]
            df.loc[a, 'k1'] = np.dot(wt, annomat[a].values) * annomat.shape[0]
            df.loc[a, 'r0'] = np.dot(pip, (1 - annomat[a]).values)
            df.loc[a, 'r1'] = np.dot(pip, annomat[a].values)
            W_new[a] = np.log(df.loc[a, 'r1'] * df.loc[a, 'k0'] / df.loc[a, 'r0'] / df.loc[a, 'k1'])
        eps = ((W_new - W_old) ** 2).sum()
        print("iteration {} with diff {}".format(ite, eps))
        if eps < 1e-2:
            break
    df = df.round(2)
    df['W'] = W_new.round(2)
    df['W_se'] = np.sqrt(1 / df['r1'] + 1 / df['r0'] - 1 / df['k1'] - 1 / df['k0']).round(2)
    df['p'] = df.apply(chi_squared, axis=1)
    return df


parser = argparse.ArgumentParser(description='Test enrichment:')
parser.add_argument('--pip', type=str, default=None, help='path to pip vectors', required=True)
parser.add_argument('--anno', type=str, default=None, help='path to annotation file', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for saving results', required=True)
parser.add_argument('--pthres', type=float, default=1e-5, help='p value cutoff for annotation enrichment')
args = parser.parse_args()

i_dict = pd.read_csv('~/utils/ukb/idx/22.rsid', sep='\t', index_col=0, header=None)[1]
anno = pd.read_csv(args.anno, sep='\t', index_col=0)
pipdf = pd.read_csv('{}/C{}.pip'.format(args.pip, 1), sep='\t', header=None, index_col=0)
annomat = anno.loc[[i_dict.get(i) for i in pipdf.index], ]
pip = pipdf[2].values
for c in range(2, 23):
    pipdf = pd.read_csv('{}/C{}.pip'.format(args.pip, c), sep='\t', header=None, index_col=0)
    annomat = pd.concat([annomat, anno.loc[[i_dict.get(i) for i in pipdf.index], ]], ignore_index=True)
    pip = np.concatenate([pip, pipdf[2].values])

df_uni = get_uni_enrich(annomat, pip)
print(df_uni)
df_uni.reset_index().to_csv("{}/{}.wsep".format(args.pip, args.prefix), sep='\t', index=False)
df_all = get_all_enrich(annomat, pip)
df_all['sigidx'] = [i for i in range(annomat.shape[1])]
print(df_all)
df_all.reset_index().to_csv("{}/{}.W1.0".format(args.pip, args.prefix), sep='\t', index=False)

sigidx = [i for i in range(annomat.shape[1]) if float(df_uni['p'][i]) < args.pthres]
if len(sigidx) == 0:
    print("None of the {} annotations is significantly enriched at p-value threshold {}.".
          format(annomat.shape[1], args.pthres))
else:
    print("{} annotations are deemed significantly enriched at {} p-value threshold and used to update priors."
          "Saving result to {}/{}.W{} file...".format(len(sigidx), args.pthres, args.pip, args.prefix, args.pthres))
    sigANNOT = annomat[annomat.columns[sigidx]]
    df = get_all_enrich(sigANNOT, pip)
    df['sigidx'] = sigidx
    df.reset_index().to_csv("{}/{}.W{}".format(args.pip, args.prefix, args.pthres), sep='\t', index=False)
