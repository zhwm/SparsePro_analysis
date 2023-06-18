import pandas as pd
import numpy as np
from scipy.special import softmax
from scipy.stats import chi2_contingency
import argparse
import os


def chi_squared(row):
    if any(row[0:4] == 0):
        return '1.00e+00'
    chi2, p, dof, expctd = chi2_contingency(np.array([[row[0], row[1]], [row[2], row[3]]]), lambda_="log-likelihood")
    return '{:.2e}'.format(p)


def get_uni_enrich(annomat, pip):
    df = pd.DataFrame(0.0, columns=['k0', 'k1', 'r0', 'r1'], index=annomat.columns)
    for a in annomat.columns:
        df.loc[a, 'k0'] = (1 - annomat[a]).sum()
        df.loc[a, 'k1'] = annomat[a].sum()
        df.loc[a, 'r0'] = np.dot(pip, (1 - annomat[a]).values)
        df.loc[a, 'r1'] = np.dot(pip, annomat[a].values)
    df['W'] = np.log(df['r1'] * df['k0'] / df['r0'] / df['k1'])
    df['W_se'] = np.sqrt(1 / df['r1'] + 1 / df['r0'] - 1 / df['k1'] - 1 / df['k0'])
    df = df.round(2)
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
parser.add_argument('--annosuf', type=str, default=None, help='annotation file suffix', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--pthres', type=float, default=1e-5, help='p value cutoff for annotation enrichment')
args = parser.parse_args()

pipdf = pd.concat([pd.read_csv('{}{}.pip'.format(args.pip, c), sep='\t', header=None, index_col=0)
                   for c in range(1, 23) if os.path.getsize('{}{}.pip'.format(args.pip, c)) > 0])
pip = pipdf[2].values
anno = pd.concat([pd.read_csv('{}/{}.{}'.format(args.anno, c, args.annosuf), sep='\t', index_col=0)
                  for c in range(1, 23)])
annomat = anno.loc[pipdf.index]
df = get_uni_enrich(annomat, pip)
df.reset_index().to_csv("{}.wsep".format(args.prefix), sep='\t', index=False)
print(df)

sigidx = [i for i in range(annomat.shape[1]) if float(df['p'][i]) < args.pthres]
if len(sigidx) == 0:
    print("None of the {} annotations is significantly enriched at p-value threshold {}.".
          format(annomat.shape[1], args.pthres))
else:
    print("{} annotations are deemed significantly enriched at {} p-value threshold and used to update priors."
          "Saving result to {}.W{} file...".format(len(sigidx), args.pthres, args.prefix, args.pthres))
    sigANNOT = annomat[annomat.columns[sigidx]]
    df = get_all_enrich(sigANNOT, pip)
    df['sigidx'] = sigidx
    df.reset_index().to_csv(os.path.join(args.save, '{}.W{}'.format(args.prefix, args.pthres)), sep='\t', index=False)
    print(df)