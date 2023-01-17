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

np.set_printoptions(precision=4, linewidth=200)

def title():
    print('**********************************************************************')
    print('* SparsePro for efficient genome-wide fine-mapping                   *')
    print('* Version 2.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')
    print()

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

class SparsePro(object):

    def __init__(self,P,K,XX,h2,var_b):
        '''initialize and set hyperparameters'''
        self.p = P
        self.k = K
        self.gamma = np.zeros((self.p,self.k))
        self.beta_mu = np.zeros((self.p,self.k))
        self.beta_prior_tau = np.tile((1.0 / var_b * np.array([k+1 for k in range(self.k)])),(self.p,1))
        self.y_tau = 1.0 / (1.0-h2) # var_Y get cancelled eventually
        self.prior_pi = np.ones((self.p,)) * (1/self.p)
        self.beta_post_tau = np.tile(XX.reshape(-1,1),(1,self.k)) * self.y_tau + self.beta_prior_tau

    def infer_q_beta(self,XX,ytX,XtX):
        '''perform variational updates'''
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:,idxall] * self.beta_mu[:,idxall]).sum(axis=1)
            self.beta_mu[:,k] = (ytX-np.dot(beta_all_k, XtX))/self.beta_post_tau[:,k] * self.y_tau # update mu
            u = -0.5*np.log(self.beta_post_tau[:,k]) + np.log(self.prior_pi.transpose()) + 0.5 * self.beta_mu[:,k]**2 * self.beta_post_tau[:,k]
            self.gamma[:,k] = softmax(u) # update gamma

    def get_elbo(self,XX,ytX,XtX):
        beta_all = (self.gamma * self.beta_mu).sum(axis=1)
        ll1 = - 2 * np.dot(beta_all,ytX)
        ll2 = ((self.gamma * self.beta_mu**2).sum(axis=1) * XX).sum()
        W = self.gamma * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(),XtX),W)
        ll3 = WtRW.sum() - np.diag(WtRW).sum()
        ll = -0.5 * self.y_tau * (ll1 + ll2 + ll3)
        betaterm1 = -0.5 * (self.beta_prior_tau * self.gamma * (self.beta_mu**2)).sum()
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1,1),(1,self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma!=0] * np.log(self.gamma[self.gamma!=0])).sum()
        eut = -0.5 * (self.gamma * np.log(self.beta_post_tau)).sum()  # extra unstandardized term
        mkl = betaterm1 + gammaterm1 - gammaterm2 + eut
        elbo = ll + mkl
        return ll,mkl,elbo

    def get_PIP(self):
        return np.max((self.gamma),axis=1)

    def update_pi(self, new_pi):
        self.prior_pi = new_pi

    def get_effect(self,sthres=0.95):
        matidx = np.argsort(-self.gamma, axis=0)
        gamma = self.gamma
        gamma[gamma<0.01]=0.0
        gammasum = gamma.sum(axis=0)
        eff={}
        eff_gamma={}
        eff_mu={}
        for k in range(self.k):
            if gammasum[k]>0.5:
                if entropy(gamma[:,k])<np.log(20):
                    for p in range(self.p):
                        if np.sum(gamma[matidx[0:p,k],k])/gammasum[k] > sthres:
                            eff[k]=matidx[0:p,k].tolist()
                            eff_gamma[k]=self.gamma[eff[k],k].round(4)
                            eff_mu[k]=self.beta_mu[eff[k],k].round(4)
                            break
        return eff,eff_gamma,eff_mu

    def train(self,XX,ytX,XtX,maxite=100,eps=0.01,verbose=False,loss=0.0):
        for ite in range(maxite):
            self.infer_q_beta(XX, ytX, XtX)
            ll, mkl, elbo = self.get_elbo(XX, ytX, XtX)
            if verbose:
                print('*'*70)
                print('Iteration-->{} . Likelihood: {:.2f} . KL: {:.2f} . ELBO: {:.2f}'.format(ite, ll, mkl, elbo))
            if abs(elbo-loss)<eps:
                break
            if ite==(maxite-1):
                print("Algorithm not converged. Please make sure matched summary statistics and LD were provided!")
            loss = elbo

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

def run_sparsepro(LD,sthres,Z,N,K,h2_hess,var_b,verbose,new_pi=None):
    XX,XtX,ytX = get_XX_XtX_ytX_z(LD,Z,N)
    model = SparsePro(len(Z),K,XX,h2_hess,var_b)
    if new_pi is not None:
        model.update_pi(new_pi)
    model.train(XX, ytX, XtX, verbose=verbose)
    mcs, eff_gamma, eff_mu = model.get_effect(sthres)
    pip_vec = model.get_PIP().round(4)
    return mcs, eff_gamma, eff_mu, pip_vec

def ukb(args):
    print("Using genome-wide fine-mapping mode with --ukb")
    if args.anno is None:
        print("Statistical fine-mapping...")
    elif args.aW is not None:
        print("Annotated fine-mapping...")
        anno = pd.read_csv(args.anno,sep='\s+',index_col=0)
        W_sig = pd.read_csv(args.aW, sep='\t')
        sigidx = W_sig['sigidx'].tolist()
        W_new = W_sig['W_sig'].values
    else:
        print("Statistical fine-mapping and estimating annotation enrichment")
        anno = pd.read_csv(args.anno,sep='\s+',index_col=0)
    z = pd.read_csv(args.zdir,sep="\s+",header=None,index_col=0) #no header, first column idx(chr.pos.ref.alt)
    print("summary statistics loaded at {}".format(time.strftime("%Y-%m-%d %H:%M")))
    ldlists=pd.read_csv(args.ukb,sep='\s+') #header with 3 columns
    pip = []
    pip_name = []
    cs = []
    cs_pip = []
    cs_eff = []
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
        if args.h2:
            h2_hess, var_b, p_chi, K=ldlists.loc[ite,['h2','varb','pval','K']]
        else:
            h2_hess, var_b, p_chi, K=get_HESS_h2_z(LD.values, z.loc[idx,1].values, args.N, gwp=args.gwp)
        if args.K is not None:
            K = args.K
        if p_chi>args.gwp:
            print('Smallest p-value in this region is larger than defined threshold, continue to next block\n')
            continue
        if args.anno is not None and args.aW is not None:
            ianno = anno.loc[idx]
            ANN = ianno.values[:,sigidx]
            new_pi= softmax(np.dot(ANN,W_new))
            mcs, eff_gamma, eff_mu, pip_vec = run_sparsepro(LD.values, args.sthres, z.loc[idx,1].values, args.N, K, h2_hess, var_b, args.verbose, new_pi)
        else:
            mcs, eff_gamma, eff_mu, pip_vec = run_sparsepro(LD.values, args.sthres, z.loc[idx,1].values, args.N, K, h2_hess, var_b, args.verbose)
        print("{} variants loaded from {} explaining {:2.2%} of trait heritability with smallest p-value {:.2e} and deemed at most {} effects  \n".format(LD.shape[1], ldfile, h2_hess, p_chi, K))
        pip.extend([pip_vec[i] for i in effidx])
        z_vec.extend([z.loc[idx,1].values[i] for i in effidx])
        pip_name.extend([idx[i] for i in effidx])
        if len(mcs)==0:
            print("No effect detected")
            print()
            continue
        print("Detected k = {}\n".format(len(mcs)))
        for e in mcs:
            if mcs[e][0] in effidx:
                mcs_idx = [idx[j] for j in mcs[e]]
                print('The {}-th effect contains effective variants:'.format(e))
                print('causal variants: {}'.format(mcs_idx))
                print('posterior inclusion probabilities: {}'.format(eff_gamma[e]))
                print('posterior causal effect size: {}'.format(eff_mu[e]))
                print()
                cs.append(mcs_idx)
                cs_pip.append(eff_gamma[e])
                cs_eff.append(eff_mu[e])
        if not args.h2:
            ldlists.at[ite,'h2'] = '{:.2e}'.format(h2_hess)
            ldlists.at[ite,'pval'] = '{:.2e}'.format(p_chi)
            ldlists.at[ite,'varb'] = '{:.2e}'.format(var_b)
            ldlists.at[ite,'K'] = '{:.2e}'.format(K)
    allPIP = pd.DataFrame({"idx":pip_name,'z':z_vec, "pip":pip})
    allPIP.to_csv(os.path.join(args.save,"{}.pip".format(args.prefix)),sep='\t',header=False,index=False)
    allcs = pd.DataFrame({"cs":['/'.join(i) for i in cs],"pip":['/'.join([str(j) for j in i]) for i in cs_pip],"beta":['/'.join([str(j) for j in i]) for i in cs_eff]})
    allcs.to_csv(os.path.join(args.save,"{}.cs".format(args.prefix)),sep='\t',header=True,index=False)
    ldlists.to_csv(os.path.join(args.save,"{}.h2".format(args.prefix)),sep='\t',header=True,index=False)
    if args.anno is not None and args.aW is None:
        annomat = anno.loc[pip_name]
        Wp=get_uni_enrich(annomat,pip)
        sigidx = [i for i in range(annomat.shape[1]) if Wp[i]<args.pthres]
        if len(sigidx)==0:
            print("None of the {} annotations is significantly enriched at p-value threshold {}.".format(annomat.shape[1], args.pthres))
        else:
            print("{} annotations are deemed significantly enriched at {} p-value threshold and used to update priors. Saving result to {}.W{} file...".format(len(sigidx),args.pthres,args.prefix,args.pthres))
            sigANNOT = annomat.values[:,sigidx]
            W_sig,W_se_sig = get_sig_enrich(sigANNOT, pip)
            df_W_sig = pd.DataFrame({'ANNO':annomat.columns[sigidx],'W_sig':W_sig, 'W_se_sig':W_se_sig, 'sigidx':sigidx})
            df_W_sig.to_csv(os.path.join(args.save,'{}.W{}'.format(args.prefix,args.pthres)),sep='\t',index=False)

def get_uni_enrich(annomat,pip):
    Wsep={}
    P = len(pip)
    K = sum(pip)
    for k in annomat.columns:
        A=(annomat[k]).sum()
        M=sum([pip[i] for i in range(len(pip)) if annomat[k][i]==1])
        if M==0:
           Wsep[k]=['{:.2e}'.format(P),'{:.2e}'.format(A),'{:.2e}'.format(K),'{:.2e}'.format(M),0.0,0.0,1.0]
           continue
        obs = np.array([[K-M,P-A-K+M],[M,A-M]])
        g, p, dof, expctd = chi2_contingency(obs, lambda_="log-likelihood")
        W = np.log(M*(P-A)/A/(K-M))
        W_se = np.sqrt(1/M + 1/(K-M) - 1/A - 1/(P-A))
        Wsep[k] = ['{:.2e}'.format(P),'{:.2e}'.format(A),'{:.2e}'.format(K),'{:.2e}'.format(M),'{:.2e}'.format(W),'{:.2e}'.format(W_se),'{:.2e}'.format(p)]
    df_Wsep = pd.DataFrame(Wsep)
    df_Wsep.index = ['Total', 'Anno', 'Causal', 'AnnoandCausal', 'W','se','p']
    df_Wsep.index.name = 'Wsep'
    df_Wsep = df_Wsep.transpose()
    print("Univariate testing finished at {}. Saving result to {}.wsep file...".format(time.strftime("%Y-%m-%d %H:%M"),args.prefix))
    df_Wsep.to_csv(os.path.join(args.save,'{}.wsep'.format(args.prefix)),sep="\t")
    return [float(i) for i in df_Wsep['p'].values]

def get_sig_enrich(sigANNOT,pip):
    P,G=sigANNOT.shape
    W = np.zeros(G)
    W_se = np.zeros(G)
    eps = 1000
    tot = sum(pip)
    for ite in range(20):
        W_old = W.copy()
        for i in range(G):
            idxall = [x for x in range(G)]
            idxall.remove(i)
            k = softmax(np.dot(sigANNOT[:,idxall],W[idxall]))
            kr = k[sigANNOT[:,i]==1].sum()
            r = np.array(pip)[np.where(sigANNOT[:,i])[0]].sum()/tot
            W_new = np.log((1-kr) * r / (1-r) / (kr))
            W[i] = W_new
            W_se_new = np.sqrt(1/(r*tot)+1/((1-r)*tot)-1/(kr*P)-1/((1-kr)*P))
            W_se[i] = W_se_new
        eps = ((W - W_old)**2).sum()
        print("iteration {} with diff {}".format(ite,eps))
        if eps < 1e-2:
            break
    return W,W_se

parser = argparse.ArgumentParser(description='SparsePro Commands:')
parser.add_argument('--ukb', type=str, default=None, help='genome-wide finemapping mode: path to LD lists')
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--LDdir', type=str, default=None, help='path to ld files', required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--anno', type=str, default=None, help='name of the annotation column in zld file')
parser.add_argument('--K', type=int, default=None, help='largest number of effect')
parser.add_argument('--gwp', type=int, default=1e-5, help='p-value threshold for heritability estimate')
parser.add_argument('--pthres', type=float, default=1e-5, help='p value threshold for enrichment')
parser.add_argument('--sthres', type=float, default=0.95, help='credible set level')
parser.add_argument('--aW', type=str, default=None, help='significant enriched file')
parser.add_argument('--h2', action="store_true", help='use previous h2 file as zld file')

args = parser.parse_args()
title()
if not os.path.exists(args.save):
    os.makedirs(args.save)

ukb(args)
