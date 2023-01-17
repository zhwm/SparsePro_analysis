import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import entropy
import os
from sklearn.metrics import precision_recall_curve
from sklearn.calibration import calibration_curve
from sklearn.metrics import auc
from sklearn.metrics import precision_score, recall_score
import seaborn as sns

import matplotlib.pylab as pylab
plt.style.use('default')
params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 10),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

color_dict={'FINEMAP':'lime',
            'SuSiE rss':'royalblue',
            'PAINTOR-':'orchid',
            'PAINTOR+':'darkmagenta',
            'SparsePro-':'orange',
            'SparsePro+':'red'}

def get_Tv(SNPlist,wdir,rep=50):
    Tv = pd.concat([SNPlist.isin(pd.read_csv(os.path.join(wdir,'C{}.causal'.format(i)),sep='\t',header=None)[0]) for i in range(1,rep+1)])
    return Tv.values

def load_finemap_PIP(wdir,rep=50):
    PIP=pd.concat([pd.read_csv(os.path.join(wdir,"C{}.snp".format(i)),sep=' ').sort_values(by=['index'])['prob'] for i in range(1,rep+1)])
    return PIP.values

def load_PAINTOR_no_PIP(wdir,rep=50):
    PIP=pd.concat([pd.read_csv(os.path.join(wdir,"PAINTOR_no/C{}.results".format(i)),sep=' ') for i in range(1,rep+1)])
    return PIP['Posterior_Prob'].values

def load_PAINTOR_anno_PIP(wdir,rep=50):
    PIP=pd.concat([pd.read_csv(os.path.join(wdir,"PAINTOR_anno/C{}.results".format(i)),sep=' ') for i in range(1,rep+1)])
    return PIP['Posterior_Prob'].values

def load_sparsepro_PIP(wdir,rep=50):
    SP_no_PIP=pd.concat([pd.read_csv(os.path.join(wdir,"sparsepro_no/C{}.txt.pip".format(i)),sep='\t',header=None) for i in range(1,rep+1)])
    return SP_no_PIP[2].values

def load_sparsepro_anno_PIP(wdir,rep=50):
    SP_anno_PIP=pd.concat([pd.read_csv(os.path.join(wdir,"sparsepro_anno/C{}.txt.pip".format(i)),sep='\t',header=None) for i in range(1,rep+1)])
    return SP_anno_PIP[2].values

def load_susie_rss_PIP(wdir,rep=50):
    PIP=pd.concat([pd.read_csv(os.path.join(wdir,"susie_rss/C{}.txt.pip".format(i)),sep='\t',header=None) for i in range(1,rep+1)])
    return PIP[2].values

def ant_plot(ax,Tv,PIP1,PIP2,name1,name2):
    '''plot PIP contrast'''
    ax.scatter(x=PIP1[[not i for i in Tv]],y=PIP2[[not i for i in Tv]],c='black',s = 100,label='Non-causal')
    ax.scatter(x=PIP1[Tv],y=PIP2[Tv],c='red',s = 100,label='Causal')
    ax.set_xlabel(name1, fontsize = 25)
    ax.set_ylabel(name2, fontsize = 25)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,-0.3))
    ax.grid()
    return ax

def ant_plot_no(wdir,K,SNPlist):
    '''PIP plot without annotations'''
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    susie_rss_PIP = load_susie_rss_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    fig, axs = plt.subplots(1,3,figsize=(30, 9))
    ant_plot(axs[0],Tv,SP_no_PIP,finemap_PIP,'SparsePro- PIP','FINEMAP PIP')
    ant_plot(axs[1],Tv,SP_no_PIP,susie_rss_PIP,'SparsePro- PIP','SuSiE PIP')
    ant_plot(axs[2],Tv,SP_no_PIP,PAINTOR_no_PIP,'SparsePro- PIP','PAINTOR- PIP')

def ant_plot_anno(wdir,K,SNPlist):
    Tv = get_Tv(SNPlist,wdir)
    SP_anno_PIP = load_sparsepro_anno_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    fig, axs = plt.subplots(1,3,figsize=(30, 9))
    ant_plot(axs[0],Tv,SP_anno_PIP,SP_no_PIP,'SparsePro+ PIP','SparsePro- PIP')
    ant_plot(axs[1],Tv,PAINTOR_anno_PIP,PAINTOR_no_PIP,'PAINTOR+ PIP','PAINTOR- PIP')
    ant_plot(axs[2],Tv,SP_anno_PIP,PAINTOR_anno_PIP,'SparsePro+ PIP','PAINTOR+ PIP')

def get_auprc(PIP,Tv):
    precision, recall, thresholds = precision_recall_curve(Tv, PIP)
    return auc(recall, precision)

def prc_plot(ax,PIP,Tv,name,color):
    precision, recall, thresholds = precision_recall_curve(Tv, PIP)
    ax.plot(recall,precision, label=name+": "+"{:.4f}".format(auc(recall, precision)),color=color)
    ax.set_ylabel('Precision')
    ax.set_xlabel('Recall')
    return ax

def prc_plot_no(wdir,K,SNPlist):
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    susie_rss_PIP = load_susie_rss_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)

    fig, ax = plt.subplots(figsize=(10,10))
    prc_plot(ax,finemap_PIP,Tv,'FINEMAP','lime')
    prc_plot(ax,susie_rss_PIP,Tv,'SuSiE','royalblue')
    prc_plot(ax,PAINTOR_no_PIP,Tv,'PAINTOR-','orchid')
    prc_plot(ax,PAINTOR_anno_PIP,Tv,'PAINTOR+','darkmagenta')
    prc_plot(ax,SP_no_PIP,Tv,'SparsePro-','orange')
    ax.legend()

def prc_plot_anno(wdir,K,SNPlist):
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    susie_rss_PIP = load_susie_rss_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    SP_anno_PIP = load_sparsepro_anno_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)

    fig, ax = plt.subplots(figsize=(10,10))
    prc_plot(ax,finemap_PIP,Tv,'FINEMAP','lime')
    prc_plot(ax,susie_rss_PIP,Tv,'SuSiE','royalblue')
    prc_plot(ax,PAINTOR_no_PIP,Tv,'PAINTOR-','orchid')
    prc_plot(ax,PAINTOR_anno_PIP,Tv,'PAINTOR+','darkmagenta')
    prc_plot(ax,SP_no_PIP,Tv,'SparsePro-','orange')
    prc_plot(ax,SP_anno_PIP,Tv,'SparsePro+','red')
    ax.legend()

def calib_plot(ax,PIP,Tv,name,color):
    cp, mpv = calibration_curve(Tv,PIP)
    ax.plot(mpv,cp,marker='.',label=name,color=color)
    ax.set_ylabel('Actual Precision')
    ax.set_xlabel('Expect Precision')
    return ax

def calibration_plot_no(wdir,K,SNPlist):
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    susie_rss_PIP = load_susie_rss_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)

    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot([0,1],[0,1],'--',color='grey')
    calib_plot(ax,finemap_PIP,Tv,'FINEMAP','lime')
    calib_plot(ax,susie_rss_PIP,Tv,'SuSiE','royalblue')
    calib_plot(ax,PAINTOR_no_PIP,Tv,'PAINTOR-','orchid')
    calib_plot(ax,PAINTOR_anno_PIP,Tv,'PAINTOR+','darkmagenta')
    calib_plot(ax,SP_no_PIP,Tv,'SparsePro-','orange')
    ax.legend()

def calibration_plot_anno(wdir,K,SNPlist):
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    susie_rss_PIP = load_susie_rss_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    SP_anno_PIP = load_sparsepro_anno_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot([0,1],[0,1],'--',color='grey')
    calib_plot(ax,finemap_PIP,Tv,'FINEMAP','lime')
    calib_plot(ax,susie_rss_PIP,Tv,'SuSiE','royalblue')
    calib_plot(ax,PAINTOR_no_PIP,Tv,'PAINTOR-','orchid')
    calib_plot(ax,PAINTOR_anno_PIP,Tv,'PAINTOR+','darkmagenta')
    calib_plot(ax,SP_no_PIP,Tv,'SparsePro-','orange')
    calib_plot(ax,SP_anno_PIP,Tv,'SparsePro+','red')
    ax.legend()

def get_Pre_Re_S_cs(wdir,folder,ite):
    df_cs = pd.read_csv(os.path.join(wdir,folder,"C{}.txt.cs".format(ite)),sep='\t',dtype={'pip':str})
    if len(df_cs)==0:
        return 0.0,0.0,0
    df_causal = pd.read_csv(os.path.join(wdir,"C{}.causal".format(ite)),sep='\t',header=None)
    cs = [j for i in df_cs['cs'].values for j in i.split('/')]
    pip = [float(j) for i in df_cs['pip'].values for j in i.split('/')]
    Pre = len(set(df_causal[0]) & set(cs)) / len(cs)
    Re = len(set(df_causal[0]) & set(cs)) / len(df_causal)
    size = len(cs)
    return Pre,Re,size

def get_all_cs(wdir,folder,rep=50):
    Pre=[0]*rep
    Re=[0]*rep
    size=[0]*rep
    for i in range(rep):
        Pre[i],Re[i],size[i] = get_Pre_Re_S_cs(wdir,folder,i+1)
    return Pre,Re,size

def get_cs_from_PIP(pip,snp,Tv,coverage):
    pip_sort_index=np.argsort(-pip)
    allPIP = sum(pip)
    if allPIP==0:
        return 0,0,0
    idx = 1
    pipsum = pip[pip_sort_index[0:idx]]
    while pipsum<coverage*allPIP:
        pipsum = pipsum + pip[pip_sort_index[idx]]
        idx = idx+1
    cs = snp[pip_sort_index[0:idx]]
    Pre = len(set(snp[Tv]) & set(cs)) / len(cs)
    Re = len(set(snp[Tv]) & set(cs)) / sum(Tv)
    size = idx
    return Pre,Re,size

def get_cs_K_W(wdir,K,W,SNPlist,rep=50):
    SP_no_cs=get_all_cs(wdir,'sparsepro_no')
    if os.path.exists(os.path.join(wdir,'sparsepro_anno/C1.txt.pip')):
        SP_anno_cs=get_all_cs(wdir,'sparsepro_anno')
    else:
        SP_anno_cs=SP_no_cs
    susie_rss_cs=get_all_cs(wdir,'susie_rss')
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)
    finemap_PIP_cs=get_all_cs_from_PIP(finemap_PIP,Tv,SNPlist.values)
    PAINTOR_no_PIP_cs=get_all_cs_from_PIP(PAINTOR_no_PIP,Tv,SNPlist.values)
    PAINTOR_anno_PIP_cs=get_all_cs_from_PIP(PAINTOR_anno_PIP,Tv,SNPlist.values)
    df_res = pd.DataFrame({'Precision':finemap_PIP_cs[0]+susie_rss_cs[0]+PAINTOR_no_PIP_cs[0]+PAINTOR_anno_PIP_cs[0]+SP_no_cs[0]+SP_anno_cs[0],
                          'Recall':finemap_PIP_cs[1]+susie_rss_cs[1]+PAINTOR_no_PIP_cs[1]+PAINTOR_anno_PIP_cs[1]+SP_no_cs[1]+SP_anno_cs[1],
                          'Size':finemap_PIP_cs[2]+susie_rss_cs[2]+PAINTOR_no_PIP_cs[2]+PAINTOR_anno_PIP_cs[2]+SP_no_cs[2]+SP_anno_cs[2],
                           'Method':['FINEMAP']*rep+['SuSiE']*rep+['PAINTOR-']*rep+['PAINTOR+']*rep+['SparsePro-']*rep+['SparsePro+']*rep
                          })
    df_res['K']=K
    df_res['W']=W

    return df_res

def get_all_cs_from_PIP(allPIP,Tv,snp,coverage=0.95,rep=50):
    allPIPs = np.split(allPIP,rep)
    Tvs = np.split(Tv,rep)
    Pre=[0]*rep
    Re=[0]*rep
    size=[0]*rep
    for i in range(rep):
        Pre[i],Re[i],size[i] = get_cs_from_PIP(allPIPs[i],snp,Tvs[i],coverage)
    return Pre,Re,size

def plot_cs_K_W(df_res):
    fig, ax = plt.subplots(3,1,figsize=(10,10))
    sns.violinplot(data=df_res,x="Method", y="Recall",palette=['lime','royalblue','orchid','darkmagenta','orange','red'],ax=ax[0])
    sns.violinplot(data=df_res,x="Method", y="Precision",palette=['lime','royalblue','orchid','darkmagenta','orange','red'],ax=ax[1])
    sns.violinplot(data=df_res,x="Method", y="Size",palette=['lime','royalblue','orchid','darkmagenta','orange','red'],ax=ax[2])
    return df_res.groupby(['Method']).mean()


CHR = 22
ST = 31000000
ED = 32000000
SNPlist = pd.read_csv('../sim/loc/{}_{}_{}/K1_W0/sparsepro_no/C1.txt.pip'.format(CHR,ST,ED),sep='\t',header=None)[0]
K = 5
W = 2

ant_plot_no('../sim/loc/{}_{}_{}/K{}_W{}'.format(CHR,ST,ED,K,W),K,SNPlist)
calibration_plot_no('../sim/loc/{}_{}_{}/K{}_W{}'.format(CHR,ST,ED,K,W),K,SNPlist)
prc_plot_no('../sim/loc/{}_{}_{}/K{}_W{}'.format(CHR,ST,ED,K,W),K,SNPlist)
get_cs_K_W('../sim/loc/{}_{}_{}/K{}_W{}'.format(CHR,ST,ED,K,W),K,W,SNPlist).groupby(['Method']).mean()


def get_auprc(PIP,Tv):
    precision, recall, thresholds = precision_recall_curve(Tv, PIP)
    return auc(recall, precision)

def get_all_auprc(wdir):
    SNPlist = pd.read_csv(os.path.join(wdir,'sparsepro_no/C1.txt.pip'),sep='\t',header=None)[0]
    Tv = get_Tv(SNPlist,wdir)
    finemap_PIP = load_finemap_PIP(wdir)
    susie_rss_PIP = load_susie_rss_PIP(wdir)
    SP_no_PIP = load_sparsepro_PIP(wdir)
    if not os.path.exists(os.path.join(wdir,'sparsepro_anno/C1.txt.pip')):
        SP_anno_PIP=SP_no_PIP
    else:
        SP_anno_PIP = load_sparsepro_anno_PIP(wdir)
    PAINTOR_no_PIP = load_PAINTOR_no_PIP(wdir)
    PAINTOR_anno_PIP = load_PAINTOR_anno_PIP(wdir)
    return [get_auprc(i,Tv) for i in [finemap_PIP,susie_rss_PIP,SP_no_PIP,SP_anno_PIP,PAINTOR_no_PIP,PAINTOR_anno_PIP]]

df_prc=[get_all_auprc('../sim/loc/{}_{}_{}/K{}_W{}'.format(CHR,ST,ED,K,W)) for CHR in [20,21,22] for K in [1,2,5,10] for W in [0,1,2]]
df_CKW=[(CHR,K,W)for CHR in [20,21,22] for K in [1,2,5,10] for W in [0,1,2]]
df_all_auprc=pd.concat([pd.DataFrame(df_CKW),pd.DataFrame(df_prc)],axis=1)
df_all_auprc.columns = ['CHR','K','W','FINEMAP','SuSiE','SparsePro-','SparsePro+','PAINTOR-','PAINTOR+']
df_all_auprc.groupby(['K','W','CHR']).mean()
df_all_auprc.to_csv('loc_auprc.csv',header=True,index=False,sep='\t')
df_all_auprc['CHR'] = df_all_auprc['CHR'].replace(20,'Locus1').replace(21,'Locus2').replace(22,'Locus3')
df_all_auprc.to_csv('loc_auprc.txt',sep='\t',header=True,index=False)



CHR = 22
ST = 31000000
ED = 32000000
SNPlist = pd.read_csv('../sim/loc/{}_{}_{}/K1_W0/SP1_no/C1.txt.pip'.format(CHR,ST,ED),sep='\t',header=None)[0]
SNPlist21 = pd.read_csv('../sim/loc/{}_{}_{}/K1_W0/SP1_no/C1.txt.pip'.format(21,ST,ED),sep='\t',header=None)[0]
SNPlist20 = pd.read_csv('../sim/loc/{}_{}_{}/K1_W0/SP1_no/C1.txt.pip'.format(20,ST,ED),sep='\t',header=None)[0]
df_cs_22=pd.concat([get_cs_K_W('../sim/loc/{}_{}_{}/K{}_W{}'.format(CHR,ST,ED,K,W),K,W,SNPlist) for K in [1,2,5,10] for W in [0,1,2]])
df_cs_21=pd.concat([get_cs_K_W('../sim/loc/{}_{}_{}/K{}_W{}'.format(21,ST,ED,K,W),K,W,SNPlist21) for K in [1,2,5,10] for W in [0,1,2]])
df_cs_20=pd.concat([get_cs_K_W('../sim/loc/{}_{}_{}/K{}_W{}'.format(20,ST,ED,K,W),K,W,SNPlist20) for K in [1,2,5,10] for W in [0,1,2]])
df_cs_22['CHR']=22
df_cs_21['CHR']=21
df_cs_20['CHR']=20
df_cs_all=pd.concat([df_cs_22,df_cs_21,df_cs_20])
df_cs_all.to_csv('loc_cs.csv',sep='\t',index=False,header=True)



df_CKWM=[(CHR,K,W,M)for CHR in [20,21,22] for K in [1,2,5,10] for W in [0,1,2] for M in ['finemap','susie_rss','PAINTOR_no','PAINTOR_anno','sparsepro_no','sparsepro_anno']]
df_t = [pd.read_csv('../sim/loc/{}_{}_{}/K{}_W{}/{}.time'.format(CHR,ST,ED,K,W,M),sep='\t',header=None).at[0,0] for CHR in [20,21,22] for K in [1,2,5,10] for W in [0,1,2] for M in ['finemap','susie_rss','PAINTOR_no','PAINTOR_anno','sparsepro_no','sparsepro_anno']]
df_time=pd.concat([pd.DataFrame(df_CKWM),pd.DataFrame(df_t)],axis=1)
df_time.columns = ['CHR','K','W','Method','time']
df_time[~((df_time['Method']=='sparsepro_anno') & (df_time['W']==0))].to_csv('loc_time.csv',sep='\t',index=False,header=True)



df_wsep_22 = pd.concat([pd.read_csv('../sim/loc/22_31000000_32000000/K{}_W{}/sparsepro_no/sparsepro_no.wsep'.format(K,W),sep='\t') for K in [1,2,5,10] for W in [0,1,2]]).reset_index()[['Unnamed: 0','W','se','p']]
df_wsep_21 = pd.concat([pd.read_csv('../sim/loc/21_31000000_32000000/K{}_W{}/sparsepro_no/sparsepro_no.wsep'.format(K,W),sep='\t') for K in [1,2,5,10] for W in [0,1,2]]).reset_index()[['Unnamed: 0','W','se','p']]
df_wsep_20 = pd.concat([pd.read_csv('../sim/loc/20_31000000_32000000/K{}_W{}/sparsepro_no/sparsepro_no.wsep'.format(K,W),sep='\t') for K in [1,2,5,10] for W in [0,1,2]]).reset_index()[['Unnamed: 0','W','se','p']]
df_KW=pd.DataFrame(np.tile([[K,W] for K in [1,2,5,10] for W in [0,1,2]],10).reshape(120,2))
pd.concat([df_wsep_22,df_KW],axis=1).to_csv('loc_sparsepro_anno_22.csv',sep='\t',header=False,index=False)
pd.concat([df_wsep_21,df_KW],axis=1).to_csv('loc_sparsepro_anno_21.csv',sep='\t',header=False,index=False)
pd.concat([df_wsep_20,df_KW],axis=1).to_csv('loc_sparsepro_anno_20.csv',sep='\t',header=False,index=False)

df_ev_22 = pd.concat([pd.read_csv('../sim/loc/22_31000000_32000000/K{}_W{}/PAINTOR_anno/Enrichment.Values'.format(K,W),sep=' ') for K in [1,2,5,10] for W in [0,1,2]]).reset_index()
df_ev_21 = pd.concat([pd.read_csv('../sim/loc/21_31000000_32000000/K{}_W{}/PAINTOR_anno/Enrichment.Values'.format(K,W),sep=' ') for K in [1,2,5,10] for W in [0,1,2]]).reset_index()
df_ev_20 = pd.concat([pd.read_csv('../sim/loc/20_31000000_32000000/K{}_W{}/PAINTOR_anno/Enrichment.Values'.format(K,W),sep=' ') for K in [1,2,5,10] for W in [0,1,2]]).reset_index()
df_KW=pd.DataFrame([[K,W] for K in [1,2,5,10] for W in [0,1,2]])
pd.concat([df_ev_22,df_KW],axis=1).to_csv('loc_paintor_anno_22.csv',sep='\t',header=False,index=False)
pd.concat([df_ev_21,df_KW],axis=1).to_csv('loc_paintor_anno_21.csv',sep='\t',header=False,index=False)
pd.concat([df_ev_20,df_KW],axis=1).to_csv('loc_paintor_anno_20.csv',sep='\t',header=False,index=False)

