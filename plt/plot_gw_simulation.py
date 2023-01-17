import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
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

i2r = pd.read_csv('../sim/gw/22/22.dict',sep='\t',header=None,index_col=0)

def get_Tv_PIP_anno(wdir,ite,i2r):
    '''load PIP'''
    df_PIP=i2r.copy()
    df_PIP['Tv']= i2r[1].isin(pd.read_csv(os.path.join(wdir,'C{}.causal'.format(ite)),sep='\t',header=None)[0])
    df_PIP['SparsePro-']=0.0
    df_sp = pd.read_csv(os.path.join(wdir,'sparsepro/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_sp.index,'SparsePro-']=df_sp[2]
    df_PIP['SuSiE']=0.0
    df_su = pd.read_csv(os.path.join(wdir,'susie_ukb/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_su.index,'SuSiE']=df_su[2]
    df_PIP['SparsePro+']=0.0
    df_asp = pd.read_csv(os.path.join(wdir,'sparsepro_anno/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_sp.index,'SparsePro+']=df_asp[2]
    df_PIP['SuSiE+PolyFun']=0.0
    df_asp = pd.read_csv(os.path.join(wdir,'susie_poly/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_sp.index,'SuSiE+PolyFun']=df_asp[2]
    return df_PIP

def get_Tv_PIP_no(wdir,ite,i2r):
    '''load PIP'''
    df_PIP=i2r.copy()
    df_PIP['Tv']= i2r[1].isin(pd.read_csv(os.path.join(wdir,'C{}.causal'.format(ite)),sep='\t',header=None)[0])
    df_PIP['SparsePro-']=0.0
    df_sp = pd.read_csv(os.path.join(wdir,'sparsepro/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_sp.index,'SparsePro-']=df_sp[2]
    df_PIP['SuSiE']=0.0
    df_su = pd.read_csv(os.path.join(wdir,'susie_ukb/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_su.index,'SuSiE']=df_su[2]
    df_PIP['SuSiE+PolyFun']=0.0
    df_asp = pd.read_csv(os.path.join(wdir,'susie_poly/C{}.pip'.format(ite)),sep='\t',header=None,index_col=0)
    df_PIP.loc[df_sp.index,'SuSiE+PolyFun']=df_asp[2]
    return df_PIP

def ant_plot(ax,Tv,PIP1,PIP2,name1,name2):
    '''plot PIP contrast'''
    ax.scatter(x=PIP1[[not i for i in Tv]],y=PIP2[[not i for i in Tv]],c='black',s = 100,label='Non-causal')
    ax.scatter(x=PIP1[Tv],y=PIP2[Tv],c='red',s = 100,label='Causal')
    ax.set_xlabel(name1, fontsize = 25)
    ax.set_ylabel(name2, fontsize = 25)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,-0.3))
    ax.grid()
    return ax

def ant_plot_gw_anno(df_PIP):
    '''ant plot for genome-wide simulation with annotation enrichment'''
    fig, axs = plt.subplots(1,2,figsize=(20, 9))
    ant_plot(axs[0],df_PIP['Tv'],df_PIP['SparsePro-'],df_PIP['SuSiE'],'SparsePro- PIP','SuSiE PIP')
    ant_plot(axs[1],df_PIP['Tv'],df_PIP['SparsePro+'],df_PIP['SuSiE+PolyFun'],'SparsePro+ PIP','SuSiE+PolyFun PIP')

def ant_plot_gw_no(df_PIP):
    '''ant plot for genome-wide simulation without annotation enrichment'''
    fig, axs = plt.subplots(1,2,figsize=(20, 9))
    ant_plot(axs[0],df_PIP['Tv'],df_PIP['SparsePro-'],df_PIP['SuSiE'],'SparsePro- PIP','SuSiE PIP')
    ant_plot(axs[1],df_PIP['Tv'],df_PIP['SuSiE+PolyFun'],df_PIP['SuSiE'],'SuSiE+PolyFun PIP','SuSiE PIP')

def prc_plot(ax,PIP,Tv,name,color):
    precision, recall, thresholds = precision_recall_curve(Tv, PIP)
    ax.plot(recall, precision, label=name+": "+"{:.4f}".format(auc(recall, precision)),color=color)
    ax.set_ylabel('Precision')
    ax.set_xlabel('Recall')
    return ax

def prc_plot_anno(wdir,i2r):
    df_PIP=pd.concat([get_Tv_PIP_anno(wdir,i+1,i2r) for i in range(22)])
    fig, ax = plt.subplots(figsize=(10,10))
    prc_plot(ax,df_PIP['SparsePro-'],df_PIP['Tv'],'SparsePro-','orange')
    prc_plot(ax,df_PIP['SparsePro+'],df_PIP['Tv'],'SparsePro+','red')
    prc_plot(ax,df_PIP['SuSiE'],df_PIP['Tv'],'SuSiE','royalblue')
    prc_plot(ax,df_PIP['SuSiE+PolyFun'],df_PIP['Tv'],'SuSiE+PolyFun','darkblue')
    ax.legend()

def prc_plot_no(wdir,i2r):
    df_PIP=pd.concat([get_Tv_PIP_no(wdir,i+1,i2r) for i in range(22)])
    fig, ax = plt.subplots(figsize=(10,10))
    prc_plot(ax,df_PIP['SparsePro-'],df_PIP['Tv'],'SparsePro-','orange')
    prc_plot(ax,df_PIP['SuSiE'],df_PIP['Tv'],'SuSiE','royalblue')
    prc_plot(ax,df_PIP['SuSiE+PolyFun'],df_PIP['Tv'],'SuSiE+PolyFun','darkblue')
    ax.legend()

def calib_plot(ax,PIP,Tv,name,color):
    cp, mpv = calibration_curve(Tv,PIP,n_bins=20)
    ax.plot(mpv,cp,marker='.',label=name,color=color)
    ax.set_ylabel('Actual Precision')
    ax.set_xlabel('Expected Precision')
    return ax

def calibration_plot_anno(wdir,i2r):
    df_PIP=pd.concat([get_Tv_PIP_anno(wdir,i+1,i2r) for i in range(22)])
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot([0,1],[0,1],'--',color='grey')
    calib_plot(ax,df_PIP['SparsePro-'],df_PIP['Tv'],'SparsePro-','orange')
    calib_plot(ax,df_PIP['SparsePro+'],df_PIP['Tv'],'SparsePro+','red')
    calib_plot(ax,df_PIP['SuSiE'],df_PIP['Tv'],'SuSiE','royalblue')
    calib_plot(ax,df_PIP['SuSiE+PolyFun'],df_PIP['Tv'],'SuSiE+PolyFun','darkblue')
    ax.legend()

def calibration_plot_no(wdir,i2r):
    df_PIP=pd.concat([get_Tv_PIP_no(wdir,i+1,i2r) for i in range(22)])
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot([0,1],[0,1],'--',color='grey')
    calib_plot(ax,df_PIP['SparsePro-'],df_PIP['Tv'],'SparsePro-','orange')
    calib_plot(ax,df_PIP['SuSiE'],df_PIP['Tv'],'SuSiE','royalblue')
    calib_plot(ax,df_PIP['SuSiE+PolyFun'],df_PIP['Tv'],'SuSiE+PolyFun','darkblue')
    ax.legend()

def get_time(wdir,folder,ite):
    return pd.read_csv(os.path.join(wdir,'C{}_{}.time'.format(ite,folder)),header=None).loc[0,0]

def get_Pre_Re_S_cs(wdir,folder,ite,i2r):
    time = get_time(wdir,folder,ite)
    df_cs = pd.read_csv(os.path.join(wdir,folder,"C{}.cs".format(ite)),sep='\t',dtype={'pip':str})
    if len(df_cs)==0:
        return 0.0,0.0,0,time
    df_causal = pd.read_csv(os.path.join(wdir,"C{}.causal".format(ite)),sep='\t',header=None,index_col=0)
    causal = set(i2r.reset_index().set_index([1]).loc[df_causal.index].values.flatten())
    cs = [j for i in df_cs['cs'].values for j in i.split('/')]
    Pre = len(causal & set(cs)) / len(cs)
    Re = len(causal & set(cs)) / len(causal)
    size = len(cs)
    return Pre,Re,size,time

def get_all_cs(wdir,folder,i2r,rep=22):
    Pre=[0]*rep
    Re=[0]*rep
    size=[0]*rep
    time=[0]*rep
    for i in range(rep):
        Pre[i],Re[i],size[i],time[i] = get_Pre_Re_S_cs(wdir,folder,i+1,i2r)
    return Pre,Re,size,time

def get_all_cs_anno(wdir,i2r,W,rep=22):
    sparsepro_cs=get_all_cs(wdir,'sparsepro',i2r)
    sparsepro_anno_cs=get_all_cs(wdir,'sparsepro_anno',i2r)
    susie_cs=get_all_cs(wdir,'susie_ukb',i2r)
    susie_poly_cs=get_all_cs(wdir,'susie_poly',i2r)
    df_res = pd.DataFrame({'Precision':sparsepro_cs[0]+sparsepro_anno_cs[0]+susie_cs[0]+susie_poly_cs[0],
                          'Recall':sparsepro_cs[1]+sparsepro_anno_cs[1]+susie_cs[1]+susie_poly_cs[1],
                          'Size':sparsepro_cs[2]+sparsepro_anno_cs[2]+susie_cs[2]+susie_poly_cs[2],
                           'Time':sparsepro_cs[3]+sparsepro_anno_cs[3]+susie_cs[3]+susie_poly_cs[3],
                           'Method':['SparsePro-']*rep+['SparsePro+']*rep+['SuSiE']*rep+['SuSiE+PolyFun']*rep
                          })
    df_res['W']=W
    return df_res

def get_all_cs_no(wdir,i2r,W=0,rep=22):
    sparsepro_cs=get_all_cs(wdir,'sparsepro',i2r)
    susie_cs=get_all_cs(wdir,'susie_ukb',i2r)
    susie_poly_cs=get_all_cs(wdir,'susie_poly',i2r)
    df_res = pd.DataFrame({'Precision':sparsepro_cs[0]+susie_cs[0]+susie_poly_cs[0],
                          'Recall':sparsepro_cs[1]+susie_cs[1]+susie_poly_cs[1],
                          'Size':sparsepro_cs[2]+susie_cs[2]+susie_poly_cs[2],
                           'Time':sparsepro_cs[3]+susie_cs[3]+susie_poly_cs[3],
                           'Method':['SparsePro-']*rep+['SuSiE']*rep + ['SuSiE+PolyFun']*rep
                          })
    df_res['W']=W
    return df_res

def plot_cs_anno(df_res):
    fig, ax = plt.subplots(4,1,figsize=(10,10))
    sns.violinplot(data=df_res,x="Method", y="Recall",palette=["orange","red","royalblue","darkblue"],ax=ax[0])
    sns.violinplot(data=df_res,x="Method", y="Precision",palette=["orange","red","royalblue","darkblue"],ax=ax[1])
    sns.violinplot(data=df_res,x="Method", y="Size",palette=["orange","red","royalblue","darkblue"],ax=ax[2])
    sns.violinplot(data=df_res,x="Method", y="Time",palette=["orange","red","royalblue","darkblue"],ax=ax[3])
    return df_res.groupby(['Method']).mean()

def plot_cs_no(df_res):
    fig, ax = plt.subplots(4,1,figsize=(10,10))
    sns.violinplot(data=df_res,x="Method", y="Recall",palette=["orange","royalblue","darkblue"],ax=ax[0])
    sns.violinplot(data=df_res,x="Method", y="Precision",palette=["orange","royalblue","darkblue"],ax=ax[1])
    sns.violinplot(data=df_res,x="Method", y="Size",palette=["orange","royalblue","darkblue"],ax=ax[2])
    sns.violinplot(data=df_res,x="Method", y="Time",palette=["orange","royalblue","darkblue"],ax=ax[3])
    return df_res.groupby(['Method']).mean()

ant_plot_gw_anno(pd.concat([get_Tv_PIP_anno('../sim/gw/{}/K{}_W{}'.format(22,100,2),i+1,i2r) for i in range(22)]))
prc_plot_anno('../sim/gw/{}/K{}_W{}'.format(22,100,2),i2r)
calibration_plot_anno('../sim/gw/{}/K{}_W{}'.format(22,100,2),i2r)
plot_cs_anno(get_all_cs_anno('../sim/gw/{}/K{}_W{}'.format(22,100,2),i2r,2,rep=22))

