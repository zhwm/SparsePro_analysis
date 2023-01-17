import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='get master file for finemap')
parser.add_argument('--save', type=str, default=None, help='path for output', required=True)
parser.add_argument('--N', type=int, default=None, help='sample size', required=True)
parser.add_argument('--rep', type=int, default=None, help='number of replication', required=True)
parser.add_argument('--ld',type=str, default=None, help='name of ld file')
parser.add_argument('--anno',type=str,default=None,help='name of annotation file')
args = parser.parse_args()

master = pd.DataFrame({'z':['C{}.z'.format(i) for i in range(1,args.rep+1)],
                       'ld':['{}'.format(args.ld)]*len(range(1,args.rep+1)),
                       'snp':['C{}.snp'.format(i) for i in range(1,args.rep+1)],
                       'config':['C{}.config'.format(i) for i in range(1,args.rep+1)],
                       'cred':['C{}.cred'.format(i) for i in range(1,args.rep+1)],
                       'log':['C{}.log'.format(i) for i in range(1,args.rep+1)],
                       'k':['C{}.k'.format(i) for i in range(1,args.rep+1)]})
master['n_samples']=args.N
master.to_csv(os.path.join(args.save,'master'),sep=';',index=False)

zld = pd.DataFrame({'z':['C{}.txt'.format(i) for i in range(1,args.rep+1)],'ld':['{}'.format(args.ld)]*args.rep,'anno':['{}'.format(args.anno)]*args.rep})
zld.to_csv(os.path.join(args.save,'zld'),sep='\t',index=False)
