import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='convert bed to annotation matrix')
parser.add_argument('--save', type=str, default=None, help='path for output', required=True)
parser.add_argument('--anno', type=str, default=None, help='large anno file', required=True)
parser.add_argument('--bim', type=str, default=None, help='bim file', required=True)
args = parser.parse_args()

anno = pd.read_csv(args.anno,sep='\t')
anno.index=anno['SNP']
bim=pd.read_csv(args.bim,sep='\t',header=None)
anno.loc[anno.index.intersection(bim[1])].to_csv(args.save, sep='\t',header=True, index=False)
