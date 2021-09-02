#!/usr/bin/env python

import pandas as pd
import os
from datetime import datetime

root = '../output/metrics/gGRCh37.75gGRCm38.78_pad162_maxRankNP1_maxCont10_maxRedSGD_ANNMethSGDe4000'
curtime = datetime.now().strftime("%Yy%mm%dd")

os.makedirs(curtime)

def wrapup(subdir):
    method = '-'.join(subdir.split('-')[:-1])

    fn = os.path.join(root, subdir, 'metrics.tsv')
    metrics = pd.read_csv(fn, sep='\t', index_col=[0,1,2], header=[0,1]).loc['HLA-MinScore']
    metrics = metrics.reset_index()
    metrics = metrics.melt(['Metric', 'Regression'], var_name=['Dataset', 'N'], value_name='value')
    metrics.to_csv(curtime + '/' + method + '.melted.tsv', sep='\t', index=False)

    fn = os.path.join(root, subdir, 'ppv.tsv')
    df = pd.read_csv(fn, sep='\t', index_col=[0,1,2,3], header=[0,1])
    metdf = df.loc['HLA-MinScore'].reset_index()
    metdf = metdf.melt(['Metric', 'Statistic', 'Regression'], var_name=['Dataset', 'N'], value_name='value')
    metdf = metdf.drop('N', axis=1).groupby(['Metric', 'Statistic', 'Regression', 'Dataset']).value.describe()
    metdf = metdf.reset_index().drop(['count', 'min', '25%', '50%', '75%', 'max'], axis=1)
    nentries = metdf.set_index('Metric').loc['nEntries'][['Dataset', 'mean']]
    nentries = nentries.drop_duplicates().rename({'mean': 'nEntries'}, axis=1).reset_index(drop=True)
    nentries.nEntries = nentries.nEntries.round().astype(int)
    nhits = metdf.set_index('Metric').loc['nHits'][['Dataset', 'mean']]
    nhits = nhits.drop_duplicates().rename({'mean': 'nHits'}, axis=1).reset_index(drop=True)
    nhits.nHits = nhits.nHits.round().astype(int)
    metdf = metdf.set_index('Metric').drop(['nEntries', 'nHits']).reset_index()
    finaldf = metdf.merge(nhits, on='Dataset').merge(nentries, on='Dataset')
    finaldf = finaldf.sort_values(['Dataset', 'Metric', 'Regression'])
    finaldf = finaldf.set_index(['Dataset', 'nEntries', 'nHits', 'Metric', 'Regression'])
    finaldf = finaldf.pivot(columns='Statistic')
    finaldf.columns = finaldf.columns.swaplevel()
    finaldf.columns = ['_'.join(x) for x in finaldf.columns.values]
    finaldf = finaldf.sort_index(axis=1).reset_index()
    finaldf = finaldf.drop('nhits_std', axis=1).rename({'nhits_mean': 'MetricHits'}, axis=1)
    finaldf.MetricHits = finaldf.MetricHits.astype(int)
    finaldf = finaldf.set_index(['Dataset', 'nEntries', 'nHits', 'Metric', 'MetricHits', 'Regression'])
    finaldf.to_csv(curtime + '/' + method + '.ppv.tsv', sep='\t', float_format='%g')


for subdir in os.listdir(root):
    wrapup(subdir)
