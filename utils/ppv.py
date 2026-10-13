#!/usr/bin/env python

import pandas as pd
import numpy as np
import os


OUT = 'output/metrics/genomeGRCh37.75genomeGRCm38.78_padding162_maxRankNP1_maxContexts10_ANNMethodSGDe4000/full/scores'


def get_ppv(df, percentage=0.01):
    a = df.sort_values(ascending=False).index.get_level_values(1).values.astype(int).cumsum()
    max_hits = a[-1]
    if max_hits == 0:
        n = 0
    elif max_hits*percentage < 1:
        n = 1
    else:
        n = round(max_hits*percentage)
    try:
        ppv = a.tolist().index(n) + 1
    except ValueError:
        assert max_hits == 0
        ppv = 0
    return ppv


def calc_reg_table(dsname, out_dir, n_replicates=12):
    out_path = os.path.join(out_dir, '%s.scores.tsv.gz' % dsname)

    df = pd.read_csv(out_path, sep='\t', index_col=[0,1,2,3,4], header=[0,1])

    return do(df.loc[dsname].loc['HLA-MinScore'], n_replicates)


def do(subdf, n_replicates=12):
    try:
        subdf = subdf.loc[:,(slice(None), 'ScoresReduced')]
        print('ScoresReduced')
    except KeyError:
        try:
            subdf = subdf.loc[:,(slice(None), 'Scores_y')]
        except KeyError:
            subdf = subdf.loc[:,(slice(None), 'Scores')]
    subdf = subdf.droplevel(level=1, axis=1)
    subdf = subdf.reset_index().drop_duplicates().set_index(['N', 'Peptide', 'Class'])

    allentries = int(subdf.shape[0]/n_replicates)
    allhits = np.mean([subdf.loc[x].index.get_level_values(1).values.sum() for x in range(1, n_replicates+1)])
    allhits = int(round(allhits))
    #print('Entries', allentries, 'Hits', allhits)

    to_use = ['CAMAP', 'BS', 'TPM', 'BS-TPM', 'BS-TPM-CAMAP']
    #resdct = {'Percentage': [], 'Metric': [], 'BS': [], 'CAMAP': [], 'TPM': [], 'BS-CAMAP': [], 'BS-TPM': [], 'BS-TPM-CAMAP': []}
    #    for params in ['BS', 'CAMAP', 'TPM', 'BS-CAMAP', 'BS-TPM', 'BS-TPM-CAMAP']:
    resdct = {x: [] for x in ['Percentage', 'Metric'] + to_use}
    #for p in [0.01] + [x/100 for x in range(5, 100, 5)] + [0.99]:
    for p in [0.01, 0.05, 0.1, 0.5]:
        nhits = np.mean([subdf.loc[x].index.get_level_values(1).values.sum()*p for x in range(1, n_replicates+1)])
        nhits = int(round(nhits))
        percentage_hits = '%.f%% (%d)' % (p*100, nhits)
        resdct['Percentage'].extend([percentage_hits]*3)
        resdct['Metric'].extend(['n', 'std', 'fpr'])
        for params in to_use:
            l = [get_ppv(subdf.loc[x][params], p) for x in range(1, n_replicates+1)]
            val = int(round(np.mean(l)))
            std = int(round(np.std(l)))
            try:
                fpr = round((val - nhits)/val*100,3)
            except ZeroDivisionError:
                if nhits == 0:
                    fpr = np.nan
                else:
                    fpr = 100
            resdct[params].extend([val, std, fpr])
            #print(
            #    p,
            #    'Allhits', allhits,
            #    'Hits', nhits, params,
            #    'N', val,
            #    'Std', std,
            #    'Percentage', round(val/allentries*100, 3),
            #    'FPR', fpr
            #)

    out_df = pd.DataFrame(resdct).set_index(['Percentage', 'Metric']).T
    out_df.columns.names = ['', '']

    return out_df


if __name__ == '__main__':
    os.makedirs('output/ppv')

    dsname = 'Human_BLCL'
    df = calc_reg_table(dsname, OUT)
    df.to_csv('output/ppv/%s.tsv' % dsname, sep='\t', float_format='%g')

    dsname = 'Human_B721'
    df = calc_reg_table(dsname, OUT)
    df.to_csv('output/ppv/%s.tsv' % dsname, sep='\t', float_format='%g')

    dsname = 'Human_PBMCs'
    df = calc_reg_table(dsname, OUT)
    df.to_csv('output/ppv/%s.tsv' % dsname, sep='\t', float_format='%g')

    dsname = 'Mouse_EL4'
    df = calc_reg_table(dsname, OUT)
    df.to_csv('output/ppv/%s.tsv' % dsname, sep='\t', float_format='%g')

    dsname = 'Mouse_CT26'
    df = calc_reg_table(dsname, OUT)
    df.to_csv('output/ppv/%s.tsv' % dsname, sep='\t', float_format='%g')
