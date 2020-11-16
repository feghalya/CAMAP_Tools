#!/usr/bin/env python3

import pickle as pkl
import gzip
from sklearn import metrics
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor

from camaptools.MLP import Model, load_models


MODELS = load_models(filters='SGD')


def calc_auc(seed, homology=10, copy='', shuffle='', dataset='BLCL'):
    if homology == 10:
        suffix = ''
        ncont = 10
    elif homology == 3:
        suffix = 'm3'
        ncont = 3
    elif homology == 0:
        suffix = 'mINF'
        ncont = 0
    else:
        sys.write.stderr('Homology value not recognized\n')
        sys.exit(1)

    if copy == 'TPM':
        suffix_2 = 'y'
        copy_dist = '_sameTPM'
    elif copy == 'BS':
        suffix_2 = 'z'
        copy_dist = '_sameBS'
    else:
        suffix_2 = ''
        copy_dist = ''

    m = MODELS[1][162]['SGD']['e4000pep9t5%s%s' % (suffix, suffix_2)][seed]

    if dataset == 'BLCL':
        dsname = 'BLCL_GRCh37.75_padding162_maxBS500%s%d_ratio5_peplen9%s%s%d_encoding-%s.pkl.gz' % (
                '_maxContexts', ncont, copy_dist, '_seed', seed+1, 'Codon' + shuffle + 'Embeddings')
        p = pkl.load(gzip.open('../output/trainDS/%s' % dsname), 'rb'))
        labels = [0]*len(p['test'][0]) + [1]*len(p['test'][1])
        scores = m.get_score(p['test'][0] + p['test'][1])[:,1].tolist()
        test_size = len(p['test'][1])
        train_size = len(p['train'][1])
        validation_size = len(p['validation'][1])
    else:
        if dataset == 'EL4' or dataset == 'CT26':
            species = 'GRCm38.78'
        else:
            species = 'GRCh37.75'
        dsname = '%s_%s_padding162_maxBS500_maxContexts0_ratio0%s%d_encoding-%s.pkl.gz' % (
                dataset, species, '_seed', seed+1, 'Codon' + shuffle + 'Embeddings')
        p = pkl.load(gzip.open('../output/trainDS/%s' % dsname), 'rb'))
        ns = p['train'][0] + p['validation'][0] + p['test'][0]
        s = p['train'][1] + p['validation'][1] + p['test'][1]
        labels = [0]*len(ns) + [1]*len(s)
        scores = m.get_score(ns + s)[:,1].tolist()
        test_size = len(s) + len(ns)
        train_size = -1
        validation_size = -1
        
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
    auc = metrics.auc(fpr, tpr)

    return auc, test_size, train_size, validation_size


def res(l):
    return tuple(zip(*[f.result() for f in l]))


def run(dataset):
    with ProcessPoolExecutor() as ex:
        c10 = [ex.submit(calc_auc, seed, 10, dataset=dataset) for seed in range(12)]
        c10s = [ex.submit(calc_auc, seed, 10, '', 'Shuffle', dataset=dataset) for seed in range(12)]

        c3 = [ex.submit(calc_auc, seed, 3, dataset=dataset) for seed in range(12)]
        c3s = [ex.submit(calc_auc, seed, 3, '', 'Shuffle', dataset=dataset) for seed in range(12)]

        c0 = [ex.submit(calc_auc, seed, 0, dataset=dataset) for seed in range(12)]
        c0s = [ex.submit(calc_auc, seed, 0, '', 'Shuffle', dataset=dataset) for seed in range(12)]

        samebs = [ex.submit(calc_auc, seed, 10, 'BS', '', dataset=dataset) for seed in range(12)]
        samebs_s = [ex.submit(calc_auc, seed, 10, 'BS', 'Shuffle', dataset=dataset) for seed in range(12)]

        sametpm = [ex.submit(calc_auc, seed, 10, 'TPM', '', dataset=dataset) for seed in range(12)]
        sametpm_s = [ex.submit(calc_auc, seed, 10, 'TPM', 'Shuffle', dataset=dataset) for seed in range(12)]

    # ** Homology **

    labels = ['C3']*12 + ['C10']*12 + ['CINF']*12
    labels = labels + [l + 'Shuffle' for l in labels]
    seeds = [x+1 for x in range(12)]*6

    outputs = [c3, c10, c0, c3s, c10s, c0s]
    aucs, test_sizes, train_sizes, validation_sizes = tuple(zip(*[res(x) for x in outputs]))
    aucs = [y for x in aucs for y in x]
    test_sizes = [y for x in test_sizes for y in x]
    train_sizes = [y for x in train_sizes for y in x]
    validation_sizes = [y for x in validation_sizes for y in x]

    df_homology = pd.DataFrame({
        "Label": labels,
        "Seed": seeds,
        "AUC": aucs,
        "TestSize": test_sizes,
        "TrainSize": train_sizes,
        "ValidationSize": validation_sizes
        }).set_index('Label')
    df_homology['Dataset'] = dataset

    # ** Binning **

    labels = ['C10']*12 + ['C10sameTPM']*12 + ['C10sameBS']*12
    labels = labels + [l + 'Shuffle' for l in labels]
    seeds = [x+1 for x in range(12)]*6

    outputs = [c10, sametpm, samebs, c10s, sametpm_s, samebs_s]
    aucs, test_sizes, train_sizes, validation_sizes = tuple(zip(*[res(x) for x in outputs]))
    aucs = [y for x in aucs for y in x]
    test_sizes = [y for x in test_sizes for y in x]
    train_sizes = [y for x in train_sizes for y in x]
    validation_sizes = [y for x in validation_sizes for y in x]

    df_binning = pd.DataFrame({
        "Label": labels,
        "Seed": seeds,
        "AUC": aucs,
        "TestSize": test_sizes,
        "TrainSize": train_sizes,
        "ValidationSize": validation_sizes
        }).set_index('Label')
    df_binning['Dataset'] = dataset

    return df_homology, df_binning


homologies = []
binnings = []
for dataset in ['BLCL', 'B721', 'PBMCs', 'EL4', 'CT26']:
    print('Running', dataset)
    df1, df2 = run(dataset)
    homologies.append(df1)
    binnings.append(df2)

homologies = pd.concat(homologies, axis=0)
binnings = pd.concat(binnings, axis=0)

homologies[homologies.Dataset == 'BLCL'].to_csv('revisions/homology.tsv', sep='\t', float_format='%g')
#^ only compare original vs shuffle
binnings[binnings.Dataset == 'BLCL'].to_csv('revisions/binning.tsv', sep='\t', float_format='%g')
#^ only compare original vs shuffle

homologies[homologies.Dataset != 'BLCL'].to_csv('revisions/homology.valid.tsv', sep='\t', float_format='%g')
#^ can compare strategies (ex: C10 vs CINF)
binnings[binnings.Dataset != 'BLCL'].to_csv('revisions/binning.valid.tsv', sep='\t', float_format='%g')
#^ can compare strategies (ex: C10 vs CINF)


# NOTES
#
# - sameBS and noBinning (C10) did not show any difference in expression.
# - sameTPM was consistently and very significantly worse than sameBS and noBinning.
# - sameTPMShuffle was also consistently worse than the other shuffle datasets but to a lesser extent
#   than what was observed with samTPM.
#
# - after testing against datasets containing all homologies versus limited to one (data not available):
#     - binning: using maxContexts0 (no limit) vs maxContexts1 didn't make any difference.
#     - homology: CINF in maxContexts1 performed better than it did in maxContexts0.
#       ==> This might indicate that a model trained with limited contexts is more robust to varying context sizes.
