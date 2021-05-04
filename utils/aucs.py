#!/usr/bin/env python3

import pickle as pkl
import gzip
from sklearn import metrics
import pandas as pd
import sys
import os
from concurrent.futures import ProcessPoolExecutor
import os

from camaptools.MLP import Model, load_models


MODELS = load_models(algorithms=['SGD', 'SGD_Shuffle', 'SGD_AminoAcid'])


def calc_auc(seed, homology=10, copy='', embeddings='original', dataset='BLCL'):
    if homology == 10:
        suffix_hom = ''
        ncont = 10
    elif homology == 3:
        suffix_hom = 'm3'
        ncont = 3
    elif homology == 0:
        suffix_hom = 'mINF'
        ncont = 0
    else:
        sys.write.stderr('Homology value not recognized\n')
        sys.exit(1)

    if copy == 'TPM':
        suffix_bin = 'y'
        copy_dist = '_sameTPM'
    elif copy == 'BS':
        suffix_bin = 'z'
        copy_dist = '_sameBS'
    else:
        suffix_bin = ''
        copy_dist = ''

    if embeddings == 'original':
        algorithm = 'SGD'
        embeddings = 'CodonEmbeddings'
    elif embeddings == 'shuffle':
        algorithm = 'SGD_Shuffle'
        embeddings = 'CodonShuffleEmbeddings'
    elif embeddings == 'aa':
        algorithm = 'SGD_AminoAcid'
        embeddings = 'AminoAcidEmbeddings'
    else:
        sys.write.stderr('Embeddings value not recognized\n')
        sys.exit(1)

    m = MODELS[1][162][algorithm]['e4000pep9t5%s%s' % (suffix_hom, suffix_bin)][seed]

    if dataset == 'BLCL':
        dsname = '%s_%s_padding162_maxBS500%s%d_ratio5_peplen9%s%s%d_encoding-%s.pkl.gz' % (
            'BLCL', 'GRCh37.75', '_maxContexts', ncont,
            copy_dist, '_seed', seed+1, embeddings
        )

        p = pkl.load(gzip.open('../output/trainDS/%s' % dsname, 'rb'))

        ns = p['test'][0]
        s = p['test'][1]

        labels = [0]*len(ns) + [1]*len(s)
        scores = m.get_score(ns + s)[:,1].tolist()

        test_size = len(p['test'][1])
        train_size = len(p['train'][1])
        validation_size = len(p['validation'][1])
    else:
        if dataset == 'EL4' or dataset == 'CT26':
            species = 'GRCm38.78'
        else:
            species = 'GRCh37.75'

        # maxContexts and ratio are set to 0 because these datasets act as validators
        dsname = '%s_%s_padding162_maxBS500_maxContexts0_ratio0%s%d_encoding-%s.pkl.gz' % (
            dataset, species, '_seed', seed+1, embeddings
        )

        p = pkl.load(gzip.open('../output/trainDS/%s' % dsname, 'rb'))

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


def run(dat):
    with ProcessPoolExecutor() as ex:
        hom = 10
        copy = ''
        c10 = [ex.submit(calc_auc, seed, hom, copy, 'original', dat) for seed in range(12)]
        c10s = [ex.submit(calc_auc, seed, hom, copy, 'shuffle', dat) for seed in range(12)]
        c10aa = [ex.submit(calc_auc, seed, hom, copy, 'aa', dat) for seed in range(12)]

        hom = 3
        copy = ''
        c3 = [ex.submit(calc_auc, seed, hom, copy, 'original', dat) for seed in range(12)]
        c3s = [ex.submit(calc_auc, seed, hom, copy, 'shuffle', dat) for seed in range(12)]

        hom = 0
        copy = ''
        c0 = [ex.submit(calc_auc, seed, hom, copy, 'original', dat) for seed in range(12)]
        c0s = [ex.submit(calc_auc, seed, hom, copy, 'shuffle', dat) for seed in range(12)]

        hom = 10
        copy = 'BS'
        samebs = [ex.submit(calc_auc, seed, hom, copy, 'original', dat) for seed in range(12)]
        samebs_s = [ex.submit(calc_auc, seed, hom, copy, 'shuffle', dat) for seed in range(12)]

        hom = 10
        copy = 'TPM'
        sametpm = [ex.submit(calc_auc, seed, hom, copy, 'original', dat) for seed in range(12)]
        sametpm_s = [ex.submit(calc_auc, seed, hom, copy, 'shuffle', dat) for seed in range(12)]

    # ** Shuffling **

    labels = ['C10']*12 + ['C10Shuffle']*12 + ['C10AminoAcid']*12
    seeds = [x+1 for x in range(12)]*3

    outputs = [c10, c10s, c10aa]
    results = [[f.result() for f in x] for x in outputs]
    results = [tuple(zip(*x)) for x in results]
    aucs, test_sizes, train_sizes, validation_sizes = tuple(zip(*results))

    aucs = [y for x in aucs for y in x]
    test_sizes = [y for x in test_sizes for y in x]
    train_sizes = [y for x in train_sizes for y in x]
    validation_sizes = [y for x in validation_sizes for y in x]

    df_shuffling = pd.DataFrame({
        "Label": labels,
        "Seed": seeds,
        "AUC": aucs,
        "TestSize": test_sizes,
        "TrainSize": train_sizes,
        "ValidationSize": validation_sizes
        }).set_index('Label')
    df_shuffling['Dataset'] = dataset

    # ** Homology **

    labels = ['C3']*12 + ['C10']*12 + ['CINF']*12
    labels = labels + [l + 'Shuffle' for l in labels]
    seeds = [x+1 for x in range(12)]*6

    outputs = [c3, c10, c0, c3s, c10s, c0s]
    results = [[f.result() for f in x] for x in outputs]
    results = [tuple(zip(*x)) for x in results]
    aucs, test_sizes, train_sizes, validation_sizes = tuple(zip(*results))
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
    df_homology['Dataset'] = dat

    # ** Binning **

    labels = ['C10']*12 + ['C10sameTPM']*12 + ['C10sameBS']*12
    labels = labels + [l + 'Shuffle' for l in labels]
    seeds = [x+1 for x in range(12)]*6

    outputs = [c10, sametpm, samebs, c10s, sametpm_s, samebs_s]
    results = [[f.result() for f in x] for x in outputs]
    results = [tuple(zip(*x)) for x in results]
    aucs, test_sizes, train_sizes, validation_sizes = tuple(zip(*results))
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
    df_binning['Dataset'] = dat

    return df_homology, df_binning, df_shuffling


homologies = []
binnings = []
shufflings = []
for dataset in ['BLCL', 'B721', 'PBMCs', 'EL4', 'CT26']:
    print('Running', dataset)
    df1, df2, df3 = run(dataset)
    homologies.append(df1)
    binnings.append(df2)
    shufflings.append(df3)

homologies = pd.concat(homologies, axis=0)
binnings = pd.concat(binnings, axis=0)
shufflings = pd.concat(shufflings, axis=0)

os.makedirs('revisions', exist_ok=True)

homologies[homologies.Dataset == 'BLCL'].to_csv('revisions/homology.tsv', sep='\t', float_format='%g')
binnings[binnings.Dataset == 'BLCL'].to_csv('revisions/binning.tsv', sep='\t', float_format='%g')
shufflings[shufflings.Dataset == 'BLCL'].to_csv('revisions/shuffling.tsv', sep='\t', float_format='%g')
#^ only compare original vs shuffle

homologies[homologies.Dataset != 'BLCL'].to_csv('revisions/homology.valid.tsv', sep='\t', float_format='%g')
binnings[binnings.Dataset != 'BLCL'].to_csv('revisions/binning.valid.tsv', sep='\t', float_format='%g')
shufflings[shufflings.Dataset != 'BLCL'].to_csv('revisions/shuffling.valid.tsv', sep='\t', float_format='%g')
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
