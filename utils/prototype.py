#!/usr/bin/env python3
"""Record of two interactive regression analyses from August 2021.

1. Random forest PPV run
   Produced output/stats/rf_ppv_aug2021.tsv, which is no longer in the tree. This is the
   only record of its parameters: Rank_NP <= 1, max 10 contexts per peptide, SGD/e4000
   scores, 12 replicates, with the BLCL train/test split taken from CAMAP's
   peptides_dataset_162. The superseding table is
   output/metrics/gGRCh37.75...maxRedSGD.../rf-full/ppv.tsv

2. TPM reduction prototype
   The prototype of what became max_reduce in camaptools: Dataset.construct_datasets(
   max_reduce=True) and _create_df's get_max(), plus the MinMaxScaler and 'Base' column
   handling. This is the only record of the sum-versus-max comparison behind that choice
   (the "maxRed" in output/metrics/gGRCh37.75...maxRedSGD..., commit 1544578).

"""

import pickle as pkl

import pandas as pd
from sklearn.preprocessing import MinMaxScaler

from camaptools.Dataset import RegressionDataset
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor
from camaptools.Regression import RegressionManager

from ppv import do


# ######################## #
# 1. Random forest PPV run #
# ######################## #

def get_blcl_peptides():
    """BLCL train/test split as used for the paper: validation folded into train."""
    train = []
    test = []
    base = '../../CAMAP/peptides_dataset_162'
    for subset in ['train', 'validation']:
        for cls in ['source', 'nonsource']:
            with open('%s/peptides.%s.%s.txt' % (base, cls, subset), 'r') as f:
                train.extend([x.strip() for x in f.readlines() if x.strip()])
    for cls in ['source', 'nonsource']:
        with open('%s/peptides.%s.test.txt' % (base, cls), 'r') as f:
            test.extend([x.strip() for x in f.readlines() if x.strip()])
    return train, test


def rf_ppv(workers=40, out='output/stats/rf_ppv_aug2021.tsv'):
    dat = RegressionDataset('GRCh37.75', 'BLCL', 162)
    dat.workers = workers
    dat.executor = EnhancedProcessPoolExecutor

    train, test = get_blcl_peptides()
    dat.construct_datasets_options(train, test)
    dat.load_peptides_options(
        max_bs_or_rank=1,
        var='Rank_NP',
        max_contexts=10,
        step='evaluateDS',
        ann_method='SGD',
        ann_params='e4000',
        max_replicates=12,
    )

    dat.load_peptides()
    dat.construct_datasets()
    dat.clear_unused()

    rem = RegressionManager(dat, workers, EnhancedProcessPoolExecutor)
    rem.initialize_trainers(['rf'])
    rem.start()
    rem.join()

    df = do(rem.results['rf-full'].df_scores.loc['Human_BLCL'].loc['HLA-MinScore'], 12)
    df.to_csv(out, sep='\t')
    return rem


# ########################## #
# 2. TPM reduction prototype #
# ########################## #

def reduce(dat, how='sum'):
    """Collapse the multiple TPM rows per peptide into one, then renormalize.

    'max' is what camaptools adopted as max_reduce; 'sum' was the alternative tested.
    """
    df_train, df_test = dat.datasets[0]

    def reduce_tpm(df):
        df = df.drop([x.replace('Base', '') for x in df.columns if x.endswith('Base')], axis=1)
        df.columns = df.columns.str.replace('Base', '')
        df = df.groupby(['Peptide', 'Class', 'BS', 'CAMAP', 'CAMAPShuff'], sort=False)
        if how == 'sum':
            df = df.sum()
        elif how == 'max':
            df = df.max()
        df = df.reset_index().set_index(['Peptide', 'Class'])
        return df

    def normalize_columns(df):
        df = df.copy()
        columns = df.columns
        df.columns = [x + 'Base' for x in columns]
        df_norm = pd.DataFrame(scaler.transform(df), columns=columns, index=df.index)
        df = pd.concat([df_norm, df], axis=1)
        df = df.reset_index().set_index('Peptide')
        return df

    df = reduce_tpm(df_train)
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(df)
    df_train = normalize_columns(df)
    df_test = normalize_columns(reduce_tpm(df_test))
    dat.datasets[0] = (df_train, df_test)


def compare_reductions(how='max', dat_pkl='dat.p', workers=3):
    """Run the regression on a reduced dataset, to compare 'sum' against 'max'."""
    dat = pkl.load(open(dat_pkl, 'rb'))
    reduce(dat, how)

    rem = RegressionManager(dat, workers, EnhancedProcessPoolExecutor)
    rem.initialize_trainers(['mlp'])
    rem.start()
    rem.join()

    return do(rem.results['mlp-full'].df_scores.loc['Human_BLCL'].loc['HLA-MinScore'], 1)
