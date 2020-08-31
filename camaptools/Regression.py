#!/usr/bin/env python3

import os
from collections import defaultdict
import numpy as np
import pandas as pd
from datetime import datetime
import json
import gzip
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegressionCV
from sklearn import metrics

from camaptools.MLP import Dataset
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor

# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.classification_report.html : Note that in binary classification, recall of the positive class is also known as “sensitivity”; recall of the negative class is “specificity”.


OUTPUT_FOLDER = "./output"


class RegressionDataset(Dataset):
    def organize_data(self):
        """
        self.load_peptides needs to be executed for all required variables to get populated
        """
        df_dictionnaries = [[], []]

        for sns in [0, 1]:
            df_dct = {
                'Peptide': [],
                'logBS': [],
                'logTPM': [],
                'ANNScore': defaultdict(list),
                'ANNShuffleScore': defaultdict(list)
                }

            for pep in self.peptides[sns]:
                p = pep[0]
                expressions = self.pepexpr[pep]
                bs_score = self.pepbs[p]
                ann_score = self.peptides[sns][pep][self.ann_method]
                ann_shuf_score = self.peptides[sns][pep][self.ann_method + '_Shuffle']
                for e in expressions:
                    df_dct['Peptide'].append(p)
                    df_dct['logBS'].append(bs_score)
                    df_dct['logTPM'].append(np.log10(e))
                    for i, v in enumerate(ann_score):
                        df_dct['ANNScore'][i].append(v)
                    for i, v in enumerate(ann_shuf_score):
                        df_dct['ANNShuffleScore'][i].append(v)

            df_dct['ANNScore'] = [df_dct['ANNScore'][k] for k in sorted(df_dct['ANNScore'])]
            df_dct['ANNShuffleScore'] = [df_dct['ANNShuffleScore'][k] for k in sorted(df_dct['ANNShuffleScore'])]

            assert len(df_dct['ANNScore']) == len(df_dct['ANNShuffleScore'])
            n_replicates = len(df_dct['ANNScore'])

            df_dct['Peptide'] = [df_dct['Peptide']]*n_replicates
            df_dct['logBS'] = [df_dct['logBS']]*n_replicates
            df_dct['logTPM'] = [df_dct['logTPM']]*n_replicates

            for r in range(n_replicates):
                subdct = {c:df_dct[c][r] for c in df_dct}
                df_dictionnaries[sns].append(subdct)

        self.df_dictionnaries = df_dictionnaries


    def construct_datasets(self):
        def split_df(df, seed=0):
            ix_train, ix_test = train_test_split(df.index.unique(), test_size=0.3, random_state=seed)
            df_train = df.loc[ix_train]
            df_test = df.loc[ix_test]
            return df_train, df_test

        self.datasets = []

        for i, (ns, s) in enumerate(zip(*self.df_dictionnaries)):
            seed = i+1

            # Merge and normalize data
            df_ns = pd.DataFrame(ns).set_index('Peptide')
            df_ns['Class'] = False
            df_s = pd.DataFrame(s).set_index('Peptide')
            df_s['Class'] = True

            df = pd.concat([df_ns, df_s], axis=0)
            for x in ['logBS', 'logTPM', 'ANNScore', 'ANNShuffleScore']:
                df[x] = df[x] - min(df[x])
                df[x] = df[x] / max(df[x])

            # Unmerge normalized data
            df_ns = df[~df.Class]
            df_s = df[df.Class]

            # Split source and non-source datasets
            df_ns_training, df_ns_test = split_df(df_ns, seed)
            df_s_training, df_s_test = split_df(df_s, seed)

            # Re-format into training and test datasets
            df_training = pd.concat([df_ns_training, df_s_training], axis=0)
            df_test = pd.concat([df_ns_test, df_s_test], axis=0)

            self.datasets.append((df_training, df_test))


class RegressionManager(object):
    def __init__(self, regression_dataset, workers, mpi=False):
        self.dataset = regression_dataset
        self.workers = workers
        self.executor = EnhancedMPIPoolExecutor if mpi else EnhancedProcessPoolExecutor
        self.trainers = []
        self.processes = []


    def initialize_trainers(self):
        for i, dat in enumerate(self.dataset.datasets):
            seed = i+1
            self.trainers.append(RegressionTrainer(dat, seed))


    def start(self):
        self.ex = self.executor(max_workers=self.workers)
        for t in self.trainers:
            p = self.ex.submit(t.run)
            self.processes.append(p)
        print(self.ex)
        self.ex.run()  # takes care of threads in main


    def join(self):
        self.trainers = [p.result() for p in self.processes]

        results_dct = defaultdict(list)
        results_max_dct = defaultdict(list)
        coef_dct = defaultdict(list)
        auc_data = defaultdict(list)
        df_scores = pd.DataFrame()

        allele = list(self.dataset.alleles)[0] if len(self.dataset.alleles) == 1 else 'HLA-MinScore'

        for i, t in enumerate(self.trainers):
            replicate = i+1

            for k in t.results_dct:
                results_dct[k].append(t.results_dct[k])
            for k in t.results_max_dct:
                results_max_dct[k].append(t.results_max_dct[k])
            for k in t.coef_dct:
                coef_dct[k].append(t.coef_dct[k])
            for k in t.auc_data:
                auc_data[k].append(t.auc_data[k])
            df_sc = t.df_scores.copy()
            df_sc.columns = pd.Index([(x[0], allele, x[1], replicate) for x in df_sc.columns])
            df_scores = pd.concat([df_scores, df_sc], axis=0)

        df_results = pd.DataFrame(results_dct).sort_index(axis=1)
        df_results.index = df_results.index + 1
        df_results.columns = pd.Index([(x[0], allele, x[1]) for x in df_results.columns])
        df_results.index = pd.Index([(self.dataset.dataset_name, x) for x in df_results.index])

        df_results_max = pd.DataFrame(results_dct).sort_index(axis=1)
        df_results_max.index = df_results_max.index + 1
        df_results_max.columns = pd.Index([(x[0], allele, x[1]) for x in df_results_max.columns])
        df_results_max.index = pd.Index([(self.dataset.dataset_name, x) for x in df_results_max.index])

        auc_data = dict(auc_data)

        dct = {}
        for prefix, val_lst in coef_dct.items():
            for rep, values in enumerate(val_lst):
                rep += 1
                minidct = {var: val for var, val in zip(values['variable'], values['coefficient'])}
                dct[(self.dataset.dataset_name, allele, prefix, rep)] = minidct
        coef_df = pd.DataFrame(dct).transpose()

        df_scores = df_scores.sort_index(axis=1)

        self.df_results = df_results
        self.df_results_max = df_results_max
        self.coef_df = coef_df
        self.auc_data = auc_data
        self.df_scores = df_scores

        self.ex.shutdown()


    def run(self, **kwargs):
        self.get_datasets()
        self.initialize_trainers(**kwargs)
        self.start()
        self.join()
        print('SUCCESS')


class RegressionTrainer(object):
    def __init__(self, dataset, seed):
        """ NOTE: Cannot be re-parallelized if using MPIPoolExecutor
        """
        self.dataset = dataset
        self.seed = seed

        self.dataset_name = '_name_'
        self.allele = '_allele_'
        self.ann_version = '_ann_version_'

        self.all_combinations = [
            ["logBS"],
            ["logTPM"],
            ["ANNShuffleScore"],
            ["ANNScore"],
            ["logBS", "logTPM"],
            ["logBS", "ANNShuffleScore"],
            ["logBS", "ANNScore"],
            ["logTPM", "ANNShuffleScore"],
            ["logTPM", "ANNScore"],
            ["logBS", "logTPM", "ANNShuffleScore"],
            ["logBS", "logTPM", "ANNScore"]
        ]


    def run(self):
        df_training, df_test = self.dataset

        regressions = []
        for x_labels in self.all_combinations:
            p = self.train_regression(x_labels, df_training, df_test, self.seed)
            regressions.append(p)

        results = {}
        results_max = {}
        coef_dct = {}
        auc_data = {}
        df_scores = pd.DataFrame()

        for r, rm, coef, auc, df_sc in regressions:
            results = {**results, **r}
            results_max = {**results_max, **rm}
            coef_dct = {**coef_dct, **coef}
            auc_data = {**auc_data, **auc}
            df_scores = pd.concat([df_scores, df_sc], axis=1)

        self.results_dct = results
        self.results_max_dct = results_max
        self.coef_dct = coef_dct
        self.auc_data = auc_data
        self.df_scores = df_scores

        return self


    @staticmethod
    def train_regression(x_labels, df_training, df_test, seed=0):
        # Which labels are used
        prefix = "-".join(x_labels)

        # Prepare return data with null values
        results = {
            ('Accuracy', prefix): np.nan,
            ('Precision', prefix): np.nan,
            ('Sensitivity', prefix): np.nan,
            ('Specificity', prefix): np.nan,
            ('AUC', prefix) : np.nan,
            ('MCC', prefix): np.nan
        }

        results_max = {
            ('Accuracy', prefix): np.nan,
            ('Precision', prefix): np.nan,
            ('Sensitivity', prefix): np.nan,
            ('Specificity', prefix): np.nan,
            ('AUC', prefix) : np.nan,
            ('MCC', prefix): np.nan
        }

        df_scores = df_test.copy()
        df_scores = df_scores.reset_index().set_index(['Peptide', 'Class'])
        df_scores = df_scores.drop([x for x in df_scores.columns if x != prefix], axis=1)
        df_scores.columns = pd.MultiIndex.from_product([df_scores.columns, ['BaseScores']])
        df_scores[(prefix, 'Scores')] = np.nan
        df_scores[(prefix, 'Predictions')] = np.nan
        df_scores[(prefix, 'ScoresReduced')] = np.nan
        df_scores[(prefix, 'PredictionsReduced')] = np.nan

        auc_dct = {'FPR': [], 'TPR': [], 'THRESHOLDS': [], 'FPR_max': [], 'TPR_max': [], 'THRESHOLDS_max': []}
        auc_dct = {prefix: auc_dct}

        coef_dct = {prefix: {'variable': x_labels, 'coefficient': [np.nan]*len(x_labels)}}

        if df_training.Class.sum() < 10 or (~df_training.Class).sum() < 10 or \
                df_test.Class.sum() < 10 or (~df_test.Class).sum() < 10:
            print('Less than 10 members in one class. Skipping...')

        else:
            # *Train*
            X, y = df_training[x_labels], df_training.Class
            clf = LogisticRegressionCV(cv=10, random_state=seed, solver='lbfgs',
                    multi_class='auto', class_weight='balanced', max_iter=100)
            #clf = SVC(gamma='auto', class_weight='balanced', probability=True)
            #clf = LogisticRegression(random_state=seed, multi_class='auto', class_weight='balanced')

            clf.fit(X, y)  # will not produce ValueError because we take care of cases where n < 10

            coefficients = clf.coef_[0]
            s = ": 10X cross validation... Regression coef %s" % coefficients
            print(prefix, s)

            coef_dct = {prefix: {'variable': x_labels, 'coefficient': list(coefficients)}}

            # *Test*
            X, y = df_test[x_labels], df_test.Class

            labels = y.to_numpy()
            scores = clf.predict_proba(X)[:,1]
            predictions = clf.predict(X)
            weights = abs(y - (y.sum() / y.shape[0]))
            assert np.array_equal(scores >= 0.5, predictions)

            # *Adjust test DF IDs*
            df_scores = df_test.copy()
            df_scores = df_scores.reset_index().set_index(['Peptide', 'Class'])
            df_scores = df_scores.drop([x for x in df_scores.columns if x != prefix], axis=1)
            df_scores.columns = pd.MultiIndex.from_product([df_scores.columns, ['BaseScores']])

            # *Append to test DF*
            df_scores[(prefix, 'Scores')] = scores
            df_scores[(prefix, 'Predictions')] = predictions

            # *Calculate metrics*

            # 1) Conventional method

            tp = (labels & predictions).sum()
            tn = (~labels & ~predictions).sum()
            fp = (~labels & predictions).sum()
            fn = (labels & ~predictions).sum()

            all_pos = labels.sum()
            all_neg = (~labels).sum()
            assert all_pos == tp + fn
            assert all_neg == tn + fp

            accuracy = clf.score(X, y, weights)
            precision, _, _, _ = metrics.precision_recall_fscore_support(
                    labels, predictions, sample_weight=weights,
                    labels=None, pos_label=1, average='binary', warn_for=())
            sensitivity = tp / (tp + fn)
            specificity = tn / (tn + fp)
            fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
            auc = metrics.auc(fpr, tpr)
            mcc = metrics.matthews_corrcoef(labels, predictions, sample_weight=weights)
            # NOTE on MCC
            #factor = len(labels[~labels])//len(labels[labels])
            #predictions_w = np.concatenate([np.repeat(predictions[labels], factor), predictions[~labels]])
            #labels_w = np.concatenate([np.repeat(labels[labels], factor), labels[~labels]])
            # ==> Calculating MCC with weights is identical to weighting predictions and labels beforehand

            results[('Accuracy', prefix)] = accuracy
            results[('Precision', prefix)] = precision
            results[('Sensitivity', prefix)] = sensitivity
            results[('Specificity', prefix)] = specificity
            results[('AUC', prefix)] = auc
            results[('MCC', prefix)] = mcc

            auc_dct =  {'FPR': list(fpr), 'TPR': list(tpr), 'THRESHOLDS': list(thresholds)}

            # 2) Keep only max score per peptide and re-calculate metrics (same strategy as reduce_tables)

            df_max = df_test[['Class']].copy()
            df_max['Scores'] = scores
            types = df_max.dtypes
            df_max = df_max.groupby("Peptide", sort=False).apply(max)
            df_max = df_max.astype(types)

            labels = df_max.Class.to_numpy()
            scores = df_max.Scores.to_numpy()
            predictions = scores >= 0.5
            df_max['Predictions'] = predictions
            weights = abs(labels - (labels.sum() / labels.shape[0]))

            tp = (labels & predictions).sum()
            tn = (~labels & ~predictions).sum()
            fp = (~labels & predictions).sum()
            fn = (labels & ~predictions).sum()

            accuracy = ((predictions == labels) * weights).sum() / weights.sum()
            precision, _, _, _ = metrics.precision_recall_fscore_support(
                    labels, predictions, sample_weight=weights,
                    labels=None, pos_label=1, average='binary', warn_for=())
            sensitivity = tp / (tp + fn)
            specificity = tn / (tn + fp)
            fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
            auc = metrics.auc(fpr, tpr)
            mcc = metrics.matthews_corrcoef(labels, predictions, sample_weight=weights)

            results_max[('Accuracy', prefix)] = accuracy
            results_max[('Precision', prefix)] = precision
            results_max[('Sensitivity', prefix)] = sensitivity
            results_max[('Specificity', prefix)] = specificity
            results_max[('AUC', prefix)] = auc
            results_max[('MCC', prefix)] = mcc

            auc_dct_max =  {'FPR_max': list(fpr), 'TPR_max': list(tpr), 'THRESHOLDS_max': list(thresholds)}

            # *Merge AUC values*
            auc_dct = {prefix: {**auc_dct, **auc_dct_max}}

            # *Append to test DF*
            df_max = df_max.reset_index().set_index(['Peptide', 'Class'])
            df_max = df_max.rename({'Scores':'ScoresReduced', 'Predictions':'PredictionsReduced'}, axis=1)
            df_max.columns = pd.MultiIndex.from_product([[prefix], df_max.columns])
            df_scores = df_scores.merge(df_max, on=['Peptide', 'Class'])

            # *Sort test DF columns*
            df_scores = df_scores.sort_index(axis=1)

        return results, results_max, coef_dct, auc_dct, df_scores
