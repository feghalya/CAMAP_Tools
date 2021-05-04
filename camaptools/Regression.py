#!/usr/bin/env python3

import os
from collections import defaultdict, namedtuple
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn import metrics
from functools import partial
import gc
import pickle as pkl
import gzip
import json
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor

# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.classification_report.html : Note that in binary classification, recall of the positive class is also known as “sensitivity”; recall of the negative class is “specificity”.


Results = namedtuple('Results', ['df_results', 'df_results_max', 'df_coef', 'df_scores', 'auc_data'])


class RegressionMetaManager(object):
    """A wrapper around RegressionManagers that also takes care of running and distributing jobs for Datasets
    """
    def __init__(self, datasets, out_dir, workers=0, executor=None):
        self.out_dir = out_dir
        self.workers = workers
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.managers = [RegressionManager(dat, self.workers, self.executor) for dat in datasets]
        for rem in self.managers:
            rem.dataset.workers = self.workers
            rem.dataset.executor = self.executor


    def set_load_peptides_options(self, *args, **kwargs):
        for rem in self.managers:
            rem.dataset.load_peptides_options(*args, **kwargs)


    def run(self):
        tex = ThreadPoolExecutor(max_workers=1)
        future = tex.submit(lambda: defaultdict(lambda: defaultdict(list)))
        for rem in self.managers:
            print(rem.name)
            self._run(rem)
            self.results = future.result()
            future = tex.submit(self._save, rem, self.results, self.out_dir)
        self.results = future.result()
        tex.shutdown()


    def join(self):
        results = self.results
        for subname in results:
            df_results = pd.concat(results[subname]['df_results'], axis=1)
            df_results_max = pd.concat(results[subname]['df_results_max'], axis=1)
            df_coef = pd.concat(results[subname]['df_coef'], axis=0)
            results[subname] = Results._make([df_results, df_results_max, df_coef, pd.DataFrame(), {}])
        self.results = dict(results)
        self.write(self.results, self.out_dir)


    @staticmethod
    def _run(rem):
        #rem.dataset.pepfiles = [x for x in rem.dataset.pepfiles if 'W9' in x]
        rem.dataset.load_peptides()
        rem.dataset.construct_datasets()
        rem.dataset.clear_unused()
        rem.initialize_trainers()
        rem.start(optimize_mem_usage=True)
        rem.join()


    @staticmethod
    def _save(rem, results, out_sub_dir):
        for subname in list(rem.results.keys()):
            res = rem.results[subname]
            out_dir = os.path.join(out_sub_dir, subname)
            out_bak_dir = os.path.join(out_sub_dir, subname, '_bak')
            out_scores_dir = os.path.join(out_sub_dir, subname, 'scores')
            out_auc_dir = os.path.join(out_sub_dir, subname, 'auc')

            os.makedirs(out_bak_dir, exist_ok=True)
            os.makedirs(out_scores_dir, exist_ok=True)
            os.makedirs(out_auc_dir, exist_ok=True)

            res.df_results.to_csv(os.path.join(out_bak_dir, rem.name + '.results.tsv'),
                    sep='\t', float_format='%g')
            res.df_results_max.to_csv(os.path.join(out_bak_dir, rem.name + '.results.max.tsv'),
                    sep='\t', float_format='%g')
            res.df_coef.to_csv(os.path.join(out_bak_dir, rem.name + '.coefficients.tsv'),
                    sep='\t', float_format='%g')
            res.df_scores.to_csv(os.path.join(out_scores_dir, rem.name + '.scores.tsv.gz'),
                    sep='\t', float_format='%g', compression='gzip')
            with gzip.open(os.path.join(out_auc_dir, rem.name + '.auc_data.json.gz'), 'wb') as f:
                f.write(str.encode(json.dumps(res.auc_data) + '\n'))

            results[subname]['df_results'].append(res.df_results)
            results[subname]['df_results_max'].append(res.df_results_max)
            results[subname]['df_coef'].append(res.df_coef)
            del res
            del rem.results[subname]
            gc.collect()
        return results


    @staticmethod
    def write(results, out_base_dir):
        for subname, res in results.items():
            out_dir = os.path.join(out_base_dir, subname)
            os.makedirs(out_dir, exist_ok=True)

            out_file_results = os.path.join(out_dir, 'metrics.tsv')
            out_file_results_max = os.path.join(out_dir, 'metrics.redmax.tsv')
            out_file_coefficients = os.path.join(out_dir, 'coefficients.tsv')

            print('Saving %s at %s' % (out_file_results, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            print('Saving %s at %s' % (out_file_results_max, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            print('Saving %s at %s' % (out_file_coefficients, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

            res.df_results.to_csv(out_file_results, sep='\t', float_format='%g', na_rep='nan')
            res.df_results_max.to_csv(out_file_results_max, sep='\t', float_format='%g', na_rep='nan')
            res.df_coef.to_csv(out_file_coefficients, sep='\t', float_format='%g', na_rep='nan')


class RegressionManager(object):
    """
    """
    def __init__(self, dataset, workers=0, executor=None):
        self.dataset = dataset
        self.workers = workers
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.allele = list(self.dataset.alleles)[0] if len(self.dataset.alleles) == 1 else 'HLA-MinScore'
        self.name = self.dataset.dataset_name


    def initialize_trainers(self):
        self.trainers = []
        subgroups = defaultdict(list)
        for i, (dat, (seed, subname)) in enumerate(zip(self.dataset.datasets, self.dataset.metadata)):
            self.trainers.append(RegressionTrainer(dat, seed))
            subgroups[subname].append(i)
        self.subgroups = dict(subgroups)


    def start(self, optimize_mem_usage=False):
        with self.executor(max_workers=self.workers, use_threads=False) as ex:
            for t in self.trainers:
                t.submit(ex)
            print(ex)

        if optimize_mem_usage:
            for t in self.trainers:
                del t.dataset
            del self.dataset.datasets
            gc.collect()


    def join(self):
        for t in self.trainers:
            t.join()

        self.results = {}

        for subname, indexes in self.subgroups.items():
            trainers = [self.trainers[i] for i in indexes]
            metadata = [self.dataset.metadata[i] for i in indexes]

        self.results[subname] = Results._make(self._join(trainers, metadata))


    def _join(self, trainers, metadata):
        results_dct = defaultdict(list)
        results_max_dct = defaultdict(list)
        coef_dct = defaultdict(list)
        auc_data = defaultdict(list)
        df_scores_list = []

        for t, (replicate, subname) in zip(trainers, metadata):
            for k in t.results_dct:
                results_dct[k].append(t.results_dct[k])
            for k in t.results_max_dct:
                results_max_dct[k].append(t.results_max_dct[k])
            for k in t.coef_dct:
                coef_dct[k].append(t.coef_dct[k])
            for k in t.auc_data:
                auc_data[k].append(t.auc_data[k])
            df_sc = t.df_scores.copy().replace(np.nan, None)
            df_sc.index = pd.Index([(self.name, self.allele, replicate, x[0], x[1])
                for x in df_sc.index], name= ['Dataset', 'Allele', 'N'] + df_sc.index.names)
            df_scores_list.append(df_sc)

        df_results = pd.DataFrame(results_dct).sort_index(axis=1)
        df_results.index = df_results.index + 1
        df_results.columns = pd.Index([(self.allele, x[0], x[1]) for x in df_results.columns])
        df_results.index = pd.Index([(self.name, x) for x in df_results.index])
        df_results = df_results.transpose()
        df_results.index.names = ['Allele', 'Metric', 'Regression']

        df_results_max = pd.DataFrame(results_max_dct).sort_index(axis=1)
        df_results_max.index = df_results_max.index + 1
        df_results_max.columns = pd.Index([(self.allele, x[0], x[1]) for x in df_results_max.columns])
        df_results_max.index = pd.Index([(self.name, x) for x in df_results_max.index])
        df_results_max = df_results_max.transpose()
        df_results_max.index.names = ['Allele', 'Metric', 'Regression']

        auc_data = {self.name + ':::' + self.allele: dict(auc_data)}

        dct = {}
        for prefix, val_lst in coef_dct.items():
            for rep, values in enumerate(val_lst):
                rep += 1
                minidct = {var: val for var, val in zip(values['variable'], values['coefficient'])}
                dct[(self.name, self.allele, prefix, rep)] = minidct
        df_coef = pd.DataFrame(dct).transpose()
        df_coef.index.names = ['Dataset', 'Allele', 'Regression', 'N']

        df_scores = pd.concat(df_scores_list, copy=False, axis=0).sort_index(axis=1)

        return df_results, df_results_max, df_coef, df_scores, auc_data


class RegressionTrainer(object):
    """
    """
    def __init__(self, dataset, seed, workers=0, executor=None):
        """ NOTE: Cannot fork if using MPIPoolExecutor (must use workers=0)
        """
        self.workers = workers
        # careful with parallelization as LogisticRegression spawns its own processes with all available CPUs
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.dataset = dataset
        self.seed = seed

        self.all_combinations = [
            ["BS"],
            ["TPM"],
            ["CAMAPShuff"],
            ["CAMAP"],
            ["BS", "TPM"],
            ["BS", "CAMAPShuff"],
            ["BS", "CAMAP"],
            ["TPM", "CAMAPShuff"],
            ["TPM", "CAMAP"],
            ["BS", "TPM", "CAMAPShuff"],
            ["BS", "TPM", "CAMAP"]
        ]


    def submit(self, ex=None):
        # expects an active executor
        ex = EnhancedProcessPoolExecutor(max_workers=0) if ex is None else ex  # use sequential pseudo-threads

        df_training, df_test = self.dataset

        self.regressions = []
        #w = self.workers + 1
        #sub_x_labels_groups = [self.all_combinations[n:n+w] for n in range(0, len(self.all_combinations), w)]
        #for sub_x_labels in sub_x_labels_groups:
        #    p = ex.submit(self._submit, sub_x_labels, df_training, df_test, self.seed, self.executor)
        #    self.regressions.append(p)
        for x_labels in self.all_combinations:
            p = ex.submit(self.train_regression, x_labels, df_training, df_test, self.seed)
            self.regressions.append(p)


    #@classmethod
    #def _submit(cls, sub_x_labels, df_training, df_test, seed, executor):
    #    workers = len(sub_x_labels) - 1
    #    sub_regressions = []
    #    with executor(max_workers=workers, use_threads=True) as sex:  # secondary ex
    #        for x_labels in sub_x_labels:
    #            p = sex.submit(cls.train_regression, x_labels, df_training, df_test, seed)
    #            sub_regressions.append(p)
    #        print(sex)
    #    return [p.result() for p in sub_regressions]


    def join(self):
        #regressions = [res for p in self.regressions for res in p.result()]
        regressions = [p.result() for p in self.regressions]

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
        df_scores = df_scores[[x for x in df_scores.columns if x == prefix or x == prefix + 'Base']].copy()
        df_scores = df_scores.rename({prefix: 'NormScores', prefix + 'Base': 'BaseScores'}, axis=1)
        df_scores.columns = pd.MultiIndex.from_product([[prefix], df_scores.columns])

        auc_dct = {'FPR': [], 'TPR': [], 'THRESHOLDS': [], 'FPR_max': [], 'TPR_max': [], 'THRESHOLDS_max': []}
        auc_dct = {prefix: auc_dct}

        coef_dct = {prefix: {'variable': x_labels, 'coefficient': [np.nan]*len(x_labels)}}

        if df_training.Class.sum() < 10 or (~df_training.Class).sum() < 10 or \
                df_test.Class.sum() < 10 or (~df_test.Class).sum() < 10:
            print('Less than 10 members in one class. Skipping...')

            df_scores[(prefix, 'Scores')] = np.nan
            df_scores[(prefix, 'Predictions')] = False
            df_scores[(prefix, 'ScoresReduced')] = np.nan
            df_scores[(prefix, 'PredictionsReduced')] = False

        else:
            # *Train*
            X, y = df_training[x_labels], df_training.Class
            clf = LogisticRegressionCV(cv=10, random_state=seed, solver='lbfgs',
                    multi_class='auto', class_weight='balanced', max_iter=1000)
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
            df_max = df_max.rename({'Scores': 'ScoresReduced', 'Predictions': 'PredictionsReduced'}, axis=1)
            df_max.columns = pd.MultiIndex.from_product([[prefix], df_max.columns])
            df_scores = df_scores.merge(df_max, on=['Peptide', 'Class'])

            # *Sort test DF columns*
            df_scores = df_scores.sort_index(axis=1)

        return results, results_max, coef_dct, auc_dct, df_scores
