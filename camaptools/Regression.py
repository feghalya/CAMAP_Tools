#!/usr/bin/env python3

import os
from collections import defaultdict, namedtuple
import numpy as np
import pandas as pd
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
    def __init__(self, datasets, classifier, out_dir, workers=0, executor=None):
        self.out_dir = out_dir
        executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.managers = [RegressionManager(dat, classifier, workers, executor) for dat in datasets]

        for rem in self.managers:
            rem.dataset.workers = workers
            rem.dataset.executor = executor


    def set_load_peptides_options(self, *args, **kwargs):
        for rem in self.managers:
            rem.dataset.load_peptides_options(*args, **kwargs)


    def run(self):
        tex = ThreadPoolExecutor(max_workers=1)
        future = tex.submit(lambda: defaultdict(lambda: defaultdict(list)))
        for rem in self.managers:
            print(rem.name)
            self._run(rem)
            future = tex.submit(self._merge, rem, future.result(), self.out_dir)
        self.results = future.result()
        tex.shutdown()


    def join(self):
        for subname, result in self.results.items():
            df_results = pd.concat(result['df_results'], axis=1)
            df_results_max = pd.concat(result['df_results_max'], axis=1)
            df_coef = pd.concat(result['df_coef'], axis=0)
            self.results[subname] = Results._make([df_results, df_results_max, df_coef, pd.DataFrame(), {}])
        self.results = dict(self.results)
        self.write(self.results, self.out_dir)


    @staticmethod
    def _run(rem):
        #rem.dataset.pepfiles = [x for x in rem.dataset.pepfiles if 'W9' in x]
        #rem.dataset.pepfiles = [x for x in rem.dataset.pepfiles if x.find('peptides_W') > -1]
        rem.dataset.load_peptides()
        rem.dataset.construct_datasets()
        rem.dataset.clear_unused()
        rem.initialize_trainers()
        rem.start()
        rem.join()


    @classmethod
    def _merge(cls, rem, results, out_sub_dir):
        for trainer_group in list(rem.results.keys()):
            res = rem.results[trainer_group]

            results[trainer_group]['df_results'].append(res.df_results)
            results[trainer_group]['df_results_max'].append(res.df_results_max)
            results[trainer_group]['df_coef'].append(res.df_coef)

            out_dir = os.path.join(out_sub_dir, trainer_group)

            cls._dump(res, out_dir, rem.name)

            # remove pointers to RegressionManager results
            del rem.results[trainer_group]

        return results


    @staticmethod
    def _dump(res, out_dir, subname):
        out_bak_dir = os.path.join(out_dir, '_bak')
        out_scores_dir = os.path.join(out_dir, 'scores')
        out_auc_dir = os.path.join(out_dir, 'auc')

        os.makedirs(out_bak_dir, exist_ok=True)
        os.makedirs(out_scores_dir, exist_ok=True)
        os.makedirs(out_auc_dir, exist_ok=True)

        out_bak = os.path.join(out_bak_dir, subname)
        out_scores = os.path.join(out_scores_dir, subname)
        out_auc = os.path.join(out_auc_dir, subname)

        def save(df, out_file, **kwargs):
            df.to_csv(out_file, sep='\t', float_format='%g', **kwargs)

        save(res.df_results, out_bak + '.results.tsv')
        save(res.df_results_max, out_bak + '.results.max.tsv')
        save(res.df_coef, out_bak + '.coefficients.tsv')
        save(res.df_scores, out_scores + '.scores.tsv.gz', compression='gzip')

        with gzip.open(out_auc + '.auc_data.json.gz', 'wb') as f:
            f.write(str.encode(json.dumps(res.auc_data) + '\n'))


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

            def save(df, out_file):
                df.to_csv(out_file, sep='\t', float_format='%g', na_rep='nan')

            save(res.df_results, out_file_results)
            save(res.df_results_max, out_file_results_max)
            save(res.df_coef, out_file_coefficients)


class RegressionManager(object):
    """
    """
    def __init__(self, dataset, classifier, workers=0, executor=None):
        self.dataset = dataset
        self.classifier = classifier
        self.workers = workers
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.allele = list(self.dataset.alleles)[0] if len(self.dataset.alleles) == 1 else 'HLA-MinScore'
        self.name = self.dataset.dataset_name


    def initialize_trainers(self):
        self.trainers = {}
        for dat, (i, subname) in zip(self.dataset.datasets, self.dataset.metadata):
            trainer = RegressionTrainer(dat, self.classifier, i)
            self.trainers[trainer] = self.classifier + '-' + subname

        self.optimize_mem_usage = False  # Can manually change to True


    def start(self):
        with self.executor(max_workers=self.workers, use_threads=False) as ex:
            for t in self.trainers:
                t.submit(ex)
            print(ex)

        if self.optimize_mem_usage:
            for t in self.trainers:
                del t.dataset
            del self.dataset.datasets
            gc.collect()


    def join(self):
        for trainer in self.trainers:
            trainer.join()

        results = {x: [] for x in set(self.trainers.values())}
        for trainer, trainer_group in self.trainers.items():
            results[trainer_group].append(trainer)

        self.results = {g: self._merge(ts) for g, ts in results.items()}


    def _merge(self, trainers):
        results_dct = defaultdict(list)
        results_max_dct = defaultdict(list)
        coef_dct = defaultdict(list)
        auc_data = defaultdict(list)
        df_scores_list = []

        for trainer in trainers:
            for k in trainer.results_dct:
                results_dct[k].append(trainer.results_dct[k])
            for k in trainer.results_max_dct:
                results_max_dct[k].append(trainer.results_max_dct[k])
            for k in trainer.coef_dct:
                coef_dct[k].append(trainer.coef_dct[k])
            for k in trainer.auc_data:
                auc_data[k].append(trainer.auc_data[k])
            df_sc = trainer.df_scores.copy().replace(np.nan, None)
            re_index = [(self.name, self.allele, trainer.seed, x[0], x[1]) for x in df_sc.index]
            re_index_names = ['Dataset', 'Allele', 'N'] + df_sc.index.names
            df_sc.index = pd.Index(re_index, name=re_index_names)
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

        results = (df_results, df_results_max, df_coef, df_scores, auc_data)

        return Results._make(results)


class RegressionTrainer(object):
    """
    """
    def __init__(self, dataset, classifier, seed): #, workers=0, executor=None):
        """ NOTE: Cannot fork if using MPIPoolExecutor (must use workers=0)
        """
        #self.workers = workers
        # careful with parallelization as LogisticRegression spawns its own processes with all available CPUs
        #self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.dataset = dataset
        self.classifier = classifier
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

        df_train, df_test = self.dataset

        self.regressions = []
        #w = self.workers + 1
        #sub_x_labels_groups = [self.all_combinations[n:n+w] for n in range(0, len(self.all_combinations), w)]
        #for sub_x_labels in sub_x_labels_groups:
        #    p = ex.submit(self._submit, sub_x_labels, df_train, df_test, self.seed, self.executor)
        #    self.regressions.append(p)
        for x_labels in self.all_combinations:
            p = ex.submit(self.train_regression, x_labels, df_train, df_test, self.classifier, self.seed)
            self.regressions.append(p)


    #@classmethod
    #def _submit(cls, sub_x_labels, df_train, df_test, seed, executor):
    #    workers = len(sub_x_labels) - 1
    #    sub_regressions = []
    #    with executor(max_workers=workers, use_threads=True) as sex:  # secondary ex
    #        for x_labels in sub_x_labels:
    #            p = sex.submit(cls.train_regression, x_labels, df_train, df_test, seed)
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
    def train_regression(x_labels, df_train, df_test, classifier, seed):
        # Which labels are used
        prefix = "-".join(x_labels)

        # Prepare return data with null values
        results = {
            ('Accuracy', prefix): np.nan,
            ('PPV', prefix): np.nan,
            ('Sensitivity', prefix): np.nan,
            ('Specificity', prefix): np.nan,
            ('AUC', prefix) : np.nan,
            ('MCC', prefix): np.nan
        }

        results_max = {
            ('Accuracy', prefix): np.nan,
            ('PPV', prefix): np.nan,
            ('Sensitivity', prefix): np.nan,
            ('Specificity', prefix): np.nan,
            ('AUC', prefix) : np.nan,
            ('MCC', prefix): np.nan
        }

        df_scores = df_test.copy()
        df_scores = df_scores.reset_index().set_index(['Peptide', 'Class'])
        df_scores = df_scores[[x for x in df_scores.columns if x == prefix or x == prefix + 'Base']].copy()
        df_scores = df_scores.rename({prefix: 'NormScores', prefix + 'Base': 'BaseScores'}, axis=1)

        auc_dct = {'FPR': [], 'TPR': [], 'THRESHOLDS': [], 'FPR_max': [], 'TPR_max': [], 'THRESHOLDS_max': []}
        auc_dct = {prefix: auc_dct}

        coef_dct = {prefix: {'variable': x_labels, 'coefficient': [np.nan]*len(x_labels)}}

        def train(x_labels, df_train, df_test, classifier, seed):
            # *Train*

            X, y = df_train[x_labels], df_train.Class

            if classifier == 'logregcv':
                from sklearn.linear_model import LogisticRegressionCV
                clf = LogisticRegressionCV(
                    cv=10,
                    random_state=seed,
                    solver='lbfgs',
                    multi_class='auto',
                    class_weight='balanced',
                    max_iter=1000,
                    n_jobs=1
                )
            elif classifier == 'logreg':
                from sklearn.linear_model import LogisticRegression
                clf = LogisticRegression(
                    random_state=seed,
                    multi_class='auto',
                    class_weight='balanced'
                )
            elif classifier == 'svc':
                from sklearn.svm import SVC
                # SVC is very slow
                clf = SVC(
                    kernel='rbf',
                    gamma='auto',
                    class_weight='balanced',
                    random_state=seed,
                    probability=True
                )
            elif classifier == 'linsvc':
                from sklearn.svm import SVC
                # still slow, even in linear mode
                # note that LinearSVC does not have an option to predict probabilities
                clf = SVC(
                    kernel='linear',
                    gamma='auto',
                    class_weight='balanced',
                    random_state=seed,
                    probability=True
                )
            elif classifier == 'rf':
                from sklearn.ensemble import RandomForestClassifier
                clf = RandomForestClassifier(
                    n_estimators=100,
                    random_state=seed,
                    class_weight='balanced',
                    n_jobs=1
                )
            elif classifier == 'knn':
                from sklearn.neighbors import KNeighborsClassifier
                clf = KNeighborsClassifier(
                    n_neighbors=101,
                    weights='distance',
                    n_jobs=1
                )
            elif classifier == 'mlp':
                from sklearn.neural_network import MLPClassifier
                clf = MLPClassifier(
                    max_iter=4000,
                    hidden_layer_sizes=(100,),
                    solver='sgd',
                    early_stopping=True
                )

            clf.fit(X, y)  # will not produce ValueError because we take care of cases where n < 10

            try:
                coefficients = clf.coef_[0]
                s = ": Regression coef %s" % coefficients
            except AttributeError:
                s = ": Coefficients not available"
                coefficients = [np.nan]*len(x_labels)

            print(prefix, s)

            # *Test*

            X, y = df_test[x_labels], df_test.Class

            weights = abs(y - (y.sum() / y.shape[0]))

            labels = y.to_numpy()
            scores = clf.predict_proba(X)[:,1]
            predictions = clf.predict(X)
            accuracy = clf.score(X, y, weights)
            calc_accuracy = (((scores >= 0.5) == labels) * weights).sum() / weights.sum()

            if np.array_equal(scores >= 0.5, predictions):
                assert np.isclose(accuracy, calc_accuracy)
            else:
                print('Scores and predictions differ, normal if not using LogReg')
                #print('Will be using calculcated accuracy for consistency')
                print(
                    'Predicated accuracy :',
                    round(accuracy, 5),
                    ' ; Calculated accuracy :',
                    round(calc_accuracy, 5)
                )

            return coefficients, labels, scores, predictions, accuracy, weights

        def calc_metrics(labels, scores, predictions):
            weights = abs(labels - (labels.sum() / labels.shape[0]))

            accuracy = ((predictions == labels) * weights).sum() / weights.sum()
            ppv = metrics.precision_score(
                labels, predictions, sample_weight=weights,
                labels=None, pos_label=1, average='binary'
            )
            npv = metrics.precision_score(
                labels, 1 - predictions, sample_weight=weights,
                labels=None, pos_label=1, average='binary'
            )
            sensitivity = metrics.recall_score(
                1 - labels, predictions, sample_weight=weights,
                labels=None, pos_label=1, average='binary'
            )
            specificity = metrics.recall_score(
                1 - labels, 1 - predictions, sample_weight=weights,
                labels=None, pos_label=1, average='binary'
            )
            mcc = metrics.matthews_corrcoef(labels, predictions, sample_weight=weights)
            fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
            auc = metrics.auc(fpr, tpr)
            auc_data = {'FPR': list(fpr), 'TPR': list(tpr), 'THRESHOLDS': list(thresholds)}

            # NOTE on MCC
            #factor = len(labels[~labels])//len(labels[labels])
            #predictions_w = np.concatenate([np.repeat(predictions[labels], factor), predictions[~labels]])
            #labels_w = np.concatenate([np.repeat(labels[labels], factor), labels[~labels]])
            # ==> Calculating MCC with weights is identical to weighting predictions and labels beforehand

            return accuracy, ppv, npv, sensitivity, specificity, mcc, auc, auc_data

        if df_train.Class.sum() < 10 or (~df_train.Class).sum() < 10 or \
                df_test.Class.sum() < 10 or (~df_test.Class).sum() < 10:
            print('Less than 10 members in one class. Skipping...')

            if df_scores.empty:
                columns = df_scores.columns.to_list() + ['Scores', 'Predictions', 'ScoresReduced', 'PredictionsReduced']
                df_scores = pd.DataFrame({'Peptide': ['X'], 'Class': [False]}).set_index(['Peptide', 'Class'])
                for col in columns:
                    df_scores[col] = np.nan

            df_scores['Scores'] = np.nan
            df_scores['Predictions'] = False
            df_scores['ScoresReduced'] = np.nan
            df_scores['PredictionsReduced'] = False

        else:
            coefficients, labels, scores, predictions, accuracy, weights = train(x_labels, df_train, df_test, classifier, seed)

            coef_dct = {prefix: {'variable': x_labels, 'coefficient': list(coefficients)}}

            # *Append to test DF*
            df_scores['Scores'] = scores
            df_scores['Predictions'] = predictions

            # *Calculate metrics*

            # 1) Conventional method

            #tp = (labels & predictions).sum()
            #tn = (~labels & ~predictions).sum()
            #fp = (~labels & predictions).sum()
            #fn = (labels & ~predictions).sum()

            #all_pos = labels.sum()
            #all_neg = (~labels).sum()
            #assert all_pos == tp + fn
            #assert all_neg == tn + fp

            accuracy, ppv, npv, sensitivity, specificity, mcc, auc, auc_data = calc_metrics(labels, scores, predictions)

            results[('Accuracy', prefix)] = accuracy
            results[('PPV', prefix)] = ppv
            results[('NPV', prefix)] = npv
            results[('Sensitivity', prefix)] = sensitivity
            results[('Specificity', prefix)] = specificity
            results[('AUC', prefix)] = auc
            results[('MCC', prefix)] = mcc

            auc_dct = auc_data

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

            #tp = (labels & predictions).sum()
            #tn = (~labels & ~predictions).sum()
            #fp = (~labels & predictions).sum()
            #fn = (labels & ~predictions).sum()

            accuracy, ppv, npv, sensitivity, specificity, mcc, auc, auc_data = calc_metrics(labels, scores, predictions)

            results_max[('Accuracy', prefix)] = accuracy
            results_max[('PPV', prefix)] = ppv
            results_max[('NPV', prefix)] = npv
            results_max[('Sensitivity', prefix)] = sensitivity
            results_max[('Specificity', prefix)] = specificity
            results_max[('AUC', prefix)] = auc
            results_max[('MCC', prefix)] = mcc

            auc_dct_max =  {k+'_max': v for k, v in auc_data.items()}

            # *Merge AUC values*
            auc_dct = {prefix: {**auc_dct, **auc_dct_max}}

            # *Append to test DF*
            df_max = df_max.reset_index().set_index(['Peptide', 'Class'])
            df_max = df_max.rename({'Scores': 'ScoresReduced', 'Predictions': 'PredictionsReduced'}, axis=1)
            df_scores = df_scores.merge(df_max, on=['Peptide', 'Class'])

        # *Sort test DF columns*
        df_scores = df_scores.sort_index(axis=1)

        # *Add prefix as level 0 column index*
        df_scores.columns = pd.MultiIndex.from_product([[prefix], df_scores.columns])

        return results, results_max, coef_dct, auc_dct, df_scores
