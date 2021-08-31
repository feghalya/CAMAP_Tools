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


#Results = namedtuple('Results', ['df_results', 'df_results_max', 'df_coef', 'df_scores', 'auc_data', 'auc_data_max'])
Results = namedtuple('Results', ['df_results', 'df_ppv', 'df_coef', 'df_scores', 'auc_data'])


class RegressionAggregator(object):
    """A wrapper around RegressionManagers that dumps results to disk
    """
    def __init__(self, regressors, out_dir, workers=0, executor=None):
        self.out_dir = out_dir
        executor = EnhancedProcessPoolExecutor if executor is None else executor
        # expects regressors to have their trainers initialized
        self.managers = regressors


    def run(self):
        tex = ThreadPoolExecutor(max_workers=1)
        known_keys = ['df_results', 'df_ppv', 'df_coef']
        future = tex.submit(lambda: defaultdict(lambda: {k: [] for k in known_keys}))
        for rem in self.managers:
            print(rem.name)
            self._run(rem)
            future = tex.submit(self._merge, rem, future.result(), self.out_dir)
        self.results = future.result()
        tex.shutdown()


    def join(self):
        for subname, result in self.results.items():
            df_results = pd.concat(result['df_results'], axis=1)
            #df_results_max = pd.concat(result['df_results_max'], axis=1)
            df_ppv = pd.concat(result['df_ppv'], axis=1)
            df_coef = pd.concat(result['df_coef'], axis=0)
            #self.results[subname] = Results._make([df_results, df_results_max, df_coef, pd.DataFrame(), {}, {}])
            self.results[subname] = Results._make([df_results, df_ppv, df_coef, pd.DataFrame(), {}])
        self.results = dict(self.results)
        self.write(self.results, self.out_dir)


    @staticmethod
    def _run(rem):
        rem.start()
        rem.join()


    @classmethod
    def _merge(cls, rem, results, out_sub_dir):
        for trainer_group in list(rem.results.keys()):
            res = rem.results[trainer_group]

            results[trainer_group]['df_results'].append(res.df_results)
            #results[trainer_group]['df_results_max'].append(res.df_results_max)
            results[trainer_group]['df_ppv'].append(res.df_ppv)
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
        save(res.df_ppv, out_bak + '.ppv.tsv')
        #save(res.df_results_max, out_bak + '.results.max.tsv')
        save(res.df_coef, out_bak + '.coefficients.tsv')
        save(res.df_scores, out_scores + '.scores.tsv.gz', compression='gzip')

        with gzip.open(out_auc + '.auc_data.json.gz', 'wb') as f:
            f.write(str.encode(json.dumps(res.auc_data) + '\n'))

        #with gzip.open(out_auc + '.auc_data.redmax.json.gz', 'wb') as f:
        #    f.write(str.encode(json.dumps(res.auc_data_max) + '\n'))


    @staticmethod
    def write(results, out_base_dir):
        for subname, res in results.items():
            out_dir = os.path.join(out_base_dir, subname)
            os.makedirs(out_dir, exist_ok=True)

            out_file_results = os.path.join(out_dir, 'metrics.tsv')
            out_file_ppv = os.path.join(out_dir, 'ppv.metrics.tsv')
            #out_file_results_max = os.path.join(out_dir, 'metrics.redmax.tsv')
            out_file_coefficients = os.path.join(out_dir, 'coefficients.tsv')

            print('Saving %s at %s' % (out_file_results, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            #print('Saving %s at %s' % (out_file_results_max, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            print('Saving %s at %s' % (out_file_ppv, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            print('Saving %s at %s' % (out_file_coefficients, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

            def save(df, out_file):
                df.to_csv(out_file, sep='\t', float_format='%g', na_rep='nan')

            save(res.df_results, out_file_results)
            #save(res.df_results_max, out_file_results_max)
            save(res.df_ppv, out_file_ppv)
            save(res.df_coef, out_file_coefficients)


class RegressionManager(object):
    """
    """
    def __init__(self, dataset, workers=0, executor=None):
        self.dataset = dataset
        self.workers = workers
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.allele = list(self.dataset.alleles)[0] if len(self.dataset.alleles) == 1 else 'HLA-MinScore'
        self.name = self.dataset.dataset_name

        self.optimize_mem_usage = False  # can be modified manually


    def initialize_trainers(self, classifiers):
        self.trainers = {}
        for dat, (i, subname) in zip(self.dataset.datasets, self.dataset.metadata):
            for classifier in classifiers:
                trainer = RegressionTrainer(dat, classifier, i)
                self.trainers[trainer] = classifier + '-' + subname

        if self.optimize_mem_usage:
            del self.dataset.datasets
            gc.collect()


    def dump_trainers(self, path_on_disk):
        for trainer, subdir in self.trainers.items():
            out_dir = os.path.join(path_on_disk, self.name, subdir)
            out_file = os.path.join(out_dir, str(trainer.seed) +'.p')
            os.makedirs(out_dir, exist_ok=True)
            with open(out_file, 'wb') as f:
                pkl.dump(trainer, f)


    def load_trainers(self, path_on_disk, classifiers):
        self.trainers = {}
        for subdir in os.listdir(os.path.join(path_on_disk, self.name)):
            if subdir.split('-')[0] in classifiers:
                for t_file in os.listdir(os.path.join(path_on_disk, self.name, subdir)):
                    assert t_file.endswith('.p')
                    t_path = os.path.join(path_on_disk, self.name, subdir, t_file)
                    if self.optimize_mem_usage:
                        trainer = RegressionTrainerLowMemory(t_path)
                    else:
                         with open(t_path, 'rb') as f:
                            trainer = pkl.load(f)
                    self.trainers[trainer] = subdir


    def start(self):
        with self.executor(max_workers=self.workers, use_threads=False) as ex:
            for t in self.trainers:
                t.submit(ex)
            print(ex)


    def join(self):
        for trainer in self.trainers:
            trainer.join()

        results = {x: [] for x in set(self.trainers.values())}
        for trainer, trainer_group in self.trainers.items():
            results[trainer_group].append(trainer)

        self.results = {g: self._merge(ts) for g, ts in results.items()}


    def _merge(self, trainers):
        results_dct = defaultdict(list)
        #results_max_dct = defaultdict(list)
        ppv_dct = defaultdict(list)
        coef_dct = defaultdict(list)
        auc_data = defaultdict(list)
        #auc_data_max = defaultdict(list)
        df_scores_list = []

        for trainer in trainers:
            for k in trainer.results_dct:
                results_dct[k].append(trainer.results_dct[k])
            #for k in trainer.results_max_dct:
            #    results_max_dct[k].append(trainer.results_max_dct[k])
            for k in trainer.ppv_dct:
                ppv_dct[k].append(trainer.ppv_dct[k])
            for k in trainer.coef_dct:
                coef_dct[k].append(trainer.coef_dct[k])
            for k in trainer.auc_data:
                auc_data[k].append(trainer.auc_data[k])
            #for k in trainer.auc_data_max:
            #    auc_data_max[k].append(trainer.auc_data_max[k])
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

        #df_results_max = pd.DataFrame(results_max_dct).sort_index(axis=1)
        #df_results_max.index = df_results_max.index + 1
        #df_results_max.columns = pd.Index([(self.allele, x[0], x[1]) for x in df_results_max.columns])
        #df_results_max.index = pd.Index([(self.name, x) for x in df_results_max.index])
        #df_results_max = df_results_max.transpose()
        #df_results_max.index.names = ['Allele', 'Metric', 'Regression']

        df_ppv = pd.DataFrame(ppv_dct).sort_index(axis=1)
        df_ppv.index = df_ppv.index + 1
        df_ppv.columns = pd.Index([(self.allele, x[0], x[1], x[2]) for x in df_ppv.columns])
        df_ppv.index = pd.Index([(self.name, x) for x in df_ppv.index])
        df_ppv = df_ppv.transpose()
        df_ppv.index.names = ['Allele', 'Metric', 'Statistic', 'Regression']

        auc_data = {self.name + ':::' + self.allele: dict(auc_data)}

        #auc_data_max = {self.name + ':::' + self.allele: dict(auc_data_max)}

        dct = {}
        for prefix, val_lst in coef_dct.items():
            for rep, values in enumerate(val_lst):
                rep += 1
                minidct = {var: val for var, val in zip(values['variable'], values['coefficient'])}
                dct[(self.name, self.allele, prefix, rep)] = minidct
        df_coef = pd.DataFrame(dct).transpose()
        df_coef.index.names = ['Dataset', 'Allele', 'Regression', 'N']

        df_scores = pd.concat(df_scores_list, copy=False, axis=0).sort_index(axis=1)

        #results = (df_results, df_results_max, df_coef, df_scores, auc_data, auc_data_max)
        results = (df_results, df_ppv, df_coef, df_scores, auc_data)

        return Results._make(results)


class RegressionTrainer(object):
    """
    """

    def __init__(self, dataset, classifier, seed):
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
        # `ex` expects an active executor
        # Cannot fork if using MPIPoolExecutor (must use workers=0)

        ex = EnhancedProcessPoolExecutor(max_workers=0) if ex is None else ex  # use sequential pseudo-threads

        df_train, df_test = self.dataset

        self.regressions = []
        for x_labels in self.all_combinations:
            p = ex.submit(self.benchmark, x_labels, df_train, df_test, self.classifier, self.seed)
            self.regressions.append(p)


    def join(self):
        regressions = [p.result() for p in self.regressions]

        results = {}
        #results_max = {}
        ppv_dct = {}
        coef_dct = {}
        auc_data = {}
        #auc_data_max = {}
        df_scores = pd.DataFrame()

        #for r, rm, coef, auc, aucm, df_sc in regressions:
        for r, ppv, coef, auc, df_sc in regressions:
            results = {**results, **r}
            #results_max = {**results_max, **rm}
            ppv_dct = {**ppv_dct, **ppv}
            coef_dct = {**coef_dct, **coef}
            auc_data = {**auc_data, **auc}
            #auc_data_max = {**auc_data_max, **aucm}
            df_scores = pd.concat([df_scores, df_sc], axis=1)

        self.results_dct = results
        #self.results_max_dct = results_max
        self.ppv_dct = ppv_dct
        self.coef_dct = coef_dct
        self.auc_data = auc_data
        #self.auc_data_max = auc_data_max
        self.df_scores = df_scores


    @classmethod
    def benchmark(cls, x_labels, df_train, df_test, classifier, seed):
        prefix = "-".join(x_labels)

        df_scores = df_test.copy()
        df_scores = df_scores.reset_index().set_index(['Peptide', 'Class'])
        df_scores = df_scores[[x for x in df_scores.columns if x == prefix or x == prefix + 'Base']].copy()
        df_scores = df_scores.rename({prefix: 'NormScores', prefix + 'Base': 'BaseScores'}, axis=1)

        tr_cls = df_train.Class
        ts_cls = df_test.Class

        if tr_cls.sum() < 10 or (~tr_cls).sum() < 10 or ts_cls.sum() < 10 or (~ts_cls).sum() < 10:
            print('Less than 10 members in one class. Skipping...')
            final_results = cls._benchmark_null(df_scores, x_labels)
        else:
            clf = cls.get_clf(classifier, seed)
            train_results = cls.train(x_labels, df_train, df_test, clf)
            final_results = cls._benchmark(df_scores, train_results)

        # *Add prefix and return*
        #results, results_max, coef_dct, auc_dct, auc_dct_max, df_scores = final_results
        results, ppv_dct, coef_dct, auc_dct, df_scores = final_results

        results = {(k, prefix): v for k, v in results.items()}
        #results_max = {(k, prefix): v for k, v in results_max.items()}
        ppv_dct = {k + (prefix,): v for k, v in ppv_dct.items()}
        coef_dct['variable'] = x_labels
        coef_dct = {prefix: coef_dct}
        auc_dct = {prefix: auc_dct}
        #auc_dct_max = {prefix: auc_dct_max}
        df_scores.columns = pd.MultiIndex.from_product([[prefix], df_scores.columns])

        #return results, results_max, coef_dct, auc_dct, auc_dct_max, df_scores
        return results, ppv_dct, coef_dct, auc_dct, df_scores


    @classmethod
    def _benchmark(cls, df_scores, train_results):
        coefficients, labels, scores, predictions, accuracy, weights = train_results

        coef_dct = {'coefficient': list(coefficients)}

        # *Append to test DF*
        #df_scores['Scores'] = scores
        #df_scores['Predictions'] = predictions

        # *Calculate metrics*

        # 1) Conventional method
        #    DEPRECATED for 2 reasons
        #       1. RedMax method is biologically more sound and returns better predictions
        #       2. Training is done on a weighted dataset where the lower class count is
        #          copied enough time to mimic the higher class count

        #results, auc_dct = cls.calc_metrics(labels, scores, predictions)
        #ppv_dct = cls.calc_ppv_metrics(labels, scores)

        # 2) Keep only max score per peptide and re-calculate metrics (same strategy as reduce_tables)

        df_max = df_scores.reset_index().set_index('Peptide')[['Class']].copy()
        df_max['Scores'] = scores
        types = df_max.dtypes
        df_max = df_max.groupby('Peptide', sort=False).apply(max)
        df_max = df_max.astype(types)

        labels = df_max.Class.to_numpy()
        scores = df_max.Scores.to_numpy()

        predictions = scores >= 0.5
        df_max['Predictions'] = predictions

        #results_max, auc_dct_max = cls.calc_metrics(labels, scores, predictions)
        #ppv_dct_max = cls.calc_ppv_metrics(labels, scores)
        results, auc_dct = cls.calc_metrics(labels, scores, predictions)
        ppv_dct = cls.calc_ppv_metrics(labels, scores)

        # *Append to test DF*
        df_max = df_max.reset_index().set_index(['Peptide', 'Class'])
        #df_max = df_max.rename(
        #    {'Scores': 'ScoresReduced', 'Predictions': 'PredictionsReduced'},
        #    axis=1
        #)

        df_scores = df_scores.merge(df_max, on=['Peptide', 'Class'])

        # *Sort test DF columns*
        df_scores = df_scores.sort_index(axis=1)

        #return results, results_max, coef_dct, auc_dct, auc_dct_max, df_scores
        return results, ppv_dct, coef_dct, auc_dct, df_scores


    @staticmethod
    def _benchmark_null(df_scores, x_labels):
        results = {
            'Accuracy': np.nan,
            'PPV': np.nan,
            'Sensitivity': np.nan,
            'Specificity': np.nan,
            'AUC' : np.nan,
            'MCC': np.nan
        }

        #results_max = {k: v for k, v in results.items()}

        coef_dct = {'coefficient': [np.nan]*len(x_labels)}

        auc_dct = {'FPR': [], 'TPR': [], 'THRESHOLDS': []}
        #auc_dct_max = {'FPR': [], 'TPR': [], 'THRESHOLDS': []}

        ppv_dct = {}
        ppv_dct[('nEntries', 'nhits')] = 0
        ppv_dct[('nHits', 'nhits')] = 0
        for p in [0.01] + [x/100 for x in range(5, 100, 5)] + [0.99]:
            p_in_str = str(int(p*100))
            stage = 'PPV%s' % ('0'*abs(len(p_in_str) - 2) + p_in_str)
            ppv_dct[(stage, 'nhits')] = 0
            ppv_dct[(stage, 'npredicted')] = 0
            ppv_dct[(stage, 'fpr')] = np.nan

        if df_scores.empty:
            cols = ['Scores', 'Predictions', 'ScoresReduced', 'PredictionsReduced']
            columns = df_scores.columns.to_list() + cols
            dummy_scores = {'Peptide': ['X'], 'Class': [False]}
            df_scores = pd.DataFrame(dummy_scores).set_index(['Peptide', 'Class'])
            for col in columns:
                df_scores[col] = np.nan

        df_scores['Scores'] = np.nan
        df_scores['Predictions'] = False
        df_scores['ScoresReduced'] = np.nan
        df_scores['PredictionsReduced'] = False

        #return results, results_max, coef_dct, auc_dct, auc_dct_max, df_scores
        return results, ppv_dct, coef_dct, auc_dct, df_scores


    @staticmethod
    def get_clf(classifier, seed):
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
                early_stopping=True,
                random_state=seed
            )
        else:
            raise

        return clf


    @staticmethod
    def train(x_labels, df_train, df_test, clf):
        # *Train*
        X, y = df_train[x_labels], df_train.Class

        clf.fit(X, y)  # will not produce ValueError because we take care of cases where n < 10

        try:
            coefficients = clf.coef_[0]
            s = ": Regression coef %s" % coefficients
        except AttributeError:
            s = ": Coefficients not available"
            coefficients = [np.nan]*len(x_labels)
        print("-".join(x_labels), s)

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
            print(
                'Predicated accuracy :',
                round(accuracy, 5),
                ' ; Calculated accuracy :',
                round(calc_accuracy, 5)
            )

        return coefficients, labels, scores, predictions, accuracy, weights


    @staticmethod
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

        # Note on MCC
        #factor = len(labels[~labels])//len(labels[labels])
        #predictions_w = np.concatenate([np.repeat(predictions[labels], factor), predictions[~labels]])
        #labels_w = np.concatenate([np.repeat(labels[labels], factor), labels[~labels]])
        # ==> Calculating MCC with weights is identical to weighting predictions and labels beforehand

        out_dct = {}
        out_dct['Accuracy'] = accuracy
        out_dct['PPV'] = ppv
        out_dct['NPV'] = npv
        out_dct['Sensitivity'] = sensitivity
        out_dct['Specificity'] = specificity
        out_dct['AUC'] = auc
        out_dct['MCC'] = mcc

        auc_dct = {'FPR': list(fpr), 'TPR': list(tpr), 'THRESHOLDS': list(thresholds)}

        return out_dct, auc_dct


    @staticmethod
    def calc_ppv_metrics(labels, scores):
        nentries = len(scores)
        nhits = sum(labels)

        def calc(percentage):
            _, sorted_labels = zip(*sorted(zip(scores, labels), key=lambda x: x[0]))
            sorted_labels = np.cumsum(sorted_labels)

            if nhits == 0:
                n = 0
            elif nhits*percentage < 1:
                n = 1
            else:
                n = round(nhits*percentage)

            try:
                npredicted = sorted_labels.tolist().index(n) + 1
            except ValueError:
                assert nhits == 0
                npredicted = 0

            try:
                fpr = round((npredicted - nhits)/npredicted*100, 3)
            except ZeroDivisionError:
                if nhits == 0:
                    fpr = np.nan
                else:
                    fpr = 100

            return n, npredicted, fpr

        out_dct = {}
        out_dct[('nEntries', 'nhits')] = nentries
        out_dct[('nHits', 'nhits')] = nhits
        for p in [0.01] + [x/100 for x in range(5, 100, 5)] + [0.99]:
            p_in_str = str(int(p*100))
            stage = 'PPV%s' % ('0'*abs(len(p_in_str) - 2) + p_in_str)
            targets, npredicted, fpr = calc(p)
            out_dct[(stage, 'nhits')] = targets
            out_dct[(stage, 'npredicted')] = npredicted
            out_dct[(stage, 'fpr')] = fpr

        return out_dct


class RegressionTrainerLowMemory(RegressionTrainer):
    def __init__(self, trainer_pkl):
        self.trainer_pkl = trainer_pkl

        classifier, seed, all_combinations = self.load_essentials()

        self.classifier = classifier
        self.seed = seed
        self.all_combinations = all_combinations

        self.dataset = self.trainer_pkl


    def load_essentials(self):
        with EnhancedProcessPoolExecutor(max_workers=1, use_threads=False) as ex:
            p = ex.submit(self._load, self.trainer_pkl)
        return p.result()


    @staticmethod
    def _load(trainer_pkl):
        with open(trainer_pkl, 'rb') as f:
            trainer = pkl.load(f)
        return trainer.classifier, trainer.seed, trainer.all_combinations


    @classmethod
    def benchmark(cls, x_labels, trainer_pkl, classifier, seed):
        with open(trainer_pkl, 'rb') as f:
            trainer = pkl.load(f)
        dataset = trainer.dataset
        return super().benchmark(x_labels, dataset, classifier, seed)
