#!/usr/bin/env python3

import os
import pickle as pkl
import pandas as pd
import shutil
import subprocess
from collections import defaultdict
import pkg_resources
from datetime import datetime
import tarfile
import gc

from camaptools.MLP import Model
from camaptools.utils import available_models
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor


OUTPUT_FOLDER = "./output"


class Peptides(object):
    def __init__(self, genome, context=162, workers=0, executor=None):
        self.output_folder = OUTPUT_FOLDER
        self.genome = genome
        assert 'GRCh' in genome or 'GRCm' in genome
        self.out_dir = os.path.join(self.output_folder, 'allPeptides_%s' % self.genome)
        self.species = 'Human' if 'GRCh' in genome else 'Mouse'

        self.context = context

        # detect pickle peptide files
        pepfiles = [os.path.join(self.out_dir, f) for f in os.listdir(self.out_dir) if f[:8] == 'peptides']
        self.pepfiles = pepfiles

        self.workers = workers
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor


    def load_models(self, filters=None):
        def load():
            return available_models(
                target = 'validation-bestMin-score.pytorch',
                context = self.context,
                filters = filters
                )

        self.model_names = {self.context: load()}

        models = {self.context: load()}
        for method, method_dct in models[self.context].items():
            for params, model_list in method_dct.items():
                method_dct[params] = [Model(m, self.context, 'cpu') for m in model_list]
        self.models = models


    def annotate(self, overwrite=False):
        self.overwrite = overwrite

        with self.executor(max_workers=self.workers, use_threads=True) as ex:
            counts = sum(list(ex.map(self._annotate, self.pepfiles)))

        pfile = os.path.join(self.out_dir, 'info.pkl')
        info = pkl.load(open(pfile, 'rb'))

        edited = False
        if 'pyTorchVersion' not in info:
            edited = True
            info['pyTorchVersion'] = pkg_resources.get_distribution("torch").version
        if 'CAMAPModels' not in info:
            edited = True
            info['modelANN'] = self.model_names
        if 'numberPeptideContexts' not in info:
            edited = True
            info['numberPeptideContexts'] = counts

        if edited:
            info['date'] = datetime.now()
            self._keep_a_copy(pfile)
            pkl.dump(info, open(pfile, 'wb'))


    def _annotate(self, pfile, chunk_size=10000):
        """
        Assumes sequences are never modified, will always run the same number of models
        on all sequence chunks.
        Assumes context size is the same for all models.

        """
        print(pfile)
        models = self.models

        def run_models():
            counter_ix, s_list = seqs.keys(), seqs.values()
            c, m, p = list(jobs)[0]
            any_model = models[c][m][p][0]
            s_embeds = any_model._get_embeds(s_list)
            for c, m, p in jobs:
                for i in counter_ix:
                    seq_scores[i][(m, c, p)] = []
                for model in models[c][m][p]:
                    scores = model.get_score(s_embeds)
                    for i, s in zip(counter_ix, scores):
                        seq_scores[i][(m, c, p)].append(float(s[1]))

        pepdict = pkl.load(open(pfile, 'rb'))

        seq_scores = {}
        seq_ids = {}
        seqs = {}
        jobs = set()
        counter = 0
        for pep, gene_dict in pepdict.items():
            gene_dict = gene_dict['genes']
            for gene_name, entries in gene_dict.items():
                for entry_i, entry in enumerate(entries):
                    seq = entry['sequenceContext']
                    if '!GA' not in seq:
                        counter += 1
                        seq_ids[counter] = (pep, gene_name, entry_i)
                        seq_scores[counter] = {}
                        seqs[counter] = seq
                        if self.overwrite:
                            del entry['CAMAPScore']
                        for context in models:
                            for method in models[context]:
                                for params in models[context][method]:
                                    if 'CAMAPScore' not in entry:
                                        entry['CAMAPScore'] = {}
                                    if method not in entry['CAMAPScore']:
                                        entry['CAMAPScore'][method] = {}
                                    if context not in entry['CAMAPScore'][method]:
                                        entry['CAMAPScore'][method][context] = {}
                                    if params not in entry['CAMAPScore'][method][context]:
                                        entry['CAMAPScore'][method][context][params] = []
                                        jobs.add((context, method, params))
                    if len(seqs) == chunk_size:
                        run_models()  # will get jobs and seqs values from the outer scope
                        jobs = set()
                        seqs = {}
                        gc.collect()
            #if counter >= 1000:
            #    break
        if jobs:
            run_models()

        # Fill dictionnary
        edited = False
        for counter_ix in seq_scores:
            pep, gene_name, entry_i = seq_ids[counter_ix]
            for (m, c, p), scores in seq_scores[counter_ix].items():
                pepdict[pep]['genes'][gene_name][entry_i]['CAMAPScore'][m][c][p] = scores
                edited = True

        if edited:
            self._keep_a_copy(pfile)
            pkl.dump(pepdict, open(pfile, 'wb'))

        return counter


    def merge_netmhc(self, overwrite=False):
        self.overwrite = overwrite
        self.columns_to_keep = ['Peptide', 'nM', 'Rank']

        netmhc_folder = 'NetMHCpan-4.0a'
        nmpan_out = os.path.join(self.out_dir, netmhc_folder)
        nmpan_pred_dir = os.path.join(nmpan_out, 'predictions')

        if not os.path.exists(nmpan_out):
            if os.path.exists(nmpan_out + '.tar.gz'):
                sys.stdout.write('Predictions folder already compressed. Uncompress if you want to rerun analysis.\n')
                return
            else:
                sys.stderr.write('ERROR: Predictions folder not found.\n')
                sys.exit(1)

        # Prepare jobs
        job_dct = defaultdict(lambda: [[], []])
        for fn_allele in os.listdir(nmpan_pred_dir):
            pinit, plen = fn_allele.split('.')[:2]
            pfile = os.path.join(self.out_dir, 'peptides_%s%s.pkl' % (pinit, plen))
            if fn_allele.split('.')[3] == 'tsv':
                job_dct[pfile][0].append(os.path.join(nmpan_pred_dir, fn_allele))
            elif fn_allele.split('.')[3] == 'NP':
                job_dct[pfile][1].append(os.path.join(nmpan_pred_dir, fn_allele))

        # Run jobs
        with self.executor(max_workers=self.workers, use_threads=False) as ex:
            allele_sets = ex.map(self._merge_netmhc, tuple(job_dct.keys()), *zip(*job_dct.values()))
            exit_code = self._tar_compress(nmpan_out)  # run while peptide files are being filled

        # Update info file
        alleles = set()
        for al in allele_sets:
            alleles.update(al)
        pfile = os.path.join(self.out_dir, 'info.pkl')
        info = pkl.load(open(pfile, 'rb'))
        edited = False
        if 'allelesNetMHC' not in info:
            info['NetMHCVersion'] = set([netmhc_folder])
            edited = True
        if 'allelesNetMHC' not in info:
            info['allelesNetMHC'] = alleles
            edited = True
        if edited:
            info['date'] = datetime.now()
            self._keep_a_copy(pfile)
            pkl.dump(info, open(pfile, 'wb'))

        # Delete NetMHC folder if compression was successful
        if exit_code:
            print('WARNING: tar compression failed for some reason')
        else:
            shutil.rmtree(nmpan_out)


    def _merge_netmhc(self, pfile, netmhc_ba_out_list, netmhc_np_out_list):
        print(pfile)
        pepdict = pkl.load(open(pfile, 'rb'))

        edited = False
        alleles = set()

        def concat(file_list):
            dfs = []
            for fn in file_list:
                with open(fn, 'r') as f:
                    allele = f.readline().strip().replace(':', '')
                    temp_df = pd.read_csv(f, sep='\t', header=0, usecols=self.columns_to_keep, index_col=0)
                temp_df.columns = pd.MultiIndex.from_product([[allele], temp_df.columns])
                dfs.append(temp_df)
            df = pd.concat(dfs, axis=1, sort=False, copy=False)
            return df

        ba_df = concat(netmhc_ba_out_list)
        np_df = concat(netmhc_np_out_list)
        np_df = np_df.drop('nM', axis=1, level=1)
        np_df.columns = np_df.columns.set_levels(np_df.columns.levels[1] + '_NP', level=1, verify_integrity=False)
        df_full = pd.concat([ba_df, np_df], axis=1, sort=False, copy=False).sort_index(axis=1)

        for pep, scores in df_full.to_dict(orient='index').items():
            if self.overwrite:
                try:
                    del pepdict[pep]['netmhcpan']
                except KeyError:
                    pass
            if 'netmhcpan' not in pepdict[pep]:
                d = defaultdict(dict)
                for allele, s in scores.items():
                    al, t = allele
                    alleles.add(al)
                    d[al][t] = s
                pepdict[pep]['netmhcpan'] = dict(d)
                edited = True

        if edited:
            self._keep_a_copy(pfile)
            pkl.dump(pepdict, open(os.path.join(pfile), 'wb'))
        else:
            print('%s: Nothing changed, keeping everything as is.' % pfile)

        return alleles


    @staticmethod
    def _keep_a_copy(pfile):
        base_dir, file_name = pfile.rsplit('/', 1)
        os.makedirs(os.path.join(base_dir, 'Backup'), exist_ok=True)
        c = 0
        while True:
            f1 = pfile
            f2 = os.path.join(base_dir, 'Backup', file_name.replace('.pkl', '.%d.pkl' % c))
            try:
                if os.path.isfile(f2 + '.gz'):
                    raise OSError
                print('%s: Moving %s to %s and replacing with updated annotations.' % (pfile, f1, f2))
                shutil.move(f1, f2)
                subprocess.call(['gzip', f2], shell=False)
                break
            except OSError:
                c += 1


    @staticmethod
    def _tar_compress(folder):
        #with tarfile.open(folder + '.tar.gz', 'w:gz') as tar:
        #    tar.add(nmpan_out, arcname=os.path.basename(nmpan_out))
        exit_code = subprocess.call(['tar', 'cf', folder + '.tar.gz', folder], shell=False)
        return exit_code
