#!/usr/bin/env python3

from sklearn.model_selection import train_test_split
import os
import pickle as pkl
from concurrent.futures import ProcessPoolExecutor as Executor, as_completed
import random
import numpy as np
import pandas as pd
import shutil
import subprocess
from collections import defaultdict
import math
import pkg_resources
from datetime import datetime
import tarfile

import camap.trainer as CAMAP

#from pyGeno.tools.UsefulFunctions import codonTable, AATable, synonymousCodonsFrequencies
from camaptools.GenomeData import codonTable, AATable, synonymousCodonsFrequencies
from camaptools.utils import available_models


OUTPUT_FOLDER = "./output"


# from ImmPred.arango.old_datasets.encodings
class CodonEmbeddings(object):
    def __init__(self, synonymousCodonsFrequencies=synonymousCodonsFrequencies):
        self.encoding = {}
        self.encodingLst = ['000']
        codonKeys = list(codonTable.keys())
        codonKeys.sort()
        self.encoding['000'] = 0
        for i, c in enumerate(codonKeys):
            self.encoding[c] = i + 1
            self.encodingLst.append(c)

        codon_frequencies = dict.fromkeys(AATable)
        for aa in codon_frequencies.keys():
            codons, cumweights = list(zip(*synonymousCodonsFrequencies[aa].items()))
            codon_frequencies[aa] = {'codons': codons, 'cumweights': cumweights}
        self.codon_frequencies = codon_frequencies

    def encode(self, seq, shuffle=False, seed=0):
        assert len(seq)%3 == 0
        X = []
        for c in range(0, len(seq), 3):
            seed += 1
            codon = seq[c:c+3]
            if shuffle:
                freqs = self.codon_frequencies[codonTable[codon]]
                random.seed(seed)
                codon = random.choices(freqs['codons'], freqs['cumweights'], k=1)[0]
            X.append(self.encoding[codon])
        return X


class Model(CAMAP.Classifier):
    def __init__(self, model, context, device='cpu'):
        self.context = context
        self.padding = int(self.context // 3)
        input_len = self.padding * 2
        super(Model, self).__init__(input_len)
        self.encoding = CodonEmbeddings()
        CAMAP.load_model(model, self, device)

    def _get_embeds(self, sequences):
        if type(sequences) == str:
            sequences = [sequences]
        seq_embeddings = np.empty((self.padding * 2, len(sequences)), dtype=int)
        for i, sequence in enumerate(sequences):
            if len(sequence) > self.context * 2:
                sequence = sequence[:self.context]+sequence[-self.context:]
            sequence = self.encoding.encode(sequence)
            seq_embeddings[:,i] = sequence
        return seq_embeddings

    def get_seq_score_from_nt(self, sequences):
        seq_embeddings = self._get_embeds(sequences)
        return self.get_score(seq_embeddings)


class Peptides(object):
    def __init__(self, genome, context=162):
        self.output_folder = OUTPUT_FOLDER
        self.genome = genome
        assert 'GRCh' in genome or 'GRCm' in genome
        self.out_dir = os.path.join(self.output_folder, 'allPeptides_%s' % self.genome)
        self.species = 'Human' if 'GRCh' in genome else 'Mouse'

        self.context = context

        # detect pickle peptide files
        pepfiles = [os.path.join(self.out_dir, f) for f in os.listdir(self.out_dir) if f[:8] == 'peptides']
        self.pepfiles = pepfiles


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


    def annotate(self, workers, executor=None, overwrite=False):
        if executor is None:
            executor = Executor
        self.overwrite = overwrite

        if workers:
            with executor(max_workers=workers) as ex:
                counts = sum(list(ex.map(self._annotate, self.pepfiles)))
        else:
            counts = sum(list(map(self._annotate, self.pepfiles)))

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
        import gc
        print(pfile)
        models = self.models

        def run_models():
            ix = seqs.keys()
            s_list = seqs.values()
            c, m, p = list(jobs)[0]
            s_embeds = models[c][m][p][0]._get_embeds(s_list)
            for c, m, p in jobs:
                for i in ix:
                    seq_scores[i][(m, c, p)] = []
                for model in models[c][m][p]:
                    scores = model.get_score(s_embeds)
                    for i, s in zip(ix, scores):
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
                        seqs[counter] = seq
                        seq_scores[counter] = {}
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
                        run_models()
                        jobs = set()
                        seqs = {}
                        gc.collect()
            #if counter >= 1000:
            #    break
        if jobs:
            run_models()

        # Fill dictionnary
        edited = False
        for ix in seq_scores:
            pep, gene_name, entry_i = seq_ids[ix]
            for (m, c, p), scores in seq_scores[ix].items():
                pepdict[pep]['genes'][gene_name][entry_i]['CAMAPScore'][m][c][p] = scores
                edited = True

        if edited:
            self._keep_a_copy(pfile)
            pkl.dump(pepdict, open(pfile, 'wb'))

        return counter


    def merge_netmhc(self, workers, executor=None, overwrite=False):
        if executor is None:
            executor = Executor
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
        if workers:
            with executor(max_workers=workers) as ex:
                allele_sets = ex.map(self._merge_netmhc, tuple(job_dct.keys()), *zip(*job_dct.values()))
                exit_code = self._tar_compress(nmpan_out)  # run while peptide files are being filled
        else:
            allele_sets = list(map(self._merge_netmhc, tuple(job_dct.keys()), *zip(*job_dct.values())))
            exit_code = self._tar_compress(nmpan_out)  # run after peptide files are filled

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


class Dataset(Peptides):
    def __init__(self, genome, cellline, context=162):
        super().__init__(genome, context)

        self.cellline = cellline

        self.dataset_name = self.species + '_' + self.cellline

        # calculate synonymous codon frequencies
        try:
            with open('output/codonCounts/%s.tsv' % self.genome, 'r') as f:
                codons = dict([s.strip().split('\t') for s in f])
            synonymousCodonsFrequencies = dict.fromkeys(AATable)
            for aa, cods in AATable.items():
                counts = [int(codons[c]) for c in cods]
                counts = [c/sum(counts) for c in counts]
                synonymousCodonsFrequencies[aa] = {c:ct for c, ct in zip(cods, counts)}
        except FileNotFoundError:
            print('output/codonCounts/%s.tsv not found, loading from pyGeno' % self.genome)
            from pyGeno.tools.UsefulFunctions import synonymousCodonsFrequencies

        self.encoding = CodonEmbeddings(synonymousCodonsFrequencies)

        assert 'GRCh' in genome or 'GRCm' in genome
        self.out_dir = os.path.join(self.output_folder, 'allPeptides_%s' % self.genome)
        self.species = 'Human' if 'GRCh' in genome else 'Mouse'

        dr = 'data/alleles'
        fn = self.species + '.tsv'
        alleles_df = pd.read_csv(os.path.join(dr, fn), sep='\t', index_col=0)
        colname = [x for x in alleles_df.columns if self.species+self.cellline in x]
        assert len(colname) == 1
        colname = colname[0]
        self.alleles = set(alleles_df[alleles_df[colname] == 1].index)

        dr = 'data/peptides'
        fn = 'detected.peptides.' + self.species + cellline
        massspec_file = [os.path.join(dr, f) for f in os.listdir(dr) if fn in f]
        assert len(massspec_file) == 1
        massspec_file = massspec_file[0]
        with open(massspec_file, 'r') as f:
            self.detected_peptides = set([l.strip() for l in f.readlines()])

        dr = 'data/expression'
        fn = 'isoforms.median.tpm.99p.' + self.species + cellline
        expression_file = [os.path.join(dr, f) for f in os.listdir(dr) if fn in f]
        assert len(expression_file) == 1
        expression_file = expression_file[0]
        with open(expression_file, 'r') as f:
            self.tpm_dct = {l.strip().split('\t')[0]:float(l.strip().split('\t')[1]) for l in f.readlines()}
        #self.tpm_dct = {x:y for x, y in self.tpm_dct.items() if float(y)>0.1}


    def load_peptides(self, max_bs=1250, max_contexts=3, workers=0,
            step='createDS', ann_method='Adam', ann_params='pep9t5'):
        self.max_bs = max_bs
        self.max_contexts = max_contexts
        self.step = step
        self.ann_method = ann_method  # only applicable if step=='evaluateDS'
        self.ann_params = ann_params  # only applicable if step=='evaluateDS'

        if workers:
            with Executor(max_workers=workers) as ex:
                out_dcts = ex.map(self._get_peptides, self.pepfiles)
        else:
            out_dcts = map(self._get_peptides, self.pepfiles)

        pepdict_nonsource = {}
        pepdict_source = {}
        peptr = {}
        pepexpr = {}
        pepbs = {}
        source_tr_dct = {k: False for k in self.tpm_dct.keys()}

        for nonsource_dct, source_dct, tdict, pit, pexp, pbs in out_dcts:
            pepdict_nonsource.update(nonsource_dct)
            pepdict_source.update(source_dct)
            peptr.update(pit)
            pepexpr.update(pexp)
            pepbs.update(pbs)
            for key, val in tdict.items():
                source_tr_dct[key] = source_tr_dct[key] or val

        self.source_transcripts = set([t for t, v in source_tr_dct.items() if v])

        # remove peptides present in source transcripts
        rem = set()
        for pep in list(pepdict_nonsource.keys()):
            p = pep[0]
            if peptr[p].intersection(self.source_transcripts):
                del pepdict_nonsource[pep]
                del pepexpr[pep]
                rem.add(p)
        for p in rem:
            del peptr[p]

        self.pepexpr = pepexpr
        self.pepbs = pepbs

        self.peptides = [pepdict_nonsource, pepdict_source]
        self.peplist_nonsource = list(pepdict_nonsource.keys())
        self.peplist_source = list(pepdict_source.keys())


    def _get_peptides(self, pfile):
        print(pfile)

        var = 'nM'
        if self.max_bs < 10:
            print('Turning nM into Rank')
            var = 'Rank'

        pdct = pkl.load(open(pfile, 'rb'))
        peptides = [{}, {}]
        # True == source peptide in transcript, else False
        tr_is_source = {k: False for k in self.tpm_dct.keys()}
        peptr = {}
        pepexpr = {}
        pepbs = {}

        for pep, pinfo in pdct.items():
            ix = 1 if pep in self.detected_peptides else 0
            try:
                mhc_scores = pinfo['netmhcpan']
            except KeyError as e:
                if 'U' not in pep:
                    print(pep)
                    raise e
                else:
                    continue
            gene_entries = pinfo['genes']
            mhc_scores = {al:mhc_scores[al] for al in mhc_scores if al in self.alleles}
            bs = min([sc[var] for sc in mhc_scores.values()])
            logbs = max([sc['1-log50k'] for sc in mhc_scores.values()])
            if bs < self.max_bs:
                contexts = []
                transcripts = []
                tr_set = set()
                for genes in gene_entries.values():
                    for g_cont in genes:
                        tr = [t for t in g_cont['transcriptID'] if t in self.tpm_dct]
                        if tr:
                            contexts.append(g_cont)
                            transcripts.append(tr)
                            tr_set.update(tr)
                            for t in tr:
                                # WARNING: Not taking into account min_BS
                                #          This would create a set of transcripts that are not source
                                #          but also not non-source, which are discarded anyway and never
                                #          appear in the negative dataset
                                tr_is_source[t] = tr_is_source[t] or (pep in self.detected_peptides)
                if len(contexts) <= self.max_contexts:
                    peptr[pep] = tr_set
                    pepbs[pep] = logbs if var == 'nM' else bs
                    for i, (g_cont, tr) in enumerate(zip(contexts, transcripts)):
                        exp = [self.tpm_dct[t] for t in tr]
                        if '!GA' not in g_cont['sequenceContext']:
                            pepexpr[(pep, i)] = exp
                            if self.step == 'createDS':
                                peptides[ix][(pep, i)] = g_cont
                            elif self.step == 'evaluateDS':
                                dct = g_cont['CAMAPScore']
                                dct = {k:dct[k] for k in dct if k.split('_')[0] == self.ann_method}
                                for key in dct:
                                    dct[key] = dct[key][self.context][self.ann_params]
                                peptides[ix][(pep, i)] = dct
                            else:
                                sys.stderr.write('ERROR: Unknown option given to step argument')
                                sys.exit(1)

        return peptides[0], peptides[1], tr_is_source, peptr, pepexpr, pepbs


    def split_dataset(self, ratio=5, same_tpm=False, seed=0):
        print('Splitting dataset')
        print(len(self.peplist_source), len(self.peplist_nonsource))

        copy_tpm_distribution = same_tpm
        #debug = False

        if copy_tpm_distribution:
            ns_exp = pd.Series(np.log2(np.array([max(self.pepexpr[pep]) for pep in self.peplist_nonsource])))
            s_exp = np.log2(np.array([max(self.pepexpr[pep]) for pep in self.peplist_source]))

            start, end = -10, 10
            l = [-np.inf] + list(range(start, end)) + [np.inf]
            bins = []
            try:
                for x in range(len(l)):
                    bins.append((l[x], l[x+1]))
            except IndexError:
                pass
            numtoix = {x: xi for xi, x in enumerate(range(start, end+1))}
            weights = [len(s_exp[(s_exp >= i) & (s_exp < j)])/len(s_exp) for i, j in bins]
            getbin = lambda x: bins[numtoix[max(start, min(math.ceil(x), end))]]
            getweight = lambda x: weights[numtoix[max(start, min(math.ceil(x), end))]]
            ns_norm = ns_exp.apply(getbin)
            ns_norm = ns_norm.groupby(ns_norm).count().to_dict()
            new_weights = [w/ns_norm[b] if w else 0 for b, w in zip(bins, weights)]
            getnewweight = lambda x: new_weights[numtoix[max(start, min(math.ceil(x), end))]]
            ns_weights = ns_exp.apply(getnewweight)

            np.random.seed(seed)
            # do this to avoid `ValueError: a must be 1-dimensional` because our items are tuples
            r = range(len(self.peplist_nonsource))
            filtered_nonsource = np.random.choice(r, len(self.peplist_source)*ratio, False, ns_weights)
            filtered_nonsource = [self.peplist_nonsource[i] for i in filtered_nonsource]

            ns_set = set(filtered_nonsource)
            unfiltered_nonsource = [pep for pep in self.peplist_nonsource if pep not in ns_set]
            keep = filtered_nonsource
            print(len(self.peplist_source), len(filtered_nonsource), len(unfiltered_nonsource))

        else:
            ratio = len(self.peplist_source)*ratio/len(self.peplist_nonsource)
            ratio = min(ratio, 1)

            if ratio == 1:
                keep = self.peplist_nonsource
            else:
                discard, keep = train_test_split(self.peplist_nonsource, test_size=ratio, random_state=seed)

        dct_pep = {
            'train': [[], []],
            'test': [[], []],
            'validation': [[], []]
        }

        dct_pep['train'][0], pre_test = train_test_split(keep, test_size=0.4, random_state=seed)
        dct_pep['validation'][0], dct_pep['test'][0] = train_test_split(pre_test, test_size=0.5, random_state=seed)

        #if debug:
        #    if copy_tpm_distribution:
        #        nr = len(dct_pep['test'][0])/len(unfiltered_nonsource)
        #        _, dct_pep['test'][0] = train_test_split(unfiltered_nonsource, test_size=nr, random_state=seed)

        dct_pep['train'][1], pre_test = train_test_split(self.peplist_source, test_size=0.4, random_state=seed)
        dct_pep['validation'][1], dct_pep['test'][1] = train_test_split(pre_test, test_size=0.5, random_state=seed)

        self.dct_pep = dct_pep


    def encode_peptides(self, seed=0, workers=0):
        """ The current implementation splits the original dataset, and then shuffles sequences.
        A direct consequence of this is when the same sequence is selected in more than 1 split,
        the resulting shuffling may or may not (most probable) be the same.
        """

        print('Encoding dataset')

        enc_dct = {
            'train': [[], []],
            'test': [[], []],
            'validation': [[], []]
        }

        shuff_enc_dct = {
            'train': [[], []],
            'test': [[], []],
            'validation': [[], []]
        }

        meta_dct = {
            'train': [[], []],
            'test': [[], []],
            'validation': [[], []]
        }

        def fill_dcts(ds, i, res):
            meta_lst, enc_lst, shuff_enc_lst = res
            meta_dct[ds][i].extend(meta_lst)
            enc_dct[ds][i].extend(enc_lst)
            shuff_enc_dct[ds][i].extend(shuff_enc_lst)

        random.seed(seed)
        iterseed = random.randint(10**6, 10**9)

        if workers:
            ex = Executor(max_workers=workers)
            chunk_size = 10000
            futures = {}
            for ds in dct_pep:
                for i in range(len(dct_pep[ds])):
                    iterseed += 1
                    for ix in range(0, len(dct_pep[ds][i]), chunk_size):
                        pep_list = dct_pep[ds][i][ix:ix+chunk_size]
                        meta_lst = [self.peptides[i][pep] for pep in pep_list]
                        future = ex.submit(self._encode, self.encoding, meta_lst, self.context, iterseed+ix)
                        futures[future] = (ds, i)
            print('Futures: %d' % len(futures))
            for future in as_completed(futures):
                (ds, i) = futures[future]
                res = future.result()
                fill_dcts(ds, i, res)
            ex.shutdown()
        else:
            for ds in dct_pep:
                for i in range(len(dct_pep[ds])):
                    iterseed += 1
                    pep_list = dct_pep[ds][i]
                    meta_lst = [self.peptides[i][pep] for pep in pep_list]
                    res = self._encode(self.encoding, meta_lst, self.context, iterseed)
                    fill_dcts(ds, i, res)

        self.meta_dct = meta_dct
        self.enc_dct = enc_dct
        self.shuff_enc_dct = shuff_enc_dct


    @staticmethod
    def _encode(encoding, meta_lst, context_len, seed):
        #print('encoding... %d' % len(meta_lst))
        embeds_lst = []
        shuffle_embeds_lst = []

        random.seed(seed)
        iterseed = random.randint(10**6, 10**9)

        for meta in meta_lst:
            iterseed += 1
            context = meta['sequenceContext']
            context = context[:context_len] + context[-context_len:]

            embeds = encoding.encode(context)
            shuffle_embeds = encoding.encode(context, shuffle=True, seed=iterseed)

            embeds_lst.append(embeds)
            shuffle_embeds_lst.append(shuffle_embeds)

        return meta_lst, embeds_lst, shuffle_embeds_lst
