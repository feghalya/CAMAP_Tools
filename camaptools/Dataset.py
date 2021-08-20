#!/usr/bin/env python3

import os
import pickle as pkl
import numpy as np
import pandas as pd
from collections import defaultdict
import math
import random
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import gc
import sys

from camaptools.GenomeData import AATable, synonymousCodonsFrequencies
from camaptools.MLP import CodonEmbeddings, CodonShuffleEmbeddings, AAEmbeddings
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor, as_completed
from camaptools.utils import ROOT, OUTPUT_FOLDER


class Dataset(object):
    def __init__(self, genome, cellline, context=162, workers=0, executor=None):
        self.workers = workers
        self.executor = EnhancedProcessPoolExecutor if executor is None else executor

        self.output_folder = OUTPUT_FOLDER
        self.genome = genome
        assert 'GRCh' in genome or 'GRCm' in genome
        self.out_dir = os.path.join(self.output_folder, 'allPeptides_%s' % self.genome)
        self.species = 'Human' if 'GRCh' in genome else 'Mouse'
        self.context = context

        # detect pickle peptide files
        fs = os.listdir(self.out_dir)
        pepfiles = [os.path.join(self.out_dir, f) for f in fs if f[:8] == 'peptides']
        self.pepfiles = pepfiles

        self.cellline = cellline
        self.dataset_name = self.species + '_' + self.cellline

        # calculate synonymous codon frequencies
        try:
            count_path = os.path.join(self.output_folder, 'codonCounts/%s.tsv' % self.genome)
            with open(count_path, 'r') as f:
                codons = dict([s.strip().split('\t') for s in f])
            synonymousCodonsFrequencies = dict.fromkeys(AATable)
            for aa, cods in AATable.items():
                counts = [int(codons[c]) for c in cods]
                counts = [c/sum(counts) for c in counts]
                synonymousCodonsFrequencies[aa] = {c:ct for c, ct in zip(cods, counts)}
        except FileNotFoundError:
            print(
                os.path.join(self.output_folder, 'codonCounts/%s.tsv' % self.genome) +\
                ' not found, loading from pyGeno'
            )
            from pyGeno.tools.UsefulFunctions import synonymousCodonsFrequencies

        self.encoding_codon = CodonEmbeddings()
        self.encoding_shuffle = CodonShuffleEmbeddings(synonymousCodonsFrequencies)
        self.encoding_aa = AAEmbeddings()

        dr = os.path.join(ROOT, 'data/alleles')
        fn = self.species + '.tsv'
        alleles_df = pd.read_csv(os.path.join(dr, fn), sep='\t', index_col=0)
        colname = [x for x in alleles_df.columns if self.species+self.cellline in x]
        assert len(colname) == 1
        colname = colname[0]
        self.alleles = set(alleles_df[alleles_df[colname] == 1].index)

        dr = os.path.join(ROOT, 'data/peptides')
        fn = 'detected.peptides.' + self.species + cellline
        massspec_file = [os.path.join(dr, f) for f in os.listdir(dr) if fn in f]
        assert len(massspec_file) == 1
        massspec_file = massspec_file[0]
        with open(massspec_file, 'r') as f:
            self.detected_peptides = set([l.strip() for l in f.readlines()])

        dr = os.path.join(ROOT, 'data/expression')
        fn = 'isoforms.median.tpm.99p.' + self.species + cellline
        expression_file = [os.path.join(dr, f) for f in os.listdir(dr) if fn in f]
        assert len(expression_file) == 1
        expression_file = expression_file[0]
        with open(expression_file, 'r') as f:
            data = [l.strip().split('\t') for l in f.readlines()]
        self.tpm_dct = {l[0]:float(l[1]) for l in data}

        count_path = os.path.join(self.output_folder, 'codonCounts/%s.trscp.tsv.gz' % self.genome)
        try:
            df = pd.read_csv(count_path, sep='\t', names=['Transcript', 'Codon', 'Count'])
            df = df[df.Codon != '!GA']
            expdct = {'Transcript': list(self.tpm_dct.keys()), 'TPM': list(self.tpm_dct.values())}
            dfexp = pd.DataFrame(expdct)
            df = df.merge(dfexp, on='Transcript', how='left')
            df['WeightedCount'] = df.TPM * df.Count
            df = df[df.WeightedCount > 0]
            self.tpm_codon_dct = df.groupby('Codon').WeightedCount.sum().to_dict()
        except FileNotFoundError:
            sys.stderr.write('WARNING: File not found. => ' + count_path + '\n')
            self.tpm_codon_dct = None


    def load_peptides_options(self,
            max_bs_or_rank=1250,
            var='nM',
            max_contexts=3,
            step='createDS',
            ann_method='Adam',
            ann_params='e500pep9t5',
            max_replicates=None,
            peptide_lengths=[8, 9, 10, 11],
            remove_neg_source=True
        ):

        self.max_bs_or_rank = max_bs_or_rank
        self.max_contexts = max_contexts if not max_contexts is None else np.inf
        self.step = step
        self.ann_method = ann_method  # only applicable if step is 'evaluateDS' or 'getDS'
        self.ann_params = ann_params  # only applicable if step is 'evaluateDS' or 'getDS'
        self.max_replicates = max_replicates
        self.remove_neg_source = remove_neg_source

        if self.step not in ['createDS', 'evaluateDS', 'getDS']:
            sys.stderr.write('WARNING: Unknown option given to step argument.\n')
            self.step = 'createDS'

        if self.max_bs_or_rank < 100 and var != 'Rank' and var != 'Rank_NP':
            print('nM < 100 --> Using Rank')
            var = 'Rank'
        elif self.max_bs_or_rank >= 100 and var != 'nM':
            print('Rank >= 100 --> Using nM')
            var = 'nM'
        else:
            assert var in ['nM', 'Rank', 'Rank_NP']
        self.var = var

        peplens = peptide_lengths if type(peptide_lengths) == list else [peptide_lengths]
        pepfiles = []
        for f in self.pepfiles:
            for peplen in peplens:
                if int(os.path.basename(f).split('_')[1][1:].replace('.pkl', '')) == peplen:
                    pepfiles.append(f)
                    break
        self.pepfiles = pepfiles


    def load_peptides(self):
        if not hasattr(self, 'var'):
            sys.stderr.write('WARNING: Load peptides options not set, using default values\n')
            self.load_peptides_options()

        # Construct datasets
        workers = min(self.workers, len(self.pepfiles))

        with self.executor(max_workers=workers, use_threads=False) as ex:
            out_dcts = ex.map(self._load_peptides, self.pepfiles)

        peptides = [{}, {}]
        source_tr_dct = {k: False for k in self.tpm_dct.keys()}

        for peps, tdict in out_dcts:
            ns, s = peps
            peptides[0].update(ns)
            peptides[1].update(s)
            for key, val in tdict.items():
                source_tr_dct[key] = source_tr_dct[key] or val

        source_transcripts = set([t for t, v in source_tr_dct.items() if v])

        # remove nonsource peptides present in source transcripts
        if self.remove_neg_source:
            key = 'expressedTranscripts'
            for pep in list(peptides[0].keys()):
                if set(peptides[0][pep][key]).intersection(source_transcripts):
                    del peptides[0][pep]

        self.peptides = peptides
        self.source_transcripts = source_transcripts


    def _load_peptides(self, pfile):
        print(pfile)

        pdct = pkl.load(open(pfile, 'rb'))
        peptides = [{}, {}]
        tr_is_source = {k: False for k in self.tpm_dct.keys()}
        #^ True == source peptide in transcript, else False

        for pep, pinfo in pdct.items():
            ix = 1 if pep in self.detected_peptides else 0

            if 'U' in pep:
                continue
            else:
                assert 'netmhcpan' in pinfo

            mhc_scores = pinfo['netmhcpan']
            mhc_scores = {al:mhc_scores[al] for al in mhc_scores if al in self.alleles}
            bs = min([sc[self.var] for sc in mhc_scores.values()])

            gene_entries = pinfo['genes']

            if bs < self.max_bs_or_rank:
                contexts = []
                transcripts = []
                for genes in gene_entries.values():
                    for g_cont in genes:
                        tr = [t for t in g_cont['transcriptID'] if t in self.tpm_dct]
                        if tr:
                            contexts.append(g_cont)
                            transcripts.append(tr)
                            for t in tr:
                                # WARNING
                                #   Not taking into account min_BS
                                #   This would create a set of transcripts that are not source
                                #   but also not non-source, which are discarded anyway and never
                                #   appear in the negative dataset
                                tr_is_source[t] = tr_is_source[t] or (pep in self.detected_peptides)

                if len(contexts) and len(contexts) <= self.max_contexts:
                    bs = 1 - np.log(bs) / np.log(50000) if self.var == 'nM' else bs
                    for i, (g_cont, tr) in enumerate(zip(contexts, transcripts)):
                        exp = [self.tpm_dct[t] for t in tr]
                        if '!GA' not in g_cont['sequenceContext']:
                            peptides[ix][(pep, i)] = {}

                            transid = g_cont['transcriptID']
                            positions = g_cont['positionInTranscript']
                            positions = [x for x, y in zip(positions, transid) if y in tr]

                            seq = g_cont['sequenceContext']
                            context_size = (len(seq) - len(pep)*3) / 2
                            assert context_size == int(context_size)
                            context_size = int(context_size)
                            seq = seq[:context_size] + seq[-context_size:]
                            gccontent = (seq.count('G') + seq.count('C')) / len(seq)

                            if self.step == 'createDS' or self.step == 'getDS':
                                for k in ['peptide', 'sequenceContext', 'geneID', 'geneName']:
                                    peptides[ix][(pep, i)][k] = g_cont[k]
                                peptides[ix][(pep, i)]['bindingAffinityMetric'] = self.var

                            peptides[ix][(pep, i)]['bindingAffinity'] = bs
                            peptides[ix][(pep, i)]['expressedTranscripts'] = tr
                            peptides[ix][(pep, i)]['transcriptsExpression'] = exp
                            peptides[ix][(pep, i)]['transcriptsPosition'] = positions
                            peptides[ix][(pep, i)]['sequenceGCContent'] = gccontent

                            if 'CAMAPScore' in g_cont and self.step == 'getDS':
                                peptides[ix][(pep, i)]['CAMAPMethod'] = self.ann_method
                                peptides[ix][(pep, i)]['CAMAPParameters'] = self.ann_params
                                peptides[ix][(pep, i)]['CAMAPContext'] = self.context

                            if 'CAMAPScore' in g_cont and self.step in ['evaluateDS', 'getDS']:
                                dct = g_cont['CAMAPScore']
                                camap_score = None
                                camap_score_shuffle = None

                                for k, v in dct.items():
                                    c, a, m = self.context, self.ann_params, self.max_replicates
                                    if k == self.ann_method:
                                        camap_score = dct[k][c][a][:m]
                                    elif k == self.ann_method + '_Shuffle':
                                        camap_score_shuffle = dct[k][c][a][:m]

                                peptides[ix][(pep, i)]['CAMAPScore'] = camap_score
                                peptides[ix][(pep, i)]['CAMAPScoreShuffle'] = camap_score_shuffle

        return peptides, tr_is_source


    def clear_unused(self):
        print('Clearing unused data')
        try:
            del self.peptides
        except AttributeError:
            print("Some attributes were not cleared because they didn't exist")
        gc.collect()


class TrainingDataset(Dataset):
    def split_dataset(self, ratio=5, same_tpm=False, same_bs=False, same_gc=False, seed=None):
        print('Splitting dataset')

        assert sum([same_tpm, same_bs, same_gc]) < 2

        if seed is None:
            seed = random.randint(0, 1e9)

        do_binning = same_tpm or same_bs or same_gc if ratio is not None else False
        #debug = False
        ratio = ratio if not ratio is None else np.inf

        peplist_nonsource, peplist_source = list(self.peptides[0].keys()), list(self.peptides[1].keys())
        print(len(peplist_source), len(peplist_nonsource))

        if do_binning:
            start, end = -10, 10

            if same_tpm:
                key = 'transcriptsExpression'
                ns_exp = [np.mean(self.peptides[0][pep][key]) for pep in peplist_nonsource]
                ns_exp = pd.Series(np.log2(np.array(ns_exp)))
                s_exp = [np.mean(self.peptides[1][pep][key]) for pep in peplist_source]
                s_exp = np.log2(np.array(s_exp))
            elif same_bs:
                key = 'bindingAffinity'
                ns_exp = pd.Series(np.array([self.peptides[0][pep][key] for pep in peplist_nonsource]))
                s_exp = np.array([self.peptides[1][pep][key] for pep in peplist_source])
                max_val = max(max(ns_exp), max(s_exp))
                ns_exp = ns_exp/max_val*30-15
                s_exp = s_exp/max_val*30-15
            elif same_gc:
                key = 'sequenceGCContent'
                ns_exp = pd.Series([self.peptides[0][pep][key] for pep in peplist_nonsource])
                s_exp = np.array([self.peptides[1][pep][key] for pep in peplist_source])
                ns_exp = ns_exp*30-15
                s_exp = s_exp*30-15

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
            r = range(len(peplist_nonsource))
            filtered_nonsource = np.random.choice(r, len(peplist_source)*ratio, False, ns_weights)
            filtered_nonsource = [peplist_nonsource[i] for i in filtered_nonsource]

            ns_set = set(filtered_nonsource)
            unfiltered_nonsource = [pep for pep in peplist_nonsource if pep not in ns_set]
            keep = filtered_nonsource
            print(len(peplist_source), len(filtered_nonsource), len(unfiltered_nonsource))

        else:
            ratio = len(peplist_source)*ratio/len(peplist_nonsource)
            ratio = min(ratio, 1)

            if ratio == 1:
                keep = peplist_nonsource
            else:
                discard, keep = train_test_split(
                    peplist_nonsource,
                    test_size=ratio,
                    random_state=seed
                )

        dct_pep = {
            'train': [[], []],
            'test': [[], []],
            'validation': [[], []]
        }

        dct_pep['train'][0], pre_test = train_test_split(
            keep,
            test_size=0.4,
            random_state=seed
        )
        dct_pep['validation'][0], dct_pep['test'][0] = train_test_split(
            pre_test,
            test_size=0.5,
            random_state=seed
        )

        #if debug:
        #    if do_binning:
        #        nr = len(dct_pep['test'][0])/len(unfiltered_nonsource)
        #        _, dct_pep['test'][0] = train_test_split(
        #            unfiltered_nonsource,
        #            test_size=nr,
        #            random_state=seed
        #        )

        dct_pep['train'][1], pre_test = train_test_split(
            peplist_source,
            test_size=0.4,
            random_state=seed
        )
        dct_pep['validation'][1], dct_pep['test'][1] = train_test_split(
            pre_test,
            test_size=0.5,
            random_state=seed
        )

        self.dct_pep = dct_pep


    def encode_peptides(self, seed=0, codon=True, shuffle=True, aa=True):
        """ The current implementation splits the original dataset, and then shuffles sequences.
        A direct consequence of this is when the same sequence is selected in more than 1 split,
        the resulting shuffling may or may not (most probable) be the same.

        """
        print('Encoding dataset')

        def generate_barebone_dct():
            return {
                'train': [[], []],
                'test': [[], []],
                'validation': [[], []]
            }

        self.meta_dct = generate_barebone_dct()
        if codon:
            self.enc_dct = generate_barebone_dct()
        if shuffle:
            self.shuff_enc_dct = generate_barebone_dct()
        if aa:
            self.aa_enc_dct = generate_barebone_dct()

        random.seed(seed)
        iterseed = random.randint(10**6, 10**9)

        ex = self.executor(max_workers=self.workers, use_threads=False)

        chunk_size = 10000
        futures = {}
        for ds in self.dct_pep:
            for i in range(len(self.dct_pep[ds])):
                iterseed += 1
                for ix in range(0, len(self.dct_pep[ds][i]), chunk_size):
                    pre_meta_lst = self.dct_pep[ds][i][ix:ix+chunk_size]
                    meta_lst = [self.peptides[i][pep] for pep in pre_meta_lst]
                    self.meta_dct[ds][i].extend(meta_lst)
                    iseed = iterseed+ix
                    if codon:
                        future = ex.submit(
                            self._encode_codon,
                            self.encoding_codon,
                            meta_lst,
                            self.context
                        )
                        futures[future] = (ds, i, 'codon')
                    if shuffle:
                        future = ex.submit(
                            self._encode_shuffle,
                            self.encoding_shuffle,
                            meta_lst,
                            self.context,
                            iseed
                        )
                        futures[future] = (ds, i, 'shuffle')
                    if aa:
                        future = ex.submit(
                            self._encode_aa,
                            self.encoding_aa,
                            meta_lst,
                            self.context
                        )
                        futures[future] = (ds, i, 'aa')

        print('Futures: %d' % len(futures))

        for future in as_completed(futures):
            (ds, i, enc_type) = futures[future]
            enc_lst = future.result()
            if enc_type == 'codon':
                self.enc_dct[ds][i].extend(enc_lst)
            elif enc_type == 'shuffle':
                self.shuff_enc_dct[ds][i].extend(enc_lst)
            elif enc_type == 'aa':
                self.aa_enc_dct[ds][i].extend(enc_lst)

        ex.shutdown()


    @staticmethod
    def _encode_codon(encoding, meta_lst, context_len):
        #print('encoding... %d' % len(meta_lst))
        embeds_lst = []
        for meta in meta_lst:
            context = meta['sequenceContext']
            context = context[:context_len] + context[-context_len:]
            embeds = encoding.encode(context)
            embeds_lst.append(embeds)
        return embeds_lst


    @staticmethod
    def _encode_shuffle(encoding, meta_lst, context_len, seed):
        embeds_lst = []
        random.seed(seed)
        iterseed = random.randint(10**6, 10**9)
        for meta in meta_lst:
            iterseed += 1
            context = meta['sequenceContext']
            context = context[:context_len] + context[-context_len:]
            embeds = encoding.encode(context, seed=iterseed)
            embeds_lst.append(embeds)
        return embeds_lst


    @staticmethod
    def _encode_aa(encoding, meta_lst, context_len):
        embeds_lst = []
        for meta in meta_lst:
            context = meta['sequenceContext']
            context = context[:context_len] + context[-context_len:]
            embeds = encoding.encode(context, codons=True)
            embeds_lst.append(embeds)
        return embeds_lst


class RegressionDataset(Dataset):
    def construct_datasets_options(self, peptides_train=[], peptides_test=[]):
        self.peptides_train = peptides_train
        self.peptides_test = peptides_test


    def construct_datasets(self):
        """Note that self.load_peptides needs to be executed for all required
        variables to get populated

        """
        try:
            n_replicates = len(self.df_dictionnaries[0])
        except AttributeError:
            self.df_dictionnaries = self.organize_data(self.peptides)
            n_replicates = len(self.df_dictionnaries[0])

        if not hasattr(self, 'peptides_train'):
            self.construct_datasets_options()

        # Construct datasets

        workers = min(self.workers, n_replicates)
        #subworkers = (self.workers - workers) // n_replicates
        #subworkers = 0
        # _process_df does nothing (for now)
        # additionaly, we cannot fork if using MPIPoolExecutor

        with self.executor(max_workers=workers, use_threads=False) as ex:
            processes = []
            for i, (nonsource, source) in enumerate(zip(*self.df_dictionnaries)):
                seed = i + 1
                p = ex.submit(
                    self._create_preprocess_split,
                    seed,
                    nonsource,
                    source,
                    self.peptides_train,
                    self.peptides_test
                )  #self.executor, subworkers)
                processes.append(p)
            print(ex)
        datasets, metadata = zip(*[p.result() for p in processes])

        self.datasets = [y for x in datasets for y in x]
        self.metadata = [y for x in metadata for y in x]


    @staticmethod
    def organize_data(all_peptides):
        df_dictionnaries = [[], []]

        for sns in [0, 1]:
            peptides = all_peptides[sns]

            df_dct = {
                'Peptide': [],
                'BS': [],
                'TPM': [],
                'CAMAP': defaultdict(list),
                'CAMAPShuff': defaultdict(list)
                }

            df_dct = defaultdict(list, df_dct)

            for pep in peptides.keys():
                dct = peptides[pep]

                expressions = dct['transcriptsExpression']
                try:
                    transcripts = dct['expressedTranscripts']
                except KeyError:
                    transcripts = [None]*len(expressions)
                try:
                    positions = dct['transcriptsPosition']
                except KeyError:
                    positions = [None]*len(expressions)

                for e, t, p in zip(expressions, transcripts, positions):
                    df_dct['Peptide'].append(pep[0])
                    df_dct['BS'].append(dct['bindingAffinity'])
                    df_dct['TPM'].append(np.log10(e))

                    for i, v in enumerate(dct['CAMAPScore']):
                        df_dct['CAMAP'][i].append(v)
                    for i, v in enumerate(dct['CAMAPScoreShuffle']):
                        df_dct['CAMAPShuff'][i].append(v)

                    #if 'expressedTranscripts' in dct:
                    #    df_dct['Transcript'].append(t)
                    #if 'transcriptsPosition' in dct:
                    #    df_dct['Position'].append(p)
                    #if 'sequenceContext' in dct:
                    #    df_dct['ContextSequence'].append(dct['sequenceContext'])
                    if 'sequenceGCContent' in dct:
                        df_dct['GCContent'].append(dct['sequenceGCContent'])

            df_dct = dict(df_dct)

            df_dct['CAMAP'] = [df_dct['CAMAP'][k] for k in sorted(df_dct['CAMAP'])]
            df_dct['CAMAPShuff'] = [df_dct['CAMAPShuff'][k] for k in sorted(df_dct['CAMAPShuff'])]

            assert len(df_dct['CAMAP']) == len(df_dct['CAMAPShuff'])
            n_replicates = len(df_dct['CAMAP'])

            for k in [x for x in df_dct.keys() if x != 'CAMAP' and x != 'CAMAPShuff']:
                df_dct[k] = [df_dct[k]] * n_replicates

            for r in range(n_replicates):
                subdct = {c:df_dct[c][r] for c in df_dct}
                df_dictionnaries[sns].append(subdct)

        if len(df_dictionnaries[0]) != len(df_dictionnaries[1]):
            if len(df_dictionnaries[0]):
                assert not len(df_dictionnaries[1])
            if not len(df_dictionnaries[0]):
                empty = 0
                full = 1
            else:
                empty = 1
                full = 0
            for l in df_dictionnaries[full]:
                df_dictionnaries[empty].append(df_dct)

        print('Peptide data organized in replicates')

        return df_dictionnaries


    @classmethod
    def _create_preprocess_split(cls, seed, nonsource, source, peptides_train, peptides_test):  #executor, workers):
        print('doing %d' % seed)
        metadata = []
        datasets = []
        if True:
        #with executor(max_workers=workers, use_threads=True) as ex:
            #procs = []
            res = cls._preprocess_df(cls._create_df(nonsource, source))
            for df, subname in res:
                dat = cls._split_and_normalize_df(df, seed, peptides_train, peptides_test)
                metadata.append((seed, subname))
                datasets.append(dat)
                #p = ex.submit(cls._split_and_normalize_df, df, seed)
                #procs.append(p)
            #print(ex)
        #datasets = [p.result() for p in procs]
        print('done %d' % seed)
        return datasets, metadata


    @staticmethod
    def _create_df(ns, s):
        df_ns = pd.DataFrame(ns).set_index('Peptide')
        df_s = pd.DataFrame(s).set_index('Peptide')
        df_ns['Class'] = False
        df_s['Class'] = True
        df_full = pd.concat([df_ns, df_s], axis=0).reset_index().set_index(['Peptide', 'Class'])
        return df_full


    @staticmethod
    def _preprocess_df(df):
        yield (df, 'full')


    @staticmethod
    def _split_and_normalize_df(df, seed, ix_train=[], ix_test=[]):
        #try:
        if True:
            ix = df.index.drop_duplicates()

            # Pre-arrange train and test entries if specified
            ix_train = list(set(ix_train))
            ix_test = list(set(ix_test))

            ix_train = [x for x in ix_train if x in df.index]
            ix_test = [x for x in ix_test if x in df.index]

            nswitch = 0
            nremove = 0
            split_ratio = 0.3/0.7
            if len(ix_train)*split_ratio > len(ix_test):
                nremove = len(ix_train) - int(len(ix_test)/split_ratio)
            elif len(ix_train)*split_ratio < len(ix_test):
                nswitch = int((len(ix_test) - len(ix_train)*split_ratio) / (1 + split_ratio))

            # NOTE: This method does not preserve class proportions
            random.seed = seed
            ix_remove = random.sample(range(len(ix_train)), nremove)
            random.seed = seed
            ix_switch = random.sample(range(len(ix_test)), nswitch)

            print('Remove:', nremove)
            print('Switch:', nswitch)

            for e in sorted(ix_remove, reverse=True):
                ix_train.pop(e)
            for e in sorted(ix_switch, reverse=True):
                ix_train.append(ix_test.pop(e))

            # If ix_train and ix_test are empty, do nothing
            try:
                ix = ix.drop(ix_train + ix_test, level='Peptide')
            except KeyError:
                pass

            # Split dataset into train and test
            classes = ix.get_level_values('Class')
            peptides = ix.get_level_values('Peptide').to_list()
            split = train_test_split(peptides, stratify=classes, test_size=0.3, random_state=seed)
            ix_train_split, ix_test_split = split
            ix_train += ix_train_split
            ix_test += ix_test_split

            df_train = df.reset_index(level='Class').loc[ix_train]
            df_train = df_train.reset_index().set_index(['Peptide', 'Class'])
            df_test = df.reset_index(level='Class').loc[ix_test]
            df_test = df_test.reset_index().set_index(['Peptide', 'Class'])

            # alternative method:
            #df_train = df.loc[(ix_train, slice(None)), :]
            #df_test = df.loc[(ix_test, slice(None)), :]
            #^ terribly inefficient

            scaler = MinMaxScaler(feature_range=(0, 1))
            scaler.fit(df_train)

            def normalize_columns(df):
                df_norm = pd.DataFrame(scaler.transform(df), columns=df.columns, index=df.index)
                df_renamed = df.rename({x: x+'Base' for x in df.columns}, axis=1)
                df = pd.concat([df_norm, df_renamed], axis=1)
                df = df.reset_index().set_index('Peptide')
                return df

            df_train = normalize_columns(df_train)
            df_test = normalize_columns(df_test)

        #except ValueError:
        #    df_train = df.copy()
        #    df_train_norm = df.copy()
        #    df_train.columns = [x + 'Base' for x in df_train.columns]
        #    df_train = pd.concat([df_train_norm, df_train], axis=1)
        #    df_train = df_train.reset_index().set_index('Peptide')
        #    df_test = pd.DataFrame(columns=df_train.reset_index().columns).set_index('Peptide')

        return df_train, df_test


    def clear_unused(self):
        super().clear_unused()
        try:
            del self.df_dictionnaries
        except AttributeError:
            print("Some attributes were not cleared because they didn't exist")
        gc.collect()
