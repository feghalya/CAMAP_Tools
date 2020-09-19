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

from camaptools.GenomeData import AATable, synonymousCodonsFrequencies
from camaptools.MLP import CodonEmbeddings
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor, as_completed


OUTPUT_FOLDER = "./output"


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
        pepfiles = [os.path.join(self.out_dir, f) for f in os.listdir(self.out_dir) if f[:8] == 'peptides']
        self.pepfiles = pepfiles

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


    def load_peptides(self, max_bs_or_rank=1250, var='nM', max_contexts=3,
            step='createDS', ann_method='Adam', ann_params='pep9t5'):

        self.max_bs_or_rank = max_bs_or_rank
        self.max_contexts = max_contexts
        self.step = step
        self.ann_method = ann_method  # only applicable if step=='evaluateDS'
        self.ann_params = ann_params  # only applicable if step=='evaluateDS'

        if self.max_bs_or_rank < 100 and var != 'Rank' and var != 'Rank_NP':
            print('nM < 100 --> Using Rank')
            var = 'Rank'
        elif self.max_bs_or_rank >= 100 and var != 'nM':
            print('Rank >= 100 --> Using nM')
            var = 'nM'
        else:
            assert var in ['nM', 'Rank', 'Rank_NP']
        self.var = var

        # Construct datasets
        workers = min(self.workers, len(self.pepfiles))

        with self.executor(max_workers=workers, use_threads=False) as ex:
            out_dcts = ex.map(self._load_peptides, self.pepfiles)

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


    def _load_peptides(self, pfile):
        print(pfile)

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
            bs_or_rank = min([sc[self.var] for sc in mhc_scores.values()])
            if bs_or_rank < self.max_bs_or_rank:
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
                    pepbs[pep] = 1 - np.log(bs_or_rank) / np.log(50000) if self.var == 'nM' else bs_or_rank
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


class TrainingDataset(Dataset):
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


    def encode_peptides(self, seed=0):
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

        if self.workers:
            ex = self.executor(max_workers=self.workers, use_threads=False)
            chunk_size = 10000
            futures = {}
            for ds in self.dct_pep:
                for i in range(len(self.dct_pep[ds])):
                    iterseed += 1
                    for ix in range(0, len(self.dct_pep[ds][i]), chunk_size):
                        pep_list = self.dct_pep[ds][i][ix:ix+chunk_size]
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
            for ds in self.dct_pep:
                for i in range(len(self.dct_pep[ds])):
                    iterseed += 1
                    pep_list = self.dct_pep[ds][i]
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


class RegressionDataset(Dataset):
    def construct_datasets(self):
        """
        self.load_peptides needs to be executed for all required variables to get populated

        """
        try:
            n_replicates = len(self.df_dictionnaries[0])
        except AttributeError:
            self._organize_data()
            n_replicates = len(self.df_dictionnaries[0])

        # Construct datasets

        workers = min(self.workers, n_replicates)
        subworkers = (self.workers - workers) // n_replicates
        subworkers = 0  # _process_df does nothing (for now)

        with self.executor(max_workers=workers, use_threads=False) as ex:
            processes = []
            for i, (ns, s) in enumerate(zip(*self.df_dictionnaries)):
                seed = i + 1
                p = ex.submit(self._create_preprocess_split, seed, ns, s, self.executor, subworkers)
                processes.append(p)
            print(ex)
        datasets, metadata = zip(*[p.result() for p in processes])

        self.datasets = [y for x in datasets for y in x]
        self.metadata = [y for x in metadata for y in x]


    def _organize_data(self):
        # Organize data

        df_dictionnaries = [[], []]

        for sns in [0, 1]:
            df_dct = {
                'Peptide': [],
                'BS': [],
                'TPM': [],
                'CAMAP': defaultdict(list),
                'CAMAPShuff': defaultdict(list)
                }

            for pep in self.peptides[sns]:
                p = pep[0]
                expressions = self.pepexpr[pep]
                bs_score = self.pepbs[p]
                ann_score = self.peptides[sns][pep][self.ann_method]
                ann_shuf_score = self.peptides[sns][pep][self.ann_method + '_Shuffle']
                for e in expressions:
                    df_dct['Peptide'].append(p)
                    df_dct['BS'].append(bs_score)
                    df_dct['TPM'].append(np.log10(e))
                    for i, v in enumerate(ann_score):
                        df_dct['CAMAP'][i].append(v)
                    for i, v in enumerate(ann_shuf_score):
                        df_dct['CAMAPShuff'][i].append(v)

            df_dct['CAMAP'] = [df_dct['CAMAP'][k] for k in sorted(df_dct['CAMAP'])]
            df_dct['CAMAPShuff'] = [df_dct['CAMAPShuff'][k] for k in sorted(df_dct['CAMAPShuff'])]

            assert len(df_dct['CAMAP']) == len(df_dct['CAMAPShuff'])
            n_replicates = len(df_dct['CAMAP'])

            df_dct['Peptide'] = [df_dct['Peptide']] * n_replicates
            df_dct['BS'] = [df_dct['BS']] * n_replicates
            df_dct['TPM'] = [df_dct['TPM']] * n_replicates

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
        self.df_dictionnaries = df_dictionnaries


    @classmethod
    def _create_preprocess_split(cls, seed, ns, s, executor, workers):
        print('doing %d' % seed)
        metadata = []
        with executor(max_workers=workers, use_threads=True) as ex:
            procs = []
            for df, subname in cls._preprocess_df(cls._create_df(ns, s)):
                metadata.append((seed, subname))
                p = ex.submit(cls._split_and_normalize_df, df, seed)
                procs.append(p)
            print(ex)
        datasets = [p.result() for p in procs]
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
    def _split_and_normalize_df(df, seed):
        try:
            ix = df.index.drop_duplicates()
            classes = ix.get_level_values('Class')
            ix_train, ix_test = train_test_split(ix, stratify=classes, test_size=0.3, random_state=seed)
            df_train = df.loc[ix_train]
            df_test = df.loc[ix_test]

            scaler = MinMaxScaler(feature_range=(0, 1)) 
            scaler.fit(df_train)

            columns = df_train.columns
            df_train.columns = [x + 'Base' for x in columns]
            df_train_norm = pd.DataFrame(scaler.transform(df_train), columns=columns, index=df_train.index)
            df_train = pd.concat([df_train_norm, df_train], axis=1)
            df_train = df_train.reset_index().set_index('Peptide')

            columns = df_test.columns
            df_test.columns = [x + 'Base' for x in columns]
            df_test_norm = pd.DataFrame(scaler.transform(df_test), columns=columns, index=df_test.index)
            df_test = pd.concat([df_test_norm, df_test], axis=1)
            df_test = df_test.reset_index().set_index('Peptide')

        except ValueError:
            df_train = df.copy()
            df_train_norm = df.copy()
            df_train.columns = [x + 'Base' for x in df_train.columns]
            df_train = pd.concat([df_train_norm, df_train], axis=1)
            df_train = df_train.reset_index().set_index('Peptide')

            df_test = pd.DataFrame(columns=df_train.reset_index().columns).set_index('Peptide')

        return df_train, df_test


    def clear_unused(self):
        print('Clearing unused data')
        try:
            del self.pepbs
            del self.pepexpr
            del self.peplist_nonsource
            del self.peplist_source
            del self.peptides
            del self.df_dictionnaries
        except AttributeError:
            print("Some attributes were not cleared because they didn't exist")

        gc.collect()
