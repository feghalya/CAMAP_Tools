#!/usr/bin/env python3

from sklearn.model_selection import train_test_split
import os
import pickle as pkl
from concurrent.futures import ProcessPoolExecutor as Executor, as_completed
import random

#from pyGeno.tools.UsefulFunctions import codonTable, AATable, synonymousCodonsFrequencies
from genomeData import codonTable, AATable, synonymousCodonsFrequencies


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


class Dataset(object):
    def __init__(self, genome, cellline):
        self.output_folder = OUTPUT_FOLDER
        self.genome = genome
        self.cellline = cellline

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
            #from pyGeno.tools.UsefulFunctions import synonymousCodonsFrequencies

        self.encoding = CodonEmbeddings(synonymousCodonsFrequencies)

        assert 'GRCh' in genome or 'GRCm' in genome
        self.out_dir = os.path.join(self.output_folder, 'allPeptides_%s' % self.genome)
        self.species = 'Human' if 'GRCh' in genome else 'Mouse'

        dr = 'data/peptides'
        fn = 'detected.peptides.' + self.species + cellline
        self.massspec_file = [os.path.join(dr, f) for f in os.listdir(dr) if fn in f]
        assert len(self.massspec_file) == 1
        self.massspec_file = self.massspec_file[0]

        dr = 'data/expression'
        fn = 'isoforms.median.tpm.99p.' + self.species + cellline
        self.expression_file = [os.path.join(dr, f) for f in os.listdir(dr) if fn in f]
        assert len(self.expression_file) == 1
        self.expression_file = self.expression_file[0]
 
        with open(self.massspec_file, 'r') as f:
            self.detected_peptides = set([l.strip() for l in f.readlines()])

        with open(self.expression_file, 'r') as f:
            self.tpm_dct = {l.strip().split('\t')[0]:float(l.strip().split('\t')[1]) for l in f.readlines()}
        #print(len(self.tpm_dct))
        #self.tpm_dct = {x:y for x, y in self.tpm_dct.items() if float(y)>0.1}
        #print(len(self.tpm_dct))

        # detect pickle peptide files
        pepfiles = [os.path.join(self.out_dir, f) for f in os.listdir(self.out_dir) if f[:8] == 'peptides']
        self.pepfiles = pepfiles


    def load_peptides(self, context_len=162, max_rank=1, max_contexts=3, workers=0):
        self.max_rank = max_rank
        self.max_contexts = max_contexts
        self.context_len = context_len

        if workers:
            with Executor(max_workers=workers) as ex:
                out_dcts = ex.map(self._get_peptides, self.pepfiles)
        else:
            out_dcts = map(self._get_peptides, self.pepfiles)

        pepdict_nonsource = {}
        pepdict_source = {}
        pep_in_tr = {}
        pepexpr = {}
        pepbs = {}
        source_tr_dct = {k: False for k in self.tpm_dct.keys()}

        for nonsource_dct, source_dct, tdict, pit, pexp, pbs in out_dcts:
            pepdict_nonsource.update(nonsource_dct)
            pepdict_source.update(source_dct)
            pep_in_tr.update(pit)
            pepexpr.update(pexp)
            pepbs.update(pbs)
            for key, val in tdict.items():
                source_tr_dct[key] = source_tr_dct[key] or val

        self.source_transcripts = set([t for t, v in source_tr_dct.items() if v])

        # remove peptides present in source transcripts
        rem = set()
        for pep in list(pepdict_nonsource.keys()):
            p = pep[0]
            if pep_in_tr[p].intersection(self.source_transcripts):
                del pepdict_nonsource[pep]
                del pepexpr[pep]
                rem.add(p)
        for p in rem:
            del pep_in_tr[p]

        self.pep_in_tr = pep_in_tr
        self.pepexpr = pepexpr
        self.pepbs = pepbs

        self.peptides = [pepdict_nonsource, pepdict_source]

        self.peplist_nonsource = list(pepdict_nonsource.keys())
        self.peplist_source = list(pepdict_source.keys())


    def _get_peptides(self, pfile):
        print(pfile)

        var = 'Rank'
        if self.max_rank > 100:
            print('Turning rank into nM')
            var = 'nM'

        pdct = pkl.load(open(pfile, 'rb'))
        peptides = [{}, {}]
        # True == source peptide in transcript, else False
        transdct = {k: False for k in self.tpm_dct.keys()}
        pep_in_tr = {}
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
            rank = min([sc[var] for sc in mhc_scores.values()])
            if rank < self.max_rank:
                contexts = []
                transcripts = []
                tr_set = set()
                for genes in gene_entries.values():
                    for g_cont in genes:
                        tr = [t for t in g_cont['transcriptID'] if t in self.tpm_dct]
                        if len(tr):
                            contexts.append(g_cont)
                            transcripts.append(tr)
                            for t in tr:
                                transdct[t] = transdct[t] or (pep in self.detected_peptides)
                                tr_set.add(t)
                if len(contexts) <= self.max_contexts:
                    pep_in_tr[pep] = tr_set
                    for i, (g_cont, tr) in enumerate(zip(contexts, transcripts)):
                        max_exp = max([self.tpm_dct[t] for t in tr])
                        if '!GA' not in g_cont['sequenceContext']: #and if max_exp >= 1:
                            peptides[ix][(pep, i)] = g_cont
                            pepexpr[(pep, i)] = max_exp
                            pepbs[(pep, i)] = rank

        return peptides[0], peptides[1], transdct, pep_in_tr, pepexpr, pepbs


    def encode_peptides(self, ratio=5, same_tpm=False, seed=0, workers=0):
        """ The current implementation splits the original dataset, and then shuffles sequences.
        A direct consequence of this is when the same sequence is selected in more than 1 split,
        the resulting shuffling may or may not (most probable) be the same.
        """

        print('Encoding dataset')
        print(len(self.peplist_source), len(self.peplist_nonsource))

        copy_tpm_distribution = same_tpm
        #filtering = False
        #debug = False

        if copy_tpm_distribution:
            import numpy as np
            import math
            import pandas as pd

            ns_exp = pd.Series(np.log2(np.array([self.pepexpr[pep] for pep in self.peplist_nonsource])))
            s_exp = np.log2(np.array([self.pepexpr[pep] for pep in self.peplist_source]))

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
            #if filtering:
            #    filtered_nonsource = list()
            #    unfiltered_nonsource = list()
            #    for pep in self.peplist_nonsource:
            #        max_exp = self.pepexpr[pep]
            #        rank = self.pepbs[pep]
            #        #if max_exp >= 1:
            #        if rank < 500:
            #            filtered_nonsource.append(pep)
            #        else:
            #            unfiltered_nonsource.append(pep)
            #    self.peplist_nonsource = filtered_nonsource

            ratio = len(self.peplist_source)*ratio/len(self.peplist_nonsource)
            ratio = min(ratio, 1)

            if ratio == 1:
                keep = self.peplist_nonsource
            else:
                discard, keep = train_test_split(self.peplist_nonsource, test_size=ratio, random_state=seed)

            #if filtering:
            #    unfiltered_nonsource.extend(discard)
            #    print(len(self.peplist_source), len(self.peplist_nonsource), len(unfiltered_nonsource))

        dct_pep = {
            'train': [[], []],
            'test': [[], []],
            'validation': [[], []]
        }

        dct_pep['train'][0], pre_test = train_test_split(keep, test_size=0.4, random_state=seed)
        dct_pep['validation'][0], dct_pep['test'][0] = train_test_split(pre_test, test_size=0.5, random_state=seed)

        #if debug:
        #    if copy_tpm_distribution: # or filtering:
        #        nr = len(dct_pep['test'][0])/len(unfiltered_nonsource)
        #        _, dct_pep['test'][0] = train_test_split(unfiltered_nonsource, test_size=nr, random_state=seed)

        dct_pep['train'][1], pre_test = train_test_split(self.peplist_source, test_size=0.4, random_state=seed)
        dct_pep['validation'][1], dct_pep['test'][1] = train_test_split(pre_test, test_size=0.5, random_state=seed)

        self.dct_pep = dct_pep

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
                        future = ex.submit(self._encode, self.encoding, meta_lst, self.context_len, iterseed+ix)
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
                    res = self._encode(self.encoding, meta_lst, self.context_len, iterseed)
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
