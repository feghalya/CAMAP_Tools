#!/usr/bin/env python3

import random
import numpy as np

import camap.trainer as CAMAP

#from pyGeno.tools.UsefulFunctions import codonTable, AATable, synonymousCodonsFrequencies
from camaptools.GenomeData import codonTable, AATable, synonymousCodonsFrequencies


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
