#!/usr/bin/env python

import camap.trainer as ca
import numpy as np
import torch
import pickle
import gzip
from lime import lime_text
from lime.lime_text import LimeTextExplainer
from tqdm import trange
from concurrent.futures import ProcessPoolExecutor
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor

import warnings
warnings.filterwarnings("ignore")


codonTable = {
    'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
    'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
    'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
    'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

    'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L',
    'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
    'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q',
    'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R',

    'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', 'ATG' : 'M',
    'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T',
    'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K',
    'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R', 'AGG' : 'R',

    'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V',
    'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
    'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E',
    'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G',
}


class CodonEmbeddings(object):
    def __init__(self):
        self.encoding = {}
        self.encodingLst = ['000']
        codonKeys = list(codonTable.keys())
        codonKeys.sort()
        self.encoding['000'] = 0
        for i, c in enumerate(codonKeys):
            self.encoding[c] = i + 1
            self.encodingLst.append(c)

    def encode(self, seq):
        #assert len(seq)%3 == 0
        X = []
        for c in range(0, len(seq), 3):
            codon = seq[c:c+3]
            X.append(self.encoding[codon])
        return X
    
    def __call__(self, seq):
        return self.encode(seq)


class SeqTokenizer:
    def __init__(self, model):
        self.last_seq = None
        self.model = model
        
    def tokenize(self, seq):
        tokens = []
        for i in range(0, len(seq), 3):
            codon = seq[i: i+3]
            tokens.append(codon)
        return tokens
    
    def encode(self, seq):
        self.last_seq = seq
        enc = CodonEmbeddings().encode(seq)
        enc = torch.tensor(enc, dtype=torch.long)
        out = self.model.embeds(enc)

        out = out.view(-1, 2*self.model.input_len).detach()
        return out

    def decode(self, *args, **kwargs):
        return self.last_seq

    def __call__(self, seq):
        return self.tokenize(seq)


def fix_mask(ref_seq, seq):
    if ref_seq == seq:
        return ref_seq

    ref_seq = ref_seq.split(" ")
    seq = seq.split(" ")
    for i in range(len(seq)):
        if seq[i] != ref_seq[i]:
            seq[i] = "000"
    return " ".join(seq)


def get_score(model):
    def do(x):
        out = model.embeds(x)
        out = out.view(-1, 2*model.input_len).detach()
        out = model.out_layer(out)
        out = model.probas(out)
        return out#[:, 1]
    return do


def get_score_encoded(model):
    def do(x):
        out = model.out_layer(x)
        out = model.probas(out)
        return out#[:, 1]
    return do


def get_seqs_score(model):
    score_fct = get_score_encoded(model)
    def do(seqs) :
        tok = SeqTokenizer(model)
        res = []
        if type(seqs) is not list:
            seqs = [seqs]
        ref_seq = seqs[0]
        for seq in seqs:
            seq = fix_mask(ref_seq, seq)
            seq = seq.replace("-", "").replace(" ", "")
            data = torch.tensor(tok.encode(seq))
            res.append(score_fct(data)[0].detach().numpy())
        res = np.array(res)
        return res
    return do


def prepare_seq(seq):
    seq = seq[:162] + "---"*9 +seq[162:]
    seq = " ".join([seq[c:c+3] for c in range(0, len(seq), 3)])
    return seq


def get_preferences(model, seq, show_middle=True):
    seq = seq.split(" ")
    zero_seq = []
    for i in range(len(seq)) :
        if seq[i] != "---":
            zero_seq.append("000")
        else:
            zero_seq.append("---")
    res = []
    for i in range(len(seq)) :
        curr_seq = list(zero_seq)
        curr_seq[i] = seq[i]
        curr_seq = " ".join(curr_seq)
        if seq[i] != "---":
            score = get_seqs_score(model)(curr_seq)[0][1]
            res.append(score)
        elif show_middle:
            res.append(0)
    return np.array(res)


def get_preferences_list(model, seq, show_middle=True):
    seq = seq.split(" ")
    zero_seq = []
    for i in range(len(seq)) :
        if seq[i] != "---":
            zero_seq.append("000")
        else:
            zero_seq.append("---")
    res = {}
    for i in range(len(seq)) :
        curr_seq = list(zero_seq)
        curr_seq[i] = seq[i]
        curr_seq = " ".join(curr_seq)
        if seq[i] not in res:
            res[seq[i]] = []
        if seq[i] != "---":
            score = get_seqs_score(model)(curr_seq)[0][1]
            res[seq[i]].append(score)
        elif show_middle:
            res[seq[i]].append(0)
    return res


def reduce_preferences(model, seq, reduction=np.mean, show_middle=True):
    ret = get_preferences_list(model, prepare_seq(seq), show_middle)
    res = { codon: reduction(values) for codon, values in ret.items() }
    return res


def run_lime(model, seq, num_samples = 5000):
    explainer = LimeTextExplainer(class_names=["Decoy", "Hit"], mask_string="000")
    exp = explainer.explain_instance(
        prepare_seq(seq),
        get_seqs_score(model),
        num_features=109,
        num_samples=num_samples
    )
    return exp


def get_lime_weights(model, seqs, num_samples=5000):
    res = {}
    print("===LIME===")
    for seq in seqs:
    #for i in trange(len(seqs)):
    #    seq = seqs[i]
        res[seq] = dict(run_lime(model, seq, num_samples).as_list())
    return res


def get_reduced_preferences(model, seqs, reduction=np.mean):
    res = {}
    print("===PREFS===")
    for seq in seqs:
    #for i in trange(len(seqs)):
    #    seq = seqs[i]
        res[seq] = reduce_preferences(model, seq, reduction)
    return res


def count_differences(lime_seq_dict, pref_seq_dict):
    eps = 1e-4
    err_count = 0
    for seq, c_values in lime_seq_dict.items():
        for codon, weight in c_values.items():
            if abs(weight) > eps and abs(pref_seq_dict[seq][codon] - 0.5) < eps:
                lime_value = weight > 0
                pref_value = pref_seq_dict[seq][codon] > 0.5
                if lime_value != pref_value:
                    #print(weight, pref_seq_dict[seq][codon])
                    err_count += 1
    return err_count


def load_data(replicate):
    ds = 'BLCL_GRCh37.75_padding162_maxBS500_maxContexts10_ratio5_peplen9'
    with gzip.open('output/trainDS/%s_seed%d-Metadata.pkl.gz' % (ds, replicate), 'rb') as f:
        meta = pickle.load(f)

    test_seqs = {}
    for subset, handle in (("Hits", 1), ("Decoys", 0)):
        test_seqs[subset] = []
        for seq in meta["test"][handle]:
            ctx = seq["sequenceContext"]
            pre = ctx[:162]
            post = ctx[-162:]
            ctx = pre + post
            test_seqs[subset].append(ctx)
    return test_seqs


def evaluate_models(mod_handle, mod):
    print("Model", mod_handle)
    test_seqs = load_data(mod_handle)
    errs = 0
    count = 0
    for subset in test_seqs:
        print("\t Subset", subset)
        data_to_use = test_seqs[subset]
        lime_weights = get_lime_weights(mod, data_to_use, 5000)
        prefs = get_reduced_preferences(mod, data_to_use, reduction = np.mean)
        count += len(prefs)
        errs += count_differences(lime_weights, prefs)
        print(errs/count)
        #print("Model", mod_handle, "Subset", subset, "err", errs/count)
    return errs, count


if __name__ == '__main__':
    models = {}
    for i in range(1, 13):
        #print("======%s======" % i)
        mod = ca.Classifier(108, False)
        mod = ca.load_model(
            "output/pytorch-sgd-nobias/p162e4000pep9t5-%s/test-bestMin-score.pytorch" % i,
            mod
        )
        models[i] = mod

    #errs, count = zip(*list(map(evaluate_models, models.keys(), models.values())))
    with EnhancedMPIPoolExecutor(max_workers=12) as ex:
        errs, count = zip(*list(ex.map(evaluate_models, models.keys(), models.values())))
    
    errs = sum(errs)
    count = sum(count)
    
    print(errs, count)
    print(errs/count)
