#!/usr/bin/env python3

import pandas as pd
import pickle as pkl
import gzip
import numpy as np
import os

from pyGeno.tools.UsefulFunctions import AATable, codonTable

from camaptools.MLP import CodonEmbeddings, Model, load_models
from camaptools.Dataset import Dataset

# ########## #
# Initialize #
# ########## #

os.makedirs('revisions', exist_ok=True)

model_names, models = load_models(162, ['SGD', 'SGD_Shuffle'])
model_names = model_names[162]
models = models[162]

codons = CodonEmbeddings().encoding
codons = sorted(codons.keys(), key=lambda x: codons[x])

# ########## #
# Embeddings #
# ########## #

model_name = model_names['SGD_Shuffle']['e4000'][6]
model = models['SGD_Shuffle']['e4000'][6]
n, p = zip(*model.embeds.weight.tolist())
df_embeds = pd.DataFrame({'Positive': p, 'Negative': n}, index=codons)
df_embeds['AA'] = [codonTable[c] if c != '000' else '0' for c in df_embeds.index]
print(df_embeds)
df_embeds.to_csv('revisions/embeddings.7.tsv', sep='\t', float_format='%g')

# ########### #
# Preferences #
# ########### #

npos = 162//3*2

codon_list = []
pos_list = []
anti_pref_list = []
pref_list = []

model = models['SGD_Shuffle']['e4000'][6]
for codon in codons:
    for pos in range(npos):
        seq = ['000']*npos
        seq[pos] = codon
        seq = ''.join(seq)
        seq = model._get_embeds(seq)
        seq = model._assert_type(seq)
        n, s = model.get_score(seq).tolist()[0]
        codon_list.append(codon)
        pos_list.append(pos)
        anti_pref_list.append(n)
        pref_list.append(s)
#predf = {'Codon': codon_list, 'Position': pos_list, 'NegativePref': anti_pref_list, 'PositivePref': pref_list}
predf = {'Position': pos_list, 'Preference': pref_list}
df_prefs = pd.DataFrame(predf, index=codon_list)
df_prefs['AA'] = [codonTable[c] if c != '000' else '0' for c in df_prefs.index]
print(df_prefs)
df_prefs.to_csv('revisions/preferences.7.tsv', sep='\t', float_format='%g')
#df_prefs.to_csv('revisions/preferences.original.7.tsv', sep='\t', float_format='%g')

# ############# #
# Distributions #
# ############# #

with gzip.open('../../CAMAP/training_datasets/OneTrue_maxBs-1250_minHLANb-1_minLogFPKM-0.0003_padding-162_multiTrans-True_nonSourceFactor-5_addFields-_multiGenes-True_sequenceAlteration-removed_sourceSize-19658_encoding-CodonShuffleEmbeddings.pkl.gz', 'rb') as f:
    trainshuff = pkl.load(f)

codon_source_dct = {codons[cod]: [0]*(162//3*2) for cod in range(65)}
codon_nonsource_dct = {codons[cod]: [0]*(162//3*2) for cod in range(65)}
for sns in trainshuff.values():
    nonsource_seqs = sns[0]
    for seq in nonsource_seqs:
        for pos, cod in enumerate(seq):
            codon_nonsource_dct[codons[cod]][pos] += 1
    source_seqs = sns[1]
    for seq in source_seqs:
        for pos, cod in enumerate(seq):
            codon_source_dct[codons[cod]][pos] += 1

aa_source_dct = {aa: [0]*(162//3*2) for aa in AATable.keys()}
for aa, cod in AATable.items():
    for pos in range(len(aa_source_dct[aa])):
        aa_source_dct[aa][pos] = sum([codon_source_dct[c][pos] for c in cod])
aa_nonsource_dct = {aa: [0]*(162//3*2) for aa in AATable.keys()}
for aa, cod in AATable.items():
    for pos in range(len(aa_nonsource_dct[aa])):
        aa_nonsource_dct[aa][pos] = sum([codon_nonsource_dct[c][pos] for c in cod])

df_cod_counts = pd.DataFrame(codon_source_dct)
df_cod_counts.index.name = 'Position'
df_aa_counts = pd.DataFrame(aa_source_dct)
df_aa_counts.index.name = 'Position'
df_cod_counts.to_csv('revisions/distributions.source.codons.tsv', sep='\t', float_format='%g')
df_aa_counts.to_csv('revisions/distributions.source.aa.tsv', sep='\t', float_format='%g')
df_cod_counts = pd.DataFrame(codon_nonsource_dct)
df_cod_counts.index.name = 'Position'
df_aa_counts = pd.DataFrame(aa_nonsource_dct)
df_aa_counts.index.name = 'Position'
df_cod_counts.to_csv('revisions/distributions.nonsource.codons.tsv', sep='\t', float_format='%g')
df_aa_counts.to_csv('revisions/distributions.nonsource.aa.tsv', sep='\t', float_format='%g')
print(df_cod_counts)
print(df_aa_counts)

# ############### #
# BS Correlations #
# ############### #

# Inefficient because it loads all pickles in memory for each allele (because we have no BS filter)

def allele(dat, al):
    peptides = []
    camap_scores = []
    bs_scores = []
    classes = []

    dat.alleles = [al]
    dat.load_peptides(max_bs_or_rank=50000, var='nM', max_contexts=10,
            step='evaluateDS', ann_method='SGD', ann_params='e4000')

    for pep in dat.peptides[1]:
        peptides.append(pep[0])
        camap_scores.append(dat.peptides[1][pep]['SGD'][6])
        bs_scores.append(dat.pepbs[pep[0]])
        classes.append('Source')
    for pep in dat.peptides[0]:
        peptides.append(pep[0])
        camap_scores.append(dat.peptides[0][pep]['SGD'][6])
        bs_scores.append(dat.pepbs[pep[0]])
        classes.append('NonSource')
    print('making df')
    df_scores = pd.DataFrame({'CAMAP': camap_scores, 'logBS': bs_scores, 'Class': classes}, index=peptides)
    df_scores['BS'] = 50000**(1-df_scores.logBS)
    print('saving df')
    df_scores.to_csv('revisions/allele_%s_scores.tsv' % al, sep='\t', float_format='%g')
    print(df_scores)
    for s in ['Source', 'NonSource']:
        outf = 'revisions/allele_%s_scores.%s.n2500.tsv' % (allele, s)
        df_scores[df_scores.Class == s].sample(2500).to_csv(outf, sep='\t', float_format='%g')
    dat.clear_unused()
    print('cleaned!')

alleles = Dataset('GRCh37.75', 'BLCL').alleles
datasets = [Dataset('GRCh37.75', 'BLCL', 162, workers=40) for _ in range(len(alleles))]
for dat, al in zip(datasets, alleles):
    allele(dat, al)

# ################ #
# TPM Correlations #
# ################ #

dat = Dataset('GRCh37.75', 'BLCL', 162, workers=40)
dat.load_peptides(max_bs_or_rank=50000, var='nM', max_contexts=10,
        step='evaluateDS', ann_method='SGD', ann_params='e4000')
peptides = []
camap_scores = []
classes = []
transcripts = []
for pep in dat.peptides[1]:
    for tr in dat.pepexpr[pep][1]:
        peptides.append(pep[0])
        camap_scores.append(dat.peptides[1][pep]['SGD'][6])
        classes.append('Source')
        transcripts.append(tr)
for pep in dat.peptides[0]:
    for tr in dat.pepexpr[pep][1]:
        peptides.append(pep[0])
        camap_scores.append(dat.peptides[0][pep]['SGD'][6])
        classes.append('NonSource')
        transcripts.append(tr)
df_scores = pd.DataFrame({'CAMAP': camap_scores, 'Transcript': transcripts, 'Class': classes}, index=peptides)
df_scores.to_csv('revisions/transcript_scores.tsv.gz', sep='\t', float_format='%g')
expression_data = '../data/expression/analysis/Project_HumanBLCLPerreault/tpm/isoforms.tpm.tsv'
df_tpm = pd.read_csv(expression_data, sep='\t', index_col=0)
for sample in df_tpm.columns:
    df = df_scores.merge(df_tpm[[sample]], how='left', left_on='Transcript', right_index=True)
    df = df.rename({sample: 'TPM'}, axis=1)
    df = df[df['TPM'] > 0]
    df['log10TPM'] = np.log10(df['TPM'])
    df.to_csv('revisions/patient_%s_scores.tsv' % sample, sep='\t', float_format='%g')
    for s in ['Source', 'NonSource']:
        outf = 'revisions/patient_%s_scores.%s.n2500.tsv' % (sample, s)
        df[df.Class == s].sample(2500).to_csv(outf, sep='\t', float_format='%g')
dat.clear_unused()
del peptides
del camap_scores
del classes
del transcripts
del df_scores
