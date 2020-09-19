#!/usr/bin/env python3

import argparse
import os
import pickle as pkl
import gzip

from camaptools.Dataset import TrainingDataset
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor


def write_to_disk(dat, ds_name):
    if not os.path.exists('output/trainDS'):
        os.makedirs('output/trainDS')

    with gzip.open('output/trainDS/%s_encoding-CodonEmbeddings.pkl.gz' % ds_name, 'wb') as f:
        pkl.dump(dat.enc_dct, f)
    with gzip.open('output/trainDS/%s_encoding-CodonShuffleEmbeddings.pkl.gz' % ds_name, 'wb') as f:
        pkl.dump(dat.shuff_enc_dct, f)
    with gzip.open('output/trainDS/%s-Metadata.pkl.gz' % ds_name, 'wb') as f:
        pkl.dump(dat.meta_dct, f)



def encode_dataset(genome, ds, peplen, context_len, max_bs_or_rank, max_contexts, ratio, same_tpm, seeds, workers, mpi):
    Executor = EnhancedMPIPoolExecutor if mpi else EnhancedProcessPoolExecutor

    dat = TrainingDataset(genome, ds, context_len, workers=workers, executor=Executor)
    dat.pepfiles = [f for f in dat.pepfiles if '%s.pkl' % peplen in f]
    dat.load_peptides(max_bs_or_rank=max_bs_or_rank, max_contexts=max_contexts)

    for seed in seeds:
        dat.split_dataset(ratio=ratio, same_tpm=same_tpm, seed=seed)
        dat.encode_peptides(seed=seed)
        ds_name = '%s_%s_padding%d_maxBS%d_maxContexts%d_ratio%d_%s%sseed%d' % (
                ds, genome, context_len, max_bs_or_rank, max_contexts,
                ratio, 'peplen%s_' % peplen if peplen else '', 'sameTPM_' if same_tpm else '', seed)
        write_to_disk(dat, ds_name)


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh38.98")
    parser.add_argument("-d", "--dataset", help="dataset [BLCL, EL4, etc.]", type=str, default="BLCL")
    parser.add_argument("-p", "--peplen", help="peptide length", type=str, default="")
    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-r", "--bs_or_rank", help="max binding score or rank (if < 25)", type=int, default=1250)
    parser.add_argument("-m", "--ncontexts", help="max contexts permitted to keep peptide", type=int, default=3)
    parser.add_argument("-t", "--ratio", help="non-source to source ratio", type=int, default=5)
    parser.add_argument("-y", "--sametpm", help="sample non-source keeping TPM proportions", action='store_true')
    parser.add_argument("-s", "--nsplits", help="number of random splittings", type=int, default=1)
    parser.add_argument("-w", "--workers", help="number of parallel workers in addition to main", type=int, default=0)
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')

    args=parser.parse_args().__dict__

    genome = args['genome']
    ds = args['dataset']
    peplen = args['peplen']
    context_len = args['context']
    max_bs_or_rank = args['bs_or_rank']
    max_contexts = args['ncontexts']
    ratio = args['ratio']
    same_tpm = args['sametpm']
    seeds = [int(s) for s in range(1, args['nsplits']+1)]

    workers = args['workers']
    mpi = args['mpi']

    encode_dataset(genome, ds, peplen, context_len, max_bs_or_rank, max_contexts, ratio, same_tpm, seeds, workers, mpi)


if __name__ == '__main__':
    main()
