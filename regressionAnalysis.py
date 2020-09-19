#!/usr/bin/env python3

import argparse
import os

from camaptools.Dataset import RegressionDataset
from camaptools.Regression import RegressionMetaManager
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor


OUTPUT_FOLDER = "./output"


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-g", "--genomes", help="genomes [GRCh37.75, etc.]", type=str, default="GRCh37.75,GRCm38.78")
    #parser.add_argument("-d", "--dataset", help="dataset [BLCL, EL4, etc.]", type=str, default="BLCL")
    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-r", "--bs_or_rank", help="max binding score (>=100) or rank (<100)", type=int, default=1250)
    parser.add_argument("-n", "--np", help="use rank score from naturally processed peptudes", action='store_true')
    parser.add_argument("-m", "--ncontexts", help="max contexts permitted to keep peptide", type=int, default=10)
    parser.add_argument("-a", "--ann_method", help="ANN method to use", type=str, default='SGD')
    parser.add_argument("-p", "--ann_parameters", help="ANN training parameters selection", type=str, default='e4000')
    parser.add_argument("-w", "--workers", help="number of parallel workers in addition to main", type=int, default=0)
    #parser.add_argument("-s", "--subworkers", help="number of parallel subworkers", type=int, default=0)
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')

    args=parser.parse_args().__dict__

    genomes = args['genomes'].split(',')
    #ds = args['dataset']
    context = args['context']
    max_bs_or_rank = args['bs_or_rank']
    var = 'Rank_NP' if args['np'] else 'nM'
    var = 'Rank' if max_bs_or_rank < 100 and var != 'Rank_NP' else var
    max_contexts = args['ncontexts']
    method = args['ann_method']
    params = args['ann_parameters']
    workers = args['workers']
    #subworkers = args['subworkers']
    mpi = args['mpi']

    print('Genomes:', genomes)
    #print('Dataset:', ds)
    print('Context:', context)
    print('Max BS/Rank:', max_bs_or_rank)
    print('BS or Rank detected:', var)
    print('Max contexts per peptide:', max_contexts)
    print('ANN method:', method)
    print('ANN parameters:', params)
    print('Workers:', workers)
    #print('Subworkers:', subworkers)
    print('MPI:', mpi)

    Executor = EnhancedMPIPoolExecutor if mpi else EnhancedProcessPoolExecutor

    out_name = 'genome%s_padding%d_max%s%d_maxContexts%d_ANNMethod%s%s' % (
            'genome'.join(genomes),
            context,
            var.replace('nM', 'BS').replace('_', ''),
            max_bs_or_rank,
            max_contexts,
            method,
            params
        )

    out_dir = os.path.join(OUTPUT_FOLDER, 'metrics', out_name)

    genome_datasets = []
    for genome in genomes:
        if 'GRCh' in genome:
            for ds in ['BLCL', 'B721', 'PBMCs']:
                genome_datasets.append((genome, ds))
        elif 'GRCm' in genome:
            for ds in ['EL4', 'CT26']:
                genome_datasets.append((genome, ds))

    datasets = [RegressionDataset(g, d, context) for g, d in genome_datasets]

    rem = RegressionMetaManager(datasets, out_dir, workers, Executor)
    rem.set_load_peptides_options(
            max_bs_or_rank=max_bs_or_rank,
            var=var,
            max_contexts=max_contexts,
            step='evaluateDS',
            ann_method=method,
            ann_params=params)
    rem.run()
    rem.join()


if __name__ == '__main__':
    main()
