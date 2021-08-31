#!/usr/bin/env python3

import argparse
import os

from camaptools.Dataset import RegressionDataset
from camaptools.Regression import RegressionManager
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor


def get_blcl_peptides():
    train = []
    test = []

    with open('../CAMAP/peptides_dataset_162/peptides.source.train.txt', 'r') as f:
        train.extend([x.strip() for x in f.readlines() if x.strip()])
    with open('../CAMAP/peptides_dataset_162/peptides.nonsource.train.txt', 'r') as f:
        train.extend([x.strip() for x in f.readlines() if x.strip()])
    with open('../CAMAP/peptides_dataset_162/peptides.source.validation.txt', 'r') as f:
        train.extend([x.strip() for x in f.readlines() if x.strip()])
    with open('../CAMAP/peptides_dataset_162/peptides.nonsource.validation.txt', 'r') as f:
        train.extend([x.strip() for x in f.readlines() if x.strip()])
    with open('../CAMAP/peptides_dataset_162/peptides.source.test.txt', 'r') as f:
        test.extend([x.strip() for x in f.readlines() if x.strip()])
    with open('../CAMAP/peptides_dataset_162/peptides.nonsource.test.txt', 'r') as f:
        test.extend([x.strip() for x in f.readlines() if x.strip()])

    return train, test


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-d", "--datasets", type=str,
        default="GRCh37.75_BLCL,GRCh37.75_B721,GRCh37.75_PBMCs,GRCm38.78_EL4,GRCm38.78_CT26",
        help="dataset [GRCh37.75_BLCL, GRCm38.78_EL4, etc.]"
    )
    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-r", "--bs_or_rank", help="max binding score (>=100) or rank (<100)", type=int, default=1250)
    parser.add_argument("-n", "--np", help="use rank score from naturally processed peptides", action='store_true')
    parser.add_argument("-m", "--ncontexts", help="max contexts permitted to keep peptide", type=int, default=10)
    parser.add_argument("-a", "--ann_method", help="ANN method to use", type=str, default='SGD')
    parser.add_argument("-p", "--ann_parameters", help="ANN training parameters selection", type=str, default='e4000')
    parser.add_argument("-x", "--max_replicates", help="ANN replicate runs to keep", type=int, default=None)
    parser.add_argument("-s", "--classifiers", type=str, default='mlp',
        help="classifiers to use for predictions, possible values are: mlp,logregcv,logreg,svc,linsvc,rf,knn [mlp]")
    parser.add_argument("-w", "--workers", help="number of parallel workers in addition to main", type=int, default=0)
    parser.add_argument("-f", "--folder", help="folder to dump trainers in", type=str, required=True)
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')
    parser.add_argument("--test", help="run in test mode", action='store_true')

    args=parser.parse_args().__dict__

    genome_datasets = [tuple(x.split('_')) for x in args['datasets'].split(',')]
    context = args['context']
    max_bs_or_rank = args['bs_or_rank']
    var = 'Rank_NP' if args['np'] else 'nM'
    var = 'Rank' if max_bs_or_rank < 100 and var != 'Rank_NP' else var
    max_contexts = args['ncontexts']
    method = args['ann_method']
    params = args['ann_parameters']
    reps = args['max_replicates']
    classifiers = args['classifiers'].split(',')
    workers = args['workers']
    folder = args['folder']
    mpi = args['mpi']
    test_run = args['test']

    genomes = sorted(list(set([x[0] for x in genome_datasets])))
    raw_datasets = [x[1] for x in genome_datasets]

    print('Genomes:', genomes)
    print('Datasets:', raw_datasets)
    print('Context:', context)
    print('Max BS/Rank:', max_bs_or_rank)
    print('BS or Rank detected:', var)
    print('Max contexts per peptide:', max_contexts)
    print('ANN method:', method)
    print('ANN parameters:', params)
    print('Max CAMAP replicates:', reps)
    print('Classifiers:', classifiers)
    print('Workers:', workers)
    #print('Subworkers:', subworkers)
    print('Folder:', folder)
    print('MPI:', mpi)
    if test_run:
        print('\n========== Running in test mode ==========\n')

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

    datasets = [RegressionDataset(g, d, context, workers, Executor) for g, d in genome_datasets]

    for i, (g, d) in enumerate(genome_datasets):
        if d == 'BLCL':
            try:
                datasets[i].construct_datasets_options(train, test)
            except NameError:
                train, test = get_blcl_peptides()
                datasets[i].construct_datasets_options(train, test)

    # necessary to avoid collisions in temp_dir 
    assert len(set(dat.dataset_name for dat in datasets)) == len(datasets)

    for dat in datasets:
        dat.load_peptides_options(
                max_bs_or_rank=max_bs_or_rank,
                var=var,
                max_contexts=max_contexts,
                step='evaluateDS',
                ann_method=method,
                ann_params=params,
                max_replicates=reps
        )
        if test_run:
            dat.pepfiles = [x for x in dat.pepfiles if 'W9' in x]
            #dat.pepfiles = [x for x in dat.pepfiles if x.find('peptides_W') > -1]
        dat.load_peptides()
        dat.construct_datasets()
        dat.clear_unused()
        rem = RegressionManager(dat)
        rem.initialize_trainers(classifiers)
        rem.dump_trainers(os.path.join(folder, out_name))


if __name__ == '__main__':
    main()
