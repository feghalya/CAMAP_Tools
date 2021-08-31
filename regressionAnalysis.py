#!/usr/bin/env python3

import argparse
import os

from camaptools.Dataset import RegressionDataset
from camaptools.Regression import RegressionAggregator, RegressionManager
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor, EnhancedMPIPoolExecutor


OUTPUT_FOLDER = "./output"


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-f", "--folder", help="location of initialized trainers", type=str, required=True)
    parser.add_argument("-d", "--datasets", type=str,
        default="GRCh37.75_BLCL,GRCh37.75_B721,GRCh37.75_PBMCs,GRCm38.78_EL4,GRCm38.78_CT26",
        help="dataset [GRCh37.75_BLCL, GRCm38.78_EL4, etc.]"
    )
    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-s", "--classifiers", type=str, default='mlp',
        help="classifiers to use for predictions, possible values are: mlp,logregcv,logreg,svc,linsvc,rf,knn [mlp]")
    parser.add_argument("-w", "--workers", help="number of parallel workers in addition to main", type=int, default=0)
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')

    args=parser.parse_args().__dict__

    folder = args['folder']
    genome_datasets = [tuple(x.split('_')) for x in args['datasets'].split(',')]
    classifiers = args['classifiers'].split(',')
    context = args['context']
    workers = args['workers']
    mpi = args['mpi']

    genomes = sorted(list(set([x[0] for x in genome_datasets])))
    raw_datasets = [x[1] for x in genome_datasets]

    print('Folder:', folder)
    print('Genomes:', genomes)
    print('Datasets:', raw_datasets)
    print('Classifiers:', classifiers)
    print('Context:', context)
    print('Workers:', workers)
    print('MPI:', mpi)

    Executor = EnhancedMPIPoolExecutor if mpi else EnhancedProcessPoolExecutor

    assert len(os.listdir(folder)) == 1
    out_name = os.listdir(folder)[0]

    out_dir = os.path.join(OUTPUT_FOLDER, 'metrics', out_name)

    datasets = [RegressionDataset(g, d, context) for g, d in genome_datasets]
    regressors = [RegressionManager(dat, workers, Executor) for dat in datasets]

    for m in regressors:
        m.optimize_mem_usage = True
        m.load_trainers(os.path.join(folder, out_name), classifiers)

    rem = RegressionAggregator(regressors, out_dir, workers, Executor)

    rem.run()
    rem.join()


if __name__ == '__main__':
    main()
