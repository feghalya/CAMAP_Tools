#!/usr/bin/env python3

import argparse

from camaptools.Regression import RegressionDataset, RegressionManager


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh37.75")
    parser.add_argument("-d", "--dataset", help="dataset [BLCL, EL4, etc.]", type=str, default="BLCL")
    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-w", "--workers", help="number of parallel workers in addition to main", type=int, default=0)
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')

    args=parser.parse_args().__dict__

    genome = args['genome']
    ds = args['dataset']
    context = args['context']
    workers = args['workers']
    mpi = args['mpi']

    print('Genome:', genome)
    print('Dataset:', ds)
    print('Context: ', context)
    print('Workers: ', workers)
    print('MPI: ', mpi)

    dat = RegressionDataset(genome, ds, context)
    dat.pepfiles = [x for x in dat.pepfiles if 'W10' in x]
    dat.load_peptides(step='evaluateDS')
    dat.organize_data()
    dat.construct_datasets()

    rem = RegressionManager(dat, workers, mpi)
    rem.initialize_trainers()
    rem.start()
    rem.join()


if __name__ == '__main__':
    main()
