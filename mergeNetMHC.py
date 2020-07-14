#!/usr/bin/env python3

import argparse

from camaptools.MLP import Peptides


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-w", "--workers", type=int, default=1, help="number of parallel workers in addition to main")
    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh37.75")
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')
    parser.add_argument("--force", "--workers", help="Overwrite previous NetMHC score values", action='store_true')

    args=parser.parse_args().__dict__

    workers = args['workers']
    genome = args['genome']

    mpi = args['mpi']
    overwrite = args['force']

    if mpi:
        from mpi4py.futures import MPIPoolExecutor
        Executor = MPIPoolExecutor
    else:
        Executor = None # use default

    print('Genome:', genome)
    print('Workers: ', workers)
    print('MPI: ', mpi)
    print('Overwrite: ', overwrite)

    peptides = Peptides(genome)
    peptides.merge_netmhc(workers=workers, executor=Executor, overwrite=overwrite)


if __name__ == '__main__':
    main()
