#!/usr/bin/env python3

import argparse

from camaptools.Peptides import Peptides
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor as EPPE, EnhancedMPIPoolExecutor as EMPE


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-w", "--workers", type=int, default=1, help="number of parallel workers in addition to main")
    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh37.75")
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')
    parser.add_argument("--force", help="Overwrite previous NetMHC score values", action='store_true')

    args=parser.parse_args().__dict__

    workers = args['workers']
    genome = args['genome']

    mpi = args['mpi']
    overwrite = args['force']

    print('Genome:', genome)
    print('Workers: ', workers)
    print('MPI: ', mpi)
    print('Overwrite: ', overwrite)

    peptides = Peptides(genome, workers=workers, executor=EMPE if mpi else EPPE)
    peptides.merge_netmhc(overwrite=overwrite)


if __name__ == '__main__':
    main()
