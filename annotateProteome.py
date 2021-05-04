#!/usr/bin/env python3

import argparse

from camaptools.Peptides import Peptides
from camaptools.EnhancedFutures import EnhancedProcessPoolExecutor as EPPE, EnhancedMPIPoolExecutor as EMPE


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-w", "--workers", help="Number of parallel workers in addition to main", type=int, default=0)
    parser.add_argument("-g", "--genome", help="Genome [GRCh37.75; GRCm38.78; etc.]", type=str, default="GRCh37.75")
    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-a", "--algorithms", help="Whitelist of algorithms (empty for all)", type=str, default="")
    parser.add_argument("-p", "--parameters", help="Whitelist of parameters (empty for all)", type=str, default="")
    parser.add_argument("--mpi", help="Parallelize using MPI", action='store_true')
    parser.add_argument("--force", help="Overwrite previous CAMAP score values", action='store_true')

    args=parser.parse_args().__dict__

    workers = args['workers']
    genome = args['genome']
    context = args['context']
    algorithms = set(args['algorithms'].split(',')) if args['algorithms'] else None
    parameters = set(args['parameters'].split(',')) if args['parameters'] else None

    mpi = args['mpi']
    overwrite = args['force']

    print('Genome:', genome)
    print('Context:', context)
    print('Algorithms:', algorithms)
    print('Parameters:', parameters)
    print('Workers: ', workers)
    print('MPI: ', mpi)
    print('Overwrite: ', overwrite)

    peptides = Peptides(genome, context=context, workers=workers, executor=EMPE if mpi else EPPE)
    #peptides.pepfiles = [x for x in peptides.pepfiles if 'W8' in x]
    peptides.annotate(overwrite=overwrite, algorithms=algorithms, parameters=parameters)


if __name__ == '__main__':
    main()
