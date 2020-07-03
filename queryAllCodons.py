#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import os
import math
import sys
import functools
import argparse


def get_all_genes(genome):
    from pyGeno.Genome import Genome
    from pyGeno.Gene import Gene
    #return [g.id for g, _ in zip(Genome(name=genome).iterGet(Gene), range(100))]  # use when debugging or testing
    return [g.id for g in Genome(name=genome).iterGet(Gene)]


def cod_in_tr(genome, gene_list):
    from pyGeno.Genome import Genome
    from pyGeno.Gene import Gene
    from pyGeno.Protein import Protein
    import pyGeno.tools.UsefulFunctions as uf
    ref = Genome(name=genome)

    def get_codons(gid):
        seen_positions = set()
        # go top down from protein to transcripts to ensure coding potential
        for protein in ref.get(Gene, id=gid)[0].iterGet(Protein):
            transcript = protein.transcript
            cds_start = len(transcript.UTR5)
            cds_end = len(transcript.UTR5) + len(transcript.cDNA)
            codpos = [x for e in transcript.exons for x in range(e.start, e.end)][cds_start:cds_end]
            tseq = transcript.cDNA
            pseq = transcript.protein.sequence
            try:
                assert '*' not in pseq[:-1] #-> see queryAllProteome to know why it is commented
            except AssertionError:
                sys.stdout.write('Misplaced stop codon in %s\n' % transcript.id)
            assert len(tseq) == len(codpos)
            while len(tseq) > 0:
                cod = tseq[:3]
                positions = tuple(codpos[:3])
                if len(cod) == 3:
                    if positions not in seen_positions:
                        yield cod
                        seen_positions.add(positions)
                tseq = tseq[3:]
                codpos = codpos[3:]

    codict = defaultdict(int)
    for gene_id in gene_list:
        for codon in get_codons(gene_id):
            codict[codon] += 1
    codict = dict(codict)  # lambda: defaultdict -> PicklingError: Can't pickle <type 'function'>

    return codict


def chunk_up(list_to_divide, div):
    n = max(1, int(math.ceil(len(list_to_divide)/float(div))))
    for i in range(0, len(list_to_divide), n):
        yield list_to_divide[i:i+n]


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-w", "--workers", type=int, default=1, help="number of parallel workers in addition to main")
    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh37.75")

    args=parser.parse_args().__dict__

    workers = args['workers']
    genome = args['genome']

    with ProcessPoolExecutor(max_workers=1) as ex:
        # run in separate process to avoid loading pyGeno in main
        all_genes = ex.submit(get_all_genes, genome).result()
    sys.stdout.write("Gene count = %d\n" % len(all_genes))

    ex = ProcessPoolExecutor(max_workers=workers)
    futures = []
    CompletedTasks = 0
    def augment_done(size, future):
        nonlocal CompletedTasks
        CompletedTasks += size
    n_jobs = 200
    n_jobs_real = 0
    for chunk in chunk_up(all_genes, n_jobs):
        n_jobs_real += 1
        future = ex.submit(cod_in_tr, genome, chunk)
        future.add_done_callback(functools.partial(augment_done, len(chunk)))
        futures.append(future)

    codon_dict = defaultdict(int)
    completed = 0
    checkpoint = 0
    for future in as_completed(futures):
        completed += 1
        for codon, count in future.result().items():
            codon_dict[codon] += count
        pcomp = completed/float(n_jobs_real)*100
        pdone = CompletedTasks/float(len(all_genes))*100
        sys.stderr.write('Added to dictionary: %.1f %%  Completed: %.1f %%\n' % (pcomp, pdone))

    ex.shutdown()

    codon_dict = dict(codon_dict)

    if not os.path.exists('output'):
        os.makedirs('output')
    if not os.path.exists('output/codonCounts'):
        os.makedirs('output/codonCounts')

    outfile = 'output/codonCounts/%s.tsv' % genome
    with open(outfile, 'w') as f:
        f.write('\n'.join(['%s\t%d' % (codon, count) for codon, count in codon_dict.items()]) + '\n')

if __name__ == '__main__':
    main()
