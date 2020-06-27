#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
import pickle as pickle
from collections import defaultdict
import os
import math
import sys
import pkg_resources
from datetime import datetime
import functools
import argparse


def get_all_genes(genome):
    from pyGeno.Genome import Genome
    from pyGeno.Gene import Gene
    #return [g.id for g, _ in zip(Genome(name=genome).iterGet(Gene), range(100))]  # use when debugging or testing
    return [g.id for g in Genome(name=genome).iterGet(Gene)]


def pep_in_prot(genome, context, gene_list):
    from pyGeno.Genome import Genome
    from pyGeno.Gene import Gene
    from pyGeno.Protein import Protein
    import pyGeno.tools.UsefulFunctions as uf
    ref = Genome(name=genome)

    def get_peptides(gid):
        for protein in ref.get(Gene, id=gid)[0].iterGet(Protein):
            pseq = protein.sequence
            tseq = protein.transcript.cDNA
            offset = 0
            accepted_pep_lengths = [8, 9, 10, 11]
            # use pseq to avoid errors caused by end stop codons in tseq
            while len(pseq) >= 116:
                pep_lengths = [x for x in accepted_pep_lengths if len(pseq) >= int(context/3)*2+x]

                contexts_pre_aa = [pseq[:int(context/3)]]*len(pep_lengths)
                contexts_pre = [tseq[:context]]*len(pep_lengths)

                peptides_aa = [pseq[int(context/3):int(context/3)+x] for x in pep_lengths]
                peptides_dna = [tseq[context:context+x*3] for x in pep_lengths]

                contexts_post_aa = [pseq[int(context/3)+x:int(context/3)+x+int(context/3)] for x in pep_lengths]
                contexts_post = [tseq[context+x*3:context+x*3+context] for x in pep_lengths]

                for cpre_aa, cpre, pep_aa, \
                    pep_dna, cpost_aa, cpost in zip(contexts_pre_aa, contexts_pre,
                                                    peptides_aa, peptides_dna,
                                                    contexts_post_aa, contexts_post):
                    context_seq = cpre+pep_dna+cpost
                    context_seq_aa = cpre_aa+pep_aa+cpost_aa
                    if not ('*' in context_seq_aa):
                        # Fails for mitochondrial genes -- rest is OK
                        #try:
                        #    assert uf.translateDNA(cpre) == cpreaa
                        #    assert uf.translateDNA(d) == p
                        #    assert uf.translateDNA(cpost) == cpostaa
                        #except AssertionError:
                        #    print cpre, cpreaa
                        #    print d, p
                        #    print cpost, cpostaa
                        #    raise AssertionError
                        yield (pep_aa, context_seq, offset+context,
                               protein.id, protein.name,
                               protein.transcript.id, protein.transcript.name,
                               protein.gene.id, protein.gene.name)
                offset += 3
                tseq = tseq[3:]
                pseq = pseq[1:]

    pepdict = defaultdict(lambda: defaultdict(list))
    for gene_id in gene_list:
        seen_contexts = set()
        for peptide in get_peptides(gene_id):
            pep, seq, pos, pid, pname, tid, tname, gid, gname = peptide
            # if by some uncanny chance, there are 2 exactly similar contexts at different
            # positions in a gene, it's considered a duplicate as well
            if seq not in seen_contexts:
                seen_contexts.add(seq)
                pepdict[pep][gid].append({
                    'peptide': pep,
                    'sequenceContext': seq,
                    'geneID': gid,
                    'geneName': gname,
                    'proteinID': [pid],
                    'transcriptID': [tid],
                    'positionInTranscript': [pos]})
            else:
                for d in pepdict[pep][gid]:
                    if d['sequenceContext'] == seq:
                        d['proteinID'].append(pid)
                        d['transcriptID'].append(tid)
                        d['positionInTranscript'].append(pos)
    pepdict = dict(pepdict)  # lambda: defaultdict -> PicklingError: Can't pickle <type 'function'>
    return pepdict


def chunk_up(list_to_divide, div):
    n = max(1, int(math.ceil(len(list_to_divide)/float(div))))
    for i in range(0, len(list_to_divide), n):
        yield list_to_divide[i:i+n]


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-c", "--context", help="mRNA context length on each side", type=int, default=162)
    parser.add_argument("-w", "--workers", type=int, default=1, help="number of parallel workers in addition to main")
    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh37.75")

    args=parser.parse_args().__dict__

    context = args['context']
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
        future = ex.submit(pep_in_prot, genome, context, chunk)
        future.add_done_callback(functools.partial(augment_done, len(chunk)))
        futures.append(future)

    peptides = defaultdict(lambda: defaultdict(dict))
    completed = 0
    checkpoint = 0
    for future in as_completed(futures):
        completed += 1
        for pep, genes_dict in future.result().items():
            k = pep[0] + str(len(pep))
            peptides[k][pep].update(genes_dict)
        pcomp = completed/float(n_jobs_real)*100
        pdone = CompletedTasks/float(len(all_genes))*100
        sys.stderr.write('Added to dictionary: %.1f %%  Completed: %.1f %%\n' % (pcomp, pdone))

    ex.shutdown()

    if not os.path.exists('output'):
        os.makedirs('output')
    if not os.path.exists('output/allPeptides_%s' % genome):
        os.makedirs('output/allPeptides_%s' % genome)

    info = {'contextsize': context,
            'genomeVersion': genome,
            'pyGenoVersion': pkg_resources.get_distribution("pyGeno").version,
            'date': datetime.now()}
    outfile = 'output/allPeptides_%s/info.pkl' % genome
    pickle.dump(info, open(outfile, 'wb'))

    peplist = []

    for key, di in peptides.items():
        peplist.extend(list(di.keys()))
        outfile = 'output/allPeptides_%s/peptides_%s.pkl' % (genome, key)
        pickle.dump(di, open(outfile, 'wb'))

    outfile = 'output/allPeptides_%s/proteome.wide.peptides.txt' % genome
    with open(outfile, 'w') as f:
        f.write('\n'.join(peplist) + '\n')

if __name__ == '__main__':
    main()
