#!/usr/bin/env python3

import pandas as pd
import os
import pickle as pkl
from functools import partial
import shutil
from collections import defaultdict
from mpi4py.futures import MPIPoolExecutor as Executor
#from concurrent.futures import ProcessPoolExecutor as Executor
import argparse
import sys
from datetime import datetime
import tarfile
import subprocess


OUTPUT_FOLDER = "./output"


def merge_netmhc_files(pfile, netmhc_out_list):
    columns_to_keep = ['Peptide', '1-log50k', 'nM', 'Rank']

    print(pfile)
    pepdict = pkl.load(open(pfile, 'rb'))

    dfs = []
    for fn in netmhc_out_list:
        with open(fn, 'r') as f:
            allele = f.readline().strip()
        temp_df = pd.read_csv(fn, sep='\t', skiprows=[0], header=0, usecols=columns_to_keep, index_col=0)
        dfs.append(temp_df)
        dfs[-1].columns = pd.MultiIndex.from_product([[allele], dfs[-1].columns])

    df = pd.concat(dfs, axis=1, sort=False, copy=False)

    edited = False
    alleles = set()
    for pep, scores in df.to_dict(orient='index').items():
        if 'netmhcpan' not in pepdict[pep]:
            d = defaultdict(dict)
            for allele, s in scores.items():
                al, t = allele
                alleles.add(al)
                d[al][t] = s
            d = dict(d)
            pepdict[pep]['netmhcpan'] = d
            edited = True

    if edited:
        keep_a_copy(pfile)
        pkl.dump(pepdict, open(os.path.join(pfile), 'wb'))
    else:
        print('$s: Nothing changed, keeping everything as is.' % pfile)

    return alleles


def keep_a_copy(pfile):
    base_dir, file_name = pfile.rsplit('/', 1)
    os.makedirs(os.path.join(base_dir, 'backup'), exist_ok=True)
    c = 0
    while True:
        f1 = pfile
        f2 = os.path.join(base_dir, 'backup', file_name.replace('.pkl', '.%d.pkl' % c))
        try:
            if os.path.isfile(f2):
                raise OSError
            print('%s: Moving %s to %s and replacing with updated annotations.' % (pfile, f1, f2))
            shutil.move(f1, f2)
            subprocess.call(['gzip', f2], shell=False)
            break
        except OSError:
            c += 1


def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-w", "--workers", type=int, default=1, help="number of parallel workers in addition to main")
    parser.add_argument("-g", "--genome", help="genome [GRCh37.75, GRCm38.78, etc.]", type=str, default="GRCh37.75")

    args=parser.parse_args().__dict__

    workers = args['workers']
    genome = args['genome']

    sys.stdout.write('Genome: %s\n' % genome)

    netmhc_folder = 'NetMHCpan-4.0a'
    out_dir = os.path.join(OUTPUT_FOLDER, 'allPeptides_%s' % genome)
    nmpan_out = os.path.join(out_dir, netmhc_folder)
    nmpan_pred = os.path.join(nmpan_out, 'predictions')

    if not os.path.exists(nmpan_out):
        if os.path.exists(nmpan_out + '.tar.gz'):
            sys.stdout.write('Predictions folder already compressed. Uncompress if you want to rerun analysis.\n')
            sys.exit()
        else:
            sys.stderr.write('ERROR: Predictions folder not found.\n')
            sys.exit(1)

    job_dct = defaultdict(list)
    for fn_allele in os.listdir(nmpan_pred):
        pinit, plen = fn_allele.split('.')[:2]
        pfile = 'peptides_%s%s.pkl' % (pinit, plen)
        job_dct[os.path.join(out_dir, pfile)].append(os.path.join(nmpan_pred, fn_allele))

    if workers:
        with Executor(max_workers=workers) as ex:
            allele_sets = ex.map(merge_netmhc_files, *zip(*job_dct.items()))
    else:
        allele_sets = map(merge_netmhc_files, *zip(*job_dct.items()))

    alleles = set()
    for al in allele_sets:
        alleles.update(al)

    pfile = os.path.join(out_dir, 'info.pkl')
    info = pkl.load(open(pfile, 'rb'))

    edited = False
    if 'allelesNetMHC' not in info:
        info['NetMHCVersion'] = set([netmhc_folder])
        edited = True
    if 'allelesNetMHC' not in info:
        info['allelesNetMHC'] = alleles
        edited = True

    if edited:
        info['date'] = datetime.now()
        keep_a_copy(pfile)
        pkl.dump(info, open(pfile, 'wb'))

    with tarfile.open(nmpan_out + '.tar.gz', 'w:gz') as tar:
        tar.add(nmpan_out, arcname=os.path.basename(nmpan_out))
    shutil.rmtree(nmpan_out)


if __name__ == '__main__':
    main()
