#!/usr/bin/env python3

import os
from collections import defaultdict
from glob import glob
from parse import parse

OUTPUT_FOLDER = "./output"

def available_models(target='validation-bestMin-score.pytorch'):
    out_dir = OUTPUT_FOLDER
    model_listing = defaultdict(lambda: defaultdict(dict))
    model_files = glob(os.path.join(out_dir, 'pytorch-*', '*', target))

    for mf in model_files:
        parse_str = os.path.join(out_dir, 'pytorch-{method}', 'p{context}e{epochs}pep{pep}t{ratio}-{n}', target)
        dct = parse(parse_str, mf)

        params = (int(dct['context']), int(dct['epochs']), dct['pep'], dct['ratio'])
        n = int(dct['n'])

        m = dct['method']
        name = m.split('_')[0].replace('sgd', 'SGD').replace('adam', 'Adam').replace('adagrad',
            'Adagrad').replace('shuffle', 'Shuffle').replace('mask', 'Mask').replace('-', '_')

        model_listing[name][params][n] = mf

    model_listing = dict(model_listing)
    for dct in model_listing.keys():
        model_listing[dct] = dict(model_listing[dct])
        for k in model_listing[dct]:
            model_listing[dct][k] = [model_listing[dct][k][v] for v in sorted(model_listing[dct][k])]

    return model_listing
