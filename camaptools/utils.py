#!/usr/bin/env python3

import os
from collections import defaultdict
from glob import glob
from parse import parse

OUTPUT_FOLDER = "./output"

def available_models(target='validation-bestMin-score.pytorch', context=162, epochs=500):
    out_dir = OUTPUT_FOLDER
    model_listing = defaultdict(lambda: defaultdict(dict))
    model_files = glob(os.path.join(out_dir, 'pytorch-*', '*', target))

    for mf in model_files:
        parse_str = os.path.join(out_dir, 'pytorch-{method}', 'p{context}e{epochs}pep{params}-{n}', target)
        dct = parse(parse_str, mf)

        if context == int(dct['context']) and epochs == int(dct['epochs']):
            params = 'pep' + dct['params']
            n = int(dct['n'])

            m = dct['method']
            name = m.split('_')[0].replace('sgd', 'SGD').replace('adam', 'Adam').replace('adagrad',
                'Adagrad').replace('shuffle', 'Shuffle').replace('mask', 'Mask').replace('-', '_')

            model_listing[name][params][n] = mf

    model_listing = dict(model_listing)
    for name in model_listing.keys():
        model_listing[name] = dict(model_listing[name])
        for k in model_listing[name]:
            model_listing[name][k] = [model_listing[name][k][v] for v in sorted(model_listing[name][k])]

    return model_listing
