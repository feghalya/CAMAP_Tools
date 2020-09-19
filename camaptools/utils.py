#!/usr/bin/env python3

import os
from collections import defaultdict
from glob import glob
from parse import parse

OUTPUT_FOLDER = "./output"

def available_models(target='validation-bestMin-score.pytorch', context=162, filters=None):
    out_dir = OUTPUT_FOLDER
    model_listing = defaultdict(lambda: defaultdict(dict))
    model_files = glob(os.path.join(out_dir, 'pytorch-*', '*', target))

    for mf in model_files:
        parse_str = os.path.join(out_dir, 'pytorch-{method}', 'p{context:d}{params}-{n}', target)
        dct = parse(parse_str, mf)

        params = dct['params']
        n = int(dct['n'])
        m = dct['method']
        c = int(dct['context'])

        if context == c:
            name = m.split('_')[0].replace('sgd', 'SGD').replace('adam', 'Adam').replace('adagrad',
                'Adagrad').replace('shuffle', 'Shuffle').replace('mask', 'Mask').replace('-', '_')
            if filters is None or name in filters:
                model_listing[name][params][n] = mf

    for name in model_listing.keys():
        for k in model_listing[name]:
            model_listing[name][k] = [model_listing[name][k][v] for v in sorted(model_listing[name][k])]
        model_listing[name] = dict(model_listing[name])
    model_listing = dict(model_listing)

    return model_listing
