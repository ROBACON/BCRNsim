#!/usr/bin/env python3

import sys, os, inspect
from pathlib import Path
import logging

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir, sim_dir = local_dir.absolute().parents[0], local_dir.absolute().parents[1]
sys.path.append(str(root_dir))
FILENAME = os.path.basename(__file__)

from pprint import pprint
import run
import json, math
from copy import deepcopy
import itertools
import numpy as np

module_logger = logging.getLogger('root')

VOLUME_ML = 1

def build(params):

    # potentially set volume
    global VOLUME_ML
    if 'volume_ml' in params.keys():
        VOLUME_ML = params['volume_ml']
        print(f'Volume has been set to {VOLUME_ML} ml.')

    events = {}
    species = {
        'A': params['init_cells'],
    }
    mappings = {
        'cell type A': ['A'],
    }


    # potentially dilute to volume less than 1ml
    for s in species.keys():
        species[s] = int( species[s] * VOLUME_ML )



    # --- reaction parameters --- #
    parameters = {
        'dup_rate':                                   (params['dup_rate'],                                  'per_min'),
    }
    reactions = {
        'res_1_duplication': {
            're': [(1, 'A')],
            'pr': [(2, 'A')],
            'kin': f"A * dup_rate"
        },
    }

    return {'species':species, 'parameters':parameters, 'reactions':reactions, 'events':events, 'mappings':mappings}


if __name__ == '__main__':
    params = parse.params(sys.argv[1])
