#!/usr/bin/env python3

#####################################################
#                                                   #
#    E.coli Donor + Recipient Antibiotic model      #
#                                                   #
#####################################################

import sys, os, inspect
from pathlib import Path
import logging

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir, sim_dir = local_dir.absolute().parents[0], local_dir.absolute().parents[1]
sys.path.append(str(root_dir))
FILENAME = os.path.basename(__file__)

import run
import json, math
from copy import deepcopy

module_logger = logging.getLogger('root')

import Donor
import Recipient_Kan

from pprint import pprint

GLOBAL_SPECIES = [
    'G',
    'AA',
    'P',
    'ADonor',
]

def __get_params_for_model(module_name: str, params):
    new_params = {}
    for p in params.keys():
        if '.' in p:
            # is a hierarchical parameter
            splitname = p.split('.')
            p_module = splitname[0]
            p_name = ".".join(splitname[1:])

            if p_module == module_name or p_module == '*':
                # is intended for this module
                new_params[p_name] = params[p]
            else:
                # intended for other module
                pass
        else:
            # not hierarchical parameter -> pass on
            new_params[p] = params[p]
    return new_params

def __rename(s, module_name, global_species=GLOBAL_SPECIES):
    if s in global_species:
        return s
    else:
        return f'{module_name}_{s}'

def __rename_list(s_list, module_name, global_species=GLOBAL_SPECIES):
    return [__rename(s, module_name, global_species) for s in s_list]

def __rename_repr_list(s_list, module_name, global_species=GLOBAL_SPECIES):
    return [(s[0], __rename(s[1], module_name, global_species)) for s in s_list]


def __rename_kinetics(kin, reactants, module_name, global_species=GLOBAL_SPECIES):
    kin_list = kin.split(' ')
    kin_list_new = []
    for term in kin_list:
        term_strip = term.strip(' ')
        reactant_species = [s[1] for s in reactants]
        if term_strip in reactant_species:
            # rename
            kin_list_new += [ __rename(term_strip, module_name, global_species) ]
        else:
            kin_list_new += [ term_strip ]
    return  " ".join(kin_list_new)


def __rename_module(build_values, module_name):
    species_new = {}
    for s in build_values['species'].keys():
        species_new[__rename(s, module_name)] = build_values['species'][s]

    parameters_new = build_values['parameters']

    reactions_new = {}
    reactions = build_values['reactions']
    for rk in reactions.keys():
        reactions_new[f'{module_name}_{rk}'] = {
            're': __rename_repr_list(reactions[rk]['re'], module_name),
            'pr': __rename_repr_list(reactions[rk]['pr'], module_name),
            'kin': __rename_kinetics(reactions[rk]['kin'], reactants=reactions[rk]['re'], module_name=module_name),
            }

    # fixme
    events_new = build_values['events']

    mappings_new = {}
    for mk in build_values['mappings'].keys():
        mapping_rule = build_values['mappings'][mk]
        if isinstance(mapping_rule, str):
            # not supperted right now
            module_logger.info(f'__rename_module: skipping mapping {mk} with formula : {mapping_rule}')
        else:
            mappings_new[__rename(mk, module_name)] = __rename_list(build_values['mappings'][mk], module_name)

    return { 'species': species_new,
        'parameters': parameters_new,
        'reactions': reactions_new,
        'events': events_new,
        'mappings': mappings_new }


def __merge_modules(params, params_to_merge):
    species_new = params['species']
    for s in params_to_merge['species'].keys():
        if s in params['species'].keys():
            if params['species'][s] == params_to_merge['species'][s]:
                module_logger.info(f'Merged: species {s}')
            else:
                module_logger.info(f'Merge conflict: species {s}')
        else:
            # merge
            species_new[s] = params_to_merge['species'][s]

    parameters_new = params['parameters']
    for s in params_to_merge['parameters'].keys():
        if s in params['parameters'].keys():
            if params['parameters'][s] == params_to_merge['parameters'][s]:
                module_logger.info(f'Merged: parameter {s}')
            else:
                module_logger.info(f'Merge conflict: parameter {s}')
        else:
            # merge
            parameters_new[s] = params_to_merge['parameters'][s]

    reactions_new = params['reactions']
    for s in params_to_merge['reactions'].keys():
        if s in params['reactions'].keys():
            module_logger.info(f'Merge conflict: reaction {s}')
        else:
            # merge
            reactions_new[s] = params_to_merge['reactions'][s]

    events_new = params['events']
    for s in params_to_merge['events'].keys():
        if s in params['events'].keys():
            module_logger.info(f'Merge conflict: event {s}')
        else:
            # merge
            events_new[s] = params_to_merge['events'][s]

    mappings_new = params['mappings']
    for s in params_to_merge['mappings'].keys():
        if s in params['mappings'].keys():
            module_logger.info(f'Merge conflict: mapping {s}')
        else:
            # merge
            mappings_new[s] = params_to_merge['mappings'][s]

    return { 'species': species_new,
             'parameters': parameters_new,
             'reactions': reactions_new,
             'events': events_new,
             'mappings': mappings_new }

def build(params):
    ret_donor = Donor.build(__get_params_for_model(module_name='Donor', params=params))
    ret_donor_new = __rename_module(ret_donor, module_name='Donor')

    ret_receiver = Recipient_Kan.build(__get_params_for_model(module_name='Recipient', params=params))
    ret_receiver_new = __rename_module(ret_receiver, module_name='Recipient')

    # get volume
    VOLUME_ML = 1.0
    if 'volume_ml' in params.keys():
        VOLUME_ML = params['volume_ml']
        print(f'DonorRecipient: volume has been set to {VOLUME_ML} ml.')


    # merge both
    ret = __merge_modules(ret_donor_new, ret_receiver_new)

    # set special inits for antibiotics
    ret['species']['ADonor']                = 0
    ret['species']['Recipient_APhage']      = 0
    ret['species']['Recipient_AReceiver']   = 0

    # set special mappings
    ret['mappings']['Phages']   = ['P']
    ret['mappings']['Ecoli']    = '[Donor_Ecoli] + [Recipient_Ecoli]'
    ret['mappings']['OD']       = '0.3*(0.01 * [Donor_dead] + [Donor_Ecoli]) + (0.01 * [Recipient_dead] + [Recipient_Ecoli])'

    # set special events
    ret['events']['add_antibiotics'] = { 'trigger': 'true',
                                         'delay': str(params['antibiotics_after_min']),
                                         'assignments': [('ADonor',              f"{params['init_antibiotics#ADonor']} * {VOLUME_ML}"  ),
                                                         ('Recipient_APhage',    f"{params['Recipient.init_antibiotics#APhage']} * {VOLUME_ML}"),
                                                         ('Recipient_AReceiver', f"{params['Recipient.init_antibiotics#AReceiver']} * {VOLUME_ML}"),
                                                        ]
                                        }
    #pprint(ret)
    return ret

if __name__ == '__main__':
    params = parse.params(sys.argv[1])
