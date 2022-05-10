#!/usr/bin/env python3

#####################################################
#                                                   #
#    E.coli with F-pili and wild-type M13           #
#                                                   #
#####################################################

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

module_logger = logging.getLogger('root')

# phages
def __pili(ecoli: str):
    """
    Returns: number of taken pili
    """
    parts = ecoli.split('_p')
    return int(parts[1].replace('_',''))

def __infectable(ecoli: str, params):
    """
    Returns: if it can be infected
    Attention: if not the case, may still be able to uptake phages on
               free pili.
    """
    return (not __infected(ecoli)) and (__pili(ecoli) < params['num_pili']) and (not __is_blocked(ecoli))

def __infected_and_can_bind_phages(ecoli: str, params):
    """
    Returns: if can bind phages
    """
    return __infected(ecoli)  # and not __is_blocked(ecoli)
#    return __infected(ecoli) and (__pili(ecoli) < params['num_pili'])

def __infected(ecoli: str):
    """
    Returns: if infected
    """
    return "_ni_" not in ecoli

def __in_early_infection_phase(ecoli: str):
    """
    Returns: if in early phase
    """
    return "_ei_" in ecoli

def __in_late_infection_phase(ecoli: str):
    """
    Returns: if in early phase
    """
    return "_li_" in ecoli


def __is_old(ecoli: str):
    """
    Returns: if old (here: has the F machinery with higher infection rate)
    """
    return "_ageo_" in ecoli


def __is_blocked(ecoli: str):
    """
    if blocked by gp3
    """
    return "_blocked" in ecoli


def __ecoli_birthdeathrates(ecoli: str, params: dict):
    """
    In: ecoli as str
    Returns: birth rates (per resource), and
             death rate
    """
    rates = {
        'birth': {
            'G': None,
            'AA': None
        },
        'death': None
    }
    # birth
    rates['birth']['G']  = f"dup_rate / {params['init_glucose']}"
    rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor / {params['init_amino_acid']}"

    # birth: plasmid load
    if __in_early_infection_phase(ecoli):
        rates['birth']['G'] += f" * {params['plasmid_penalty_factor']}"
        rates['birth']['AA'] += f" * {params['plasmid_penalty_factor']}"

    elif __in_late_infection_phase(ecoli):
        rates['birth']['G'] += f" * {params['plasmid_penalty_factor_late']}"
        rates['birth']['AA'] += f" * {params['plasmid_penalty_factor_late']}"

    # death
    rates['death'] = 'death_rate'

    return rates


def build(params):

    ATTACHED_OUTFLOW_FACTOR = params['ATTACHED_OUTFLOW_FACTOR']

    if params['debug']:
        if params['num_pili'] == 1:
            module_logger.info(f"num pili currently set to 1")
        if params['init_phages'] != 0:
            module_logger.info(f"using phages")

    events = {}
    species = {
        'G':     params['init_glucose'],
        'AA':    params['init_amino_acid'],
        'P':     params['init_phages'],
        'dead':  0,
        'gp3':   params['init_gp3'],
        'other': params['init_other'],
    }

    # --- species & mappings --- #

    # Various E. coli:
    # keep each __identifier__ unique, since the code uses string replacement for reactions
    mappings = {
        'Resources':['G','AA'],
        'Phages':['P'],
        'Ecoli':[],
        'Ecoli uninfected':[],
        'Ecoli infected':[],
        # infection
        'Ecoli early infection':[],
        'Ecoli late infection':[],
        # gp3
        'Ecoli blocked':[],
        'Ecoli nonblocked':[],
        # age
        'Ecoli young':[],
        'Ecoli old':[],
        'OD': '0.01 * [dead] + [Ecoli]',
        # sums
        'allphages': '[P] + [Ecoli infected]',
    }
    for i in range(params['num_pili']+1):
        mappings[f'Ecoli {i} pili'] = []

    for infct in ['ni','ei','li']:
        for i in range(params['num_pili']+1):
            pili_name = f'p{i}'
            for gp3 in ['blocked', 'nonblocked']:
                for age in ['agey', 'ageo']: # young and old
                    s = f'E_{infct}_{gp3}_{age}_{pili_name}'
                    # init species
                    init_species_young =    f'E_ni_nonblocked_agey_p0'
                    init_species_old =      f'E_ni_nonblocked_ageo_p0'
                    if s == init_species_young:
                        species[s] = int(params['init_ecoli'] * (1-params['init_ratio_old']))
                    elif s == init_species_old:
                        species[s] = int(params['init_ecoli'] * params['init_ratio_old'])
                    else:
                        species[s] = 0
                    # mapping
                    mappings['Ecoli'] += [s]
                    if not __infected(s):
                        mappings['Ecoli uninfected'] += [s]
                    else:
                        mappings['Ecoli infected'] += [s]
                        if __in_early_infection_phase(s):
                            mappings['Ecoli early infection'] += [s]
                        elif __in_late_infection_phase(s):
                            mappings['Ecoli late infection'] += [s]

                    if __is_blocked(s):
                        mappings['Ecoli blocked'] += [s]
                    else:
                        mappings['Ecoli nonblocked'] += [s]

                    if __is_old(s):
                        mappings['Ecoli old'] += [s]
                    else:
                        mappings['Ecoli young'] += [s]

                    mappings[f'Ecoli {i} pili'] += [s]


    # not preinfected
    moi = 'na' if params['init_ecoli'] == 0 else params['init_phages']/params['init_ecoli']
    module_logger.info(f"MOI = {moi}")


    # --- reaction parameters --- #
    parameters = {
        'gp3_block_rate':                       (params['viscosity_factor']*params['gp3_block_rate'],   'per_min'),
        'gp3_unblock_rate':                     (params['gp3_unblock_rate'],                            'per_min'),

        'dup_rate':                             (params['dup_rate'],                                    'per_min'),
        'other_dup_rate':                       (params['other_dup_rate'],                              'per_min'),
        'amino_acid_penalty_factor':            (params['amino_acid_penalty_factor'],                   'dimensionless'),
        'death_rate':                           (params['death_rate'],                                  'per_min'),

        'inf_rate':                             (params['viscosity_factor']*params['inf_rate'],         'per_min'),
        'inf_factor_young':                     (params['inf_factor_young'],                            'dimensionless'),
        'young_to_old_rate':                    (params['young_to_old_rate'],                           'per_min'),
        'uninf_rate':                           (params['uninf_rate'],                                  'per_min'),
        'early_to_late_inf_rate':               (params['early_to_late_inf_rate'],                      'per_min'),

        'sec_rate':                             (params['sec_rate'],                                    'per_min'),
        'gp3_sec_rate':                         (params['gp3_sec_rate'],                                'per_min'),
        'outflow_rate':                         (params['outflow_rate'],                                'per_min'),
        'late_sec_penalty':                     (params['late_sec_penalty'],                            'dimensionless'),
        
        'G_inflow_rate':                        (params['G_inflow_rate'],                               'per_min'),
        'AA_inflow_rate':                       (params['AA_inflow_rate'],                              'per_min'),
        'E_ni_nonblocked_agey_p0_inflow_rate':  (params['E_ni_nonblocked_agey_p0_inflow_rate'],         'per_min'),
        'other_inflow_rate':                    (params['other_inflow_rate'],                           'per_min'),
    }
    reactions = {}

    # RESOURCE + Ecoli + other INFLOW
    for res in ['G', 'AA', 'E_ni_nonblocked_agey_p0', 'other']:
        reactions[f'resource_{res}_inflow'] = {
            're': [ ],
            'pr': [(1, res)],
            'kin': f"{res}_inflow_rate"
        }

    # PHAGE DECAY
    reactions[f'decay_P'] = {
        're': [(1, 'P')],
        'pr': [ ],
        'kin': f"P * {params['phage_decay_rate']}"
    }

    # GP3, PHAGE, RESOURCE OUTFLOW
    for thing in ['gp3', 'P', 'G', 'AA']:
        reactions[f'outflow_{thing}'] = {
                're': [(1, thing)],
                'pr': [],
                'kin': f'{thing} * outflow_rate'
                }

    # go over all species
    for s in species:

        # OTHER SPECIES
        if s == 'other':
            # other duplicate
            other_rate = {}
            other_rate['G']  = f"other_dup_rate / {params['init_glucose']}"
            other_rate['AA'] = f"other_dup_rate * amino_acid_penalty_factor / {params['init_amino_acid']}"
            for res in ['G','AA']:
                # birth per resource
                reactions[f'birth_other_{res}'] = {
                    're': [(1, 'other'), (1, res)],
                    'pr': [(1, 'other'), (1, 'other')],
                    'kin': f"other * {res} * {other_rate[res]}"
                }

            # outflow
            reactions[f'outflow_{s}'] = {
                're': [(1, s)],
                'pr': [],
                'kin': f'{ATTACHED_OUTFLOW_FACTOR} * {s} * outflow_rate'
                }

            # DEATH
            reactions[f'death_{s}'] = {
                're': [(1, s)],
                'pr': [(1, 'dead')],
                'kin': f"{s} * death_rate"
            }

        # ECOLI REACTIONS
        if 'E_' in s:
            ecoli = s
            birthdeathrates = __ecoli_birthdeathrates(ecoli, params=params)

            # BIRTH

            # pili, blocked
            parent_pili = __pili(ecoli)
            if params['new_pili_on_duplication']:
                l_child = ecoli.replace('_p'+str(parent_pili),'_p'+str(int(parent_pili/2)))
                r_child_nonblocked = ecoli.replace('_blocked','_nonblocked').replace('_p'+str(parent_pili),'_p'+str(int((parent_pili+1)/2)))
                r_child = r_child_nonblocked
            else:
                l_child = r_child = ecoli

            # age: always old -> young
            l_child = l_child.replace('_ageo','_agey')
            r_child = r_child.replace('_ageo','_agey')

            # birth per resource
            for res in ['G','AA']:
                reactions[f'birth_{ecoli}_{res}'] = {
                    're': [(1, ecoli), (1, res)],
                    'pr': [(1, l_child), (1, r_child)],
                    'kin': f"{ecoli} * {res} * {birthdeathrates['birth'][res]}"
                }

                # # also other duplicate
                # other_rate = {}
                # other_rate['G']  = f"{params['other_dup_rate']} / {params['init_glucose']}"
                # other_rate['AA'] = f"{params['other_dup_rate']} * amino_acid_penalty_factor / {params['init_amino_acid']}"
                
                # reactions[f'birth_other_{res}'] = {
                #     're': [(1, 'other'), (1, res)],
                #     'pr': [(1, 'other'), (1, 'other')],
                #     'kin': f"other * {res} * {other_rate[res]}"
                # }

            # DEATH
            reactions[f'death_{ecoli}'] = {
                're': [(1, ecoli)],
                'pr': [(1, 'dead')],
                'kin': f"{ecoli} * {birthdeathrates['death']}"
            }

            # INFECTION
            # kinetics
            if __is_old(ecoli):
                inf_kin = f'inf_rate * P * {ecoli}'
            else:
                inf_kin = f'inf_rate * inf_factor_young * P * {ecoli}'

            # clogging
            if params['clogging_infection']:
                inf_kin += f" * {1 - parent_pili / params['num_pili']}"

            # reactions        
            if __infectable(ecoli, params):
                # reaction
                reactions[f'infect_{ecoli}'] = {
                    're': [(1, ecoli), (1, 'P'), ],
                    'pr': [(1, ecoli.replace('_ni','_ei').replace(f"_p{parent_pili}",f"_p{parent_pili+1}").replace('_nr','_ar'))],
                    'kin' : inf_kin,
                }

            if __infected_and_can_bind_phages(ecoli, params=params):
                # uptake of further phages on pili
                reactions[f'already_infected_take_phage_{ecoli}'] = {
                    're': [(1, ecoli), (1, 'P'), ],
                    'pr': [(1, ecoli.replace(f'_p{parent_pili}',f'_p{min(parent_pili+1,params["num_pili"])}'))],
                    'kin' : inf_kin,
                }

            if __infected(ecoli):
                # spontaneous uninfected
                reactions[f'spont_uninf_{ecoli}'] = {
                    're': [(1, ecoli)],
                    'pr': [(1,ecoli.replace('_ei','_ni').replace('_li','_ni').replace('_ar','_nr'))],
                    'kin' : f'uninf_rate * {ecoli}',
                }

            if __in_early_infection_phase(ecoli):
                # early to late infection:
                reactions[f'early_to_late_inf_{ecoli}'] = {
                    're': [(1, ecoli)],
                    'pr': [(1, ecoli.replace('_ei_','_li_'))],
                    'kin' : f'early_to_late_inf_rate * {ecoli}',
                }


            # BLOCKING
            pili = __pili(ecoli)
            if __is_blocked(ecoli):
                reactions[f'unblocking_{ecoli}'] = {
                    're': [(1, ecoli)],
                    'pr': [(1, ecoli.replace('_blocked', '_nonblocked').replace(f'_p{pili}',f'_p{max(pili-1,0)}') )],
                    'kin' : f'gp3_unblock_rate * {ecoli}',
                }

            if not __is_blocked(ecoli):
                reactions[f'blocking_{ecoli}'] = {
                    're': [(1, ecoli), (1, 'gp3')],
                    'pr': [(1, ecoli.replace('_nonblocked','_blocked').replace(f'_p{pili}',f'_p{min(pili+1,params["num_pili"])}') )],
                    'kin' : f'gp3_block_rate * gp3 * {ecoli}',
                }

            # AGE
            if not __is_old(ecoli):
                # is young
                reactions[f'grow_up_{ecoli}'] = {
                    're': [(1, ecoli) ],
                    'pr': [(1, ecoli.replace('_agey','_ageo'))],
                    'kin' : f"{ecoli} * young_to_old_rate",
                }

            # PHAGE AND GP3 SECRETION
            if __infected(ecoli):
                if __in_early_infection_phase(ecoli):
                    for res in ['G','AA']:
                            # birth per resource
                            reactions[f'secrete_phage_{ecoli}_{res}'] = {
                                're': [(1, ecoli), (1, res)],
                                'pr': [(1, ecoli), (1, res), (1,'P'), ],
                                'kin': f"sec_rate * {ecoli} * {res} * {birthdeathrates['birth'][res]} / dup_rate",
                            }
                            # birth per resource
                            reactions[f'secrete_gp3_{ecoli}_{res}'] = {
                                're': [(1, ecoli), (1, res)],
                                'pr': [(1, ecoli), (1, res), (1,'gp3'), ],
                                'kin': f"gp3_sec_rate * {ecoli} * {res} * {birthdeathrates['birth'][res]} / dup_rate",
                            }
                else:
                    for res in ['G','AA']:
                            # birth per resource
                            reactions[f'secrete_phage_{ecoli}_{res}'] = {
                                're': [(1, ecoli), (1, res)],
                                'pr': [(1, ecoli), (1, res), (1,'P'), ],
                                'kin': f"late_sec_penalty * sec_rate * {ecoli} * {res} * {birthdeathrates['birth'][res]} / dup_rate",
                            }
                            # birth per resource
                            reactions[f'secrete_gp3_{ecoli}_{res}'] = {
                                're': [(1, ecoli), (1, res)],
                                'pr': [(1, ecoli), (1, res), (1,'gp3'), ],
                                'kin': f"late_sec_penalty * gp3_sec_rate * {ecoli} * {res} * {birthdeathrates['birth'][res]} / dup_rate",
                            }

            # ECOLI OUTFLOW
            reactions[f'outflow_{ecoli}'] = {
                    're': [(1, ecoli)],
                    'pr': [],
                    'kin': f'{ATTACHED_OUTFLOW_FACTOR} * {ecoli} * outflow_rate'
                    }

    # hack (not used):
    # events to prevent values <= eps
    # if getting negative values and >0 values
    # which result in nonterminating simultions
    if False:
        eps = 0.1
        days = 41
        for d in [60*24*i for i in range(days)]:
            for s in species.keys():
                trigger = f'{s} <= {eps}'
                assignments = [ (s, '0') ]
                events[f'eps_{s}_{d}'] = { 'trigger': trigger,
                                           'delay': str(d),
                                           'assignments': assignments }


    return {'species':species, 'parameters':parameters, 'reactions':reactions, 'events':events, 'mappings':mappings}


if __name__ == '__main__':
    params = parse.params(sys.argv[1])
