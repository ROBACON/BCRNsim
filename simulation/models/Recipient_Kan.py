#!/usr/bin/env python3

#####################################################
#                                                   #
#    E.coli Recipient Antibiotic model              #
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
import numpy as np

module_logger = logging.getLogger('root')

ANTIBIOTIC_TYPES = ['APhage', 'AReceiver', 'ADonor']
VOLUME_ML = 1


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
    # return __infected(ecoli)  # THIS WORKS
    return __infected(ecoli) and (not __is_blocked(ecoli))  # NEW
    # return __infected(ecoli) and (__pili(ecoli) < params['num_pili'])

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


# antibiotics
def __resistant(a_type: str, ecoli: str):
    """
    Returns: if resistant a_type
    """
    if a_type == 'APhage':
        return '_nr' not in ecoli
    elif a_type == 'ADonor':
        return True
    else:
        return False

def __can_uptake_antibiotic(a_type: str, ecoli: str, params):
    """
    Returns: if can uptake antibiotic
    """
    assert a_type in ANTIBIOTIC_TYPES, "unknown a_type"
#    return True
    return (__antibiotic_slots(a_type, ecoli) < params[f'antibiotic_slots'])

def __antibiotic(a_type: str, ecoli: str, antibiotic_threshold=0):
    """
    Returns: if antibiotic inside
    """
    assert a_type in ANTIBIOTIC_TYPES, "unknown a_type"
    return (__antibiotic_slots(a_type, ecoli) > antibiotic_threshold)


def __any_antibiotic(ecoli: str, antibiotic_threshold=0):
    """
    Returns: if any antibiotic inside
    """
    return any([ __antibiotic(a_type, ecoli, antibiotic_threshold) for a_type in ['APhage', 'AReceiver'] ])


def __antibiotic_slots(a_type: str, ecoli: str):
    """
    Returns: number of antibiotc slots taken
    """
    parts = ecoli.split(f'_{a_type}')
    parts2 = parts[1].split('_')
    return int(parts2[0])

def __is_resistant_to_all_antibiotics_inside(ecoli: str):
    return all([ not __antibiotic(a_type, ecoli) or __resistant(a_type, ecoli) for a_type in ANTIBIOTIC_TYPES ])


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
    global VOLUME_ML
    rates = {
        'birth': {
            'G': None,
            'AA': None
        },
        'death': None
    }
    init_glucose = int( params['init_glucose'] * VOLUME_ML )
    init_amino_acid = int( params['init_amino_acid'] * VOLUME_ML )

    # birth
    if __any_antibiotic(ecoli, antibiotic_threshold=0):
        if __is_resistant_to_all_antibiotics_inside(ecoli):
            rates['birth']['G']  = f"dup_rate * antibiotic_resistant_birth_penalty_factor / {init_glucose}"
            rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor * antibiotic_resistant_birth_penalty_factor / {init_amino_acid}"
        else:
            rates['birth']['G']  = f"dup_rate * antibiotic_birth_penalty_factor / {init_glucose}"
            rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor * antibiotic_birth_penalty_factor / {init_amino_acid}"

    else:
        # no antibiotic inside
        rates['birth']['G']  = f"dup_rate / {init_glucose}"
        rates['birth']['AA'] = f"dup_rate * amino_acid_penalty_factor / {init_amino_acid}"

    # birth: plasmid load
    if __in_early_infection_phase(ecoli):
        rates['birth']['G'] += f" * {params['plasmid_penalty_factor']}"
        rates['birth']['AA'] += f" * {params['plasmid_penalty_factor']}"

    elif __in_late_infection_phase(ecoli):
        rates['birth']['G'] += f" * {params['plasmid_penalty_factor_late']}"
        rates['birth']['AA'] += f" * {params['plasmid_penalty_factor_late']}"

    # death
    if __any_antibiotic(ecoli, antibiotic_threshold=params['antibiotic_threshold_slots']):
        if __is_resistant_to_all_antibiotics_inside(ecoli):
            rates['death'] = f"antibiotic_resistant_death_rate"

        else:
            rates['death'] = f"antibiotic_death_rate"

    else:
        rates['death'] = 'death_rate'

    return rates


def build(params):

    # potentially set volume
    global VOLUME_ML
    if 'volume_ml' in params.keys():
        VOLUME_ML = params['volume_ml']
        print(f'Receiver: volume has been set to {VOLUME_ML} ml.')

    assert (params['antibiotic_threshold_slots'] < params['antibiotic_slots']), "antibiotic threshold is too large"

    if params['debug']:
        if params['num_pili'] == 1:
            module_logger.info(f"num pili currently set to 1")
        if params['init_phages'] != 0:
            module_logger.info(f"using phages")

    events = {}
    species = {
        'G': params['init_glucose'],
        'AA': params['init_amino_acid'],
        'P': params['init_phages'],
        'APhage': params['init_antibiotics#APhage'],
        'AReceiver': params['init_antibiotics#AReceiver'],
        'ADonor': params['init_antibiotics#ADonor'],
        'dead': 0,
        'gp3': params['init_gp3'],
        'inf_count': 0,  # how many msgs where received
    }

    # --- species & mappings --- #

    # Various E. coli:
    # keep each __identifier__ unique, since the code uses string replacement for reactions
    mappings = {
        'Resources':['G','AA'],
        'Phages':['P'],
        'Antibiotics APhage':['APhage'],
        'Ecoli':[],
        'Ecoli uninfected':[],
        'Ecoli infected':[],
        # resitant
        'Ecoli not resistant to APhage':[],
        'Ecoli resistant to APhage':[],
        # infection
        'Ecoli early infection':[],
        'Ecoli late infection':[],
        # antibiotics
        'Ecoli not antibiotic':[],
        'Ecoli antibiotic':[],
        # gp3
        'Ecoli blocked':[],
        'Ecoli nonblocked':[],
        # age
        'Ecoli young':[],
        'Ecoli old':[],
        'OD': '0.01 * [dead] + [Ecoli]'
    }
    for i in range(params['num_pili']+1):
        mappings[f'Ecoli {i} pili'] = []

    for infct in ['ni','ei','li']:
        for aslots in itertools.product(range(params['antibiotic_slots']+1), repeat=len(ANTIBIOTIC_TYPES)):
            # aslots = stple that holds the taken slots for both antibiotics
            antib_name = '_'.join([f'{ANTIBIOTIC_TYPES[i]}{aslots[i]}' for i in range(len(ANTIBIOTIC_TYPES))])
            for resistance in ['nr', 'ar']:
                for i in range(params['num_pili']+1):
                    pili_name = f'p{i}'
                    for gp3 in ['blocked', 'nonblocked']:
                        for age in ['agey', 'ageo']: # young and old
                            s = f'E_{infct}_{antib_name}_{resistance}_{gp3}_{age}_{pili_name}'
                            
                            # init species
                            init_species_young =    f'E_ni_APhage0_AReceiver0_ADonor0_nr_nonblocked_agey_p0'
                            init_species_old =      f'E_ni_APhage0_AReceiver0_ADonor0_nr_nonblocked_ageo_p0'
                            if s == init_species_young:
                                species[s] = int(params['init_ecoli'] * (1-params['init_ratio_old']))
                            elif s == init_species_old:
                                species[s] = int(params['init_ecoli'] * params['init_ratio_old'])
                            else:
                                species[s] = 0
                            
                            # mappings
                            mappings['Ecoli'] += [s]
                            if not __infected(s):
                                mappings['Ecoli uninfected'] += [s]
                            else:
                                mappings['Ecoli infected'] += [s]
                                if __in_early_infection_phase(s):
                                    mappings['Ecoli early infection'] += [s]
                                elif __in_late_infection_phase(s):
                                    mappings['Ecoli late infection'] += [s]

                            if __any_antibiotic(s):
                                mappings['Ecoli antibiotic'] += [s]
                            else:
                                mappings['Ecoli not antibiotic'] += [s]

                            if f'Ecoli {antib_name} antibiotics' in mappings.keys():
                                mappings[f'Ecoli {antib_name} antibiotics'] += [s]
                            else:
                                mappings[f'Ecoli {antib_name} antibiotics'] = [s]

                            if '_nr' in s:
                                mappings["Ecoli not resistant to APhage"] += [s]
                            else:
                                mappings["Ecoli resistant to APhage"] += [s]

                            if __is_blocked(s):
                                mappings['Ecoli blocked'] += [s]
                            else:
                                mappings['Ecoli nonblocked'] += [s]

                            if __is_old(s):
                                mappings['Ecoli old'] += [s]
                            else:
                                mappings['Ecoli young'] += [s]

                            mappings[f'Ecoli {i} pili'] += [s]


    # potentially dilute to volume less than 1ml
    for s in species.keys():
        species[s] = int( species[s] * VOLUME_ML )


    # --- pre-infection --- #
    if params['pre_inf']:
        module_logger.info(f"Pre infecting cells via simulation with frozen cells...")
        params_pre = deepcopy(params)

        # -- do not recurse
        params_pre['pre_inf'] = False
        # -- do not save and plot
        params_pre['save_data'] = False
        params_pre['plot'] = False
        # -- just a single run
        params_pre['repetitions'] = 1
        # -- duration
        params_pre['sim_stepsize'] = 1
        params_pre['duration'] = params['pre_inf_duration']

        # -- no antibiotics
        for  a_type in ANTIBIOTIC_TYPES:
            params_pre[f'init_antibiotics#{a_type}'] = 0

        # -- gp3 is already in medium
        #params_pre['init_gp3'] = 0

        # freezing?
        if params['pre_inf_freeze']:
            for key in ['dup_rate',
                'death_rate',
                'early_to_late_inf_rate',
                'young_to_old_rate']:
                params_pre[key] = 0

        # simulate
        # not rescale for preinfection - this will be done in the end by 'run'
        data = run.simulate(params_pre, verbose=False, plot_run=False, rescale_to_1ml=False)

        # take result as initial values
        for s in species.keys():
            species[s] = np.mean([ run[-1] for run in data[s]['runs'] ])

        # but:
        # take antibiotics from parameter file
        for a_type in ANTIBIOTIC_TYPES:
            species[a_type] = params[f'init_antibiotics#{a_type}']
        # take gp3 from parameter file
        #species['gp3'] = params['init_gp3']

        # some statistics
        num_infected_after, total_E_after = 0, 0
        for s in species.keys():
            if 'E_' in s:
                total_E_after += species[s]
                if __infected(s):
                    num_infected_after += species[s]

        if params['init_phages'] > 0:

            module_logger.warning(f"Preinfection: # infected = {num_infected_after}")

            growth_percent = total_E_after / params['init_ecoli']
            module_logger.info(f"Preinfection: factor growth = {growth_percent:.1f}")

    else:
        # not preinfected
        moi = 'na' if params['init_ecoli'] == 0 else params['init_phages']/params['init_ecoli']
        module_logger.info(f"MOI = {moi}")


    # --- reaction parameters --- #
    parameters = {
        # gp3
        'gp3_block_rate':                               (params['gp3_block_rate'],                              'per_min'),
        'gp3_unblock_rate':                             (params['gp3_unblock_rate'],                            'per_min'),
        # tune with no phages, no antibiotics
        'dup_rate':                                     (params['dup_rate'],                                    'per_min'),
        'antibiotic_birth_penalty_factor':              (params['antibiotic_birth_penalty_factor'],             'dimensionless'),
        'antibiotic_resistant_birth_penalty_factor':    (params['antibiotic_resistant_birth_penalty_factor'],   'dimensionless'),
        'amino_acid_penalty_factor':                    (params['amino_acid_penalty_factor'],                   'dimensionless'),
        'death_rate':                                   (params['death_rate'],                                  'per_min'),
        'antibiotic_death_rate':                        (params['antibiotic_death_rate'],                       'per_min'),
        'antibiotic_resistant_death_rate':              (params['antibiotic_resistant_death_rate'],             'per_min'),
        # tune with phages, but no antibiotics
        'inf_rate':                                     (params['inf_rate'],                                    'per_min'),
        'inf_factor_young':                             (params['inf_factor_young'],                            'dimensionless'),
        'inf_factor_antibiotic':                        (params['inf_factor_antibiotic'],                       'per_min'),
        'young_to_old_rate':                            (params['young_to_old_rate'],                           'per_min'),
        'uninf_rate':                                   (params['uninf_rate'],                                  'per_min'),
        'early_to_late_inf_rate':                       (params['early_to_late_inf_rate'],                      'per_min'),
        'phage_decay_rate':                             (params['phage_decay_rate'],                            'per_min'),
    }
    reactions = {}
    for s in species:

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

            for res in ['G','AA']:
                # birth per resource
                if  __antibiotic('AReceiver', ecoli, antibiotic_threshold=0) and not __resistant('AReceiver', ecoli):
                    # birth leads to death with some probability 'antibiotic_death_on_duplication_prob'
                    # Remark: we do not consume food on this
                    reactions[f'antibiotic_death_on_birth_AReceiver_{ecoli}_{res}'] = {
                        're': [(1, ecoli), (1, res)],
                        'pr': [(1, 'dead'), (1, res)],
                        'kin': f"{ecoli} * {res} * {params['antibiotic_death_on_duplication_prob#AReceiver']} * {birthdeathrates['birth'][res]}"
                    }

                elif  __antibiotic('APhage', ecoli, antibiotic_threshold=0) and not __resistant('APhage', ecoli):
                    # birth leads to death with some probability 'antibiotic_death_on_duplication_prob'
                    # Remark: we do not consume food on this
                    reactions[f'antibiotic_death_on_birth_APhage_{ecoli}_{res}'] = {
                        're': [(1, ecoli), (1, res)],
                        'pr': [(1, 'dead'), (1, res)],
                        'kin': f"{ecoli} * {res} * {params['antibiotic_death_on_duplication_prob#APhage']} * {birthdeathrates['birth'][res]}"
                    }

                else:
                    reactions[f'birth_{ecoli}_{res}'] = {
                        're': [(1, ecoli), (1, res)],
                        'pr': [(1, l_child), (1, r_child)],
                        'kin': f"{ecoli} * {res} * {birthdeathrates['birth'][res]}"
                    }

            # DEATH
            reactions[f'death_{ecoli}'] = {
                're': [(1, ecoli)],
                'pr': [(1, 'dead')],
                'kin': f"{ecoli} * {birthdeathrates['death']}"
            }
            reactions[f'decay_P'] = {
                're': [(1, 'P')],
                'pr': [],
                'kin': f"P * phage_decay_rate"
            }

            # INFECTION
            # kinetics
            if __is_old(ecoli):
                inf_kin = f'inf_rate * P * {ecoli} / {VOLUME_ML}'
            else:
                inf_kin = f'inf_rate * inf_factor_young * P * {ecoli} / {VOLUME_ML}'

            # clogging
            if params['clogging_infection']:
                inf_kin += f" * {1 - parent_pili / params['num_pili']}"

            # infection reactions
            if __infectable(ecoli, params):
                assert parent_pili+1 <= params['num_pili'], "Check: at least 1 pili free"
                
                # reaction
                if __antibiotic('APhage', ecoli, antibiotic_threshold=0):
                    # already antibiotic with 'APhage'
                    # -> use a potentially smaller infection rate
                    inf_kin = f'inf_rate * inf_factor_antibiotic * P * {ecoli} / {VOLUME_ML}'
            
                reactions[f'infect_{ecoli}'] = {
                    're': [(1, ecoli), (1, 'P'), ],
                    'pr': [(1, ecoli.replace('_ni','_ei').replace(f"_p{parent_pili}",f"_p{parent_pili+1}").replace('_nr','_ar')),
                           (1, 'inf_count'),
                          ],
                    'kin' : inf_kin,
                }

            # further takeup of phages
            if __infected_and_can_bind_phages(ecoli, params=params):
                # uptake of further phages on pili
                reactions[f'already_infected_take_phage_{ecoli}'] = {
                    're': [(1, ecoli), (1, 'P'), ],
                    'pr': [(1, ecoli.replace(f'_p{parent_pili}',f'_p{min(parent_pili+1,params["num_pili"])}'))],
                    'kin' : inf_kin,
                }

            reactions['dead_takeup_phages'] = {
              're': [(1, 'dead'), (1, 'P')],
              'pr': [(1, 'dead')],
              'kin': f'dead * P * inf_rate / {VOLUME_ML}'
            }

            if __infected(ecoli):
                # spontaneous uninfected
                reactions[f'spont_uninf_{ecoli}'] = {
                    're': [(1, ecoli)],
                    'pr': [(1, ecoli.replace('_ei','_ni').replace('_li','_ni').replace('_ar','_nr'))],
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
                  'pr': [(1, ecoli.replace('_blocked','_nonblocked').replace(f'_p{pili}',f'_p{max(0,pili-1)}') )],
                  'kin': f'{ecoli} * gp3_unblock_rate',
                }
            else:
                # kinetics
                if __is_old(ecoli):
                    gp3_kin = f'gp3_block_rate * gp3 * {ecoli} / {VOLUME_ML}'
                else:
                    gp3_kin = f'gp3_block_rate * inf_factor_young * gp3 * {ecoli} / {VOLUME_ML}'

                reactions[f'blocking_{ecoli}'] = {
                    're': [(1, ecoli), (1, 'gp3')],
                    'pr': [(1, ecoli.replace('_nonblocked','_blocked').replace(f'_p{pili}',f'_p{min(pili+1,params["num_pili"])}') )],
                    'kin' : gp3_kin,
                }

            # ANTIBIOTICS
            for a_type in ANTIBIOTIC_TYPES:
                if __can_uptake_antibiotic(a_type, ecoli=ecoli, params=params):
                    num_slots_taken = __antibiotic_slots(a_type, ecoli)
                    reactions[f'antibiotic_exposure_{a_type}_{ecoli}'] = {
                        're': [(1, ecoli),(1, a_type)],
                        'pr': [(1, ecoli.replace(f'_{a_type}{num_slots_taken}',f'_{a_type}{min(num_slots_taken+1,params["antibiotic_slots"])}'))],
                        'kin': f"{a_type} * {ecoli} * {params[f'antibiotic_uptake_rate#{a_type}']} / {VOLUME_ML}"
                    }

            if __antibiotic('APhage', ecoli=ecoli) and (not __infected(ecoli)):
                reactions[f'antibiotic_spont_res_{ecoli}'] = {
                    're': [(1, ecoli) ],
                    'pr': [(1, ecoli.replace('_nr','_ar'))],
                    'kin' : f"{ecoli} * {params['antibiotic_spont_res']}",
                }

            # AGE
            if not __is_old(ecoli):
                # is young
                reactions[f'grow_up_{ecoli}'] = {
                    're': [(1, ecoli) ],
                    'pr': [(1, ecoli.replace('_agey','_ageo'))],
                    'kin' : f"{ecoli} * young_to_old_rate",
                }

    # potentially stop simulation early by setting all species to 0
    # this is used to speed up long stochastic simulations with large species numbers
    if 'stop_simulation_at_resources_consumed' in params.keys():
        module_logger.warning(f'Stopping simulation earlier when resources have been consumed: {params["stop_simulation_at_resources_consumed"]}.')
        
        max_raise = params['stop_simulation_at_resources_consumed']
        
        th = sum( [species[s] for s in mappings['Resources']] ) - max_raise
        trigger = '(' + ' + '.join( [f'{s}' for s in  mappings['Resources']] ) + f') <= {th}'
        assignments = [ (s, '0') for s in species.keys() ]
        events['stop'] = { 'trigger': trigger,
                            'delay': '0',
                            'assignments': assignments }


    return {'species':species, 'parameters':parameters, 'reactions':reactions, 'events':events, 'mappings':mappings}


if __name__ == '__main__':
    params = parse.params(sys.argv[1])
