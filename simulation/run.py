#!/usr/bin/env python3

#################################################
#                                               #
#   Primary script for running a simulation     #
#                                               #
#################################################

import sys, os, inspect
from pathlib import Path
import argparse
import logging

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir, sim_dir = local_dir.absolute().parents[0], local_dir.absolute().parents[1]
sys.path.append(str(root_dir))
sys.path.append(str(sim_dir))

from models import builder
from lib import parse, pybash, util, confidence_intervals, stats, init_solver, log, rewrite_cps
import math
import copy
import pickle
from time import time
import optimize, plot
from pprint import pprint
import contextlib
import joblib
from tqdm import tqdm
from joblib import Parallel, delayed
import tempfile
import shutil

import numpy as np

module_logger = logging.getLogger('root')


# path progress bar
@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


############################ PRIMARY SIMULATION FUNCTIONS ####################################

def simulate(params, verbose=False, plot_run=True, parallel=True, rescale_to_1ml=True):

    if params['simulation_method'] in ['deterministic','Deterministic','LSODA']:
        if 'deterministic_reps' in params.keys() and params['deterministic_reps']:
            data = simulate_many(copy.deepcopy(params), verbose=verbose, parallel=parallel)
        else:
            data = simulate_one(params, verbose=verbose, reformat=True)
    elif params['repetitions'] > 1:
        data = simulate_many(copy.deepcopy(params), verbose=verbose, parallel=parallel)
    else:
        data = simulate_one(params, verbose=verbose)

    # potentially print state
    if 'print_final_state' not in params.keys():
        params['print_final_state'] = False
    if 'print_initial_state' not in params.keys():
        params['print_initial_state'] = False
    if 'print_header' not in params.keys():
        params['print_header'] = False

    # potentially rescale data to 1ml
    if rescale_to_1ml:
        if 'volume_ml' in params.keys():
            VOLUME_ML = params['volume_ml']
            for s in data.keys():
                if s != 'Time':
                    # real species
                    # print(f'rescale from: {data[s]}')
                    runs = data[s]['runs']
                    new_runs = [ [ value / VOLUME_ML for value in run] for run in runs ]
                    data[s]['runs'] = new_runs
                    # print(f'rescale to  : {data[s]}')

    # potentially save data
    if params['save_data']:
        module_logger.debug("Saving data (reason: parameter <save_data>)")

        plot_params = None
        if plot_run and params['plot']:
            # load and also include plot params
            plot_params = plot.get_plot_params(plot_json_filename=params['plot_json'])
        pickled = { 'data': data,
                    'params': params,
                    'plot_params': plot_params }

        # save it
        util.pickle_it(pickled, os.path.join(params['output_dir'], 'simulation.pickle'))

    else:
        module_logger.debug("NOT saving data (reason: parameter <save_data>)")

    if plot_run and params['plot']:
        # why 2 toggles for plot? One is if calling many times, such as during optimize, the other is for user control
        plot.plot_data_fromjsonfile(plot_json_filename=params['plot_json'], data=data)

    module_logger.info('simulate: Done')

    return data

def simulate_one(params, verbose=True, reformat=True):

    if not os.path.isdir(params['output_dir']):
        module_logger.info("Creating desired output directory: %s..." %(params['output_dir']))
        os.makedirs(params['output_dir'], exist_ok=True)

    if 'verbose' not in params.keys():
        params['verbose'] = verbose

    with tempfile.TemporaryDirectory() as tempdir:
        tempdir += '/'

        temp_params = copy.deepcopy(params)
        temp_params['cps_dir'] = tempdir
        temp_params['output_dir'] = tempdir
        temp_params['output_file'] = os.path.join(tempdir, os.path.basename(params['output_file']))

        # execute model py
        mappings = init(temp_params)
        module_logger.info("Built Copasi model in %s, now running it..." %(temp_params['cps_dir']+'model.cps'))

        shutil.copyfile(os.path.join(temp_params['cps_dir'],'model.cps'), os.path.join(params['cps_dir'],'model.cps'))

        data = run_solver(temp_params, mappings, verbose=verbose)

        shutil.copyfile(temp_params['output_file'], params['output_file'])


        if data is None or len(data['Time']) < params['duration']/params['sim_stepsize']-1: # copasi error
            module_logger.error("encountered copasi error.") # copasi error, should also be found in the sim itself
            return None

        module_logger.info("Finished Copasi run, saving and plotting...")

        if reformat:    # to have 1 common format btwn stoch and deterministic
            data = reformat_deterministic(params,data)

        if not params['save_data']: # will be only for deterministic runs
            os.remove(params['output_file'])

        return data


def simulate_many(params, verbose=False, parallel=False, JOBS=-1):

    def __single_run(packed):
        pa, i = packed
        instance_params = copy.deepcopy(pa)
        instance_params['output_file'] = orig_output.replace('.tsd',f'_rep{i}.tsd')
        instance_params['save_data'] = True # need to for multiprocessing
        # this currently does not work in parallel, since
        # copasise writes to the same .tsd file
        return simulate_one(instance_params, verbose=False, reformat=False)

    # typically for stochastic sims, where merge=True
    module_logger.info("Running the same simulation many times to obtain average and variance...")

    data = []
    orig_output = params['output_file']

    items = list(range(params['repetitions']))

    if parallel:
        if verbose:
            with tqdm_joblib(tqdm(desc="Simulate repetitions", total=len(items))) as progress_bar:
                data = Parallel(n_jobs=JOBS)(delayed(__single_run)((params,i)) for i in items)
        else:
            data = Parallel(n_jobs=JOBS)(delayed(__single_run)((params,i)) for i in items)
    else:
        for i in items:
            data += [ __single_run((params,i)) ]

    module_logger.info("Finished running the simulations. Merging them ...")
    merged_data = merge(params, data)

    return merged_data


#def sim_multi(multi_params,orig_params, target_data, opt=False):
#    # runs multiple sims, which have different initial conditions (ex. during opt or multi plot)
#    # as opposed to simulate_many() which merges multiple identical stochastic sims
#
#    # multi_params should come from opt_params or plot_params (json)
#    # whereas orig_params should come from a single-run params (txt)
#
#    sim_data = []
#    for datapoint in multi_params['dataset']:
#        if orig_params['debug']:
#            assert(len(datapoint['sim'])==len(datapoint['empirical']))
#            if opt:
#                #point in opt is to measure dist btwn sim and data, whereas plotting might do just one or the other
#                assert(len(datapoint['init_params'])==len(datapoint['init_vals']))
#
#        for j in range(len(datapoint['init_vals'])):
#            if isinstance(datapoint['init_vals'][j],str) and datapoint['init_vals'][j].lower() in ["time0","t0","zero"]:
#                #print('correctly setting',datapoint['init_params'][j],'=',target_data[datapoint["empirical"][j]][0])
#                if not isinstance(target_data[datapoint["empirical"][j]],dict):
#                    orig_params[datapoint['init_params'][j]] = target_data[datapoint["empirical"][j]][0]
#                else:
#                    orig_params[datapoint['init_params'][j]] = target_data[datapoint["empirical"][j]]['val'][0]
#            else:
#                orig_params[datapoint['init_params'][j]] = datapoint["init_vals"][j]
#
#        if datapoint['sim'] != []:
#            sim_data += [simulate(orig_params, verbose=False, plot_run=False)]
#            if sim_data[-1] is None:
#                if opt:
#                    optimize.console_log(multi_params,"ERROR optimize.copasi_with_targets(): Copasi error, see above.",verbose_lvl=0)
#                    params_used = [str(target_param['name']) + '=' + str(orig_params[target_param['name']]) for target_param in multi_params['target_params']]
#                    optimize.console_log(multi_params,'During attempt with params:' + str(params_used) + '\n',verbose_lvl=0)
#                else:
#                    module_logger.error("run.multi_sim: Copasi error, see above.")
#                return None
#    return sim_data

################################## HELPER FUNCTIONS ####################################

def merge(params, data):
    """
    a basic merge of the data with possible resampling in time to save less data

    params: dict of parameters
    data:   Array of single runs: data[0], data[1], ...
    """

    assert ( data[0]['Time'] == data[1]['Time'] ) #should be set at exact same times
    
    N = params['repetitions']

    # generate 'Time'
    if params['stats_stepsize'] is None or params['stats_stepsize'] == 0:
        # take time points from first run
        merged_data = {'Time': data[0]['Time']}
        time_points = data[0]['Time']
    else:
        # construct time points by resampling time
        merged_time, time_points =[],[]
        for t in range(len(data[0]['Time'])):
            if t*params['sim_stepsize'] % params['stats_stepsize'] == 0:
                merged_time += [data[0]['Time'][t]] 
                time_points += [t]
        merged_data = {'Time': merged_time} 

    for key in data[0].keys():
        if key not in ['Time']:
            merged_data[key] = {'runs':  [ data[i][key] for i in range(N) ],
                               }

    return merged_data
    

def multi_sim_for_plot(params):
    mapped_plots = plot.load_plot_params(params)
    plot.init_mpl(params)
    pickle_data=[]

    for title in mapped_plots.keys():
        if 'multi' in mapped_plots[title]:
            if mapped_plots[title]['ON'] and mapped_plots[title]['multi']:
                multi_params = mapped_plots[title]
                ext_data = multi_params["external_data"]
                if 'time_units' in ext_data.keys():
                    time_units = ext_data['time_units']
                else:
                    time_units = 'minutes'
                if 'pickle' not in mapped_plots[title].keys() or not mapped_plots[title]['pickle']:
                    target_data = parse.read_csv_target_file(ext_data['file'], time_units=time_units, reformat=False, transpose=False)
                else:
                    with open(ext_data['file'], 'rb') as f:
                        target_data = pickle.load(f)

                sim_data = sim_multi(multi_params, params, target_data, opt=False)

                pickle_data+=[{'sims':sim_data,'target_data':target_data,'params':params, 'title':title}]

                plot.multi(multi_params, params, target_data, sim_data,title=title)
    if params['save_data']:
        module_logger.info("Saving multi run data")
        util.pickle_it(pickle_data, os.path.join(params['output_dir'], params['timestamp'] + '_multi.pickle'))


def init(params):
    species, mappings = builder.build(params)

    if params['solver'].lower() in ['copasi','either']:
        # the cps and sbml file names must match those in models/builder.py and ../bin/rewrite_cps.py
        err, num_err = True, 0
        while err:
            err = pybash.run(command=['CopasiSE','-i',os.path.join(params['cps_dir'],'model.sbml'),'-s',os.path.join(params['cps_dir'],'pre.cps'),'--nologo'],logfile=params['solver_log'])

            if num_err > params['max_solver_failures']:
                module_logger.error("run.init(): solver returned an error max allowed times.")
                return None

        rewrite_cps.light_version(params, species)

    if params['solver'].lower() in ['ibiosim','either']:
        init_solver.ibiosim(params)

    #else:
    #   print("ERROR unknown solver parameter:",params['solver'])
    #   assert(False)

    return mappings


def run_solver(params, mappings, verbose=True):

    err, num_err = True, 0
    orig_solver = params['solver']
    while err:
        if num_err >= params['max_solver_failures']:
            if params['solver'].lower() == 'either':
                module_logger.error("run_solver: copasi failed, using iBioSim for this run.")
                os.remove(os.path.join(params['output_dir'],params['output_file']))#del the partially completed copasi file
                params['solver'] = 'ibiosim'
                num_err=0
            else:
                module_logger.error("run_solver: solver returned an error max allowed times.")
                return None

        if params['solver'].lower() in ['copasi','either']:
            err = pybash.run(command=['CopasiSE',os.path.join(params['cps_dir'],'model.cps'),'--nologo',],logfile=params['solver_log'], cwd=params['cps_dir'])
            if err and num_err == 0:
                # was first try
                module_logger.warning(f"retrying with stochastic solver {params['cps_dir']}")
                rewrite_cps.cps_set_stochmethod(params, method_to="stochastic")
                err = pybash.run(command=['CopasiSE',os.path.join(params['cps_dir'],'model.cps'),'--nologo',],logfile=params['solver_log'], cwd=params['cps_dir'])

        elif params['solver'].lower() == 'ibiosim':
            #TODO autofind reb2sac, ex based on OS?

            model_sbml = os.path.join(root_dir, params['cps_dir'], 'model.sbml')
            properties = os.path.join(root_dir, params['cps_dir'], 'ibiosim_props.txt')
            err = pybash.run(command=['reb2sac.exe','--target.encoding=rkf45','--reb2sac.properties.file=' + properties,model_sbml],logfile=params['solver_log'])
            # original bash: $REB2SAC --target.encoding=rkf45 --reb2sac.properties.file=$PROP $SBML 1>&2

            os.rename(os.path.join(params['output_dir'],'run-1.tsd'), os.path.join(params['output_dir'],params['output_file']))
            os.remove(os.path.join(params['output_dir'],'statistics.txt'))
            os.remove(os.path.join('./term-time.txt'))
        else:
            module_logger.error(f"unknown solver parameter: {params['solver']}")
            exit(1)

        num_err += 1

    data = parse.tsd(params)

    if data == [] or data is None:
        module_logger.error("solver simulation error inferred: I got an empty data set.")
        return None

    if mappings is not None:
        data = remap_species(data, params, mappings)

    params['solver'] = orig_solver
    return data


def reformat_deterministic(params, data):
    # to have 1 common format btwn stoch and deterministic
    re_data = {}
    for k in data.keys():
        if k == 'Time':
            # time and no species -> just take it
            re_data['Time'] = data['Time']
        
        else:
            # not time, but species
            re_k = k.replace('[','').replace(']','') # remove [] used for concentrations by the solver
            re_data[re_k] = {'runs': [] }
            re_data[re_k]['runs'] = [ data[k] ] # only single run

    return re_data


def remap_species(data, params, mapping):
    """
    Takes the simulated species (data) and add ones defined by
    mapping (mapping).

    By default, a mapping is a list of species in which case the
    sum is taken: ['a', 'b', ...]

    A different behavior can be specified by a string, however:
    '[a] * 2 + 7 - [b]'.
    MIND: dont forget the '[]'.
    """

    mapped_data = {'Time': data['Time']}
    T = len(data['Time'])

    # copy over all unmapped ones
    for k in data.keys():
        mapped_data[k]= data[k]

    # 1st pass with sum mappings
    for group in mapping.keys():
        the_mapping = mapping[group]
        try:
            # check if is a list -> sum
            if type(the_mapping) is list:
                species = the_mapping
                mapped_data[group] = [sum([data[species[i]][t] for i in range(len(species))]) for t in range(T)]

        except Exception as e:
            module_logger.error(f'run: remap_species: error when remapping "{the_mapping}".')
            module_logger.error(e)
            exit(1)
    # save as data 
    data = mapped_data

    # 2nd pass with function mappings
    for group in mapping.keys():
        the_mapping = mapping[group]
        try:
            # check if is str -> treat as formula
            if type(the_mapping) is str:
                # rewrite: [S] -> data['S'][t]
                the_mapping_new = the_mapping.replace("[", "data['").replace("]", "'][t]")
                # map with local scope
                mapped_data[group] = []
                for t in range(T):
                    locs = locals()
                    mapped_data[group] += [ eval(the_mapping_new,locs) ]

        except Exception as e:
            module_logger.error(f'run: remap_species: error when remapping "{the_mapping}".')
            module_logger.error(e)
            exit(1)

    return mapped_data

########################################################################################

if __name__ == "__main__":
    module_logger = log.setup_logger('root')

    parser = argparse.ArgumentParser(description='Simulation framework for CRNs.')
    parser.add_argument('PARAM_FILE',
                        nargs=1,
                        type=str,
                        help='the parameter file')
    parser.add_argument('--multi',
                        action=argparse.BooleanOptionalAction,
                        help='run multiple simulations')
    parser.add_argument('--sequential',
                        action=argparse.BooleanOptionalAction,
                        help='run simulations sequentially')
    parser.add_argument("-log", "--log", 
                        default="warning",
                        help=(
                            "Provide logging level. "
                            "Example --log debug', default='warning'")
                        )
    levels = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warn': logging.WARNING,
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG
    }
    
    # get arguments
    args = parser.parse_args()

    # set debug level
    level = levels.get(args.log.lower())
    if level is None:
        module_logger.error(
            f"log level given: {options.log}"
            f" -- must be one of: {' | '.join(levels.keys())}")
        exit(1)
    module_logger.setLevel(level)

    # check params file
    param_file = args.PARAM_FILE[0]
    if not os.path.isfile(param_file):
        module_logger.error(f"Cannot find param file {param_file}. Check its path.")
        exit(1)

    module_logger.info(f"Starting simulation with param file: {param_file}")
    params = parse.params(param_file)
    if args.multi:
        #TODO: add parallel here
        multi_sim_for_plot(params)
    else:
        simulate(params, parallel=(not args.sequential), verbose=True)
