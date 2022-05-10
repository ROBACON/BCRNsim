#!/usr/bin/env python3

#####################################################
#                                                    #
#    Low-level parser for reading in a tab-delimited #
#    parameter file and for reading out a            #
#    tab-delimited data file                         #    
#                                                    #
#####################################################

import sys, os, inspect
from pathlib import Path

this_file_path = Path(inspect.getfile(inspect.currentframe()))
this_dir = this_file_path.absolute().parents[0]
LIGHT_DIR = this_file_path.absolute().parents[0]
root_dir = this_file_path.absolute().parents[0]
sys.path.append(str(this_dir))
sys.path.append(str(root_dir))

from util import bool, timestamp, rng
from lib.TSDReader import readtsd, writetxt
import math
import pandas as pd


def params(param_file):
    # EXPECT A TAB DELITED FILE: name\value\dtype
    if not os.path.isfile(param_file):
        print("ERROR: can't find parameter file at: ", param_file)
        assert(False) # unknown file
    with open(param_file) as pfile:
        lines = pfile.readlines()
        pdict = {} # PARAMS
        for line in lines[1:]: #skip header line
            if not line.isspace():
                if line[0] != '#': #used as comment
                    param = line.strip().split('\t')
                    # handle removing multiple tabs
                    param = [p for p in param if p != '']
                    assert (len(param) in [3,4]), f"ERROR: Cannot parse the following line in {param_file}:\n{line}"
                    if param[1][0] == '[': #ie is list
                        param[1] = param[1].replace('[','').replace(']','')
                        param[1] = param[1].split(',')

                        is_list = True
                    else:
                        is_list = False


                    if param[2] in ['str','string']:
                        if not is_list:    val = param[1]
                        else: val = [p for p in param[1]]
                    elif param[2] == 'int':
                        if not is_list: val = int(float(param[1]))
                        else: val = [int(p) for p in param[1]]
                    elif param[2] == 'float':
                        if not is_list: val = float(param[1])
                        else: val = [float(p) for p in param[1]]
                    elif param[2] == 'bool':
                        if not is_list: val = bool(param[1])
                        else: val = [bool(p) for p in param[1]]
                    elif param[2] == 'eval':
                        if not is_list: val = eval(param[1])
                        else: val = [eval(p) for p in param[1]]
                    elif param[2] == 'exp':
                        pieces = param[1].split('e')
                        val = float(pieces[0]) * math.pow(10,float(pieces[1]))
                    elif param[2] == 'csvcolumn':
                        pieces = param[1].split('[')
                        assert ( len(pieces) == 2 )
                        assert ( pieces[1][-1] == ']' )
                        filename = pieces[1][:-1]
                        colname = pieces[0]
                        # reead csv
                        df = pd.read_csv(filename,
                            sep=',',
                            comment='#')
                        val = df[colname][0]
                    else:
                        print("\nERROR: unknown datatype for parameter %s. Check tab-delimiters too.\n" %(param[0]))
                        assert(False) #unknown val
                    pdict[param[0]] = val


    pdict['timestamp'] = timestamp()

    for dir_param in ['output_dir','cps_dir']:
        pdict[dir_param] = os.path.join(LIGHT_DIR, pdict[dir_param])

    # appending output_dir to output_file is necessary before writing full path of output_dir
    # why? because the copasi script will dump into a local path from run.py
    pdict['output_file'] = pdict['output_dir'] + pdict['output_file'].replace('TIMESTAMP',pdict['timestamp'])

    if 'solver_log' not in pdict.keys():
        pdict['solver_log'] = None #TODO: should clean this, defaults to std.out is 'solver_log' is None

    return pdict


def tsd(params, root_light_dir=True):
    """
    Time is assumed to be FLOAT, all other also assumed to be FLOAT
    
    Returns:
    ...

    """
    data = None

    tsd_file = params['output_file']
    if root_light_dir:
        tsd_file = os.path.join(LIGHT_DIR, tsd_file)

    if params['solver'].lower() in ['copasi','either']:
        with open(tsd_file, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines)):
                line = lines[i]
                line = line.replace('\n','')
                parts = line.split('\t')
                if i==0:
                    data = {part.replace('[','').replace(']',''):[] for part in parts}
                    header = [part.replace('[','').replace(']','') for part in parts]
                else:
                    for j in range(len(parts)):
                        part = float(parts[j])
                        data[header[j]] += [part]

    elif params['solver'].lower() == 'ibiosim':
        with open(tsd_file, 'r') as file:
            header, data = readtsd(tsd_file, epsilon=0.0)
            data['Time'] = data['time']

    return data



def read_csv_target_file(target_file, transpose=False, reformat=True, time_units='minutes'):
    """
    If reformat returns data as:
    { 'Time': [t0, t1, ..]
      species0 : {'runs': [run0, run1, run2, ..]},
      species1 : {'runs': [run0, run1, run2, ..]}, ... }
    """

    if time_units in ['hour', 'hours']:
        time_multiplier = 60
    elif time_units in ['minutes','min']:
        time_multiplier = 1
    else:
        print("\nparse: read_csv_target_file: unrecognized time_units: ",time_units,'\n')
        assert(False)

    if not transpose:
        with open(target_file) as f:
            lines = f.readlines()
            lines = [lines[i].replace('\n','') for i in rng(lines) if not lines[i][0] == '#']
            header = lines[0].split(',')
            lines = lines[1:]
            data = {k: [] for k in header}
            for i in range(len(lines)):
                line = lines[i].split(',')
                for j in range(len(header)):
                    if header[j] in ['Time','time']:
                        data[ 'Time' ] += [ float(line[j]) * time_multiplier ]
                    else:
                        data[ header[j] ] += [ float(line[j]) ]
    else:
        assert False, "fix me"
        """
        with open(target_file) as f:
            lines = f.readlines()
            lines = [lines[i].replace('\n','') for i in rng(lines) if not lines[i][0] == '#']
            data = {}
            for i in range(1,len(lines)):
                line = lines[i].split(',')
                line = [line[0]] + [float(line[i]) for i in range(1,len(line))]
                if line[0] in ['Time','time']:
                    line = [line[0]] + [line[i]*time_multiplier for i in range(1,len(line))]
                data[line[0]] = line[1:]
        """

    if reformat:
        re_data={ 'Time': data['Time'] }
        for key in data.keys():
            if key not in ['Time','time']:
                re_data[key] = { 'runs': [data[key]] }
        data = re_data

    return data

######################################################################

