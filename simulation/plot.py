#!/usr/bin/env python3

#####################################################
#                                                   #
#   Collection of plotting function                 #
#   Originally for bacteria/simulation/light        #
#   But could easily be adapted for other uses      #                                               
#                                                   #
#####################################################

"""

Expects data to be in the following format:

'runs': [run0, run1, run2, ..]

In case of external data and deterministic simulations this is just
added as:

'runs': [run0]

with a possible entry also for

'Time': [t0, t1, ...]

The Time entry is if the data (expecially external data) has its own timebase.

For plotting this is converted into lines and shades
via the following steps:

0) precompute different runs
1) statistics: runs -> lines & shades
2) plot these

"""


import logging
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.lines import Line2D
rcParams['font.family'] = 'serif'

import sys, os, inspect, pprint, json
from pathlib import Path
import pprint

from copy import deepcopy

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir = local_dir.absolute().parents[0]
# sys.path.append(str(root_dir))

from lib.util import rng, none, true_in_dict
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np, pickle, json, copy
from lib import parse, log, stats
from PIL import Image

# -- begin style
import figure_styles.nature_small
from figure_styles.nature_small import *
# -- end style

module_logger = logging.getLogger('root')


LINESTYLES = ['-',':','--','-.']
MARKERS = ['o','v','s','D','X']
COLORS = ['#009999', '#cc0066','red','blue','green','purple','cyan','orange','magenta','olive','grey','black']


# global lines
legend = None
lines = []
lined = {} # legend line -> line


# path json encoder
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)



######################## INTERACTION ###########################################

def __on_pick(event):
    global lined
    legline = event.artist
    origline = lined[legline]
    visible = not origline.get_visible()

    origline.set_visible(visible)
    legline.set_alpha(1.0 if visible else 0.2)
    ax = plt.gca()
    ax.figure.canvas.draw()

def start_picking():
    global lines, lined, legend
    for legline, origline in zip(legend.get_lines(), lines):
        legline.set_picker(True)
        lined[legline] = origline
    ax = plt.gca()
    ax.figure.canvas.mpl_connect('pick_event', __on_pick)


########################## MERGING PLOTS #####################################

def merge_images(titles, width_figs=5, output_file='out.png'):
    """
    Merges the files in titles into a single figure
    """
    images = [ Image.open(x) for x in titles ]
    widths, heights = zip(*(i.size for i in images))
    W = max(widths)
    H = max(heights)

    total_width = width_figs * W
    height_figs = len(images) // width_figs + (0 if len(images) % width_figs == 0 else 1)
    total_height = height_figs * H
    new_im = Image.new('RGB', (total_width, total_height))

    for i, im in enumerate(images):
        x_offset = (i % width_figs) * W
        y_offset = (i // width_figs) * H
        new_im.paste(im, (x_offset, y_offset))
    new_im.save(output_file)
    module_logger.info(f"written: {output_file}")


########################## HELPER FUNCTIONS #####################################

def init_mpl(mpl_use=None, save_fig=True):
    """
    Initializes the plot backend
    """
    try:
        if mpl_use is not None:
            if mpl_use == 'Agg':
                matplotlib.use('Agg')
            elif mpl_use == 'TkAgg':
                matplotlib.use('TkAgg')
            elif mpl_use == 'WX':
                matplotlib.use('WX')
            elif mpl_use == 'QTAgg':
                matplotlib.use('QTAgg')
            elif mpl_use == 'QT4Agg':
                matplotlib.use('QT4Agg')

        elif save_fig:
            matplotlib.use('Agg')

        else:
            matplotlib.use('TkAgg')

    except ModuleNotFoundError:
        # fall back to MacOSX lib
        matplotlib.use('MACOSX')


def get_plot_params(plot_json_filename):
    """
    In:
      plot_json_filename: json file name
    
    Returns:
      plot parameter dictionary
    """
    with open(plot_json_filename,'r') as file:
        try:
            json_data = json.load(file)
        except json.decoder.JSONDecodeError as e:
            module_logger.error(f'Error while decoding json file "{plot_json_filename}".')
            module_logger.error(f'{e}')
            module_logger.error(f'File content:\n{file.read()}')
            exit(1)

        if 'default' not in json_data.keys():
            module_logger.error('Please set a <default> parameter in the plot params file.')
            exit(1)

        if ('plots' not in json_data.keys()):
            module_logger.error('Please set a <plots> parameter in the plot params file.')
            exit(1)

    return json_data


def apply_default_params(plots, default_params):
    """
    Returns: plots filled up with default params
    """
    plots_new = deepcopy(plots)

    for plotname, plot in plots.items():
        # a plot
        for default_key in default_params.keys():
            if default_key not in plot.keys():
                # not specified so add
                plots_new[plotname][default_key] = default_params[default_key]

    return plots_new


def getColor(col_parameter: str) -> str:
    """
    Takes a color parameter from the plot file and
    returns a valid color str for matplotlib.

    Can resolve colors defined in a figure style
    with the keyword 'STYLE:'.

    examples:
    getColor('red') = 'red'
    getColor('#aaffbb') = '#aaffbb'
    getColor('STYLE:ECOLI_PINK') = '#'
    """
    if col_parameter[:len('STYLE:')] == 'STYLE:':
        # try to resolve with color style
        COLOR_VAR = col_parameter[len('STYLE:'):]
        candidates = [v for v in globals() if type(v) == str]

        if COLOR_VAR in candidates:
            return eval(COLOR_VAR)
        else:
            module_logger.error(f'Problem when trying to resolve the color parameter {col_parameter}')
            module_logger.error(f'Could not find a color with name {COLOR_VAR}')
            module_logger.error(f'Candiates are {candidates}')
            exit(1)

    else:
        return col_parameter


########################################################################################

def from_pickle(pickle_file, multiplot=False):

    with open(pickle_file, 'rb') as pickle_file_handle:
        payload = pickle.load(pickle_file_handle)

        data = payload['data']
        plot_params = payload['plot_params']
        
        # pprint.pprint(data)
        # pprint.pprint(plot_params)

        plot_data(data=data, plot_params=plot_params)


########################## DEFAULT SIMULATION PLOT #####################################

def plot_combine_runs(plot_params, data_dict):

    merged_plot_params = deepcopy(plot_params)
    for plot_name in plot_params['plots'].keys():
        merged_plot_params['plots'][plot_name]['species'] = []

    timebase = list(data_dict.values())[0]['Time']
    merged_data = {'Time': timebase}
    for name, data in data_dict.items():
        if not (data['Time'] == timebase):
            module_logger.error('plot_combine_runs: different timebases')
            raise Exception('fatal error')

        for species_name, species_data in data.items():
            merged_data[f'{name}: {species_name}'] = species_data

        for plot_name in plot_params['plots'].keys():
            # for a plot do:
            for species_name in plot_params['plots'][plot_name]['species']:
                merged_plot_params['plots'][plot_name]['species'] += [f'{name}: {species_name}']

            for style_name in plot_params['plots'][plot_name]['styles'].keys():
                if style_name == 'default':
                    # just copy
                    merged_plot_params['plots'][plot_name]['styles'][style_name] = plot_params['plots'][plot_name]['styles'][style_name]
                else:
                    merged_plot_params['plots'][plot_name]['styles'][f'{name}: {style_name}'] = plot_params['plots'][plot_name]['styles'][style_name]

    # print(merged_data)
    plot_data(data=merged_data, plot_params=merged_plot_params)



def plot_data_fromjsonfile(data, plot_json_filename):
    """
    The main plot function if started from json file.
    """
    plot_params = get_plot_params(plot_json_filename=plot_json_filename)
    plot_data(data=data, plot_params=plot_params)


def plot_data(data, plot_params):
    """
    The main plot function if started with plot parameters as dictionary.
    """
    if 'default' in plot_params.keys():
        default_params = plot_params['default']
    else:
        default_params = {}

    if 'plots' in plot_params.keys():
        plots = plot_params['plots']
    else:
        plots = {}

    # apply default parameters to the plots dict
    plots_dict = apply_default_params(plots=plots, default_params=default_params)

    try:
        save_fig = default_params['save_fig']
    except KeyError as e:
        module_logger.error(e)
        module_logger.error('In the plot parameters there must be a section <default> with a bool <save_fig>.')
        module_logger.error(f'I only found: {list(default_params.keys())}')
        exit(1)

    try:
        output_dir = default_params['output_dir']
    except KeyError as e:
        module_logger.error(e)
        module_logger.error('In the plot parameters there must be a section <default> with a str <output_dir>.')
        module_logger.error(f'I only found: {list(default_params.keys())}')
        exit(1)

    # check if plots to do
    if plots_dict == {}:
        # nothing to do
        return

    # setup plotting backend
    if 'mpl_use' in plot_params.keys():
        mpl_use = plot_params['mpl_use']
    else:
        mpl_use = None
    init_mpl(mpl_use=mpl_use, save_fig=save_fig)

    # loop over plots
    for title in plots_dict.keys():
        current_plot_dict = plots_dict[title]

        if current_plot_dict['ON']:
            # this plot is activated

            # create new data that will be used for plotting
            to_plot_data = { 'Time': data['Time'] }

            # handle <convert_to_time_unit>
            if 'convert_to_time_unit' in current_plot_dict.keys():
                if current_plot_dict['convert_to_time_unit'] == 's':
                    factor = 60
                elif current_plot_dict['convert_to_time_unit'] == 'min':
                    factor = 1
                elif current_plot_dict['convert_to_time_unit'] == 'h':
                    factor = 1/60
                elif current_plot_dict['convert_to_time_unit'] == 'd':
                    factor = 1/60 * 1/24
                else:
                    module_logger.error(f'Unknown <convert_to_time_unit> value: {convert_to_time_unit}')
                    exit(1)
            else:
                factor = 1
            # convert time
            to_plot_data['Time'] = [ t*factor for t in to_plot_data['Time'] ]

            # add internal data
            if 'species' not in current_plot_dict.keys():
                module_logger.error(f'While plotting the plot "{title}" I did not find <species>.')
                module_logger.error(f'Hint: If you do not need external data, set "species": []')
                exit(1)

            for s in current_plot_dict['species']:
                to_plot_data[s] = {'runs': data[s]['runs'] }

            # add external data
            if 'external_data' not in current_plot_dict.keys():
                module_logger.error(f'While plotting the plot "{title}" I did not find <external_data>.')
                module_logger.error(f'Hint: If you do not need external data, set "external_data": []')
                exit(1)

            for target in current_plot_dict['external_data']:
                if 'time_units' in target.keys():
                    time_units = target['time_units']
                else:
                    time_units = 'minutes'

                try:
                    target_data = parse.read_csv_target_file(target['file'], time_units=time_units)
                    for s in target['species']:
                        # add them with their own timebase

                        # look if should be called differently:
                        if 'with_file_prefix' in target.keys() and target['with_file_prefix']:
                            base = os.path.basename(target['file'])
                            callme = f"{base}_{s}"
                        else:
                            callme = s

                        # add runs
                        to_plot_data[callme] = { 'runs': target_data[s]['runs'],
                                                 'Time': [ t*factor for t in target_data['Time'] ] }

                except KeyError as e:
                    module_logger.error(f"plot_data: KeyError for {e}")
                    module_logger.error(f'Did not find plot <species> parameter "{s}" in target_data with keys={list(target_data.keys())}')
                    exit(1)


            # handle <multi>
            if 'multi' in current_plot_dict.keys():
                if current_plot_dict['multi']:
                    # do not plot it
                    break

            # handle <filename_type>
            if 'filename_type' in current_plot_dict.keys():
                filename_type = current_plot_dict['filename_type']
                module_logger.info(f'Using file extension(s): {filename_type}')
            else:
                filename_type = 'png'

            # handle <statistics>
            if 'statistics' in current_plot_dict.keys():
                stat =  current_plot_dict['statistics']
            else:
                stat = 'avg'

            # handle <write_data_into>
            # potentially write data into a json file
            if 'write_data_into' in current_plot_dict.keys():
                data_filename = current_plot_dict['write_data_into']
                with open(data_filename, 'w') as outfile:
                    module_logger.info(f'Writing data as requested by write_data_into. Dumping into: {data_filename}')
                    json.dump(to_plot_data, outfile, indent=2)

            # compute statistics:
            # data -> lines and shades
            lines_and_shades = stats.get_statistics(to_plot_data, statistics=stat)

            # set styles
            styles = current_plot_dict['styles']

            # set axis labels
            if 'xlabel' in current_plot_dict.keys():
                xlabel = current_plot_dict['xlabel']
            else:
                # use default
                xlabel = 'time (min)'

            if 'ylabel' in current_plot_dict.keys():
                ylabel = current_plot_dict['ylabel']
            else:
                # use defaults
                ylabel='species concentration [ml$^{-1}$]'

            # set figure size
            if 'figsize' in current_plot_dict.keys():
                figsize = current_plot_dict['figsize']
            else:
                # use default
                figsize = [8, 6]

            # set xlim
            if 'xlim' in current_plot_dict.keys():
                xlim = current_plot_dict['xlim']
            else:
                xlim = None

            # set xticks
            if 'xticks' in current_plot_dict.keys():
                xticks = current_plot_dict['xticks']
            else:
                xticks = None

            # set ylim
            if 'ylim' in current_plot_dict.keys():
                ylim = current_plot_dict['ylim']
            else:
                ylim = None

             # set yticks
            if 'yticks' in current_plot_dict.keys():
                yticks = current_plot_dict['yticks']
            else:
                yticks = None

            # handle <no_xaxis> and <no_yaxis>
            if 'no_xaxis' in current_plot_dict.keys():
                no_xaxis = current_plot_dict['no_xaxis']
            else:
                no_xaxis = False
            if 'no_yaxis' in current_plot_dict.keys():
                no_yaxis = current_plot_dict['no_yaxis']
            else:
                no_yaxis = False

            # create legend labels dictionary
            legend_labels = {}
            for s in lines_and_shades.keys():
                if s in current_plot_dict["styles"].keys():
                    # there is a style for it
                    if "name" in current_plot_dict["styles"][s]:
                        # use this name instead
                        legend_labels[s] = current_plot_dict["styles"][s]["name"]
                    else:
                        legend_labels[s] = s
                else:
                    legend_labels[s] = s

            # handle <no_legend>
            if 'no_legend' in current_plot_dict.keys():
                no_legend = current_plot_dict['no_legend']
            else:
                no_legend = False

            # handle <legend.bbox_to_anchor>
            if 'legend.bbox_to_anchor' in current_plot_dict.keys():
                legend_bbox_to_anchor = tuple(current_plot_dict['legend.bbox_to_anchor'])
            else:
                legend_bbox_to_anchor = None

            # handle <legend.ncol>
            if 'legend.ncol' in current_plot_dict.keys():
                legend_ncol = current_plot_dict['legend.ncol']
            else:
                legend_ncol = 1

            # handle <no_title>
            if 'no_title' in current_plot_dict.keys():
                no_title = current_plot_dict['no_title']
            else:
                no_title = False

            # handle <no_grid>
            if 'no_grid' in current_plot_dict.keys():
                no_grid = current_plot_dict['no_grid']
            else:
                no_grid = False

            # handle <title>
            if 'title' in current_plot_dict.keys():
                use_as_title = current_plot_dict['title']
            else:
                use_as_title = title

            # plot it
            default_plot(lines_and_shades,
                title=use_as_title,
                output_dir=output_dir,
                filename_type=filename_type,
                filename_base=title,
                styles=styles,
                figsize=figsize,
                save_fig=save_fig,
                ax=None,
                logscale=current_plot_dict['logscale'],
                xlabel=xlabel,
                ylabel=ylabel,
                legend_title=None,
                legend_labels=legend_labels,
                alpha=1.0,
                xlim=xlim,
                ylim=ylim,
                xticks=xticks,
                yticks=yticks,
                err_area=True,
                no_legend=no_legend,
                legend_bbox_to_anchor=legend_bbox_to_anchor,
                legend_ncol=legend_ncol,
                no_title=no_title,
                no_grid=no_grid,
                no_xaxis=no_xaxis,
                no_yaxis=no_yaxis)





########################## QUICK DEFAULT PLOTS #####################################

def default_plot(lines_and_shades,
    styles,
    figsize,
    title,
    output_dir,
    save_fig,
    filename_type='png',
    filename_base=None,
    xlim=None,
    ylim=None,
    ax=None,
    logscale=None,
    xlabel=None,
    xticks=None,
    yticks=None,
    ylabel=None,
    legend_title=None,
    legend_labels=None,
    alpha=1.0,
    err_area=True,
    no_legend=False,
    legend_bbox_to_anchor=None,
    legend_ncol=1,
    no_title=False,
    no_grid=False,
    no_xaxis=False,
    no_yaxis=False):

    # create new figure if does not exist
    if ax is None:
        plt.figure(1, figsize, dpi=1200)

    for i, current_species in enumerate(lines_and_shades.keys()):

        # check if special style exists for this species

        # standard style
        plot_type = 'line'
        color = COLORS[i%len(COLORS)]
        linestyle = LINESTYLES[i%len(LINESTYLES)]
        linewidth = 1
        marker = MARKERS[i%len(MARKERS)]
        dotsize = 10

        # update style with default
        if "default" in styles.keys():
            current_style = styles["default"]
            if 'plot_type' in current_style.keys():
                plot_type = current_style['plot_type']
            if 'color' in current_style.keys():
                color =  getColor(current_style['color'])
            if 'linestyle' in current_style.keys():
                linestyle = current_style['linestyle']
            if 'linewidth' in current_style.keys():
                linewidth = current_style['linewidth']
            if 'marker' in current_style.keys():
                marker = current_style['marker']
            if 'dotsize' in current_style.keys():
                dotsize = current_style['dotsize']

        # update style with particular style for this species
        if current_species in styles.keys():
            current_style = styles[current_species]
            if 'plot_type' in current_style.keys():
                plot_type = current_style['plot_type']
            if 'color' in current_style.keys():
                color = getColor(current_style['color'])
            if 'linestyle' in current_style.keys():
                linestyle = current_style['linestyle']
            if 'linewidth' in current_style.keys():
                linewidth = current_style['linewidth']
            if 'marker' in current_style.keys():
                marker = current_style['marker']
            if 'dotsize' in current_style.keys():
                dotsize = current_style['dotsize']

        plot_lines_and_shades(lines_and_shades[current_species],
            legend_label=legend_labels[current_species],
            plot_type=plot_type,
            linestyle=linestyle,
            marker=marker,
            color=color,
            ax=ax,
            index=i,
            logscale=logscale,
            linewidth=linewidth,
            dotsize=dotsize,
            err_area=err_area,
            alpha=alpha)

    default_finish(title=title,
        output_dir=output_dir,
        filename_type=filename_type,
        filename_base=filename_base,
        save_fig=save_fig,
        ax=ax,
        legend_title=legend_title,
        xlabel=xlabel,
        ylabel=ylabel,
        xticks=xticks,
        yticks=yticks,
        xlim=xlim,
        ylim=ylim,
        no_legend=no_legend,
        legend_bbox_to_anchor=legend_bbox_to_anchor,
        legend_ncol=legend_ncol,
        no_title=no_title,
        no_grid=no_grid,
        no_xaxis=no_xaxis,
        no_yaxis=no_yaxis)


def plot_lines_and_shades(lines_and_shades_for_this_species,
    legend_label,
    plot_type,
    linestyle,
    marker,
    color,
    ax=None,
    index=0,
    logscale='Y',
    alpha=1.0,
    shade_alpha=.2,
    linewidth=2,
    dotsize=10,
    err_area=True,
    scatter_with_edges=False):
    """
    The function that that the main plotting of lines and shades
    """
    las = lines_and_shades_for_this_species
    
    # checks
    if 'x' not in las.keys():
        module_logger.error('plot_lines_and_shades: <lines_and_shades_for_this_species> requires an <x> entry.')
        exit(1)
    if 'lines' not in las.keys():
        module_logger.error('plot_lines_and_shades: <lines_and_shades_for_this_species> requires a <lines> entry.')
        exit(1)
    if 'shades' not in las.keys():
        module_logger.error('plot_lines_and_shades: <lines_and_shades_for_this_species> requires a <shades> entry.')
        exit(1)

    # set axis
    if ax is None:
        ax = plt.gca()

    # plot lines
    first = True
    for line in las['lines']:
        plot_line(las['x'], line,
            legend_label=legend_label if first else '_nolegend_',
            plot_type=plot_type,
            ax=ax,
            alpha=alpha,
            linewidth=linewidth,
            index=index,
            logscale=logscale,
            linestyle=linestyle,
            dotsize=dotsize,
            marker=marker,
            color=color,
            scatter_with_edges=scatter_with_edges)
        first = False

    # plot shades
    x = np.array(las['x'])
    for shade in las['shades']:
        lower_y = shade[0]
        upper_y = shade[1]
        ax.fill_between(x, lower_y, upper_y,
            alpha=shade_alpha,
            color=color,
            label='_nolegend_')
    
    # set axis scales
    if logscale in ['Y','both','XY']:
        ax.set_yscale('log')
    if logscale in ['X','both','XY']:
        ax.set_xscale('log')


def plot_line(X, Y,
    legend_label,
    plot_type,
    color,
    marker,
    linestyle,
    ax=None,
    alpha=.9,
    linewidth=8,
    dotsize=10,
    index=0,
    logscale=None,
    scatter_with_edges=False):
    """
    Plots a single line
    """

    if plot_type not in ['scatter','both','line']:
        module_logger.error(f'The <plot_type> parameter was set to "{plot_type}".')
        module_logger.error('Currently only supported: scatter, line, both.')
        exit(1)

    if plot_type in ['line','both']:
        line = None
        if ax is None:
            if none(logscale):
                line, = plt.plot(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
            elif logscale=='Y':
                line, = plt.semilogy(X,Y,alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
            elif logscale=='X':
                line, = plt.semilogx(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
            elif logscale in ['both','XY']:
                line, = plt.loglog(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
        else:
            if none(logscale):
                line, = ax.plot(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
            elif logscale=='Y':
                line, = ax.semilogy(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
            elif logscale=='X':
                line, = ax.semilogx(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)
            elif logscale in ['both','XY']:
                line, = ax.loglog(X,Y, alpha=alpha, linewidth=linewidth, linestyle=linestyle, color=color, label=legend_label)

        if line is not None:
            global lines
            lines += [line]

    if plot_type in ['scatter','both']:
        if ax is None:
            ax = plt.gca()
        if scatter_with_edges:
            ax.scatter(X,Y, s=dotsize, marker=marker, alpha=alpha, color=color, linewidth=linewidth, edgecolors='black', label=legend_label)
        else:
            ax.scatter(X,Y, s=dotsize, marker=marker, alpha=alpha, color=color, label=legend_label)

        if logscale in ['Y','both','XY']:
            ax.set_yscale('log')
        if logscale in ['X','both','XY']:
            ax.set_xscale('log')


def default_finish(
    title,
    output_dir,
    filename_base=None,
    filename_type='png',
    save_fig=False,
    xlim=None,
    ylim=None,
    ax=None,
    xlabel=None,
    ylabel=None,
    xticks=None,
    yticks=None,
    legend_title=None,
    sep_save_title=None,
    no_legend=False,
    legend_bbox_to_anchor=None,
    legend_ncol=1,
    no_title=False,
    no_grid=False,
    no_xaxis=False,
    no_yaxis=False):
    """
    adjusts appearance like font sizes and
    potentially saves the figure to a file
    """

    global legend

    if ax is None:

        # xlim
        if xlim is not None:
            plt.gca().set_xlim(left=xlim[0], right=xlim[1])

        # ylim
        if ylim is not None:
            plt.gca().set_ylim(bottom=ylim[0], top=ylim[1])

        # grid
        if not no_grid:
            plt.grid(b=True, alpha=.4, which='major')
            plt.grid(b=True, alpha=.3, which='minor')

        # legend
        if not no_legend:
            # placement of legend (legend.bbox_to_anchor)
            if legend_bbox_to_anchor is not None:
                legend = plt.legend(title=legend_title,
                                    frameon=False,
                                    bbox_to_anchor=legend_bbox_to_anchor,
                                    ncol=legend_ncol,
                                    loc='upper left')
            else:
                legend = plt.legend(title=legend_title,
                                    frameon=False)
        else:
            module_logger.info(f'No legend as requested by <no_legend>')
            legend = None

        # x/y-labels
        plt.xlabel(xlabel)
        if xticks is None:
            plt.xticks()
        else:
            plt.gca().set_xticks(xticks)

        plt.ylabel(ylabel)
        
        if yticks is None:
            plt.yticks()
        else:
            plt.gca().set_yticks(yticks)

        # potentially remove ticks and labels on axis
        if no_xaxis:
            plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False) # labels along the bottom edge are off
        if no_yaxis:
            ax = plt.gca()
            ax.yaxis.set_ticklabels([])
            # plt.tick_params(
            #     axis='y',
            #     which='both',
            #     left=False,
            #     right=False,
            #     labelleft=False)

        # set title
        if not no_title:
            plt.title(title)

        if save_fig:
            # set filename base
            if filename_base is None:
                filename_base = title

            # save (potentially several filenames)
            if type(filename_type) == str:
                filename_type_list = [ filename_type ]
            elif type(filename_type) == list:
                filename_type_list = filename_type
            else:
                module_logger.error(f'Expecting either a str or list of str as filename_type, but got this: {filename_type}')
                exit(1)

            for filename_type in filename_type_list:
                filename = output_dir + f'{filename_base}.{filename_type}'
                try:
                    plt.savefig(filename, bbox_inches='tight')
                    module_logger.info(f'Figure saved at: {filename}')

                except FileNotFoundError as e:
                    module_logger.error(f'FileNotFoundError: {e}')
                    module_logger.error(f'I could not write to the directory of the file "{filename}".')
                    module_logger.error('Hint: create the directory or change the <output_dir> parameter.')
                    exit(1)

        else:
            if legend is not None:
                start_picking()
            plt.show()

        plt.clf()
        plt.close()

    else: #subplot

        # xlim
        if xlim is not None:
            ax.set_xlim(left=xlim[0], right=xlim[1])

        # ylim
        if ylim is not None:
            ax.set_ylim(bottom=ylim[0], top=ylim[1])

        # grid
        ax.grid(b=True, alpha=.4, which='major')
        ax.grid(b=True, alpha=.1, which='minor')

        # legend
        lgd = ax.legend(title=legend_title)
        lgd.get_title().set_fontsize()

        # x/y labels
        ax.set_xlabel(xlabel)
        ax.tick_params(axis='both', which='major')
        if xticks is not None:
            ax.set_xticks(xticks[0])
            ax.set_xticklabels(xticks[1])
        ax.set_ylabel(ylabel)

        # title
        if title is not None:
            ax.set_title(title)




# ##################################################################################################

if __name__ == "__main__":
    module_logger = log.setup_logger('root')

    parser = argparse.ArgumentParser(description='Plotting of executions.')
    parser.add_argument('PICKLE_FILE',
        nargs=1,
        type=str,
        help='the pickle file')
    parser.add_argument('--multi',
        action=argparse.BooleanOptionalAction,
        help='run multiple simulations')
    parser.add_argument('--opt',
        action=argparse.BooleanOptionalAction,
        help='optimize')
    parser.add_argument('--printplotparams',
        action=argparse.BooleanOptionalAction,
        help='also prints the plot parameters')
    parser.add_argument('--printparams',
        action=argparse.BooleanOptionalAction,
        help='also prints the simulation parameters')
    parser.add_argument('--printdata',
        action=argparse.BooleanOptionalAction,
        help='also prints the simulation parameters')
    args = parser.parse_args()
    pickle_file = args.PICKLE_FILE[0]

    if args.printplotparams:
        with open(pickle_file, 'rb') as pickle_file_handle:
            payload = pickle.load(pickle_file_handle)
            plot_params = payload['plot_params']
            print('-'*10 + ' plot parameters ' + '-'*10)
            json.dump(plot_params, sys.stdout, indent=4, cls=NumpyEncoder)
            print()

    if args.printparams:
        with open(pickle_file, 'rb') as pickle_file_handle:
            payload = pickle.load(pickle_file_handle)
            params = payload['params']
            print('-'*10 + ' parameters ' + '-'*10)
            pprint.pprint(params)
            print()

    if args.printdata:
        with open(pickle_file, 'rb') as pickle_file_handle:
            payload = pickle.load(pickle_file_handle)
            data = payload['data']
            print('-'*10 + ' data ' + '-'*10)
            pprint.pprint(data)
            print()

    if args.opt and args.multi:
        module_logger.error("Flags --opt and --multi cannot be used together.")
        exit(1)

    if not os.path.isfile(pickle_file):
        module_logger.error("Cannot find file to plot (1st arg), check its path.")
        exit(1)

    if args.opt:
        from_pickle_opt(pickle_file)
        exit(0)

    if args.multi:
        from_pickle(pickle_file, multiplot=True)
    else:
        from_pickle(pickle_file, multiplot=False)
