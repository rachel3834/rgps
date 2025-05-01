from os import path, getcwd
from sys import path as pythonpath
pythonpath.append(path.join(getcwd(), '..'))
from astropy.table import Table, Column, vstack
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import config_utils
import regions
import visualization_utils
import numpy as np
from astropy import units as u
from mw_plot import MWSkyMap


def plot_metric_optic_heatmap(sim_config, metric_results, strategy_name, category,
                              metric_name, metric_label, file_path):
    """
    Function to plot a heatmap of a set of metric results, comparing the all optics for a category of
    science cases.
    """

    # Downselect the metric results table for the entries corresponding to the requested strategy
    idx = metric_results['Survey_strategy'] == strategy_name
    metric_filter = metric_results[idx]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 20))

    # Use only those science cases where metric results were calculated for this filter
    # This excludes science cases that did not request a given filter.
    science_cases = list(set(metric_filter['Science_case'].data))
    science_cases.sort()

    data = np.zeros((len(science_cases), len(sim_config['OPTICAL_COMPONENTS'])))

    # Select all metric results for the current filter
    for ioptic, optic in enumerate(sim_config['OPTICAL_COMPONENTS']):
        jdx = np.where(metric_filter['Optic'] == optic)[0]

        for iscience, science_name in enumerate(science_cases):
            kdx = np.where(metric_filter[jdx]['Science_case'] == science_name)[0]

            # Not all science cases request all filters, so it is possible to have no results here
            if len(metric_filter[jdx][kdx]) > 0:
                data[iscience, ioptic] = metric_filter[jdx][kdx][metric_name].data
            else:
                data[iscience, ioptic] = np.nan

    # The plotted grids for a heat map have to account for Python's indexing
    xgrid = np.arange(0, len(sim_config['OPTICAL_COMPONENTS']) + 1, 1)
    ygrid = np.arange(0, len(science_cases) + 1, 1)

    # Plot normalised metric data
    norm = mpl.colors.Normalize(0.0, 100.0)
    ax.pcolormesh(xgrid, ygrid, data, cmap="magma", norm=norm)

    # Label axes
    ax.set_frame_on(False)

    ax.set_xticks(xgrid[0:-1] + 0.5)
    ax.set_yticks(ygrid[0:-1] + 0.5)
    ax.set_ylabel('Science case', fontsize=20)
    ax.set_xlabel('Optic', fontsize=20)
    ax.set_title(strategy_name + ' survey strategy and science category ' + category, fontsize=20)
    ax.set_xticklabels(sim_config['OPTICAL_COMPONENTS'], rotation=45.0, horizontalalignment='right', fontsize=20)
    ax.set_yticklabels(science_cases, fontsize=20, horizontalalignment='right')

    cb = fig.colorbar(
        ScalarMappable(norm=norm, cmap="magma"),
        ax=ax,  # Pass the new axis
        orientation="vertical")
    cb.set_label(metric_label, fontsize=20)
    cb.ax.tick_params(labelsize=18)

    plt.tight_layout()
    plt.savefig(file_path)
    plt.close(fig)


def plot_metric_multifilter(results, metric_name, sim_config, metric_label, science_name='all'):
    """
    Function to plot metric values for a given science use case over multiple survey designs, for multiple filters.

    Input
        selected_results  Table   Metric results table selecting the output for a specific science case only.
    """

    # Select the results for the specified science case
    if science_name != 'all':
        idx = np.where(results['Science_case'] == science_name)
        selected_results = results[idx]
    else:
        selected_results = results

    # Identify which filters were requested for this science case
    filter_set = list(set(selected_results['Optic'].data))
    filter_set.sort()

    # Make a list of the survey strategies included
    surveys_options = list(set(selected_results['Survey_strategy'].data))
    survey_options.sort()

    # Plot the metric results
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

    for optic in filter_set:

        fdx = selected_results['Optic'] == optic
        metric_filter = selected_results[fdx]

        data = []
        for survey_name in survey_options:
            jdx = np.where(metric_filter['Survey_strategy'] == survey_name)[0]
            data.append(metric_filter[jdx][metric_name][0])

        ax.plot(range(0, len(survey_options), 1), data,
                color=sim_config['PLOT_COLORS'][optic], marker=sim_config['PLOT_SYMBOLS'][optic],
                label=optic)
        ax.plot(range(0, len(survey_options), 1), data, color=sim_config['PLOT_COLORS'][optic], linestyle='-')

    ax.set_xlabel('Survey strategy', fontsize=20)
    ax.set_ylabel(metric_label, fontsize=20)
    ax.set_title(metric_name + ' for ' + science_name + ' catalog', fontsize=20)
    ax.set_xticks(range(0, len(survey_options), 1))
    ax.set_xticklabels(survey_options, rotation=45.0, horizontalalignment='right', fontsize=20)
    yticks = ax.get_yticks()
    yticklabels = ax.get_yticklabels()
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=20, horizontalalignment='right')
    ax.legend()
    ax.grid()

    plt.tight_layout()
    plt.savefig(path.join(sim_config['root_dir'], 'metric_results', 'm3_results.png'))