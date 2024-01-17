import sys, os, h5py
import numpy as np
from matplotlib import pyplot as plt
from pycbc.catalog import Catalog
from pycbc.catalog import Merger
from pprint import pprint
from pycbc.types.timeseries import *
from . import gw_cwt

os.environ['LAL_DATA_PATH'] = '/scratch/lalsimulation'

fig_width_pt = 1020.
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width/golden_mean       # height in inches
fig_size = [fig_width,fig_height]

fontsize = 16
legendfontsize = 14

plot_params={'text.usetex': False,
        'axes.labelsize': fontsize,
        'font.size': fontsize,
        'legend.fontsize': legendfontsize,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'figure.figsize': fig_size,
        'font.weight': 'normal'
       }


import pylab
pylab.rcParams.update(plot_params)
pylab.rcParams['axes.linewidth'] = 1
pylab.rc('axes', linewidth=1)


def touchbox(ax):
    ax.minorticks_on()
    ax.tick_params('both', length=5, width=1, which='major')
    ax.tick_params('both', length=3.5, width=1, which='minor')
    ax.spines['top'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    return

# mkdir function to create directory for saving
def mkdir(dir_name):
    if os.path.exists(dir_name):
        return False
    temp_path = ''
    for i in range(len(dir_name)):
        if dir_name[i] == "/":
            try:
                os.mkdir(temp_path)
            except:
                pass
        temp_path += dir_name[i]
    return True
            
# mergerDict - get data and return it in a dictionary
# takes in a list of the catalogs from which to retrieve events, and a list of detectors to get information from.
# returns a dictionary mapping:
# 'name' (str) - the name of the event
# 'merger_time' (TimeSeries object) - the time of the merger
# 'strain_info' - a dictionary mapping names of detectors to the whitened data on the merger
# 'catalog' - the name of the catalog where the event is located
def mergerDict(catalogs, ifos):
    print('Retrieving data...')
    numDetections = 0
    data = {}
    print('Loading events...')
    for c in catalogs:
        print(f'Retrieving data from {c}...')
        for e in Catalog(c): # Catalog was imported from pycbc.catalog (see preamble)
            m = Merger(e, c)
            numDetections += 1
            strain_info = {}
            for ifo in ifos:
                try:
                    # add whitened data to strain_info dictionary for event
                    strain_info[ifo] = m.strain(ifo).whiten(4,4)
                    try:
                        e = e[:e.index("-")]
                    except:
                        pass
                except:
                    pass
            data[e] = {'name': e, 'merger_time': m.time, 'strain_info': strain_info, 'catalog': c}
    print(f'Done! retrieved data for {numDetections} detections.')
    return data

# function to create merger dictionary for single event
def oneEventMergerDict(eventName, ifos):
    print('Retrieving data...')
    data = {}
    strain_info = {}
    found = False
    for c in ['GWTC-1-confident', 'GWTC-2.1-confident', 'GWTC-3-confident']:
        if found:
            break
        for e in Catalog(c).names:
            if eventName in e:
                if eventName != e:
                    print(f'Full name for event {eventName} not entered. Set it to {e}.')
                    eventName = e
                catalog = c
                found = True
    m = Merger(eventName, catalog)
    for ifo in ifos:
        try:
            # add whitened data to strain_info dictionary for event
            strain_info[ifo] = m.strain(ifo).whiten(4,4)
            if "-" in e:
                e = e[:e.index("-")]
        except:
            pass
    data[e] = {'name': eventName, 'merger_time': m.time, 'strain_info': strain_info, 'catalog': c}
    print('Done retrieving data!')
    return data

# retrieves data from specified catalogs or event and adds to an h5py file cwt_data located at the directory specified
# directory is your current working directory if not specified
# returns the h5file path for use with running CWT

def get_data(catalogs=None, ifos=['H1', 'L1', 'V1'], path='', event=None):
    if path:
        if not mkdir(path):
            sys.exit('Path already exists, danger of overwriting existing data.')
        if path[-1] == '/':
            h5file = path + 'cwt_data'
        else:
            h5file = path + '/' + 'cwt_data'
    else:
        h5file = 'cwt_data'
    if os.path.exists(h5file):
        sys.exit('h5py file already exists. Make sure to create a new one to avoid overwriting existing data! Try specifying a new path.')
    if not catalogs:
        if not event:
            sys.exit('did not specify catalog or event.')
        else:
            data = oneEventMergerDict(event, ifos)
    else:
        data = mergerDict(catalogs, ifos)
    mkh5(data, h5file)
    if path:
        print(f'Done! Your h5py file is located at {path}')
    else:
        print('Done! Your h5py file is located in this directory and is called cwt_data.')
    return h5file

# save whitened data (timeseries.data), times = s_data.sample_times (timestamps for raw strains), timeslice, 
# timeslice.data, timeslice.sample_times
# in the format: {h5py file name}/event/{detector/strain_info/sliced/times etc, other keys}/data
def mkh5(data, path):
    print('Starting to make h5py file')
    try:
        h5f = h5py.File(path, 'w')
    except:
        sys.exit('h5py file already exists. Make sure to create a new one to avoid overwriting existing data!')
    for event in data:
        event = data[event]
        time = event['merger_time']
        name = event['name']
        for key in event:
            if type(event[key]) == list:
                h5f.create_dataset(f'{name}/{key}', data=event[key])
            elif key == 'strain_info':
                for detector in event[key]:
                    path = f'{name}/strain_info/{detector}'
                    timeseries = event[key][detector]
                    odata = timeseries.data #unsliced data
                    otimes = timeseries.sample_times #unsliced data times
                    timeslice = timeseries.time_slice(time - 1, time + 0.5) #slicing it
                    timeslice_times = timeslice.sample_times #sliced data times
                    timeslice_data = timeslice.data #sliced data
                    h5f.create_dataset(f'{path}/unsliced/data', data = odata)
                    h5f.create_dataset(f'{path}/unsliced/times', data = otimes)
                    h5f.create_dataset(f'{path}/sliced/data', data = timeslice_data)
                    h5f.create_dataset(f'{path}/sliced/times', data = timeslice_times)
            else:
                h5f.create_dataset(f'{name}/{key}', data = [event[key]]) # if string, it will be byte and must be decoded
                # use .decode('utf-8')
    h5f.close()
    print('h5py file complete')

# cwt - to run it
# takes in an h5py file path for the event, returned by mkh5() function, and parameters used for build_cwt_dev.
# takes in optional list of events to run CWT on if not all desired
# TODO add better documentation
# returns a dictionary mapping each event to a list of tuples for each detector with the detector name, data with CWT run 
# on it, the specific slice of data on which CWT was run, the times corresponding to the data in a TimeSeries object, the 
# top of the scale range, and a tuple with Q and df.
def run_it(h5file, Q=6.0, chirp_rate=0.0, f_range=(10.0, 500.0), freq_spacing='Log', n_conv=400, df=None, da=None, f_list=None, Norm=True, events=None):
    print('Starting CWT')
    returnDict = {}
    file = h5py.File(h5file, 'r')
    for event in file.keys():
        if events and event not in events:
            continue
        addList = []
        mergerDict = file[event]
        time = mergerDict['merger_time']
        name = mergerDict['name'][0].decode('utf-8')
        catalog = mergerDict['catalog'][0].decode('utf-8')
        returnDict[name] = {'catalog': catalog}
        for ifo in file[event]['strain_info'].keys():
            print(f'running CWT on {name}, detector {ifo}...')
            signal = np.array(file[event]['strain_info'][ifo]['sliced']['data'])
            times = np.array(file[event]['strain_info'][ifo]['sliced']['times'])
            cwt_result = gw_cwt.build_cwt(signal, times, Q=Q, chirp_rate=chirp_rate, f_range=f_range, freq_spacing=freq_spacing, n_conv=n_conv, df=df, da=da, f_list=f_list, Norm=Norm)
            print(f'CWT done for {name}, detector {ifo}')
            times = cwt_result['times']
            addList.append((ifo, cwt_result, signal, times, (Q, chirp_rate)))
        returnDict[name]['cwt_result'] = addList
    print('Done running cwt, saving data...')
    file.close()
    return returnDict

# save map, frequencies, scales
def cwt_savedata(returnDict, h5file):
    file = h5py.File(h5file, 'a')
    for event in returnDict:
        for tup in returnDict[event]['cwt_result']:
            result = tup[1]
            detector = tup[0]
            path = f'{event}/cwt_result/{detector}'
            file.create_dataset(f'{path}/map', data = result['map'])
            file.create_dataset(f'{path}/frequencies', data = result['frequencies'])
            file.create_dataset(f'{path}/scales', data = result['scales'])
    file.close()
    print('Data saved')
    
# Runs cwt and saves to an h5py file.
def run_cwt(h5file, Q=6.0, chirp_rate=0.0, path='', events=None, f_range=(10.0, 500.0), freq_spacing='Log', n_conv=400, df=None, da=None, f_list=None, Norm=True):
    cwt_dict = run_it(h5file, Q=Q, f_range=f_range, freq_spacing=freq_spacing, n_conv=n_conv, df=df, da=da, f_list=f_list, Norm=Norm, events=events)
    cwt_savedata(cwt_dict, h5file)
    print(f'Done! Data saved to {h5file}')

# function to plot CWT results
def plot_cwt(h5file, savefig=True, savefigpath='', freq_to_plot=(10, 300)):
    f = h5py.File(h5file, 'r')
    if savefigpath:
        if savefigpath[-1] != "/":
            savefigpath += "/"
        if not mkdir(savefigpath):
            sys.exit(f'Error with creating directory. Path already exists! Path: {savefigpath}')
    for event in f.keys():
        if 'cwt_result' not in f[event].keys():
            continue
        for detector in f[event]['cwt_result'].keys():
            cwt_result = f[event]['cwt_result'][detector]
            wfreqs = np.array(cwt_result['frequencies'])
            fmap = np.array(cwt_result['map'])
            scales = np.array(cwt_result['scales'])
            times = np.array(f[event]['strain_info'][detector]['sliced']['times'])
            signal = np.array(f[event]['strain_info'][detector]['sliced']['data'])
            merger_time = float(np.array(f[event]['merger_time'][0]))
            fig_width = 8.0 #in
            golden_ratio = 1.61803398875
            fig_height = fig_width / golden_ratio
            plt.close('all')
            fig, ax = plt.subplots(figsize=(fig_width, fig_height),sharex=True, ncols=1, nrows=1)
            fig.patch.set_facecolor('white')
            # shifting times so merger is at t = 0 seconds
            times = times - merger_time
            plot_domain = (-1, 0.5)
            #
            # Spectrograms
            #
            spec = ax.pcolormesh(times, wfreqs, fmap, rasterized=False, vmin=0, vmax=1, shading='auto', cmap='viridis')
            ax.set_ylabel('Frequency [Hz]', fontsize=10)
            ymin, ymax = freq_to_plot
            ax.set_ylim(ymin, ymax)
            ax.tick_params(axis='y',labelsize=10)
            ax.set_xlabel('Time from merger [s]', fontsize=10)
            ax.set_xlim(plot_domain)
            ax.tick_params(axis='x',labelsize=10)
            touchbox(ax)
            clb = fig.colorbar(spec, ax=ax, shrink=0.95, pad=0.05)
            clb.ax.set_title(r'$T\left(t, f\right)$', fontsize=10, pad=5)
            clb.ax.tick_params(labelsize=10)
            if savefig:
                print('Saving figure')
                mkdir(savefigpath + detector + "/")
                plt.savefig(f'{savefigpath}{detector}/{event}.png', format='png', bbox_inches='tight')
                print(f'Saved figure for merger {event} at {detector}.')
            else:
                print(f'Done with merger {event} at {detector}.')
                plt.show()
    print('Done!')
    f.close()
    
