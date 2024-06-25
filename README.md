# gw_cwt
This package is used to create time-frequency maps using the continuous wavelet transform, and is tuned for analysis of gravitational wave data. If used for analysis, please cite the accompanying paper: https://arxiv.org/abs/2402.16533.

## Installation

To install `gw_cwt`, you can follow these steps:

1. **Clone the Git Repository**:
   - Open your terminal or command prompt.
   - Navigate to the directory where you want to install `gw_cwt`.
   - Clone the Git repository by running the following command:
     ```
     git clone https://github.com/chadhenshaw/gw_cwt.git
     ```
   - This will create a local copy of the `gw_cwt` package on your machine.

2. **Install the Package**:
   - Navigate into the `gw_cwt` directory:
     ```
     cd gw_cwt
     ```
   - Install the package using `pip` with the following command:
     ```
     pip install .
     ```
   - The `pip install .` command will install the package locally from the cloned repository.

3. **Verify Installation**:
   - You can verify that `gw_cwt` is installed correctly by importing it in Python:
     ```python
     import gw_cwt
     ```
     
## Usage and Features

### Primary Use: Continuous Wavelet Transform for Gravitational Wave Data

The primary purpose of the `gw_cwt` package is to create time-frequency maps of gravitational wave data using the continuous wavelet transform (CWT). This is achieved through the `build_cwt` function, which enables users to analyze gravitational wave signals using various options. See the paper (link here) and Examples notebook for more.

```python
import gw_cwt

# Example usage of build_cwt
timefreq_map = gw_cwt.build_cwt(timeseries, timestamps, Q=8.0, chirp_rate=0.0, f_range=(10.0, 500.0), freq_spacing='Log', n_conv=400, df=None, da=None, f_list=None, Norm=True)

# timeseries: Input time series data for analysis.
# timestamps: Timestamps corresponding to the time series data.
# Q: Quality factor, defines the mother frequency and relative resolution.
# chirp_rate: Chirp rate for the Chirplet wavelet. The default is zero (Morlet-Gabor wavelet).
# f_range: Frequency range (min, max) for the transform.
# freq_spacing: Method for spacing frequencies ('Log' or 'Linear').
# n_conv: Number of convolutions - i.e. frequencies within f_range.
# df: Frequency increment for linear spacing (optional).
# da: Scale increment for log spacing (optional).
# f_list: Custom list of frequencies (optional).
# Norm: Flag to normalize the time-frequency map (optional).

```

A secondary feature is the `cwt_catalog` script, which provides a set of functions to retrieve data from catalogs of gravitational wave events, run CWT on the data, and save the results in an HDF5 file. Additionally, it provides functions for plotting CWT scalograms.

```python

import gw_cwt.cwt_catalog

# Example usage of cwt_catalog functions

#Retrieves data from specified catalogs or a single event and adds it to an HDF5 file named cwt_data.
cwt_catalog.get_data(catalogs=None, ifos=['H1', 'L1', 'V1'], path='', event=None)
#catalogs: A list of catalogs from which to retrieve events.
#ifos: A list of detectors to get information from.
#path: The directory where the HDF5 file will be saved (default is current working directory).
#event: The name of a specific event to retrieve.

# wrapper for running build_cwt on catalog events. See above for further option details.
cwt_catalog.run_cwt(h5file, Q=6.0, chirp_rate=0.0, path='', events=None, f_range=(10.0, 500.0), freq_spacing='Log', n_conv=400, df=None, da=None, f_list=None, Norm=True)
# h5file: path to h5 file containing event data retrieved with get_data.
# path: The directory where the CWT results will be saved (default is current working directory).
# events: A list of specific events to run CWT on.

# Function for plotting CWT scalograms. 
cwt_catalog.plot_cwt(h5file, savefig=True, savefigpath='', freq_to_plot=(10, 300))
# h5file: path to h5 file containing event data retrieved with get_data and evaluated with run_cwt.
# savefigpath: The directory where the figures will be saved (default is current working directory).
# freq_to_plot: Tuple specifying the frequency range to plot.

```


## Dependencies

- Python 3.6 or higher
- NumPy
- SciPy
- Matplotlib
- h5py
- PyCBC

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Authors

- Chad Henshaw
- Megan Arogeti
- Alice Heranval
- Laura Cadonati

## Contact

For questions, suggestions, or support, you can reach out to us at [cgh3@gatech.edu] or visit our [GitHub Profile](https://github.com/chadhenshaw).
