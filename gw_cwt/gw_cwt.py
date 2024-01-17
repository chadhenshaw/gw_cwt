# Continuous wavelet transform analysis for gravitational waves data
# Modification of pyCWT code (https://github.com/Unidata/pyCWT)
# Author: Chad Henshaw (cgh3@gatech.edu)
# Last modification: 12/20/23

import numpy as np
from scipy.fftpack import fft, ifft, fftshift
import warnings

__all__ = ['cwt', 'Chirplet']


class Chirplet:
    """
    Represents a Chirplet (a modified Morlet-Gabor wavelet) used for continuous wavelet transform.

    Attributes:
        sampf (float): Sampling frequency of the signal.
        scales (ndarray): Array of scales used in the wavelet transform.
        len_signal (int): Length of the input signal.
        name (str): Name of the wavelet, set to 'Chirplet'.
        len_wavelet (int): Total length of the wavelet, considering zero padding.
        fc (float): Characteristic frequency of the wavelet.
        chirp_rate (float): Chirp rate of the wavelet.
        coefs (ndarray): Computed coefficients of the dilated wavelet.
    """

    def __init__(self, len_signal=None, pad_to=None, scales=None, sampf=1, f0=0.849, d=0.00):
        """
        Initializes the Chirplet wavelet with specified parameters.

        Parameters:
            len_signal (int): Length of the signal to be analyzed.
            pad_to (int, optional): Length to which the signal is zero-padded.
            scales (ndarray): Scales at which the wavelet transform is computed.
            sampf (float): Sampling frequency of the signal.
            f0 (float): Characteristic frequency of the Chirplet.
            d (float): Chirp rate.
        """
        
        self.sampf = sampf
        self.scales = scales
        self.len_signal = len_signal
        self.name = 'Chirplet'

        # set total length of wavelet to account for zero padding
        if pad_to is None:
            self.len_wavelet = len_signal
        else:
            self.len_wavelet = pad_to

        # define characteristic frequency
        self.fc = f0
        
        # define chirp rate
        self.chirp_rate = d

        # compute coefficients for the dilated wavelet
        self.coefs = self.get_coefs()

    def get_coefs(self):
        """
        Calculates the coefficients for the mother Chirplet.

        Returns:
            ndarray: Computed wavelet coefficients.
        """

        # Create array containing values used to evaluate the wavelet function
        xi=np.arange(-self.len_wavelet / 2., self.len_wavelet / 2.)

        # find mother wavelet coefficients at each scale
        xsd = xi / (self.scales[:,np.newaxis])
        
        # calculate
        amp = np.power(np.pi,-0.25)
        sine = np.exp(complex(1j) * 2. * np.pi * (self.fc * xsd + 0.5*self.chirp_rate*xsd**2))
        gauss = np.exp(-np.power(xsd, 2) / 2.)
        cp = amp*sine*gauss

        self.coefs = cp

        return cp

class Wavelet(object):
    """
    Represents the result of a continuous wavelet transform.

    Attributes:
        coefs (ndarray): Coefficients of the wavelet transform.
        _pad_coefs (ndarray, optional): Coefficients corresponding to the padded part of the signal.
        motherwavelet (Chirplet): Instance of the Chirplet wavelet used for the transform.
        weighting_function (function): Function used to weight the wavelet transform.
        _signal_dtype (dtype): Data type of the original signal.
    """
    
    def __init__(self, wt, wavelet, weighting_function, signal_dtype, deep_copy=True):
        """
        Initializes the Wavelet object with the transform results and wavelet information.

        Parameters:
            wt (ndarray): Wavelet transform results.
            wavelet (Chirplet): Chirplet wavelet used in the transform.
            weighting_function (function): Function to weight the wavelet transform.
            signal_dtype (dtype): Data type of the original signal.
            deep_copy (bool): Whether to deep copy the wavelet instance or not.
        """

        from copy import deepcopy
        self.coefs = wt[:,0:wavelet.len_signal]

        if wavelet.len_signal !=  wavelet.len_wavelet:
            self._pad_coefs = wt[:,wavelet.len_signal:]
        else:
            self._pad_coefs = None
        if deep_copy:
            self.motherwavelet = deepcopy(wavelet)
        else:
            self.motherwavelet = wavelet

        self.weighting_function = weighting_function
        self._signal_dtype = signal_dtype

def cwt(x, wavelet, weighting_function=lambda x: x**(-0.5), deep_copy=True):
    
    """
    Performs the continuous wavelet transform on a given signal.

    Parameters:
        x (ndarray): Input signal.
        wavelet (Chirplet): Wavelet to be used for the transform.
        weighting_function (function): Weighting function for the transform.
        deep_copy (bool): Indicates whether to deep copy the wavelet instance.

    Returns:
        Wavelet: Object containing the results of the wavelet transform.
    """
    
    signal_dtype = x.dtype

    if len(x) < wavelet.len_wavelet:
        n = len(x)
        x = np.resize(x, (wavelet.len_wavelet,))
        x[n:] = 0

    # Transform the signal and mother wavelet into the Fourier domain
    xf=fft(x)
    mwf=fft(wavelet.coefs.conj(), axis=1)

    # Convolve (multiply in Fourier space)
    wt_tmp=ifft(mwf*xf[np.newaxis,:], axis=1)

    # shift output from ifft and multiply by weighting function
    wt = fftshift(wt_tmp,axes=[1]) * weighting_function(wavelet.scales[:, np.newaxis])

    # if mother wavelet and signal are real, only keep real part of transform
    wt=wt.astype(np.lib.common_type(wavelet.coefs, x))

    return Wavelet(wt,wavelet,weighting_function,signal_dtype,deep_copy)

def build_cwt(timeseries, timestamps, Q=8.0, chirp_rate=0.0, f_range=(10.0, 500.0), freq_spacing='Log', n_conv=400, df=None, da=None, f_list=None, Norm=True):
    
    """
    Builds a continuous wavelet transform time-frequency map for a given time series.

    Parameters:
        timeseries (ndarray): Time series data to be analyzed.
        timestamps (ndarray): Timestamps corresponding to the time series data.
        Q (float): Quality factor defining the mother frequency.
        chirp_rate (float): Chirp rate for the Chirplet wavelet.
        f_range (tuple): Frequency range for the transform (min, max).
        freq_spacing (str): Method for spacing frequencies ('Log' or 'Linear').
        n_conv (int): Number of convolutions in frequency spacing.
        df (float, optional): Frequency increment for linear spacing.
        da (float, optional): Scale increment for linear spacing.
        f_list (list, optional): Custom list of frequencies.
        Norm (bool): Flag to normalize the time-frequency map.

    Returns:
        dict: Dictionary containing the results of the wavelet transform.
    """

    delta_t = np.diff(timestamps)[0]
    sample_rate = 1./delta_t
    mother_freq= Q/(2*np.sqrt(2)*np.pi) # fiduciary frequency of mother wavelet
    
    # protection against infinite scale, and low frequency error
    if float(f_range[0]) < 0.4:
        warnings.warn("Wavelets below f=0.4 Hz are less accurate! Setting lower limit to 0.4 Hz.")
        f_range = (0.4, f_range[-1])
    # can't go above the nyquist frequency
    if float(f_range[-1]) > float(sample_rate/2):
        warnings.warn("Nyquist limit of %.1f exceeded, setting upper frequency range limit to the Nyquist frequency." % float(sample_rate/2))
        frange = (f_range[0], float(sample_rate/2))
    
    if freq_spacing == 'Log':
        # convert frequency to wavelet scale
        # calc max_scale for f_min
        max_scale = mother_freq*sample_rate / f_range[0]
        # calc min_scale for f_max
        min_scale = mother_freq*sample_rate / f_range[-1]
        if da is not None:
            # Assign linear scale spacing
            print('Forcing scale spacing to da = %.4f' % da)
            scales = np.arange(min_scale, max_scale, da)
        else:
            # space by number of convolutions
            scales = np.linspace(min_scale, max_scale, n_conv)
            
    elif freq_spacing == 'Linear':
        if df is not None:
            # Assign linear frequency spacing
            print('Forcing freq spacing to df = %.4f' % df)
            freqs = np.arange(f_range[0], f_range[-1], df)
            scales = mother_freq*sample_rate / freqs
        else:
            # space by number of convolutions
            freqs = np.linspace(f_range[0], f_range[-1], n_conv)
            scales = mother_freq*sample_rate / freqs
            
    elif f_list is not None:
        # Read in custom list of frequencies
        print('Using custom frequency distribution')
        freqs = np.array(f_list)
        scales = mother_freq*sample_rate / freqs

    # Construct the mother wavelet
    mother_wavelet = Chirplet(len_signal=len(timeseries), scales=scales,
        sampf=delta_t, f0=mother_freq, d=chirp_rate)    

    # Compute the CWT
    wavelet = cwt(timeseries, mother_wavelet)

    # Take the absolute values of coefficients - returns magnitude
    tfmap = np.abs(wavelet.coefs)

    # Normalise
    if Norm:
        tfmap /= max(map(max,abs(wavelet.coefs)))

    # Determine frequency at each scale
    freqs = sample_rate * wavelet.motherwavelet.fc \
            / wavelet.motherwavelet.scales

    # Return a dictionary
    timefreq = dict()
    timefreq['analysed_data'] = timeseries
    timefreq['map'] = tfmap # note that this is just normalized magnitudes of wavelet coefs
    timefreq['times'] = timestamps
    timefreq['frequencies'] = freqs
    timefreq['scales'] = scales
    timefreq['mother_wavelet'] = mother_wavelet
    timefreq['image_shape'] = np.shape(tfmap)
    timefreq['mother_wavelet_scales'] = wavelet.motherwavelet.scales
    timefreq['raw_wavelet_coefs'] = wavelet.coefs
    timefreq['wavelet_phase'] = np.arctan2(np.imag(wavelet.coefs),np.real(wavelet.coefs))

    return timefreq


# *******************************************************************************

