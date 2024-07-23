import numpy as np
from scipy.fftpack import dct, idct
from scipy import linalg


def get_kc(CI_data: np.ndarray):
    """
    Determine spectral cutoff kc for the low-pass filter.
    
    :param CIdata: FRB spectrum in the discrete cosine domain, absolute values for all delta DM
    :type CIdata: :class:`np.ndarray`
    :return: Spectral cutoff `kc` for use in filtering the supplied FRB.
    :rtype: int
    """
    cumulative_window_size = 5
    noise_margin = int(CI_data.shape[1]/2) #assuming all of the higher frequency data is pure noise
    CI_data_transposed = np.transpose(np.abs(CI_data))
    CI_data_maxima = np.max(CI_data_transposed, axis = 1)
    kc = None

    #find the level where the noise flattens out
    noise_top = np.mean(CI_data_maxima[noise_margin:])

    #find where the signal dips down into the noise
    for i in range(cumulative_window_size,len(CI_data_transposed)):
        rolling_average=np.mean(CI_data_maxima[i-cumulative_window_size:i+1])
        if rolling_average <= noise_top:
            kc = i
            break
    
    if kc is not None:
        print(f"Found kc of: {kc}")
    else:
        raise RuntimeError('Could not find a value for kc')

    return kc


def uncertainty_calc(detrended_noise: np.ndarray, LPF_data: np.ndarray, filter_diag: np.ndarray) -> np.ndarray:
    """
    Calculate uncertainty from the detrended noise and smoothed data.
    
    :param detrended_noise: The I(DM, t) array of the detrended noise (`I_data - I_smooth`)
    :type detrended_noise: :class:`np.ndarray`
    :param LPF_data: The discrete cosine space smoothed data.
    :type LPF_data: :class:`np.ndarray`
    :param filter_diag: The low pass filter as a diagonal matrix.
    :type filter_diag: :class:`np.ndarray`
    :return: Array containing the uncertainty at each delta DM
    :rtype: :class:`np.ndarray`
    """
    CI_dbl_filtered=filter_diag@LPF_data #double smoothing
    norm_CI_dbl_filtered=linalg.norm(CI_dbl_filtered, axis=0) #compute norm

    C_delta_I = dct(detrended_noise, norm='ortho') #DCT the DeltaI
    C_delta_I_filtered = filter_diag@np.transpose(C_delta_I) #pass DCT data through the combined filter
    norm_C_delta_I_filtered = linalg.norm(C_delta_I_filtered, axis=0) #compute norm
    
    uncertainty = norm_C_delta_I_filtered/norm_CI_dbl_filtered
    return uncertainty


def get_ranges_above_max(max_SP: float, adjusted_SPs: np.ndarray) -> list:
    """
    Find the start and end index for every unbroken sequence of values in adjusted_SPs above max_SP.
    
    :param max_SP: The largest structure parameter prior to scaling
    :type max_SP: float
    :param adjusted_SPs: Calculated structure parameters scaled up by their relative uncertainty.
    :type adjusted_SPs: :class:`np.ndarray`
    :return: List of lists where possible_max_ranges[i][0] represents the 
    starting index of an uncertainty range, and possible_max_ranges[i][1] represents the end.
    :rtype: list
    """
    possible_max_ranges = []
    continuous = False
    for idx, SP in enumerate(adjusted_SPs):
        if not continuous and SP >= max_SP:
            continuous = True
            possible_max_ranges.append([idx])
        elif continuous and SP < max_SP:
            continuous = False
            possible_max_ranges[-1].append(idx)
    return possible_max_ranges