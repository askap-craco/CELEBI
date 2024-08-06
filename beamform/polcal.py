##################################################
# Author:   Tyson Dial                           #
# Email:    tdial@swin.edu.au                    #
# Date (created):     11/11/2023                 #
# Date (updated):     05/08/2024                 #
##################################################
# Polarisation Calibration [polcal.py]           #
#                                                #
# This script solves for free parameters         #
# that describe the rotational offset and        #
# polarisation leakage of the observation.       #
#                                                #
# This uses a given polarisation calibrator      #
# (for example the VELA pulsar) with well        #
# known polarisation properties to solve for.    #
##################################################

## Imports 
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from copy import deepcopy
import bilby
from bilby.core.utils.io import check_directory_exists_and_if_not_mkdir
import yaml
import inspect


## import basic libraries
import argparse
from os import path, mkdir
import shutil


## constants
c = 2.997924538e8 # Speed of light [m/s]


## plotting properties
rc('font', **{'size' : 16})  # change figure axes font size
default_col = plt.rcParams['axes.prop_cycle'].by_key()['color']

# GLOBAL container for holding data for plotting purposes
class GLOBAL_():
    pass


##-------------------------##
# additional util functions #
##-------------------------##

def calc_ratio(I, X, Ierr = None, Xerr = None):
    """
    Calculate Stokes Ratio X/I 
    """
    XIerr = None

    # calc XI
    XI = X/I

    # calc error?
    if Ierr is not None and Xerr is not None:
        XIerr = np.sqrt((Xerr/I)**2 + (Ierr*X/I**2)**2)

    return XI, XIerr


def rotate_stokes(A, B, angle):
    """
    Apply rotation between 2 stokes parameters, also applies rotation to stokes 
    errors.

    ##== inputs ==##
    A:                  array 1 to rotate 
    B:                  array 2 to rotate
    angle:              angle to rotate A and B by

    ##== outputs ==##
    X:                  rotated A array
    Y:                  rotated B array
        
    """

    # apply rotation between A and B (and errors)
    X = A*np.cos(angle) + B*np.sin(angle)
    Y = -A*np.sin(angle) + B*np.cos(angle)

    return X, Y



def get_args():
    """
    Get arguments passed during script call


    ##==== outputs ====##
    args:               Arguments for POLCAL.py script

    """

    parser = argparse.ArgumentParser(
        description = "Fit for pol cal solutions"
    )

    ## data arguments
    parser.add_argument("-i", help = "Stokes I dynspec", type = str)
    parser.add_argument("-q", help = "Stokes Q dynspec", type = str)
    parser.add_argument("-u", help = "Stokes U dynspec", type = str)
    parser.add_argument("-v", help = "Stokes V dynsepc", type = str)
    parser.add_argument("--l_model", help = "L/I coefficient", type = str)
    parser.add_argument("--v_model", help = "V/I coefficient", type = str)


    ## data reduction arguments
    parser.add_argument("--peak_w", help = "Width of window [in bins] to average data [boost S/N]", type = int, default = 0)
    parser.add_argument("--rms_w", help = "number of bins to estimate rms", type = int, default = 40)
    parser.add_argument("--tN", help = "Factor to scrunch in time", type = int, default = 100)
    parser.add_argument("--fN", help = "Factor to scrunch in frequency", type = int, default = 1)
    parser.add_argument("--RFIguard", help = "Total phase of pulse period between pulse peak and RFI windows",
                        type = float, default = 0.1)
    parser.add_argument("--chanflag", help = "string of channels to flag", type = str, default = None)


    ## calibrator arguments
    parser.add_argument("--pa0", help = "PA at reference frequency", type = float, default = None)
    parser.add_argument("--f0", help = "reference frequency for rm fitting", type = float, default = None)

    ## observation arguments
    parser.add_argument("--cfreq", help = "Central frequency [MHz]", type = float)
    parser.add_argument("--bw", help = "Bandwidth [MHz]", type = float)


    ## sampler arguments
    parser.add_argument("--cpus", help = "Number of cpus to dedicate to sampling, i.e. number of walkers.",
                        type = int, default = 1)
    parser.add_argument("--live", help = "Number of walker iterations?", type = int, default = 1000)
    parser.add_argument("--redo", help = "Redo Sampling", action = "store_true")
    parser.add_argument("--priors", help = "File to YAML file containing priors", default = "pol_priors.yaml")
    parser.add_argument("--ellipse", help = "Also sample possible ellipticity angle between X and Y.", action = "store_true")


    ## output arguments
    parser.add_argument("--odir", help = "output dir", type = str, default = None)

    parser.add_argument("--ofile", help = "Name out polcal solutions .txt file", type = str, 
                        default = "polcal_solutions.txt")
    parser.add_argument("--save_data", help = "Will save measured and corrected stokes data to .npy file",
                        action = "store_true")

    args = parser.parse_args()

    ## processing some arguments
    if args.odir is None:
        args.odir = ""

    return args









def load_data(args):
    """
    Load in Stokes I, Q, U & V data along with 
    estimates on L/I and V/I coefficients.

    ##== inputs ==##
    args:           Arguments of POLCAL.py script


    ##== outputs ==##
    stk:            dict of stokes dynspecs [IQUV]
    freqs:          array of frequencies
    l_model:        L/I model for calibrator
    v_model:        V/I model for calibrator
    
    """

    ## load in folded pulsar stokes data
    stk = {}


    # for now zap top band (works for low and mid band frbs, not sure about upper band frbs)
    # (DC component)
    stk['I'] = np.load(args.i, mmap_mode = "r")[1:]
    stk['Q'] = np.load(args.q, mmap_mode = "r")[1:]
    stk['U'] = np.load(args.u, mmap_mode = "r")[1:]
    stk['V'] = np.load(args.v, mmap_mode = "r")[1:]


    ## Get Frequencies 
    freqs = np.linspace(args.cfreq - args.bw/2 + 0.5, args.cfreq + args.bw/2 - 0.5,
                        int(args.bw))[:-1]
    

    ## models  
    l_model = np.polyval(np.loadtxt(args.l_model), freqs)
    v_model = np.polyval(np.loadtxt(args.v_model), freqs) 


    return stk, freqs, l_model, v_model












def get_spectra(args, stk, freqs, l_model, v_model):
    """
    Get spectra from folded stokes data

    Process:
        1. scrunch in time and frequency
        2. window pulse profile and off-pulse (rms) region to produce spectra


    ##== inputs ==##
    args:           arguments
    stk:            dict of stokes parameters
    freqs:          array of frequencies


    ##== outputs ==##
    stk_s:          noise subtracted stokes ratios (i.e. U/I)
                    ["I", "Ierr", "U", "Uerr", "V", "Verr", "Q", "Qerr"]
    freqs:          scrunched array of frequencies 

    
    """
    print("[POLCAL]: Folding/Scrunching Stokes...")

    ##=================##
    ## local functions ##
    ##=================##

    def scrunch(x, N, axis = 0):
        """
        scrunch data in either f or t

        ##==== inputs ====##
        x:              data
        N:              scrunch factor
        axis:           axis to scrunch in

        ##==== outputs ====## 
        x:              scrunched data

        """
        if axis == -1:
            # time
            t_new = (x.shape[-1]//N) * N
            shape_new = list(x.shape)
            shape_new[-1] = t_new // N
            shape_new.append(N)

            return np.nanmean(x[...,:t_new].reshape(shape_new), axis = -1)
        
        elif axis == 0:
            # frequency
            f_new = (x.shape[0]//N) * N
            shape_new = list(x.shape)[::-1]
            shape_new[-1] = f_new // N
            shape_new.append(N)

            return np.nanmean(x.T[...,:f_new].reshape(shape_new), axis = -1).T

        else:
            print("invalid scrunching axis")
            return x
        

    
    def flag_chan(f, chanflag):
        """
        Mask channels based on string in format...
        chanflag = "a:b,c,d" -> channels a to b, channel c and channel d will be masked

        NOTE: frequency channels are assumed to be at 1MHz resolution to begin with, and in descending order

        ##====== inputs ======##
        f:              freq array
        chanflag:       str for channel masking

        ##====== outputs ======##
        flagged_chans:             indicies for channel masking 

        """

        flagged_chans = []

        # split flag string into segments using the ',' delimeter
        chan_segments = chanflag.split(',')
        for i, segment in enumerate(chan_segments):
            # check if a channel range was given based on if ':' is present
            if ":" in segment:
                # start and end of range will be in ascending order, so just reverse it
                seg_start = float(segment.split(':')[0])
                seg_end = float(segment.split(':')[1])

                # check if within bounds
                if (seg_end < f[0] - 0.5) or (seg_start > f[-1] + 0.5):
                    # outside bounds
                    continue
                # snap to bounds of bandwidth
                if (seg_start < f[0] - 0.5):
                    seg_start = f[0] - 0.5

                if (seg_end > f[-1] + 0.5):
                    seg_end = f[-1] + 0.5

                segid_start = int(seg_start - (f[0] - 0.5))
                segid_end = int(seg_end - (f[0] - 0.5))

                flagged_chans += list(range(segid_start, segid_end + 1))

            else:
                seg_chan = float(segment)
                if (seg_chan < f[0] - 0.5) or (seg_chan > f[-1] + 0.5):
                    # outside bounds
                    continue
                flagged_chans += [int(seg_chan - (f[0] - 0.5))]

        
        return flagged_chans

                 

        

    ## plot data for diagnostics ##
    PLOTDATA.ds = {}
    PLOTDATA.r_raw_spec = {}
    PLOTDATA.r_spec = {}
    PLOTDATA.off_pulse = {}



    ##=======================##
    ## fold and make spectra ##
    ##=======================##


    ## create new dict for folded and scrunched stokes dynspecs
    stk_s = {}
    peak_w = args.peak_w
    rms_w  = args.rms_w
    tN     = args.tN
    fN     = args.fN

    ## channel flagging
    if (args.chanflag is None) or (args.chanflag == '') or (args.chanflag == "''"):
        chan_flag = []
    else:
        chan_flag = flag_chan(freqs, args.chanflag)

    ## scrunch freq array to match stokes_spec resolution
    freqs = scrunch(freqs, fN, -1)      
    l_model = scrunch(l_model, fN, -1)    # L/I model
    v_model = scrunch(v_model, fN, -1)    # V/I model


    # iterate over "IQUV" items
    for S in "IQUV":
        print(f"Averaging Stokes {S}...")

        ## scrunch dynamic spectra
        ds_r = scrunch(stk[S], tN, -1)        # Time

        ## mask (flag) dynspec in freq
        ds_r = ds_r[::-1]
        ds_r[chan_flag] = np.nan
        
        # flip band
        ds_r = scrunch(ds_r, fN, 0)     # Frequency
        
        ## get peak of mean pulse in time (using power dynspec)
        if S == "I":
            peak = np.argmax(np.nanmean(ds_r, axis = 0))
            rfi_guard = int(ds_r.shape[1] * args.RFIguard)

        stk_s[S] = np.nanmean(ds_r[:,peak-peak_w:peak+peak_w+1], axis = 1)

        # get offpulse spectra rms by taking bins of either windows of off-pulse
        # on either side
        off_pulse1 = ds_r[:,peak-rfi_guard-rms_w:peak-rfi_guard+1]
        off_pulse2 = ds_r[:,peak+rfi_guard:peak+rfi_guard+rms_w+1]
        off_pulseT = np.concatenate((off_pulse1, off_pulse2), axis = 1)

        # Error accounts for sampling in time along with ratio of off-pulse region to on-pulse region
        # in time.
        off_rms    = np.nanstd(off_pulseT, axis = 1)/(off_pulseT.shape[1]**0.5) * ((2*rms_w + 1)/(peak_w + 1))**0.5
        # account for added noise when doing RFI subtraction
        stk_s[f"{S}err"] = off_rms * (1 + 1/(4*rms_w + 2)**0.5)

        # RFI subtraction (it may be the case that Baseline correction is)
        # being done over a much larger range than what could be allowed when
        # approximating the RFI to be constant. So this step may be added just
        # in case. 
        stk_s[S] -= np.nanmean(off_pulseT, axis = 1)

        # get rid of flagged stuff (np.nan)
        freq_mask = ~np.isnan(stk_s[S])
        stk_s[S] = stk_s[S][freq_mask]
        stk_s[f"{S}err"] = stk_s[f"{S}err"][freq_mask]


        ##=================##
        ## For Diagnostics ##
        ##=================##

        ## plot data for diagnostics ##
        PLOTDATA.r_raw_spec[S] = stk_s[S].copy()
        PLOTDATA.r_raw_spec[f"{S}err"] = stk_s[f"{S}err"].copy()
        PLOTDATA.off_pulse[S] = np.nanmean(off_pulseT, axis = 1)[freq_mask]
        PLOTDATA.ds[S] = ds_r.copy()
        PLOTDATA.peak = peak
        PLOTDATA.rfi_guard = rfi_guard
        PLOTDATA.freq_mask = freq_mask.copy()
        PLOTDATA.freqs = freqs.copy()


    ## mask freqs
    freqs = freqs[freq_mask]
    l_model = l_model[freq_mask]
    v_model = v_model[freq_mask]

    ## calculate Stokes Ratios
    for S in "QUV":
        stk_s[S], stk_s[f"{S}err"] = calc_ratio(stk_s['I'], stk_s[S],
                                        stk_s['Ierr'], stk_s[f"{S}err"])
        PLOTDATA.r_spec = deepcopy(stk_s)

    return stk_s, freqs, l_model, v_model



















def get_solutions(args, stk_s, freqs, l_model, v_model):
    """
    Derive solutions for PAF rotational Offset (PA offset) and Polarisation
    leakage.

    Process:
        1. Set up Sampler
        2. Run sampler and retrieve solutions

    ##==== inputs ====##
    args:           arguments for POLCAL script
    stk_s:          stokes spectra
    freqs:          array of frequencies [MHz]
    l_model:        polcal model of linear fraction
    v_model:        polcal model of circular fraction

    """

    print("[POLCAL]: Deriving Polarisation calibrator Solutions...")



    ##=============================================##
    ## POL LEAKAGE CAL LIKELIHOOD CLASS + FUNCTION ##
    ##=============================================##

    ## Likelyhood function for Stokes calibration
    class pol_likelihood(bilby.Likelihood):
        def __init__(self, f, Q, U, V, function, Qerr, Uerr, Verr):
            """
            Likelihood function for Stokes calibration correction

            ##== inputs ==##
            f:                  Frequencies
            Q, U, V:            Stokes spectra
            Qerr, Uerr, Verr:   Stokes spectra rms
            function:           Function to evaluate model
            
            """
            # f data
            self.f = f
            self.N = f.size

            # stk data
            self.Q = Q
            self.U = U
            self.V = V

            # errors
            # should work regardless of being an array of errors or single value
            self.Qsigma = Qerr
            self.Usigma = Uerr
            self.Vsigma = Verr

            # function
            self.function = function

            # These lines of code infer parameters from provided function
            parameters = inspect.getfullargspec(function).args[1:]
            super().__init__(parameters = dict.fromkeys(parameters))
            self.parameters = dict.fromkeys(parameters)

            self.function_keys = self.parameters.keys()

        def log_likelihood(self):
            """
            Log Likelihood, adding Q, U and V likelihoods
            
            """
            # calculate likelihoods of Q, U and V
            model_parameters = {k: self.parameters[k] for k in self.function_keys}
            Qm, Um, Vm = self.function(self.f, **model_parameters)

            # assume we have error for each frequency
            def lhood(x, mu, sig):
                return -0.5 * np.sum(np.log(2*np.pi*sig**2) + ((x-mu)/sig)**2)
            
            # calculate total
            lhood_total = (lhood(self.Q, Qm, self.Qsigma) 
                            + lhood(self.U, Um, self.Usigma) 
                            + lhood(self.V, Vm, self.Vsigma))


            return lhood_total




    ##=========================##
    ## STOKES FITTING FUNCTION ##
    ##=========================##

    def fit_stk_ellipse(f, psi, rm, tau, phi, Lfrac, Vfrac, alpha):
        """
        Fit for predicted Stk Q, U and V based on pol cal solutions

        ##== inputs ==##
        f:              Frequencies (MHz)
        psi:            rotation offset in [X,Y] PAF basis
        rm:             rotation measure
        tau:            time delay in U and V (ms)
        phi:            phase delay in U and V
        Lfrac:          Linear polarisation fraction free parameter 
        Vfrac:          Circular polarisation fraction free parameter
        alpha:          elliptical angle for X,Y coupling


        ##== outputs ==##
        Qi:           Predicted stokes Q
        Ui:           Predicted stokes U
        Vi:           Predicted stokes V
        
        """

        # combined pol leakage angle
        theta = 2 * np.pi * f * tau + phi
        
        # Faraday rotation and PAF rotational offset
        PA = rm * c**2 / 1e12 * (1/f**2 - 1/args.f0**2) + args.pa0 + psi

        # base vela stokes
        Qi = Lfrac * l_model * np.cos(2*PA)
        Ui = Lfrac * l_model * np.sin(2*PA)
        Vi = Vfrac * v_model

        # apply X,Y elliptical pol leakage (between Q and V)
        Qi, Vi = rotate_stokes(Qi, Vi, 2 * alpha)

        # apply Polarisation leakage (between U and V)
        Ui, Vi = rotate_stokes(Ui, Vi, theta)

        return Qi, Ui, Vi

    def fit_stk(f, psi, rm, tau, phi, Lfrac, Vfrac):
        """
        Fit for predicted Stk Q, U and V based on pol cal solutions

        ##== inputs ==##
        f:              Frequencies (MHz)
        psi:            rotation offset in [X,Y] PAF basis
        rm:             rotation measure
        tau:            time delay in U and V (ms)
        phi:            phase delay in U and V
        Lfrac:          Linear polarisation fraction free parameter 
        Vfrac:          Circular polarisation fraction free parameter
        alpha:          elliptical angle for X,Y coupling


        ##== outputs ==##
        Qi:           Predicted stokes Q
        Ui:           Predicted stokes U
        Vi:           Predicted stokes V
        
        """

        # combined pol leakage angle
        theta = 2 * np.pi * f * tau + phi
        
        # Faraday rotation and PAF rotational offset
        PA = rm * c**2 / 1e12 * (1/f**2 - 1/args.f0**2) + args.pa0 + psi

        # base vela stokes
        Qi = Lfrac * l_model * np.cos(2*PA)
        Ui = Lfrac * l_model * np.sin(2*PA)
        Vi = Vfrac * v_model

        # apply Polarisation leakage (between U and V)
        Ui, Vi = rotate_stokes(Ui, Vi, theta)

        return Qi, Ui, Vi



    
    ##===============================##
    ## SET UP AND RUN NESTED SAMPLER ##
    ##===============================##

    # stk func to sample
    if args.ellipse:
        stk_func = fit_stk_ellipse
    else:
        stk_func = fit_stk

    # load priors
    with open(args.priors, "r") as f:
        p = yaml.safe_load(f)

        if not args.ellipse:
            del p['alpha']
        
        print("\nFitting the following parameters:")
        for key in p.keys():
            print(f"{key}:  {p[key]}")
        print("\n")

        if args.ellipse:
            print("\033[36m In addition, the ellipticity angle between X and Y will also be modelled. \033[0m")


    # convert to Uniform (Gaussian) prior object (Bilby)
    priors = {}
    for key in p.keys():
        priors[key] = bilby.core.prior.Uniform(*p[key], key)


    # check if saved sampler is being scrapped, if so remake directory 
    # for bilby output
    sampler_dir = path.join(args.odir, "polcal_sampler", "")
    check_directory_exists_and_if_not_mkdir(sampler_dir)
    if args.redo:
        shutil.rmtree(sampler_dir)

        mkdir(sampler_dir)


    # create likelihood instance
    likelihood = pol_likelihood(freqs, stk_s['Q'], stk_s['U'], stk_s['V'], stk_func, 
                                stk_s['Qerr'], stk_s['Uerr'], stk_s['Verr'])

    # run sampler -> by default runs "dynesty" sampler
    result = bilby.run_sampler(likelihood = likelihood, priors = priors, outdir = sampler_dir,
                               label = "polcal", npool = args.cpus, nlive = args.live)

    
    # get posterior
    posterior = {}
    posterior_err = {}

    for key in priors.keys():
        # get bestfit value
        par = result.get_one_dimensional_median_and_error_bar(key)
        posterior[key] = par.median
        posterior_err[key] = (abs(par.plus) + abs(par.minus))/2
        print(f"{key}: {par.median}    +{par.plus}    -{par.minus}")


    # get prediction for measures stokes spectra
    stk_pred = {}
    stk_pred['Q'], stk_pred['U'], stk_pred['V'] = stk_func(freqs, **posterior)

    # add alpha if not fitting ellipticity
    if not args.ellipse:
        posterior['alpha'] = 0.0
        posterior_err['alpha'] = 0.0



    ##=======================##
    ## SAVE RESULTS TO FILES ##
    ##=======================##

    # save solutions to .txt file
    with open(path.join(args.odir, args.ofile), "w") as file:
        for key in posterior.keys():
            file.write(f"{key}:{posterior[key]}:{posterior_err[key]}\n")

        if not args.ellipse:
            file.write("ELLIPTICITY:0\n")
        else:
            file.write("ELLIPTICITY:1\n")


    ## plot data for diagnostics ##
    PLOTDATA.posterior = deepcopy(posterior)
    PLOTDATA.stk_pred  = deepcopy(stk_pred)

    # make corner plot
    result.plot_corner()


    return









    





def plot_diagnostics(args, stk_s, freqs, l_model, v_model):
    """
    Plot diagnostics
        1. Stokes dynamic spectra showing windowed pulse
        2. Stokes spectra
        3. PA fit
        4. Pol leakage fit
        5. measured/corrected data VS calibrator model 

    ##==== inputs ====##
    args:               arguments for POLCAL.py
    stk_s:              stokes spectra              
    freqs:              array of frequencies [MHz]
    l_model:            polcal model of linear fraction
    v_model:            polcal model of circular fraction

    
    """
    ##=====================##
    ## Auxillary functions ##
    ##=====================##

    def save_stokes(freqs, stk, filename):
        """
        Save stokes data to file

        ##==== inputs ====##
        freqs:          frequency array [MHz]
        stk:            Stokes spectra
        filename:       Filename of file to save to

        """
        with open(filename, "wb") as f:
            for S in "IQUV":
                np.save(f, stk[S])
            
            np.save(f, freqs)

        return 


    
    print("[POLCAL]: Plotting Diagnostics...")
  
    ##===========================##
    ## 1. Stokes dynamic spectra ##
    ##===========================##

    peak = PLOTDATA.peak                # peak of pulse profile (I)
    peak_w = args.peak_w                # width of pulse profile window
    rms_w = args.rms_w                  # width of rms off-pulse window
    rfi_guard = PLOTDATA.rfi_guard      # rfi guard

    for i,S in enumerate("IQUVIQUV"):

        #ds data
        ds = PLOTDATA.ds[S]
        spec = PLOTDATA.r_raw_spec[S]

        filename_ext = ""
        if i > 3:
            ds[PLOTDATA.freq_mask] -= PLOTDATA.off_pulse[S][:,None]
            spec -= PLOTDATA.off_pulse[S]
            filename_ext = "_RFIsub"

        fig = plt.figure(figsize = (10, 10))
        ax1 = fig.add_axes([0.10, 0.07, 0.74, 0.74])  # dynspec
        ax2 = fig.add_axes([0.10, 0.81, 0.74, 0.1])  # time profile
        ax3 = fig.add_axes([0.84, 0.07, 0.10, 0.74])  # on-pulse spectra

        # plot dynamic spectrum
        im_xlim = [.0, ds.shape[1] - 1]
        im_ylim = [PLOTDATA.freqs[0], PLOTDATA.freqs[-1]]
        ax1.imshow(X = ds[::-1], aspect = 'auto', extent = [*im_xlim, *im_ylim])
        # region markers
        ax1.plot([peak - peak_w, peak - peak_w], im_ylim, 'r--')
        ax1.plot([peak + peak_w, peak + peak_w], im_ylim, 'r--')
        ax1.plot([peak - rms_w - rfi_guard, peak - rms_w - rfi_guard], im_ylim, 'm--')
        ax1.plot([peak - rfi_guard, peak - rfi_guard], im_ylim, 'm--')
        ax1.plot([peak + rfi_guard, peak + rfi_guard], im_ylim, 'm--')
        ax1.plot([peak + rms_w + rfi_guard, peak + rms_w + rfi_guard], im_ylim, 'm--')
        # axes labels
        ax1.set(xlim = im_xlim, ylim = im_ylim, xlabel = "Time [samples]", ylabel = "Frequency [MHz]")

        # plot time series
        ax2.plot(np.linspace(*im_xlim, ds.shape[1]), np.nanmean(ds, axis = 0), 'k')
        t_ylim = ax2.get_ylim()
        # region markers
        ax2.plot([peak - peak_w, peak - peak_w], t_ylim, 'r--')
        ax2.plot([peak + peak_w, peak + peak_w], t_ylim, 'r--')
        ax2.plot([peak - rms_w - rfi_guard, peak - rms_w - rfi_guard], t_ylim, 'm--')
        ax2.plot([peak - rfi_guard, peak - rfi_guard], t_ylim, 'm--')
        ax2.plot([peak + rfi_guard, peak + rfi_guard], t_ylim, 'm--')
        ax2.plot([peak + rms_w + rfi_guard, peak + rms_w + rfi_guard], t_ylim, 'm--')
        # axes labels
        ax2.set(xlim = im_xlim, ylim = t_ylim, ylabel = "Flux Density (arb.)")
        ax2.get_xaxis().set_visible(False)

        # plot frequency spectrum
        ax3.plot(spec, freqs, 'r')
        # axes labels
        ax3.set(ylim = im_ylim, xlabel = "Flux Density (arb.)")
        ax3.get_yaxis().set_visible(False)

        plt.suptitle(f"Stokes [{S}] dynamic Spectrum")

        plt.savefig(path.join(args.odir,f"polcal_{S}{filename_ext}_ds.png"))



    ##=========================##
    ## 2. Stokes ratio spectra ##
    ##=========================##

    fig, AX = plt.subplots(2, 2, figsize = (16, 10))
    AX = AX.flatten()

    # plot
    PLOTDATA.r_spec['L'] = np.sqrt(PLOTDATA.r_spec['Q']**2 + 
                                       PLOTDATA.r_spec['U']**2)
    PLOTDATA.r_spec["Lerr"] = np.sqrt(PLOTDATA.r_spec['Q']**2*PLOTDATA.r_spec['Qerr']**2 + 
                                      PLOTDATA.r_spec['U']**2*PLOTDATA.r_spec['Uerr']**2)/PLOTDATA.r_spec['L']
    for i,S in enumerate("LQUV"):
        # plot
        AX[i].plot(freqs, PLOTDATA.r_spec[S], default_col[i])
        AX[i].fill_between(freqs, PLOTDATA.r_spec[S]-PLOTDATA.r_spec[f"{S}err"],
                         PLOTDATA.r_spec[S] + PLOTDATA.r_spec[f"{S}err"], color = default_col[i], alpha = 0.4)
        AX[i].set(title = f"{S}/I")

    # axes labels
    AX[0].set(ylabel = "S/I")
    AX[2].set(xlabel = "Frequency [MHz]", ylabel = "S/I")
    AX[3].set(xlabel = "Frequency [MHz]")
    
    fig.tight_layout()

    plt.savefig(path.join(args.odir,f"polcal_stokes_ratios.png"))



    ##=========================##
    ## 2.2. Stokes raw spectra ##
    ##=========================##

    fig, AX = plt.subplots(2, 3, figsize = (16, 10))
    AX = AX.flatten()

    # plot
    PLOTDATA.r_raw_spec['L'] = np.sqrt(PLOTDATA.r_raw_spec['Q']**2 + 
                                       PLOTDATA.r_raw_spec['U']**2)
    PLOTDATA.r_raw_spec["Lerr"] = np.sqrt(PLOTDATA.r_raw_spec['Q']**2*PLOTDATA.r_raw_spec['Qerr']**2 + 
                                      PLOTDATA.r_raw_spec['U']**2*PLOTDATA.r_raw_spec['Uerr']**2)/PLOTDATA.r_raw_spec['L']
    for i,S in enumerate("ILQUV"):
        # plot
        AX[i].plot(freqs, PLOTDATA.r_raw_spec[S], default_col[i])
        AX[i].fill_between(freqs, PLOTDATA.r_raw_spec[S]-PLOTDATA.r_raw_spec[f"{S}err"],
                         PLOTDATA.r_raw_spec[S] + PLOTDATA.r_raw_spec[f"{S}err"], color = default_col[i], alpha = 0.4)
        AX[i].set(title = f"{S}")

    # axes labels
    AX[0].set(ylabel = "Flux Density (arb.)")
    AX[3].set(xlabel = "Frequency [MHz]", ylabel = "Flux Density (arb.)")
    AX[4].set(xlabel = "Frequency [MHz]")
    AX[5].set_axis_off()
    
    fig.tight_layout()

    plt.savefig(path.join(args.odir,f"polcal_stokes_raw.png"))


    ##=========================##
    ## 2.5. Stokes time series ##
    ##=========================##

    fig, ax = plt.subplots(1, 1, figsize = (10,10))

    for i, S in enumerate("IQUV"):
        
        #ds data
        ds = PLOTDATA.ds[S]
        
        # plot
        ax.plot(np.linspace(0.0, 1.0, ds.shape[1]), np.mean(ds, axis = 0), 
                default_col[i], label = S)

    ax.legend()
    ax.set(xlabel = "Phase (Time units)", ylabel = "Flux Density (arb.)",
               title = "Time series stokes spectra")
        
    plt.savefig(path.join(args.odir, f"polcal_stokes_t.png"))




    ##======================================================##
    ## Making models, correcting data for further modelling ##
    ##======================================================##

    # angles needed for correction
    params = PLOTDATA.posterior
    leak_ang    = 2 * np.pi * freqs * params['tau'] + params['phi']    # pol leakage
    psi_ang     = 2 * params['psi']                                    # offset rotation
    faraday_ang = 2 * (params['rm'] * c**2 / 1e12 * 
                    (1/freqs**2 - 1/args.f0**2) + args.pa0)            # faraday rotation
    ellipse_ang  = 2 * params['alpha']                                 # elliptical angle
    
    # stokes model data
    stk_m = {}
    stk_m['Q'] = l_model * np.ones(freqs.size)
    stk_m['U'] = np.zeros(freqs.size)
    stk_m['V'] = v_model * np.ones(freqs.size)
    stk_m['L'] = l_model * np.ones(freqs.size)
    stk_m['fQ'] = l_model * np.cos(faraday_ang + psi_ang) # with faraday rotation
    stk_m['fU'] = l_model * np.sin(faraday_ang + psi_ang) # with faraday rotation

    # correct data
    stk_corr = {}
    stk_corr['I'] = stk_s['I'].copy()
    stk_corr['U'], stk_corr['V'] = rotate_stokes(stk_s['U'], stk_s['V'], -leak_ang)                     # correct leakage
    stk_corr['Q'], stk_corr['V'] = rotate_stokes(stk_s['Q'], stk_corr['V'], -ellipse_ang)               # correct ellipticity
    stk_corr['Q'], stk_corr['U'] = rotate_stokes(stk_corr['Q'], stk_corr['U'], psi_ang)                 # correct rotation
    stk_corr['Q'], stk_corr['U'] = rotate_stokes(stk_corr['Q'], stk_corr['U'], faraday_ang)             # correct faraday
    stk_corr['L'] = np.sqrt(stk_corr['Q']**2 + stk_corr['U']**2)        # Linear pol fraction (L/I)

    # make seperate spectra with pol leakage and rotation offset added
    # back in, this allowing us to remove faraday rotation for diagnostic
    # purposes.
    stk_meas = {}
    stk_meas['Q'], stk_meas['U'] = rotate_stokes(stk_corr['Q'],stk_corr['U'], -psi_ang)  
    stk_meas['Q'], stk_meas['V'] = rotate_stokes(stk_meas['Q'], stk_corr['V'], ellipse_ang)          
    stk_meas['U'], stk_meas['V'] = rotate_stokes(stk_meas['U'], stk_meas['V'], leak_ang)
    stk_meas['L'] = np.sqrt(stk_meas['Q']**2 + stk_meas['U']**2)        # Linear pol fraction (L/I)

    # make seperate spectra for predicted data
    stk_pred = PLOTDATA.stk_pred
    stk_pred['L'] = np.sqrt(stk_pred['Q']**2 + stk_pred['U']**2)

    # Linear pol Fraction for plotting 
    stk_s['L'] = np.sqrt(stk_s['Q']**2 + stk_s['U']**2)

    # # make dataset with only leakage
    U_leak, V_leak = rotate_stokes(stk_s['U'], stk_s['V'], -leak_ang)
    Q_leak, V_leak = rotate_stokes(stk_s['Q'], V_leak, -ellipse_ang)

    Q_ellipse = Q_leak.copy()
    V_ellipse = V_leak.copy()

    U_leak, V_leak = rotate_stokes(U_leak, V_leak, leak_ang)

    # make dataset with only ellipticity
    Q_ellipse, V_ellipse = rotate_stokes(Q_ellipse, V_ellipse, ellipse_ang)

    # create dataset with rotation offset removed for PA diagnostics
    Q_PA, U_PA = rotate_stokes(stk_corr['Q'], stk_corr['U'], -faraday_ang)


    ##===================##
    ## 3. PA fitted plot ##
    ##===================##

    plt.figure(figsize = (10, 10))
    PA = 0.5 * np.arctan2(U_PA, Q_PA)
    # plot
    plt.scatter(freqs, PA%np.pi, marker = '+', s = 20, c = 'k', label = "Measured (offset corrected, de-faraday)")
    plt.plot(freqs, 0.5*faraday_ang%np.pi, 'r--', label = "Model")
    # axes labels
    plt.title(f"PA = 0.09 * [{params['rm']:.2f}]$(1/f^{{2}} - 1/[{args.f0:.2f}]^{{2}}) + [{args.pa0:.2f}]$")
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("PA [rad]")

    plt.legend()

    plt.savefig(path.join(args.odir,f"polcal_PA_fit.png"))




    ##========================##
    ## 4. leakage fitted plot ##
    ##========================##

    plt.figure(figsize = (10, 10))
    # TODO: Will i need to wrap this?
    theta_pred = np.arctan2(params['Lfrac']*stk_m['fU'], params['Vfrac']*stk_m['V'])
    theta_meas = np.arctan2(U_leak, V_leak)
    theta_res = theta_meas-theta_pred
    # plot
    plt.scatter(freqs, theta_res, c = 'k', s = 20, marker = '+', label = "residual $\\theta$")
    plt.plot(freqs, leak_ang, 'r--', label = "Model $\\theta$")

    # axes labels
    plt.title(f"$\\theta = 2\\pi$f[{params['tau']:.2g}] + [{params['phi']:.2f}]")
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("$\\theta$ [rad]")
    plt.ylim([-np.pi, np.pi])

    plt.legend()

    plt.savefig(path.join(args.odir,f"polcal_leakage_fit.png"))




    ##==========================##
    ## 4.5. ellipse fitted plot ##
    ##==========================##

    if args.ellipse:

        plt.figure(figsize = (10, 10))
        # TODO: Will i need to wrap this?
        theta_pred = np.arctan2(params['Lfrac']*stk_m['fQ'], params['Vfrac']*stk_m['V'])
        theta_meas = np.arctan2(Q_ellipse, V_ellipse)
        theta_res = theta_meas - theta_pred
        # plot
        plt.scatter(freqs, theta_res, c = 'k', s = 20, marker = '+', label = "residual $\\theta$")
        plt.plot(freqs, ellipse_ang*np.ones(freqs.size), 'r--', label = "Model $\\theta$")

        # axes labels
        plt.title(f"$\\theta = 2\\pi$f[{params['tau']:.2g}] + [{params['phi']:.2f}]")
        plt.xlabel("Frequency [MHz]")
        plt.ylabel("$\\theta$ [rad]")
        plt.ylim([-np.pi, np.pi])

        plt.legend()

        plt.savefig(path.join(args.odir,f"polcal_ellipse_fit.png"))




    ##==========================================##
    ## 5. measured/corrected data VS model data ##
    ##==========================================##

    fig = plt.figure(figsize = (16, 14))
    ax1 = fig.add_axes([0.1, 0.54, 0.4, 0.4])   # for measured data plot
    ax2 = fig.add_axes([0.57, 0.54, 0.4, 0.4])  # for corrected data plot
    ax3 = fig.add_axes([0.3, 0.07, 0.4, 0.4])   # for predicted data plot
    # plot
    ax3.scatter([],[],marker = '+', s = 20, c = 'k', label = "Data")
    ax3.plot([],[], 'k', label = "Calibrator model")

    for i,S in enumerate("QUVL"):
        lfrac = 1.0
        if S in "QUL":   # account for L/I fraction free parameter
            lfrac = params['Lfrac']
        elif S == "V":
            lfrac = params['Vfrac']
        # plot de-faraday measured data
        ax1.plot(freqs, stk_meas[S], linestyle = "--", color = default_col[i])
        ax1.plot(freqs, lfrac * stk_m[S], color = default_col[i])

        # plot de-faraday corrected data
        ax2.plot(freqs, stk_corr[S], linestyle = "--", color= default_col[i])
        ax2.plot(freqs, lfrac * stk_m[S], color = default_col[i])

        # plot predicted data
        ax3.plot(freqs, stk_s[S], linestyle = "--", color = default_col[i])
        ax3.plot(freqs, stk_pred[S], color = default_col[i], label = S)


    # axes labels
    lim1 = ax1.get_ylim()
    lim2 = ax2.get_ylim()
    lim_new = [min([lim1[0],lim2[0]]), max([lim1[1],lim2[1]])]
    ax1.set_ylim(lim_new)
    ax2.set_ylim(lim_new)
    ax1.set(xlabel = "Frequency [MHz]", ylabel = "Stokes/I", title = "Measured (de-faraday)")
    ax2.set(xlabel = "Frequency [MHz]", title = "Corrected (de-faraday)")
    ax3.set(title = "Predicted")

    ax3.legend(loc = "center left", bbox_to_anchor=(1.15, 0.5))

    plt.savefig(path.join(args.odir,f"polcal_stokes_corrected.png"))


    ##=========================================================##
    ## 6. measured/corrected data VS model data, seperate axes ##
    ##=========================================================##

    fig, AX = plt.subplots(2, 2, figsize = (10,10))
    AX = AX.flatten()

    for i, S in enumerate("QUVL"):
        AX[i].set_title(f"S = {S}")
        AX[i].plot(freqs, stk_s[S], linestyle = "--", color = default_col[i])
        AX[i].plot(freqs, stk_pred[S], color = 'k')

        if i == 0:
            AX[i].plot([], [], linestyle = "--", color = 'k', label = "data")
            AX[i].plot([], [], color = 'k', label = "model")
            AX[i].legend()
    
    AX[0].set_ylabel("S/I")
    AX[2].set_ylabel("S/I")
    AX[2].set_xlabel("Freq [MHz]")
    AX[3].set_xlabel("Freq [MHz]")

    plt.savefig(path.join(args.odir, "polcal_stokes_predicted.png"))





    ##======================##
    ## optional - save data ##
    ##======================##

    # save data if enabled
    if args.save_data:
        print("Saving Stokes data...")
        save_stokes(freqs, stk_s, path.join(args.odir, "measured_polcal_spec.npy"))
        stk_corr['Q'], stk_corr['U'] = rotate_stokes(stk_corr['Q'], stk_corr['U'], -faraday_ang)
        save_stokes(freqs, stk_corr, path.join(args.odir, "corrected_polcal_spec.npy"))
        for S in "IQUV":
            np.save(path.join(args.odir, f"polcal_{S}_ds.npy"),PLOTDATA.ds[S])


    return







#### ENTRY POINT ####
if __name__ == "__main__":
    """
    Main code block 
    """

    PLOTDATA = GLOBAL_()

    ## get arguments of script call
    args = get_args()


    ## load in data
    stk, freqs, l_model, v_model = load_data(args)


    ## process Stokes data (fold, scrunch)
    stk_s, freqs, l_model, v_model = get_spectra(args, stk, freqs,
                                                 l_model, v_model)


    ## channel zapping?


    ## Derive Calibrator Solutions
    get_solutions(args, stk_s, freqs, l_model, v_model)


    ## diagnostics
    plot_diagnostics(args, stk_s, freqs, l_model, v_model)


    print("[POLCAL]: Completed.")
    # END OF SCRIPT.