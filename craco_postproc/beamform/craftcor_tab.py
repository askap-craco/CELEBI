#!/usr/bin/env python3
"""
Tied-array beamforming vcraft files, based on "craftcor.py".

Copyright (C) CSIRO 2017
"""
#import pylab
from re import M
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import vcraft
import glob
from calc11 import ResultsFile
#from corruvfits import CorrUvFitsFile
from astropy.coordinates import SkyCoord
import multiprocessing
import signal
import warnings
from scipy.interpolate import interp1d
from joblib import Parallel, delayed
from timeit import default_timer as timer

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

CLIGHT=299792458.0

def print_delay(xx):
    xxang = np.angle(xx)
    punwrap = np.unwrap(xxang)
    f = np.arange(len(punwrap)) - len(punwrap)/2
    gradient, phase = np.polyfit(f, punwrap, 1)
    '''
    pylab.figure(20)
    pylab.plot(f, xxang)
    pylab.plot(f, punwrap, 'x')
    pylab.plot(f, np.polyval((gradient, phase), f))
    pylab.show()
    '''
    delay = gradient/2./np.pi*len(punwrap)
    delayns = delay*27./32.*1e3*(54./len(punwrap))
    print(('Unwrapped phase = {} rad = {} deg, gradient={} rad per channel, delay={}samples ={} ns nsamp={}' \
        .format(phase, np.degrees(phase), gradient, delay, delay*27./32.*1e3, len(punwrap))))

    return (delay, np.degrees(phase))

class FitsOut(object):
    def __init__(self, fname, corr):
        ants = self.corr.ants

    def put_product(self, corr, a1, a2, xx):
        pass

class PlotOut(object):
    def __init__(self, corr):
        self.stuff = []
        self.corr = corr
        self.delayout = open(corr.values.outfile.replace('.fits','')+'.delays.parset', 'w')

    def put_product(self, a1, a2, xxp):
        if self.corr.values.show:
            if a1 != a2:
                xx= xxp[:,0]
                self.stuff.append(xx)
                self.last_xx = (a1, a2, xx)
                if len(self.stuff) >0:
                    self.plot_stuff(a1, a2)

    def plot_stuff(self, a1, a2):

        '''(a1, a2, xx) = self.last_xx
        fig, (ax1, ax2, ax3, ax4, ax5,ax6) = pylab.subplots(6,1)
        xxang = np.angle(xx)
        print_delay(xx)
        ax1.plot(abs(xx))
        ax2.plot(np.degrees(xxang), 'o')
        nf = self.corr.nfine_out_per_coarse
        for i in xrange(self.corr.ncoarse_chan):
            ax2.axvline(nf*i, c='r')
            ax1.axvline(nf*i, c='r')
            print 'PLOTSTUFF', nf, i, nf*i, nf*(i+1), self.corr.ncoarse_chan, xx.shape
            print_delay(xx[nf*i:nf*(i+1)])

        lag = np.fft.fftshift(abs(np.fft.fft(xx)))
        lagidx = np.argmax(lag)
        lagoff = lagidx - len(xx)/2.0
        lagns= lagoff/self.corr.full_bw*1e3
        lagsn = lag[lagidx]/np.std(lag[0:lagidx-100])
        lagx = np.arange(len(lag)) - len(lag)/2.0
        ax3.plot(lagx/self.corr.full_bw*1e3, lag,label='lag')
        ax3.axvline(0, c='r')
        ax3.set_xlabel('Delay (ns)')
        angxx = np.angle(np.array(self.stuff))
        ax4.imshow(angxx, aspect='auto')
        phase_vt = np.degrees(angxx[:, 10:50].mean(axis=1))
        t = np.arange(len(phase_vt))
        rate, offset = np.polyfit(t, phase_vt, 1)
        fit_std = (np.polyval((rate, offset), t) - phase_vt).std()
        ax5.plot(phase_vt)
        ax6.imshow(abs(np.array(self.stuff)), aspect='auto')
        print 'Phase rate={} offset={} std={} deg lagoff={}samples = {}ns lag S/N={}'.format(rate, offset, fit_std, lagoff, lagns, lagsn)
        fig.suptitle('a1 {} a2 {}'.format(a1.antname, a2.antname))
        curr_delay = self.corr.get_fixed_delay_usec(a2.antno)*1e3
        if a1 != a2:
            self.delayout.write('#{}={} S/N={} delay={}ns\n'.format(a1.antname, a2.antname, lagsn, lagns))
            self.delayout.write('common.antenna.ant{}.delay={}ns\n'.format(a2.antno, curr_delay+lagns))
            self.delayout.flush()

        pylab.show()'''
        return # added by hyerin 

    def finish(self):
        stuff = np.array(self.stuff)
        #pylab.imshow(np.angle(stuff))
        #pylab.show()


class AntennaSource(object):
    def __init__(self, vfile):
        self.vfile = vfile
        self.antname = self.vfile.hdr['ANT'][0].lower()
        self.antno = int(self.vfile.hdr['ANTENNA_NO'][0])
        self.mjdstart = self.vfile.start_mjd
        self.trigger_frame = self.vfile.start_frameid
        self.hdr = self.vfile.hdr
        self.init_geom_delay_us = None
        self.all_geom_delays = []
        self.all_mjds = []
        self.pol = self.vfile.pol.lower()
        print(('antenna {} {}'.format(self.antname, self.vfile.freqconfig)))

    def do_f(self, corr):
        self.frparams = FringeRotParams(corr, self)
        # calculate sample start
        framediff_samp = corr.refant.trigger_frame - self.trigger_frame
        framediff_us = framediff_samp / corr.fs
        (geom_delay_us, geom_delay_rate_us) = corr.get_geometric_delay_delayrate_us(self)
        self.all_geom_delays.append(geom_delay_us)
        self.all_mjds.append(corr.curr_mjd_mid)

        #logging.debug'ALL GEOM DELAYS', self.all_geom_delays, type(geom_delay_us)
        #print 'ALL MJDs', self.all_mjds
        # test with 0 delay rate
        #geom_delay_us = self.all_geom_delays[0]
        
        geom_delay_samp = geom_delay_us * corr.fs
        fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
        fixed_delay_samp = fixed_delay_us*corr.fs
        total_delay_samp = framediff_samp
        whole_delay = int(np.round(total_delay_samp))
        total_delay_us = total_delay_samp / corr.fs
        whole_delay_us = whole_delay / corr.fs

        frac_delay_samp = (total_delay_samp - whole_delay)
        frac_delay_us = frac_delay_samp *corr.fs

        # get data
        nsamp = corr.nint*corr.nfft
        logging.debug('F %s sample delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d phase=%f deg',
            self.antname, framediff_samp, geom_delay_samp, fixed_delay_samp,
            total_delay_samp, whole_delay, frac_delay_samp, nsamp,
            360.*frac_delay_us*corr.f0)
        logging.debug('F %s us delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d phase=%f deg rate=%e',
                self.antname, framediff_us, geom_delay_us, fixed_delay_us,
                total_delay_us, whole_delay_us, frac_delay_us, nsamp,
                      360.*frac_delay_us*corr.f0, geom_delay_rate_us)

        sampoff = corr.curr_samp*corr.nfft + whole_delay
        rawd = self.vfile.read(sampoff, nsamp)
        assert rawd.shape == (nsamp, corr.ncoarse_chan), 'Unexpected shape from vfile: {} expected ({},{})'.format(rawd.shape, nsamp, corr.ncoarse_chan)
        self.data = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in), dtype=np.complex64)
        d1 = self.data
        nfine = corr.nfft - 2*corr.nguard_chan

        for c in range(corr.ncoarse_chan):
            cfreq = corr.freqs[c]
            freqs = (np.arange(nfine, dtype=np.float) - float(nfine)/2.0)*corr.fine_chanbw
            if corr.sideband == -1:
                freqs = -freqs
                
            x1 = rawd[:, c].reshape(-1, corr.nfft)
            xf1 = np.fft.fft(x1, axis=1)
            xf1 = np.fft.fftshift(xf1, axes=1)

            xfguard = xf1[:, corr.nguard_chan:corr.nguard_chan+nfine:] # scale because oterhwise it overflows
            np.save('xfguard_c{}.npy' % c, xfguard)
            delta_t = -fixed_delay_us + geom_delay_us
            phases = delta_t * (freqs + cfreq)
            logging.debug('PHASOR %s[%s] chan=%s freq=%sfixed=%f us geom=%f us delta_t %s us coff*fixed = %f deg coff*geom = %f deg',
                          self.antname, self.ia, c, cfreq, fixed_delay_us, geom_delay_us, delta_t, cfreq*fixed_delay_us*360., cfreq*geom_delay_us*360.)


            # If you plot the phases you're about to correct, after adding a artificial
            # 1 sample delay ad tryig to get rid of it with a phase ramp, it becaomes
            # blatetly clear what you should do
            phasor = np.exp(np.pi*2j*phases, dtype=np.complex64)
            '''
            pylab.figure(10)
            pylab.plot(np.angle(phasor[0, :]))
            pylab.plot(np.angle(phasor[-1:, :]))
            '''


            xfguard *= phasor
            # slice out only useful channels
            fcstart = c*nfine
            fcend = (c+1)*nfine
            self.data[:, fcstart:fcend, 0] = xfguard

    def do_f_tab(self, corr, iant):
        start = timer()
        self.frparams = FringeRotParams(corr, self)
        # calculate sample start
        framediff_samp = corr.refant.trigger_frame - self.trigger_frame
        (geom_delay_us, geom_delay_rate_us) = corr.get_geometric_delay_delayrate_us(self)
        #samp_us = 1.0/336.0
        #geom_delay_rate_us *= samp_us   # set the units of geom_delay_rate_us to us/us instead of us/sample
        self.all_geom_delays.append(geom_delay_us)
        self.all_mjds.append(corr.curr_mjd_mid)

        #logging.debug'ALL GEOM DELAYS', self.all_geom_delays, type(geom_delay_us)
        #print 'ALL MJDs', self.all_mjds
        # test with 0 delay rate
        #geom_delay_us = self.all_geom_delays[0]
        
        fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
        total_delay_samp = framediff_samp
        whole_delay = int(np.round(total_delay_samp))

        data_out = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in), dtype=np.complex64)
        nfine = corr.nfft - 2*corr.nguard_chan

        nsamp = corr.nint*corr.nfft

        # time-dependent geometric delays
        # np.linspace(0, 1, nsamp) == time in units of integrations
        geom_delays_us = geom_delay_us + geom_delay_rate_us * np.linspace(0, 1, nsamp) - fixed_delay_us
        np.save('TEMP_geom_delays_us.npy', geom_delays_us)
        del geom_delays_us
        geom_delays_us = np.load("TEMP_geom_delays_us.npy", mmap_mode="r")

        sampoff = whole_delay + corr.abs_delay
        print(("antenna #: ",iant, self.antname))
        frameid = self.vfile.start_frameid+sampoff
        print(('FRAMEID: '+str(frameid)+', remainder from 32: '+str(frameid % 32)))
        # To avoid iPFB fractional delay, set FRAMEID such that the remainder is 0

        rawd = self.vfile.read(sampoff, nsamp)  # TODO HERE'S THE DATA
        np.save("TEMP_rawd.npy", rawd)
        del rawd
        rawd = np.load("TEMP_rawd.npy", mmap_mode="r")

        assert rawd.shape == (nsamp, corr.ncoarse_chan), 'Unexpected shape from vfile: {} expected ({},{})'.format(rawd.shape, nsamp, corr.ncoarse_chan)

        # Fine channel frequencies from -0.5 MHz to 0.5 MHz
        freqs = (np.arange(nfine, dtype=float) - float(nfine)/2.0)*corr.fine_chanbw
        # TODO: what is sideband?
        if corr.sideband == -1:
            freqs = -freqs

        # Save and reload to improve memory efficiency
        np.save("TEMP_freqs.npy", freqs)
        del freqs
        freqs = np.load("TEMP_freqs.npy", mmap_mode="r")
        print(f"do_f_tab (stage 1): {timer()-start} s")
        start = timer()
        def process_chan(i, c):
            # Some things here are amalgamated to attempt to improve memory efficiency
            # Channel frequency
            cfreq = corr.freqs[c]

            '''
            rawd's shape: (nsamp, corr.ncoarse_chan)
            nsamp = input.i * (64 * input.n)
            '''
            # x1 = rawd[:, c].reshape(-1, corr.nfft)

            # Fringe rotation for Earth's rotation
            # turn_fringe = cfreq * geom_delays_us

            # phasor_fringe = np.exp(np.pi * 2j * cfreq * geom_delays_us, dtype=np.complex64)

            x1 = rawd[:, c].reshape(-1, corr.nfft) * np.exp(np.pi * 2j * cfreq * geom_delays_us, dtype=np.complex64)

            # corr.nguard_chan = 5 * input.n = number of fine channels on each side to cut off
            # xfguard is xf1 with the ends trimmed off
            # xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)
            # xfguard_f = xf1[:, corr.nguard_chan:corr.nguard_chan+nfine:] # scale because oterhwise it overflows
            xfguard_f = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)[:, corr.nguard_chan:corr.nguard_chan+nfine:]

            del x1

            #logging.debug('PHASOR %s[%s] chan=%s freq=%sfixed=%f us geom=%f us delta_t %s us coff*fixed = %f deg coff*geom = %f deg',
            #             self.antname, self.ia, c, cfreq, fixed_delay_us, geom_delay_us, delta_t, cfreq*fixed_delay_us*360., cfreq*geom_delay_us*360.)

            #   If you plot the phases you're about to correct, after adding a artificial
            #   1 sample delay ad tryig to get rid of it with a phase ramp, it becaomes
            #   blatetly clear what you should do

            # Fractional sample phases
            # turn_frac = freqs * np.mean(geom_delays_us)

            # phasors to rotate the data with constant amplitude = 1
            phasor = np.exp(np.pi * 2j * freqs * np.mean(geom_delays_us), dtype=np.complex64)

            # get absolute frequencies in gigahertz
            freq_ghz = (cfreq+freqs)/1e3

            # get the calibration solutions and apply them to the phasors
            mir_cor = corr.mir.get_solution(iant, 0, freq_ghz)
            del freq_ghz

            #if mir_cor[0] == 0: # if correction is 0, flag data
            #    phasor *= 0
            #else:
            phasor /= mir_cor

            xfguard_f *= phasor
            del phasor

            # select the channels for this coarse channel
            fcstart = c*nfine
            fcend = (c+1)*nfine

            '''
            RECAP
            xfguard is a "dynamic" spectrum of the current coarse channel in
            only the fine channels not trimmed (oversampled PFB).
            xfguard.shape == (input.i, 64 * input.n)
            '''

            '''
            Slot xfguard (a trimmed spectrum for this coarse channel) into
            the corresponding slice of the fine channels
            '''
            data_out[:, fcstart:fcend, 0] = xfguard_f
            del xfguard_f
        
        Parallel(n_jobs=4, require="sharedmem")(delayed(process_chan)(i, c) for i, c in enumerate(range(corr.ncoarse_chan)))
        print(f"do_f_tab (stage 2): {timer()-start} s")
        return data_out

            
    def get_data(self, chan, pol):
        return self.data[:, chan, pol]


class FringeRotParams(object):
    cols = ('U (m)', 'V (m)', 'W (m)', 'DELAY (us)')

    def __init__(self, corr, ant):
        print("frdata_mid keys:", corr.frdata_mid.keys())
        mid_data = corr.frdata_mid[ant.antname]
        self.u,self.v,self.w,self.delay = list(map(float, [mid_data[c] for c in FringeRotParams.cols]))
        self.delay_start = float(corr.frdata_start[ant.antname]['DELAY (us)'])
        self.delay_end = float(corr.frdata_end[ant.antname]['DELAY (us)'])
        self.delay_rate = (self.delay_end - self.delay_start)/float(corr.nint)
        self.ant = ant
        self.corr = corr

    def __str__(self):
        s = 'FR {} uvw=({},{},{}) m = {} us'.format(self.ant.antname, self.u, self.v, self.w, self.delay)
        return s

    __repr__ = __str__

class Correlator(object):
    def __init__(self, ants, sources, values, abs_delay=0):
        self.running = True
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)
        self.ants = ants
        self.values = values
        self.pool = None
        if self.values.num_threads > 1:
            self.pool = multiprocessing.Pool(processes=values.num_threads)

        self.parse_parset()

        for ia, a in enumerate(self.ants):
            a.ia = ia
            a.antpos = self.get_ant_location(a.antno)

        refantname = self.parset['cp.ingest.tasks.FringeRotationTask.params.refant'].lower()
        self.abs_delay = abs_delay
        #self.refant = filter(lambda a:a.antname == refantname, ants)[0]
        self.refant = ants[0]
        self.calcresults = ResultsFile(values.calcfile)
        self.dutc = 0
        self.mjd0 = self.refant.mjdstart + self.dutc/86400.0
        self.frame0 = self.refant.trigger_frame
        self.nint = values.nint
        self.nfft = 64*values.fft_size
        self.nguard_chan = 5*values.fft_size
        self.oversamp = 32./27.
        self.fs = self.oversamp # samples per microsecnd
        self.ncoarse_chan = len(self.refant.vfile.freqs)
        self.sideband = -1
        self.coarse_chanbw = 1.0
        self.nfine_per_coarse = self.nfft - 2*self.nguard_chan
        self.nfine_chan = self.ncoarse_chan*self.nfine_per_coarse
        self.fine_chanbw = self.coarse_chanbw / float(self.nfine_per_coarse)
        self.full_bw = self.fine_chanbw * self.nfine_chan
        self.fscrunch = values.fscrunch
        assert self.fscrunch >= 1
        assert self.nfine_per_coarse % self.fscrunch == 0, 'Fsrunch must yield an integer number of fine channels per coarse channel'
        self.nfine_out_per_coarse = self.nfine_per_coarse / self.fscrunch
        self.nfine_out_chan = self.nfine_out_per_coarse*self.ncoarse_chan
        self.out_chanbw = self.coarse_chanbw / float(self.nfine_out_per_coarse)
        self.npol_in = 1
        self.npol_out = 1
        #self.f0 = self.ants[0].vfile.freqs.mean() # centre frequency for fringe rotation
        self.f0 = self.ants[0].vfile.freqs[0]
        self.freqs = self.ants[0].vfile.freqs
        self.fmid = self.freqs.mean()
        self.inttime_secs = float(self.nint*self.nfft)/(self.fs*1e6)
        self.inttime_days = self.inttime_secs/86400.
        self.curr_intno = 0
        self.curr_samp = self.curr_intno*self.nint + 1000
        self.prodout = PlotOut(self)
        self.calcmjd()
        self.get_fr_data()
        self.pol = self.ants[0].pol
        self.parse_mir()
        #self.fileout = CorrUvFitsFile(values.outfile, self.fmid, self.sideband*self.out_chanbw, \
                                      #self.nfine_out_chan, self.npol_out, self.mjd0, sources, ants, self.sideband)

        logging.debug('F0 %f FINE CHANNEL %f kHz num=%d freqs=%s', self.f0, self.fine_chanbw*1e3, self.nfine_chan, self.freqs)

    def exit_gracefully(self, signum, frame):
        self.running = False

    def parse_parset(self):
        self.parset = {}

        # open the fcm file
        with open(self.values.parset, 'r') as f:
            for line in f:
                if '=' not in line or line.startswith('#'):
                    continue

                name, value = line.strip().split('=')
                name = name.strip()
                value = value.strip()
                self.parset[name] = value

    def parse_mir(self):
        #self.mir = None
        #if self.values.mirsolutions is not None or self.values.aips_c is not None:
        self.mir = MiriadGainSolutions(self.values.mirsolutions,self.values, self.values.aips_c, self.pol, self.freqs)
                

    def get_ant_location(self, antno):
        key = 'common.antenna.ant{}.location.itrf'.format(antno)
        value = self.parset[key]
        location = list(map(float, value.replace('[','').replace(']','').split(',')))
        return location

    def get_fixed_delay_usec(self, antno):
        key = 'common.antenna.ant{}.delay'.format(antno)
        value = self.parset[key]
        delayns =  float(value.replace('ns',''))
        delayus = delayns/1e3

        return delayus


    def get_uvw(self, ant1, ant2):
        fr1 = FringeRotParams(self, ant1)
        fr2 = FringeRotParams(self, ant2)
        uvw = np.array([fr1.u - fr2.u, fr1.v - fr2.v, fr1.w - fr2.w])/CLIGHT

        return uvw

    def get_geometric_delay_delayrate_us(self, ant):
        fr1 = FringeRotParams(self, ant)
        fr2 = FringeRotParams(self, self.refant)

        # TODO: There is a discrepancy here, below comment says fr1 is ref ant, but above suggests fr2 is?
        # fr1: reference antenna
        # Account for effects of Earth's rotation
        #delay = fr1.delay - fr2.delay
        #delayrate = fr1.delay_rate - fr2.delay_rate
        delay = fr1.delay_start - fr2.delay_start
        delayrate = fr1.delay_rate - fr2.delay_rate

        with open('delays/{}_ant_delays.dat'.format(ant.antno), 'w') as f:
            f.write('#field fr1({}) fr2({})\n'.format(ant, self.refant))
            f.write('delay_start {} {}\n'.format(fr1.delay_start, fr2.delay_start))
            f.write('delay {} {}\n'.format(fr1.delay, fr2.delay))
            f.write('delay_end {} {}\n'.format(fr1.delay_end, fr2.delay_end))
            f.write('delay_rate {} {}\n'.format(fr1.delay_rate, fr2.delay_rate))

        return (delay, delayrate)

    def calcmjd(self):
        i = float(self.curr_intno)
        abs_delay_days = float(self.abs_delay)/86400./(self.fs*1e6)
        self.curr_mjd_start = self.mjd0 + self.inttime_days*(i + 0.0) + abs_delay_days
        self.curr_mjd_mid = self.mjd0 + self.inttime_days*(i + 0.5) + abs_delay_days
        self.curr_mjd_end = self.mjd0 + self.inttime_days*(i + 1.0) + abs_delay_days

    def next_integration(self):
        #self.curr_intno +=
        self.curr_intno += 1
        self.curr_samp += self.nint
        self.calcmjd()
        self.get_fr_data()

    def get_calc_results(self, mjd):
        #res = self.calcresults.scans[0].eval_src0_poly_delta(mjd, self.refant.antname.lower())
        res = self.calcresults.scans[0].eval_src0_poly(mjd)  # Calcresults is a ResultsFile

        return res

    def get_fr_data(self):
        self.frdata_start = self.get_calc_results(self.curr_mjd_start)
        self.frdata_mid = self.get_calc_results(self.curr_mjd_mid)
        self.frdata_end = self.get_calc_results(self.curr_mjd_end)


    def do_f(self):
        for iant, ant in enumerate(self.ants):
            if not self.running:
                raise KeyboardInterrupt()

            ant.do_f(self)

    def do_x(self):
        nant = len(self.ants)
        for ia1 in range(nant):
            for ia2 in range(ia1, nant):
                if not self.running:
                    raise KeyboardInterrupt()
                    
                a1 = self.ants[ia1]
                a2 = self.ants[ia2]
                self.do_x_corr(a1, a2)

    def do_x_corr(self, a1, a2):
        npolout = self.npol_out
        xx = np.empty([self.nfine_chan, npolout], dtype=np.complex64)
        #np.seterr(all='raise')
        rfidelay = self.values.rfidelay
        for p1 in range(self.npol_in):
            for p2 in range(self.npol_in):
                d1 = a1.data[:, :, p1]
                d2 = a2.data[:, :, p2]
                pout = p2 + p1*self.npol_in
                ntimesi = d1.shape[0] - rfidelay
                ntimes = float(ntimesi)

                try:
                    for c in range(self.nfine_chan):
                        #xx[:,pout] = (d1 * np.conj(d2)).mean(axis=0)
                        # vdot conjugates the first argument
                        # this is equivalent to (d1 * conj(d2)).mean(axis=0)
                        if self.sideband == -1:
                            xx[c, pout] = np.vdot(d2[:ntimesi, c], d1[rfidelay:, c])/ntimes
                        else:# take complex conjugate if inverted
                            xx[c, pout] = np.vdot(d1[:ntimesi, c], d2[rfidelay:, c])/ntimes
                            
                except Exception as e:
                    print(('Error', e))
                    import ipdb
                    ipdb.set_trace()


        if self.fscrunch > 1:
            chans = np.arange(self.nfine_chan, step=self.fscrunch)
            xx = np.add.reduceat(xx, chans, axis=0)/float(self.fscrunch)

        assert xx.shape == (self.nfine_out_chan, self.npol_out)
            
        self.put_product(a1, a2, xx)


    def do_tab(self, an=None):
        # Tied-array beamforming
        
        nsamp = self.nint
        nchan = self.ncoarse_chan*self.nfine_per_coarse
            
        sum_aligned = np.zeros((nsamp, nchan, self.npol_in), dtype=np.complex64)
        
        if an == None: # add all antennas
            print(('## Summing up all '+str(len(self.ants))+' antennas'))
            if self.values.tab:
                print('coherent sum')
            else:
                print('incoherent sum')
            for iant, ant in enumerate(self.ants):
                if not self.running:
                    raise KeyboardInterrupt()
            
                temp = ant.do_f_tab(self,iant)
                
                if self.values.tab:
                    sum_aligned += temp
                else:
                    sum_aligned += np.power(abs(temp),2)

            return sum_aligned
        else:
            print(('## Operate on only antenna #: '+str(an)))
            ant = self.ants[an]
            iant = an
            start = timer()
            temp = ant.do_f_tab(self,iant)
            print(f"do_f_tab (total): {timer()-start} s")
            return temp
    
    def put_product(self, a1, a2, xx):
        self.prodout.put_product(a1, a2, xx)
        uvw = self.get_uvw(a1, a2)
        #self.fileout.put_data(uvw, self.curr_mjd_mid, a1.ia, a2.ia,
            #self.inttime_secs, xx)



def parse_delays(values):
    delayfile = values.calcfile.replace('.im','.hwdelays')
    if os.path.exists(delayfile)==False:
        delayfile = values.hwfile
    #print(delayfile)
    delays = {}
    if delayfile is not None and os.path.exists(delayfile):
        with open(delayfile, 'rU') as dfile:
            for line in dfile:
                bits = line.split()
                if not line.startswith('#') and len(bits) == 2:
                    raw = -int(bits[1])
                    if raw % 8 !=0: # if it is not a multiple of 8, round
                        new = int(8 * round(float(raw)/8))
                        print(('hwdelay ',raw,' rounded to ',new))
                        delays[bits[0].strip()] = new
                    else:
                        delays[bits[0].strip()] = raw

        logging.info('Loaded %s delays from %s', len(delays), delayfile)
    else:
        logging.info('No delays loaded. %s does not exist', delayfile)
    

    return delays

def parse_gpplt(fin):
    ''' 
    Parse a miriad gpplt exported log file

    :returns: tuple(x, values) where x is an array of strings
    containing x values (dependant on the type of file)
    and values is a (len(x), Nant) numpy array
    '''
    with open(fin, 'rU') as f:

        all_values = [] 
        curr_values = []
        x = []
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue

            # X value is in the first 16 columns
            # if empty, it's a continuation of previous antennas
            sz = 14
            xfield = line[0:sz].strip()
            bits = list(map(float, line[sz:].split()))
            if xfield == '':
                curr_values.extend(bits)
            else:
                x.append(xfield)
                curr_values = []
                all_values.append(curr_values)
                curr_values.extend(bits)

        v = np.array(all_values)
        # make 2D
        if v.ndim == 1:
            v = v[np.newaxis, :]

        assert v.ndim == 2

    return np.array(x), v

class MiriadGainSolutions(object):
    def __init__(self, file_root, values, bp_c_root=None, pol=None, freqs=None):
        '''Loads gpplt exported bandpass and gain calibration solutions.
        Expects 4 files at the following names, produced by miriad gpplt
        with the given options

        Limitations: Currently does no time interpolation - just uses first time
        $file_root.gains.real - yaxis=real
        $file_root.gains.imag - yaxis=imag
        $file_root.bandpass.real - options=bandpass, yaxis=real
        $file_root.bandpass.imag - options=bandpass, yaxis=imag
        
        '''
        

        #print(file_root)
        #print(bp_c_root)
        if bp_c_root == None and file_root == None:
            print("No bandpass solutions are given")
            # temporary
            g_real = np.full((1,36),1,dtype=np.complex64)
            g_imag = np.full((1,36),0,dtype=np.complex64)
            self.bp_real = None
        elif bp_c_root == None:
            print('Using MIR bandpass solutions')
            times1, g_real = parse_gpplt(file_root+'.gains.real')
            times2, g_imag = parse_gpplt(file_root+'.gains.imag')
            if 1: # temp
                g = g_real + 1j * g_imag
                g = 1/np.conj(g)
                g_real = np.real(g)
                g_imag = np.imag(g)

            assert all(times1 == times2), 'Times in gains real/imag dont match'
            assert g_real.shape == g_imag.shape, 'Unequal shapes of gain files'

            freqs1, bp_real = parse_gpplt(file_root+'.bandpass.real')
            freqs2, bp_imag = parse_gpplt(file_root+'.bandpass.imag')
            assert all(freqs1 == freqs2), 'Freqs in bandpass real/imag dont match'

            assert bp_real.shape == bp_imag.shape, 'Unequal shapes of bandpass files'
            nant = g_real.shape[1]
            assert bp_real.shape[1] == nant, 'Unequal number of antennas in gain and bandpass files'
            self.nant = nant
            self.freqs = np.array(freqs1).astype(np.float) # convert to float
            self.times = times1 # TODO: Parse times.
            if len(times1) > 1:
                warnings.warn('MiriadSolution can only handle 1 time step')
            self.bp_real = bp_real
            bp_imag *= -1  # -1 for 181112, temp
            self.bp_imag = bp_imag
            self.bp_real_interp = [interp1d(self.freqs, bp_real[:, iant], fill_value=(self.bp_real[0,iant], self.bp_real[-1,iant]), bounds_error=False) for iant in range(nant)]
            self.bp_imag_interp = [interp1d(self.freqs, bp_imag[:, iant], fill_value=(self.bp_imag[0,iant], self.bp_imag[-1,iant]), bounds_error=False) for iant in range(nant)]
            self.bp_coeff = None
        else:
            print('Using AIPS bandpass solutions')
            if "polyfit_coeff" in bp_c_root: # AIPS polyfit coefficient
                self.bp_coeff = np.load(bp_c_root)
            else:
                from parse_aips import aipscor
                nant = None
                nfreq = None
                with open(bp_c_root,'r') as fl:
                    for line in fl:
                        if 'NAXIS2' in line:
                            nant = int(line.split()[2])
                            print(f"nant = {nant}")
                            if '190608' in bp_c_root:
                              nant -=2 # exclude ak31 and ak32 from calibration solutions (from Adam)
                        if 'TFDIM11' in line:
                            nfreq = int(line.split()[2])
                if nant is None or nfreq is None:
                    print('WARNING! nant or nfreq not assigned while parsing AIPS bandpass')
                print(f"nant = {nant}")
                fmax = freqs[0]+0.5 # in MHz
                bw = len(freqs) # in MHz
                self.freqs = (-np.arange(float(nfreq))/nfreq*bw+fmax-float(bw)/nfreq/2)/1e3 # reassign freqs in GHz
                self.bp_real = np.full((nfreq,nant),np.nan,dtype=np.complex64)
                self.bp_imag = np.full((nfreq,nant),np.nan,dtype=np.complex64)
                g_real = np.full((1,nant),np.nan,dtype=np.complex64)
                g_imag = np.full((1,nant),np.nan,dtype=np.complex64)

                # look for a README and get fring, selfcal filenames
                drcal = os.path.dirname(bp_c_root)
                import glob
                readme = glob.glob("README*")
                if len(readme) == 1:
                    with open(readme[0],'r') as fl:
                        for line in fl:
                            if "delays" in line and ".sn.txt" in line:
                                fring_f = line.split()[0]
                            if "selfcal" in line and ".sn.txt" in line:
                                sc_f = line.split()[0]
                                print(('FOUND SELFCAL FILE: {}'.format(sc_f)))
                else:
                    print('No or multiple readme file exists for AIPS')
                    fring_f = bp_c_root.replace(bp_c_root.split('/')[-1],"delays.sn.txt")#("bandpasses.bp.txt","delays.sn.txt")
                    sc_f = bp_c_root.replace(bp_c_root.split('/')[-1],"selfcal.sn.txt")
                    #("bandpasses.bp.txt","selfcal.sn.txt")
                    
                specific_frb = None
                if '190608' in bp_c_root:
                  specific_frb = '190608'
                aips_cor = aipscor(fring_f,sc_f,bp_c_root)
                for iant in range(nant):
                    bp = aips_cor.get_phase_bandpass(iant,pol)
                    bp = np.fliplr([bp])[0] # decreasing order
                    

                    # fring delay
                    delta_t_fring_ns = aips_cor.get_delay_fring(iant,pol)*1e9 
                    phases = delta_t_fring_ns * self.freqs
                    phases -= phases[int(len(phases)/2)] # TODO! READ THE REFERENCE FREQUENCY AND SET TO THAT REFERENCE
                    #print(np.shape(phases),np.shape(bp))
                    bp *= np.exp(np.pi*2j*phases, dtype=np.complex64)
                    phase_offset = values.polcal_offset
                    polfring_delay = values.polcal_delay
                    polphases = 2*np.pi*polfring_delay*1e9*self.freqs + phase_offset
                    bp *= np.exp(polphases*1j, dtype=np.complex64)

                    try:
                        g = aips_cor.get_phase_fring(iant,pol)*aips_cor.get_phase_selfcal(iant,pol)
                        g = 1/g # inverse of gain
                    except Exception as e:
                        print(e)
                        g = 0
                    # g = np.conj(g)
                    bp = np.conj(bp)
                    if '180924' in bp_c_root and iant == 15:
                        print('FRB 180924 antenna 15 for AIPS signal negated. fixing this...')
                        bp *= -1
                    if '190102' in bp_c_root:
                        print('FRB 190102 has bandpass magnitude inversed. fixing this...')
                        bp = 1/np.conj(bp)
                    self.bp_real[:,iant] = np.real(bp)
                    self.bp_imag[:,iant] = np.imag(bp)
                    g_real[0,iant] = np.real(g)
                    g_imag[0,iant] = np.imag(g)
                #self.freqs = (-np.arange(2016)/2016.*336+1488.5-1/12.)/1e3 # reassign freqs in GHz
                #self.bp_real, self.bp_imag = parse_aips_bp(bp_c_root, pol)
                print(f"nant = {nant}")
                self.bp_real_interp = [interp1d(self.freqs, self.bp_real[:, iant], fill_value=(self.bp_real[0,iant], self.bp_real[-1,iant]), bounds_error=False) for iant in range(nant)]
                self.bp_imag_interp = [interp1d(self.freqs, self.bp_imag[:, iant], fill_value=(self.bp_imag[0,iant], self.bp_imag[-1,iant]), bounds_error=False) for iant in range(nant)]
                self.bp_coeff = None

        self.g_real = g_real
        self.g_imag = g_imag
        # arrays indexed by antenna index


    def get_solution(self, iant, time, freq_ghz):
        '''
        Get solution including time and bandpass
        iant - antenna index
        time - some version of time. Ignored for now
        freq_ghz - frequency float in Ghz
        '''
        if self.bp_real is None: # this means no bandpass/gain solution was passed
            bp_value = np.array([1])
        elif self.bp_coeff is not None: # Use AIPS polyfit coefficient
            bp_fit = np.poly1d(self.bp_coeff[iant,0,:])+1j*np.poly1d(self.bp_coeff[iant,1,:])
            bp_value = bp_fit(freq_ghz*1e3)
        else: # AIPS polyfit coefficient doesn't exist. Use Miriad/AIPS bandpass interpolation
            print(iant, len(self.bp_real_interp))
            f_real = self.bp_real_interp[iant](freq_ghz)
            f_imag = self.bp_imag_interp[iant](freq_ghz)
            bp_value = f_real + 1j*f_imag

        g_value = self.g_real[0,iant] + 1j*self.g_imag[0,iant]
        total_value = bp_value * g_value
        
        return total_value

    def plot(self):
        '''fig, ax = pylab.subplots(3,3)
        ax = ax.flatten()
        for i in xrange(min(9, self.nant)):
            freq_ghz = self.freqs
            sol = self.get_solution(i, 0, freq_ghz)
            ax[i].plot(freq_ghz, abs(sol))
            ax2 = ax[i].twinx()
            ax2.plot(freq_ghz, np.degrees(np.angle(sol)), 'ro')
            ax[i].set_xlabel('Freq(GHz)')
           '''
        return #added by hyerin

def load_sources(calcfile):
    calc_input = calcfile.replace('.im','.calc')
    d = {}
    for line in open(calc_input, 'r'):
        if len(line) == 0 or line.startswith('#'):
            continue
        bits = line.split(':')
        if len(bits) != 2:
            continue
        
        k,v = bits

        d[k.strip()] = v.strip()

    assert d['NUM SOURCES'] == '1'
    name = d['SOURCE 0 NAME']
    # ra/dec in radians
    ra = float(d['SOURCE 0 RA'])
    dec = float(d['SOURCE 0 DEC'])
    pos = SkyCoord(ra, dec, unit=('rad','rad'), frame='icrs')
    sources = [{'name':name,'ra':pos.ra.deg,'dec':pos.dec.deg}]

    return sources



def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--data', type=str, required=True, help='Directory containing data to process')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--outfile', help='Output fits/.npy file', default='corr.fits')
    parser.add_argument('-c','--channel', type=int, help='Channel to plot', default=0)
    parser.add_argument('-n','--fft-size', type=int, help='Multiple of 64 channels to make channels- default=1', default=1)
    parser.add_argument('-t','--num-threads', type=int, help='Number of threads to run with', default=1)
    parser.add_argument('--calcfile', help='Calc file for fringe rotation')
    parser.add_argument('-w','--hwfile', help='Hw delay file')
    parser.add_argument('-p','--parset', help='Parset for delays')
    parser.add_argument('--show', help='Show plot', action='store_true', default=False)
    parser.add_argument('-i','--nint', help='Number of fine spectra to average', type=int, default=128)
    parser.add_argument('-f','--fscrunch', help='Frequency average by this factor', default=1, type=int)
    parser.add_argument('--rfidelay', type=int, help='Delay in fine samples to add to second component to make an RFI data set', default=0)
    parser.add_argument('--mirsolutions', help='Root file name for miriad gain solutions')
    parser.add_argument('--aips_c', help='AIPS banpass polynomial fit coeffs',default=None)
    parser.add_argument('--an', type=int, help='Specific antenna', default=None)
    parser.add_argument('--offset', type=int, help='FFT offset to add', default=0)
    parser.add_argument('--tab', help='Do tied-array beamforming', action='store_true', default=False)
    parser.add_argument('--ics', help='Do incoherent beamforming', action='store_true', default=False)
    parser.add_argument('--pol', type=str, help='Polarisation to process (x or y)')
    parser.add_argument('--polcal_delay', type=float, default=0, help='Polarisation calibration delay (seconds)')
    parser.add_argument('--polcal_offset', type=float, default=0, help='Polarisation calibration offset')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # get vcraft files from data directory - copied from processTimeStep.py
    antennadirs = sorted(glob.glob(values.data + '/ak*'))
    vcraftfiles = []
    for a in antennadirs:
        beamdirs = sorted(glob.glob(a + "/*"))

        if values.pol == 'x':
            vcraftfiles += sorted(glob.glob(beamdirs[0] + "/*[ac]*vcraft"))
        elif values.pol == 'y':
            vcraftfiles += sorted(glob.glob(beamdirs[1] + "/*[ac]*vcraft"))
        else:
            print((args.pol + ' is not a valid polarisation! Must be x or y'))
            sys.exit(1)

    start = timer()
    sources = load_sources(values.calcfile)
    print(f"load_sources: {timer()-start} s")
    #solutions = None
    #import pdb
    #pdb.set_trace()
    if 0:#values.mirsolutions is not None:
        solutions = MiriadGainSolutions(values.mirsolutions, values, values.aips_c, 'y')
        if 0:
            solutions.plot()
            pylab.show()
        
    # hacking delays
    start = timer()
    delaymap = parse_delays(values)
    antennas = [AntennaSource(mux) for mux in vcraft.mux_by_antenna(vcraftfiles, delaymap)]
    print(f"Parse antennas: {timer()-start} s")
    print(('NUMBER OF ANTENNAS',len(antennas)))
    #antennas = [AntennaSource(vcraft.VcraftFile(f)) for f in values.files]
    given_offset = values.offset
    start = timer()
    corr = Correlator(antennas, sources, values, abs_delay=given_offset)
    print(f"setup Correlator: {timer()-start} s")
    try:
        if values.tab or values.ics:
            import time
            t0 = time.time()
            if values.tab and values.ics:
                print('Confused...you chose both tab and ics. PERFORMING TAB')
            elif values.tab:
                print('PERFORMING TIED-ARRAY BEAMFORMING')
            else:
                print('PERFORMING INCOHERENT BEAMFORMING')
            temp = corr.do_tab(values.an)
            fn = values.outfile
            print(('saving output to '+fn))
            np.save(fn,temp)
        else:
            print('PERFORMING CORRELATION')
            while(corr.running):
                #fq = corr.fileout.fq_table()
                #an = corr.fileout.an_table(corr.ants)
                #su = corr.fileout.su_table(sources)
                corr.do_f()
                corr.do_x()
                corr.next_integration()
    finally:
        if values.tab:
            print(('craftcor.py running time: '+str(time.time()-t0)))
            print('done')
        else:
            corr.fileout.close()




if __name__ == '__main__':
    _main()
