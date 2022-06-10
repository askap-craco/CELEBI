#!/usr/bin/env python
"""
Header utility classes

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import os

from collections import namedtuple, OrderedDict

class DadaHeader(OrderedDict):
    def __init__(self, add_comment=True):
        OrderedDict.__init__(self)
        self.add_comment = add_comment

    def add_card(self, name, value, comment):
        self[name] = (value, comment)

    def __iadd__(self, bits):
        self.add_card(*bits)
        return self

    def get_value(self, item, default=None):
        if item in list(self.keys()):
            return self[item][0]
        else:
            return default

    def set_value(self, item, value):
        comment = self.get_comment(item)
        self[item] = (value, comment)
        return self[item]

    def get_comment(self, item):
        if item in list(self.keys()):
            comment = self[item][1]
        else:
            comment = ''
        return comment
        
    def reset_hdr_size(self, block_size=4096):
        s = self.str(check=False)
        nblocks = len(s)/block_size + 1
        assert nblocks >= 1
        hsize = nblocks*block_size
        assert hsize % block_size == 0
        self.set_value('HDR_SIZE', hsize)
        return hsize

    def str(self, check=True):
        '''Convert to string, without checks'''

        s = ''
        for name, (value, comment) in self.items():
            if self.add_comment and comment is not None:
                s += '%s %s # %s\n' % (name, value, comment)
            else:
                s += '%s %s\n' % (name, value)

        if check:
            hsize = int(self['HDR_SIZE'][0])
            assert len(s) <= hsize, 'Header size is %s but HDR_SIZE=%s' % (len(s), hsize)

        return s

    
    def tofile(self, fout, add_zeros=True):
        self.reset_hdr_size()
        s = str(self)
        fout.write(s)
        hsize = self.get_value('HDR_SIZE')
        nzeros = hsize - len(s)
        assert nzeros >= 0
        if add_zeros:
            fout.write('\x00'*nzeros)

        return hsize

    @staticmethod
    def _fromfile(filename, hdr_size):
        d = DadaHeader()
        with open(filename, 'rU') as fin:
            first_lines = fin.read(hdr_size)
            for iline, line in enumerate(first_lines.split('\n')):
                if ' ' not in line:
                    continue
                
                bits = line.split(None, 1)
                if len(bits) != 2:
                    continue
                    
                name, rest = bits
                valuecomment = rest.split('#', 1)
                if len(valuecomment) >= 2:
                    value, comment = valuecomment
                    comment = comment.strip()
                else:
                    value = valuecomment[0]
                    comment = None
                d[name] = (value.strip(), comment)

        return d


    @staticmethod
    def fromfile(filename, init_size=4096):
        d = DadaHeader._fromfile(filename, init_size)
        if 'HDR_SIZE' in list(d.keys()):
            hdr_size = int(d.get_value('HDR_SIZE'))
        else:
            hdr_size = os.path.getsize(filename)

        if hdr_size != init_size:
            d = DadaHeader._fromfile(filename, hdr_size)

        d['HDR_SIZE'] = (hdr_size, 'Set from file')
        return d

    def __str__(self):
        return self.str(check=True)
        
    __repr__ = __str__


class CraftHeader(DadaHeader):
    def __init__(self):
        DadaHeader.__init__(self)

def default_header(source, obsid, sbid, scanid, captureid, freq, bw, npol, nbeam, nchan, nbit, tsamp, ra, dec, ras, decs, beam_pol, beam_id, antno, sum_weights, int_cycles, int_time, bfaddr, chanmaps, freqmaps, file_size=1<<32):
    hdr = Craftheader()
    hdr += ('HDR_VERSION', 1.0, 'Version of the header')
    hdr += ('HDR_SIZE', 4096, 'Size of header')
    hdr += ('SOURCE', source, 'Observing source')
    hdr += ('OBS_ID', obsid, 'OBS_ID - see PSRDADA??')
    hdr += ('SBID',sbid, 'ASKAP Schedblock id')
    hdr += ('SCANID', scanid, 'ASKAP scan ID')
    hdr += ('CAPTUREID', captureid, 'Capture ID - within a scan')
    hdr += ('FREQ', freq, 'Observing frequency start (MHz)')
    hdr += ('BW', bw, 'Channel bandwidth (MHz)')
    hdr += ('NPOL', npol, 'Number of polarisations')
    hdr += ('NBEAM', nbeam, 'Number of beams')
    hdr += ('NCHAN', nchan, 'Number of channels')
    hdr += ('NBIT', nbit, 'Number of bits')
    hdr += ('DTYPE', '<f4', 'Data type (see numpy dtype for meaning)')
    hdr += ('DORDER', 'TFBP', 'Data ordering. Last is fastest. T=Time, B=Beam, P=Polarsation, F=Freq')
    hdr += ('TSAMP', tsamp, 'Sampling time')
    hdr += ('RA', ra, 'Right Ascension of pointing (J2000 - degrees)')
    hdr += ('DEC', dec, 'Declination of pointing (J2000 - degrees)')
    hdr += ('BEAM_RA', ','.join(map(str, ras)), 'Right Ascensions of all beams (J2000 - degrees)')
    hdr += ('BEAM_DEC', ','.join(map(str, decs)), 'Declinations of all beams (J2000 - degrees)')
    hdr += ('BEAM_POL', ','.join(map(str, beam_pol)), 'Polarisation mnemonic for beams. e.g. XX or YY')
    hdr += ('BEAM_ID', ','.join(map(str, beam_id)), 'Beam numbers. 0 based')
    hdr += ('TELESCOPE', 'ASKAP', 'Telescope name')
    hdr += ('INSTRUMENT', 'CRAFT', 'Name of instrument')
    hdr += ('ANTNO', antno, 'Antenna number')
    hdr += ('SUM_WEIGHTS', ','.join(map(str, sum_weights)), 'Weights of each antenna going into an incoherent sum')
    hdr += ('ANTNAME', 'AK%02d'%antno, 'Antenna name')
    hdr += ('INT_CYCLES', int_cycles, 'Number of integrations per beamformer packet')
    hdr += ('INT_TIME', int_time, 'Number of samples per integration')
    hdr += ('OBS_OVERLAP', 0, 'The amount by which neighbouring fiels overlap - See PSRDADA')
    hdr += ('OBS_OFFSET', 0, 'The number of byte from the start of the observation - See PSRDADA')
    #hdr += ('RESOLUTION', resolution, '???? The number of bytes in a single packet/chunk? - see dada_dbdisk.c???')
    #hdr += ('BYTES_PER_SECOND', bytes_per_second, 'Number of bytes per second (useful for dada_dbdisk)')
    hdr += ('FILE_SIZE', file_size, 'Number bytes to write to a file before opening another one')


    hdr += ('NUM_BEAMFORMERS', len(bfaddr), 'Number of beamformers')
    #chanmap = np.arange(nchan)

    for ibf, (bf, chanmap, freqmap) in enumerate(zip(bfaddr, chanmaps, freqmaps)):
        #freqmap = chanmap*bw + fstart + nchan*bw*ibf
        hdr += ('BEAMFORMER%d_ADDR'%ibf, bf, 'Address of beamformer %d' % ibf)
        hdr += ('BEAMFORMER%d_CHANMAP'%ibf, ','.join(map(str, chanmap)), 'Channel map for beamformerr %d' % ibf)
        hdr += ('BEAMFORMER%d_FREQMAP'%ibf, ','.join(map(str, freqmap)), 'Frequency map for beamformerr %d' % ibf)


    return hdr

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Captures CRAFT autocorrelation spectra and writes to SIGPROC file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-a', '--ant', help='Antenna number', type=int, default=-1)
    parser.add_argument('-c','--card', help='Card number', type=int, default=-1)
    parser.add_argument('-o','--output', help='Output file name')
    parser.add_argument('-i','--int-time', help='Integration time (1-65535)', default=1000, type=int)
    parser.add_argument('-s','--int-cycles', help='Number of cycles to combine (1-7)', default=3, type=int)
    parser.add_argument('-b','--beam-number', help='Beam number to save voltages for (0-35). Unspecified for all beams', default=None)
    parser.add_argument('-m','--buffer-mode', help='Buffer number of bits', choices=[16, 8, 4, 1], type=int, default=16)
    parser.add_argument('bfaddr', nargs='+', help='Beamformer addresses')
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    int_time = int(values.int_time)
    assert int_time >=1 and int_time<= 65535, 'Invalid integration time {}'.format(int_time)
    
    int_cycles = int(values.int_cycles)
    assert int_cycles >=1 and int_cycles <= 7, 'Invalid int cycles {}'.format(int_cycles)

    samp_rate = 32./27.*1e6 # samples/sec
    tsamp = float(int_time)/samp_rate

    bufmode_map = {16:0, 8:1, 4:2, 1:3}
    assert values.buffer_mode in list(bufmode_map.keys()), 'Invalid buffer mode'
    bufmode = bufmode_map[values.buffer_mode]

    if values.beam_number is None:
        beam_number = -1
    else:
        beam_number = int(values.beam_number)
        assert beam_number >= 0 and beam_number < 36, 'Invalid beam number to save {}'.format(beam_number)
        bufmode += 4 # This tells the firmware to save only 1 beam
    
    antno = int(values.ant)
    cardno = int(values.card)
    bfaddr = values.bfaddr
    ra = 15.0
    dec = -30.0
    fstart = 1200.0
    nbit = 32
    nchan = 48*len(bfaddr)
    nbeams = 72
    ras = np.linspace(ra, ra+5, 36)

    # Repeat same RA for each 2 beams
    ras = np.vstack((ras, ras)).T.flatten()

    decs = np.linspace(dec, dec+5, 36)

    # repeat same dec for each 2 beams
    decs = np.vstack((decs, decs)).T.flatten()
    
    pols = ['XX','YY']
    beam_pol = pols*(nbeams/2)
    beam_id = np.arange(nbeams)

    nant = 12
    bw = -1.0

    # only a single antenna in the sum
    sum_mask = 1<<antno
    sum_weights = np.zeros(nant+1)
    sum_weights[antno] = 1.
    
    bytes_per_second = int((nbit/8)*nchan*nbeams/tsamp)
    resolution = (nbit/8)*nchan*nbeams
    
    for (h, v, c) in hdr:
        print(h, v)


if __name__ == '__main__':
    _main()
