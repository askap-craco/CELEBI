# Convert Stokes IQUV dynamic spectra from numpy to filterbank format

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import matplotlib.pyplot as plt
import numpy as np
from your import Your
from your.formats.filwriter import make_sigproc_object

def _main():
    args = get_args()

    data = np.load(args.infile[0]).T

    nsamps, nchans = data.shape
    tsamp = args.tsamp / 1e6    # us -> s

    fmax = args.f0 + args.bw/2

    sigproc_object = make_sigproc_object(
        rawdatafile = args.outfile,
        source_name = args.source,
        nchans = nchans,
        foff = args.df,
        fch1 = fmax,
        tsamp = tsamp,
        tstart = args.tstart,
        src_raj = args.ra,
        src_dej = args.dec,
        nbits = 32
    )

    sigproc_object.write_header(args.outfile)
    print(data.shape)
    sigproc_object.append_spectra(data, args.outfile)


def get_args() -> Namespace:
    parser = ArgumentParser(
        description="Convert dynamic spectrum from numpy format to "\
            "filterbank format. Assumes the shape of the numpy array "\
            "is (nchans, nsamps), with the first channel being the "\
            "highest in frequency, and the data type being float32.",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "infile",
        type=str,
        nargs=1,
        help="Input .npy file to convert"
    )
    parser.add_argument(
        "-s", "--source",
        type=str,
        required=True,
        help="Source name"
    )
    parser.add_argument(
        "--tsamp",
        type=float,
        required=True,
        help="Sampling time/time resolution of dynamic spectra in us"
    )
    parser.add_argument(
        "--tstart",
        type=float,
        required=True,
        help="Start time of data in MJD"
    )
    parser.add_argument(
        "--f0",
        type=float,
        required=True,
        help="Centre frequency of dynamic spectrum in MHz"
    )
    parser.add_argument(
        "-b", "--bw",
        type=float,
        default=336,
        help="Bandwidth of dynamic spectrum in MHz"
    )
    parser.add_argument(
        "--df",
        type=float,
        default=-1,
        help="Channel step size in MHz. Assumes first channel is " \
             "highest frequency, so should be negative."
    )
    parser.add_argument(
        "-r", "--ra",
        type=float,
        required=True,
        help="Right ascension of source as HHMMSS.SS"
    )
    parser.add_argument(
        "-d", "--dec",
        type=float,
        required=True,
        help="Declination of source as DDMMSS.SS"
    )
    parser.add_argument(
        "-o", "--outfile",
        type=str,
        required=True,
        help="Output .fil file"
    )
    return parser.parse_args()


if __name__ == "__main__":
    _main()