# Convert Stokes IQUV dynamic spectra from numpy to filterbank format

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import matplotlib.pyplot as plt
import numpy as np
from your import Your
from your.formats.filwriter import make_sigproc_object

def _main():
    args = get_args()

    data = [np.load(f).T for f in args.infiles]

    file0 = args.infiles[0]
    nsamps, nchans = data[0].shape
    for i in range(1, len(data)):
        file = args.infiles[i]
        assert data[0].shape == data[i].shape, \
            f"MISMATCH IN FILE SHAPES:\n" \
            f"\t{file0}\t{data[0].shape}\n" \
            f"\t{file}\t{data[i].shape}"

    data = np.array(data)           # (pol, t, f)
    data = np.swapaxes(data, 0, 1)  # (t, pol, f)

    nifs = data.shape[1]

    assert (nsamps, nifs, nchans) == data.shape
    print(data.shape)

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
        nbits = 32,
        nifs = nifs,
    )

    sigproc_object.write_header(args.outfile)
    sigproc_object.append_spectra(data, args.outfile)

    # read in to verify
    your_object = Your(args.outfile)
    print(your_object.your_header)
    read_data = np.array([
        your_object.get_data(nstart=0, nsamp=nsamps, pol=i) for i in range(nifs)
    ])
    print(read_data.shape)  # (pol, t, f)

    fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(10, 10))
    for i in range(nifs):
        axs[i][0].imshow(
            data[:,i].T,
            aspect='auto',
            interpolation='none'
        )
        axs[i][0].set_xlim(975, 1050)
        axs[i][1].imshow(
            read_data[i].T,
            aspect='auto',
            interpolation='none'
        )
        axs[i][1].set_xlim(975, 1050)
    plt.tight_layout()
    plt.savefig('test.png')



def get_args() -> Namespace:
    parser = ArgumentParser(
        description="Convert dynamic spectrum from numpy format to "\
            "filterbank format. Assumes the shape of the numpy array "\
            "is (nchans, nsamps), with the first channel being the "\
            "highest in frequency, and the data type being float32.",
        formatter_class=ArgumentDefaultsHelpFormatter
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
    parser.add_argument(
        "infiles",
        type=str,
        nargs='+',
        help="Input .npy files to convert"
    )

    return parser.parse_args()


if __name__ == "__main__":
    _main()