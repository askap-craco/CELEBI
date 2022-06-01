from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np


def _main():
    args = get_args()
    spec = np.load(args.f)
    dedispersed_spec = dedisperse(spec, args.DM, args.f0, args.bw)
    np.save(args.o, dedispersed_spec)


def get_args() -> ArgumentParser:
    """Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Coherently dedisperses given fine spectrum",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", help="Fine spectrum file to dedisperse")
    parser.add_argument(
        "--DM",
        type=float,
        help="Dispersion measure to dedisperse to in pc/cm3",
    )
    parser.add_argument(
        "--f0", type=float, help="Central frequency of observation in MHz"
    )
    parser.add_argument(
        "--bw", type=float, help="Bandwidth of observation in MHz", default=336
    )
    parser.add_argument(
        "-o", help="Output file to save dedispersed spectrum to"
    )
    return parser.parse_args()


def get_freqs(f0: float, bw: float, nchan: int) -> np.ndarray:
    """Create array of frequencies.

    The returned array is the central frequency of `nchan` channels
    centred on `f0` with a bandwidth of `bw`.

    :param f0: Central frequency (arb. units, must be same as `bw`)
    :type f0: float
    :param bw: Bandwidth (arb. units, must be same as `f0`)
    :type bw: float
    :param nchan: Number of channels
    :type nchan: int
    :return: Central frequencies of `nchan` channels centred on `f0`
        over a bandwidth `bw`
    :rtype: :class:`np.ndarray`
    """
    fmin = f0 - bw / 2
    fmax = f0 + bw / 2

    chan_width = bw / nchan

    freqs = np.linspace(fmax, fmin, nchan, endpoint=False) + chan_width / 2

    return freqs


def dedisperse(
    spec: np.ndarray, DM: float, f0: float, bw: float
) -> np.ndarray:
    """
    Coherently dedisperse the given complex spectrum.

    Coherent dedispersion is performed by applying the inverse of the
    transfer function that acts on radiation as it travels through a
    charged medium. This is detailed in Lorimer & Kramer's Handbook of
    Pulsar Astronomy (2005, Cambridge University Press).

    In practice, this is a frequency-dependent rotation of the complex
    spectrum. None of the amplitudes are altered.

    :param spec: Complex 1D-spectrum in a single polarisation
    :type spec: :class:`np.ndarray`
    :param DM: Dispersion measure to dedisperse to (pc/cm3)
    :type DM: float
    :param f0: Central frequency of the spectrum (MHz)
    :type f0: float
    :param bw: Bandwidth of the spectrum (MHz)
    :type bw: float
    :return: Coherently dedispersed complex spectrum
    :rtype: :class:`np.ndarray`
    """
    nchan = spec.shape[0]

    """
    This value of k_DM is not the most precise available. It is used 
    because to alter the commonly-used value would make pulsar timing 
    very difficult. Also, to quote Hobbs, Edwards, and Manchester 2006:
        ...ions and magnetic fields introduce a rather uncertain 
        correction of the order of a part in 10^5 (Spitzer 1962), 
        comparable to the uncertainty in some measured DM values...
    """
    k_DM = 2.41e-4

    freqs = get_freqs(f0, bw, nchan)

    dedisp_phases = np.exp(
        2j * np.pi * DM / k_DM * ((freqs - f0) ** 2 / f0 ** 2 / freqs * 1e6)
    )

    spec *= dedisp_phases

    return spec


if __name__ == "__main__":
    _main()
