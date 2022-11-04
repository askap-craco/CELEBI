import numpy as np
from argparse import ArgumentParser

from requests import head

def _main():
    args = parse_args()

    ds = normalise(np.load(args.ds))
    t = np.load(args.t)
    freqs = get_freqs(args.f0, args.bw, ds.shape[0])

    #parse snoopy
    cand = np.loadtxt(args.snoopy, delimiter=" ", comments="#")
    DM_cand = cand[5]

    DMs = np.arange(
        DM_cand - args.DMrange, DM_cand + args.DMrange, args.DMstep
    )

    widths = np.arange(1, 11)
    
    SN_peak, t_peak, DM_peak, width_peak = incoh_search(ds, freqs, DMs, widths)

    print("Original -> refined:")
    print(f"DM:\t{DM_cand} -> {DM_peak}")
    print(f"mjd:\t{cand[7]} -> {t[t_peak]}")
    print(f"width:\t{cand[3]} -> {width_peak}")

    new_cand = cand.copy()
    new_cand[3] = width_peak
    new_cand[5] = DM_peak
    new_cand[7] = t[t_peak]

    np.savetxt(args.o, new_cand, delimiter=" ", 
        header="# S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, latency_ms"
    )


def parse_args():
    parser = ArgumentParser(
        description="Searches an ICS dynamic spectrum for an FRB and "\
                    "edit snoopy file to reflect optimal values. DM "\
                    "range for search is centred on snoopy candidate's "\
                    "DM."
    )

    parser.add_argument(
        "--ds",
        type=str,
        help="Dynamic spectrum to search filename"
    )
    parser.add_argument(
        "--DMrange",
        type=float,
        default=10.0,
        help="How far above and below the snoopy candidate's DM to search"
    )
    parser.add_argument(
        "--DMstep",
        type=float,
        default=0.01,
        help="DM search step size"
    )
    parser.add_argument(
        "-s", "--snoopy",
        type=str,
        help="Snoopy candidate filename"
    )
    parser.add_argument(
        "-f", "--f0",
        type=float,
        help="Dynamic spectrum central frequency"
    )
    parser.add_argument(
        "-t",
        type=str,
        help="Time (MJD) axis filename"
    )
    parser.add_argument(
        "-o",
        type=str,
        help="Refined candidate output filename"
    )
    parser.add_argument(
        "-b", "--bw",
        type=float,
        default=336.0,
        help="Dynamic spectrum bandwidth"
    )
    parser.add_argument()


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


def normalise(ds):
    """
    Normalise a dynamic spectrum along frequency channels to cancel out RFI
    """
    F, T = ds.shape
    out_ds = ds.copy()
    norm = lambda x: (x - np.mean(x)) / np.std(x)
    for f in range(F):
        out_ds[f] = norm(ds[f])
    return out_ds


def incoh_dedisp(ds, freqs, DM):
    """
    Incoherently dedisperse a dynamic spectrum. Assumes 1 ms time resolution
    """
    ds_dd = ds.copy()

    # reference frequency 1 MHz higher than lowest channel for consistency with polyco
    f_ref = np.min(freqs) + 1

    # index to roll for dedispersion
    k_dm = 1 / 2.41e-4
    idt = lambda DM, f: -int(k_dm * (f ** (-2) - f_ref ** (-2)) * DM * 1000)

    for i, f in enumerate(freqs[::-1]):
        ds_dd[i] = np.roll(ds[i], idt(DM, f))
    return ds_dd


def incoh_search(ds, freqs, DMs, widths):
    SNs = np.zeros((ds.shape[1], DMs.shape[0], widths.shape[0]))
    for d, DM in enumerate(DMs):
        ds_dd = incoh_dedisp(ds, freqs, DM)
        prof = np.sum(ds_dd, axis=0)
        for w, width in enumerate(widths):
            smth_prof = np.convolve(prof, np.ones(width), mode="same")
            SNs[:, d, w] = smth_prof / np.std(smth_prof)

    t_peak, d_peak, w_peak = np.unravel_index(np.argmax(SNs), SNs.shape)
    SN_peak = SNs[t_peak, d_peak, w_peak]
    return SN_peak, t_peak, DMs[d_peak], widths[w_peak]



if __name__ == "__main__":
    _main()