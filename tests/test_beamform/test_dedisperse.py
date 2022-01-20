#!/usr/bin/env python

from argparse import Namespace

import numpy as np

from craco_postproc.beamform import dedisperse


def test_get_args():
    default_namespace = Namespace(
        f=None,
        DM=None,
        f0=None,
        bw=336,
        o=None,
    )

    assert dedisperse.get_args() == default_namespace


def test_get_freqs():
    # Case 1: 8 channels, centred on 1000 MHz, bandwidth 8 MHz
    f0 = 1000
    bw = 8
    nchan = 8
    freqs = np.array(
        [996.5, 997.5, 998.5, 999.5, 1000.5, 1001.5, 1002.5, 1003.5]
    )

    assert np.all(dedisperse.get_freqs(f0, bw, nchan) == freqs)
