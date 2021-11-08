import os
import sys
from datetime import datetime
from pathlib import Path

import espressolib
import pandas as pd


def parse_eops(eops_file):
    col_names = [
        "jd",
        "xpole",
        "ypole",
        "ut1_tai",
        "xpole_err",
        "ypole_err",
        "ut1_tai_err",
    ]
    return pd.read_csv(
        eops_file,
        delim_whitespace=True,
        comment="#",
        skiprows=1,
        names=col_names,
        usecols=list(range(7)),
    )


class EOPHandler:
    def __init__(self, eops_file=None, leapsec_file=None):
        if eops_file is None:
            eops_file = os.getenv("DIFX_EOPS")

        if leapsec_file is None:
            leapsec_file = os.getenv("DIFX_UT1LS")

        assert eops_file is not None, "You must set $DIFX_EOPS or pass a file"
        assert (
            leapsec_file is not None
        ), "You must set $DIFX_UT1LS or pass a file"

        self._eops_file = eops_file
        self._leapsec_file = leapsec_file

        self._eops = parse_eops(self._eops_file)
        self.logger.info(f"Reading EOP data from {self._eops_file}")
        self._leapsec = parse_eops(self._leapsec_file)
        self.logger.info(f"Reading leap second data from {self._leapsec_file}")

    @property
    def eops_file(self):
        return self._eops_file

    @property
    def leapsec_file(self):
        return self._leapsec_file

    @property
    def eop_update_time(self):
        ctime = datetime.fromtimestamp(Path(self.eops_file).stat().st_ctime)
        return ctime.isoformat(sep=" ", timespec="seconds")

    @property
    def eops(self):
        return self._eops

    @property
    def v2d_header(self):
        return f"# EOPS from {self.eops_file} last updated at {self.eop_update_time}"

    @property
    def v2d_format(self):
        return "EOP {mjd:d} {{ xPole={xp:f} yPole={yp:f} tai_utc={tai_utc:d} ut1_utc={ut1_utc:f} }}"

    def eops(self):
        mjd_jd = 2400000.5
        # get target MJD
        target_mjd = args[0]
        # convert to MJD if necessary
        target_mjd = espressolib.convertdate(target_mjd, "mjd")
        target_jd = round(target_mjd) + mjd_jd

        # dates before June 1979 not valid (earliest EOPs)
        if target_jd < 2444055.5:
            raise Exception("Date too early. No valid EOPs before July 1979")

    def get_leapsec(leapsec_page, target_jd):
        """parse the leap seconds page"""

        tai_utc = None
        for line in leapsec_page:
            linedate = float(line[17:27])
            if linedate > target_jd:
                break
            else:
                tai_utc = int(float(line[38:49]))
        return tai_utc


# parse the eop page
neop = -1
nlines = None
for nlines, line in enumerate(eop_page):
    # skip the first line, which isn't commented
    if nlines == 0:
        continue
    if not line:
        continue
    # strip comments
    if line[0] == "#":
        continue
    # split the line on whitespace and convert to floats
    eop_fields = [float(field) for field in line.split()]
    # print an EOP line if we're within 3 days of the target day
    if abs(eop_fields[0] - target_jd) < 3:
        neop += 1
        tai_utc = get_leapsec(leapsec_page, eop_fields[0])
        if tai_utc is None:
            raise Exception("Leap seconds not found! Check your UT1LS file")
        xpole = eop_fields[1] / 10.0
        ypole = eop_fields[2] / 10.0
        ut1_utc = tai_utc + eop_fields[3] / 1.0e6
        eopdate = int(eop_fields[0] - mjd_jd)
        print(eopformat.format(eopdate, xpole, ypole, tai_utc, ut1_utc, neop))

sys.stderr.write(f"Processed {nlines:d} lines\n")
