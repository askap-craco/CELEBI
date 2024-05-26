import os
import pandas as pd
import numpy as np
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from astropy.coordinates import SkyCoord as sc
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus


def _main():
    parser = ArgumentParser(
        description="Lookup source position in RACS",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-l", "--localracspath", default="", 
                        help="Use this local RACS catalog instead of CASDA")
    parser.add_argument("--racsrasystematic", type=int, default=0.0, 
                        help="Add this value in quadrature to the final RA uncertainty")
    parser.add_argument("--racsdecsystematic", type=int, default=0.0,
                        help="Add this value in quadrature to the final Dec uncertainty")
    parser.add_argument("-o", required=True, help="Output RACS positions file")
    parser.add_argument("-a", required=True, help="Output ASKAP positions file")
    parser.add_argument("-n", required=True, help="Output names file")
    parser.add_argument("-r", required=True, help="Output region file")
    parser.add_argument("-j", required=True, help="Output file with list of jmfits")
    parser.add_argument("files", nargs="+", help="Source stats files")
    args = parser.parse_args()

    print(args.files)

    askap_file = open(args.a, "a")
    pos_file = open(args.o, "a")
    name_file = open(args.n, "a")
    region_file = open(args.r, "a")
    jmlist_file = open(args.j, "a")

    names = []

    if args.localracspath != "":
        print("Reading local catalogue at ", args.localracspath)
        if not os.path.exists(args.localracspath):
            print("This path does not exist! Aborting.")
            sys.exit()
        cat = pd.read_csv(args.localracspath)
        print("Local catalog read successfully")
    else:
        print("Opening casdatap")
        casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
        print("casdatap open")

    for f in args.files:
        print(f)
        coord = Coord(f)

        if coord.sn < 7:
            continue

        if args.localracspath != "":
            t1 = RACS_lookup_local(coord.ra_hms, coord.dec_dms, cat)
        else:
            t1 = RACS_lookup1(coord.ra_hms, coord.dec_dms, casdatap)

        print(t1)

        table = t1
        '''
        if len(t1) > 0:
            table = t1
        else:
            # fall back on less robust query
            # e.g. if source is near the edge of the RACS footprint
            table = RACS_lookup2(coord.ra_hms, coord.dec_dms, casdatap)
        '''
        if len(table) == 0:
            continue

        # if more than one source: skip
        if len(table) > 1:
            continue
        
        source = table[0]

        print(source)

        # highest precision position is with the deg coords
        try:
            source_sc = sc(
                source["ra"], source["dec"], unit="deg"
            )
        except KeyError:
            source_sc = sc(
                source["ra_deg_cont"], source["dec_deg_cont"], unit="deg"
            )
        ra_hms, dec_dms = source_sc.to_string("hmsdms").split()

        askap_file.write(
            writestr(
                coord.ra_hms, coord.ra_err, coord.dec_dms, coord.dec_err
            )
        )

        try:
            pos_file.write(
                writestr(
                    ra_hms, max(float(source["e_ra"]), 0.01), dec_dms, max(float(source["e_dec"]), 0.01)
                )
            )
        except KeyError:
            pos_file.write(
                writestr(
                    ra_hms, max(float(source["ra_err"]), 0.01), dec_dms, max(float(source["dec_err"]), 0.01)
                )
            )
        
        try:
            name = source["source_id"]
        except KeyError:
            name = source["component_name"]

        name_file.write(name + "\n")

        region_str = f'{coord.ra_hms}, {coord.dec_dms}, {coord.ra_err}", {coord.dec_err}", 0'
        region_file.write(
            f'fk5;ellipse({region_str.replace("h", ":").replace("d", ":").replace("m", ":").replace("s", ":")}) # text="{name} - S/N={coord.sn:.2f}"\n'
        )
        jmlist_file.write(f)
        jmlist_file.write("\n")

    askap_file.close()
    pos_file.close()
    name_file.close()
    jmlist_file.close()


class Coord:
    def __init__(self, stats):
        fields = {}
        f = open(stats)
        for line in f:
            line = line[:-1].split(": ")  # [:-1] to trim newline
            fields[line[0]] = line[1].strip()  # strip extra whitespace
        f.close()

        self.ra_hms = fields["Actual RA"]  # hms
        self.ra_err = max(float(fields["Est. RA error (mas)"]) / 1e3, 0.01)  # arcseconds
        self.dec_dms = fields["Actual Dec"]  # dms
        self.dec_err = max(float(fields["Est. Dec error (mas)"]) / 1e3, 0.01) # arcseconds
        self.sn = float(fields["S/N"])


def RACS_lookup1(ra_hms, dec_dms, casdatap):
    c = sc(ra_hms, dec_dms, unit="hour,deg")
    print(f"Looking up {c}")
    job1 = casdatap.launch_job_async(
        f"SELECT * FROM AS110.racs_dr1_gaussians_galacticcut_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {c.ra.value},{c.dec.value},0.0014))",
    )
    return job1.get_results()


def RACS_lookup2(ra_hms, dec_dms, casdatap):
    c = sc(ra_hms, dec_dms, unit="hour,deg")
    print(f"Looking up {c}")
    job2 = casdatap.launch_job_async(
        f"SELECT * FROM casda.continuum_component where 1=CONTAINS(POINT('ICRS', ra_deg_cont, dec_deg_cont),CIRCLE('ICRS',{c.ra.value},{c.dec.value},0.0014)) and project_id = 23",
    )
    return job2.get_results()

def RACS_lookup_local(ra_hms, dec_dms, cat, radius=0.0014):
    """
    Query the RACS local catalog based on the given coordinates and radius.

    Parameters:
    ra (float): Right Ascension in degrees.
    dec (float): Declination in degrees.
    radius (float, optional): Search radius in degrees. Default is 5 arcsec.
    catpath (str, optional): Path to the catalog file. Default is an empty string.

    Returns:
    astropy.Table: Subset of the catalog containing sources within the specified radius.
    """
    ra = sc(ra_hms, dec_dms, unit="hour,deg").ra.deg
    dec = sc(ra_hms, dec_dms, unit="hour,deg").dec.deg
    
    cat_ra, cat_dec = np.array(cat['ra']), np.array(cat['dec'])
    
    phi1 = ra * np.pi / 180
    theta1 = dec * np.pi / 180
    phi2 = cat_ra * np.pi / 180
    theta2 = cat_dec * np.pi / 180
    
    cos_sep_radian = np.sin(theta1) * np.sin(theta2) + np.cos(theta1) * np.cos(theta2) * np.cos(phi1-phi2)
    
    sep = np.arccos(cos_sep_radian) * 180 / np.pi
    select_bool = sep < radius
    
    return Table.from_pandas(cat.iloc[select_bool])


def writestr(ra_hms, ra_err, dec_dms, dec_err):
    return f"{ra_hms},{ra_err},{dec_dms},{dec_err}\n"


if __name__ == "__main__":
    _main()
