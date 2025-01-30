import sys
import os
import pandas as pd
import numpy as np
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from astropy.coordinates import SkyCoord as sc
from astropy.table import Table
from astropy.table import vstack
from astropy import units as u

#from astroquery.utils.tap.core import TapPlus
#from astroquery.vizier import Vizier


def _main():
    parser = ArgumentParser(
        description="Lookup source position in RACS",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--referencecatalog", default="RACS", type=str,
                        help="Reference catalog to use [RACS or VLASS currently supported]")
    parser.add_argument("--localracsgausspath", default="",
                        help="Use this local RACS catalog (gaussian components) instead of CASDA")
    parser.add_argument("--localracssourcepath", default="",
                        help="Use this local RACS catalog (source components) instead of CASDA")
    parser.add_argument("--localvlasspath", default="",
                        help="Use this local VLASS catalog instead of Vizier")
    parser.add_argument("--racsrasystematic", type=int, default=0.0,
                        help="Add this value in quadrature to the final RA uncertainty")
    parser.add_argument("--racsdecsystematic", type=int, default=0.0,
                        help="Add this value in quadrature to the final Dec uncertainty")
    parser.add_argument("-o", required=True, help="Output reference positions file")
    parser.add_argument("-a", required=True, help="Output ASKAP positions file")
    parser.add_argument("-n", required=True, help="Output names file")
    parser.add_argument("-r", required=True, help="Output region file")
    parser.add_argument("-j", required=True, help="Output file with list of jmfits")
    parser.add_argument("files", nargs="+", help="Source stats files")
    args = parser.parse_args()

    print(args.files)

    askap_file = open(args.a, "w")
    pos_file = open(args.o, "w")
    name_file = open(args.n, "w")
    region_file = open(args.r, "w")
    jmlist_file = open(args.j, "w")

    names = []
    
    localgausscatalogue = None
    localsourcecatalogue = None
    localvlasscatalogue = None
    if args.referencecatalog != "RACS" and args.referencecatalog != "VLASS":
        print("args.referencecatalog was set to {0}, must be either RACS or VLASS".format(args.referencecatalog))
        sys.exit()
    if args.localracsgausspath != "":
        print("Reading local catalogue at ", args.localracsgausspath)
        if not os.path.exists(args.localracsgausspath):
            print("This path does not exist! Aborting.")
            sys.exit()
        localgausscatalogue = pd.read_csv(args.localracsgausspath)
        print("Local gaussians catalog read successfully")
    if args.localracssourcepath != "":
        print("Reading local catalogue at ", args.localracssourcepath)
        if not os.path.exists(args.localracssourcepath):
            print("This path does not exist! Aborting.")
            sys.exit()
        localsourcecatalogue = pd.read_csv(args.localracssourcepath)
        print("Local sources catalog read successfully")
    if args.referencecatalog == "RACS" and args.localracsgausspath == "" and args.localracssourcepath == "":
        print("Opening casdatap")
        casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
        print("casdatap open")
    if args.localvlasspath != "":
        print("Reading local catalogue at ", args.localvlasspath)
        if not os.path.exists(args.localvlasspath):
            print("This path does not exist! Aborting.")
            sys.exit()
        localvlasscatalogue = pd.read_csv(args.localvlasspath)

    for f in args.files:
        print(f)
        coord = Coord(f)

        if coord.sn < 7:
            continue

        racs_avail = True # Flag to indicate whether RACS is available in the region
        if float(coord.dec_dms.split(':')[0]) > 29:
            racs_avail = False
        if args.referencecatalog == "VLASS": # Going to use VLASS as primary:
            # First check the declination
            if float(coord.dec_dms.split(':')[0]) < -45:
                print("VLASS is not available below -45 degrees, can't use!")
                print("You should probably just use --referencecatalog=RACS")
                sys.exit()
            t1 = VLASS_lookup_local(coord.ra_hms, coord.dec_dms, localvlasscatalogue)
        else: # Going to use RACS
            if args.localracssourcepath == "": # No local source catalogue specified
                if args.localracsgausspath != "": # But there is a local gaussians catalogue - use that
                    t1 = RACS_lookup_local(coord.ra_hms, coord.dec_dms, localgausscatalogue)
                else: # No local catalogues supplied at all - use CASDA
                    t1 = RACS_lookup1(coord.ra_hms, coord.dec_dms, casdatap)
            else: # A local source catalogue was supplied - use that as the primary
                t1 = RACS_lookup_local(coord.ra_hms, coord.dec_dms, localsourcecatalogue)

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
        
        # Filter out sources that are likely unacceptably resolved, if we had both gauss and source catalogues
        if args.referencecatalog == "VLASS":
            if not VLASS_compact(table):
                print("Skipping source {0} that is not compact in VLASS".format(f))
                continue
        else: # Using RACS
            if args.localracssourcepath != "" and args.localracsgausspath != "":
                if RACS_exclude_resolved(localgausscatalogue, table):
                    continue
            # If we don't have both source and local gaussians, can't do any compactness exclusion
            # An improvement would be to use VLASS if available
        
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
            name = source["Component_name"]

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
    region_file.close()
    jmlist_file.close()
    
    indices = find_duplicates(args.n)
    
    remove_duplicates(args.n, indices)
    remove_duplicates(args.a, indices)
    remove_duplicates(args.o, indices)
    remove_duplicates(args.r, indices)
    remove_duplicates(args.j, indices)


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

def VLASS_lookup(ra_hms, dec_dms, radius=0.0014):
    """
    Performs a VLASS lookup for a given celestial position.

    Args:
        ra_hms (str): Right Ascension in hours, minutes, and seconds (e.g., '12h34m56s').
        dec_dms (str): Declination in degrees, minutes, and seconds (e.g., '12d34m56s').
        radius (float): Search radius in degrees (default is 0.0014 degrees).

    Returns:
        job (Table): A table containing the VLASS lookup results
    """
    try:
        c = sc(ra_hms, dec_dms, unit="hour,deg")
        print(f"Looking up {c}")
        v = Vizier(columns=['CompName', 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'Ftot', 'e_Ftot', 'Fpeak', 'e_Fpeak', 'Maj', 'e_Maj', 'Min', 'e_Min', 'PA', 'e_PA'], row_limit=-1)
        job = v.query_region(c, radius=radius*u.degree, catalog='VLASS')[0]
    except:
        job = Table(names=['CompName', 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'Ftot', 'e_Ftot', 'Fpeak', 'e_Fpeak', 'Maj', 'e_Maj', 'Min', 'e_Min', 'PA', 'e_PA'])
    
    job.rename_column('CompName', 'Component_name')
    job.rename_column('RAJ2000', 'ra')
    job.rename_column('DEJ2000', 'dec')
    job.rename_column('e_RAJ2000', 'ra_err')
    job.rename_column('e_DEJ2000', 'dec_err')
    
    # Add 0.25 arcsec to all rows in the dec column to compensate for the known VLASS offset
    job['dec'] -= 0.25 / 3600

    return job

def VLASS_lookup_local(ra_hms, dec_dms, cat, radius=0.0014):
    """
    Perform a VLASS lookup based on the given coordinates and catalog.

    Parameters:
    ra_hms (str): Right Ascension in hours, minutes, and seconds.
    dec_dms (str): Declination in degrees, minutes, and seconds.
    cat (pandas.DataFrame): Catalog containing RA and DEC columns.
    radius (float, optional): Search radius in degrees. Default is 0.0014 degrees.

    Returns:
    astropy.table.Table: Subset of the input catalog containing selected rows.
    """
    ra = sc(ra_hms, dec_dms, unit="hour,deg").ra.deg
    dec = sc(ra_hms, dec_dms, unit="hour,deg").dec.deg
    
    cat_ra, cat_dec = np.array(cat['RA']), np.array(cat['DEC'])
    
    phi1 = ra * np.pi / 180
    theta1 = dec * np.pi / 180
    phi2 = cat_ra * np.pi / 180
    theta2 = cat_dec * np.pi / 180
    
    cos_sep_radian = np.sin(theta1) * np.sin(theta2) + np.cos(theta1) * np.cos(theta2) * np.cos(phi1-phi2)
    
    sep = np.arccos(cos_sep_radian) * 180 / np.pi
    select_bool = sep < radius
    
    vlass_table = Table.from_pandas(cat.iloc[select_bool])
    
    vlass_table.rename_column('RA', 'ra')
    vlass_table.rename_column('DEC', 'dec')
    vlass_table.rename_column('E_RA', 'ra_err')
    vlass_table.rename_column('E_DEC', 'dec_err')
    
    return vlass_table

def VLASS_compact(table):
    if (table[0]['Total_flux'] / table[0]['Peak_flux']) < 1.5:
        return True
    else:
        return False

def RACS_lookup1(ra_hms, dec_dms, casdatap):
    c = sc(ra_hms, dec_dms, unit="hour,deg")
    print(f"Looking up {c}")
    job1 = casdatap.launch_job_async(
        f"SELECT * FROM AS110.racs_dr1_gaussians_galacticcut_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {c.ra.value},{c.dec.value},0.0014))",
    )
    job2 = casdatap.launch_job_async(
        f"SELECT * FROM AS110.racs_dr1_gaussians_galacticregion_v2021_08_v02 WHERE 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {c.ra.value},{c.dec.value},0.0014))",
    )
    return vstack([job1.get_results(), job2.get_results()])


def RACS_lookup2(ra_hms, dec_dms, casdatap):
    c = sc(ra_hms, dec_dms, unit="hour,deg")
    print(f"Looking up {c}")
    job2 = casdatap.launch_job_async(
        f"SELECT * FROM casda.continuum_component where 1=CONTAINS(POINT('ICRS', ra_deg_cont, dec_deg_cont),CIRCLE('ICRS',{c.ra.value},{c.dec.value},0.0014)) and project_id = 23",
    )
    return job2.get_results()

def RACS_lookup_local(ra_hms, dec_dms, cat, radius=0.0014):
    """
    Perform a local RACS lookup based on the given coordinates and catalog.

    Parameters:
    ra_hms (str): Right Ascension in hours, minutes, and seconds (e.g., '12:34:56').
    dec_dms (str): Declination in degrees, minutes, and seconds (e.g., '+12:34:56').
    cat (pandas.DataFrame): Catalog containing 'ra' and 'dec' columns.
    radius (float): Search radius in degrees. Default is 0.0014 degrees.

    Returns:
    astropy.table.Table: Subset of the input catalog containing rows within the search radius.
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

def RACS_exclude_resolved(cat, source_table):
    """
    Compare the Gaussian parameters of a catalog entry with a source table entry.
    
    Parameters:
    - cat: pandas DataFrame
        The catalog containing Gaussian parameters.
    - source_table: pandas DataFrame
        The source table containing source information.
    
    Returns:
    - bool
        True if the Gaussian parameters meet the specified conditions, False otherwise.
    """
    cat_gaus = cat[(cat['source_id'] == source_table['source_id'][0])]
    
    source_coords = sc(source_table['ra'], source_table['dec'], unit='deg')
    gaus_coords = sc(cat_gaus['ra'], cat_gaus['dec'], unit='deg')
    
    sep = source_coords.separation(gaus_coords).arcsec
    cat_gaus['separation'] = sep
    
    print(cat_gaus)
    cat_test = cat_gaus[(cat_gaus['total_flux_gaussian'] > 0.2 * cat_gaus['total_flux_source']) & (cat_gaus['separation'] > 10)]
    
    if len(cat_test) > 0:
        print(cat_test['total_flux_gaussian'], cat_test['total_flux_source'], cat_test['separation'])
        return True # Reject this source, it has a gaussian component that indicates it is unacceptably resolved
    else:
        return False # No need to reject this source

def writestr(ra_hms, ra_err, dec_dms, dec_err):
    return f"{ra_hms},{ra_err},{dec_dms},{dec_err}\n"

def find_duplicates(filename):
    """
    Find duplicate lines in a file and return the indices of the lines to keep.

    Args:
        filename (str): The path to the file to be processed.

    Returns:
        list: A list of indices representing the lines to keep.
    """
    lines_to_keep = set()
    index_to_keep = []

    with open(filename, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if line not in lines_to_keep:
                index_to_keep.append(index)
                lines_to_keep.add(line)
    
    return index_to_keep

def remove_duplicates(filename, index_to_keep):
    """
    Removes duplicate lines from a file.

    Args:
        filename (str): The path to the file.
        index_to_keep (list): List of indices to keep.

    Returns:
        None
    """
    lines_to_keep = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if index in index_to_keep:
                lines_to_keep.append(line)
    
    with open(filename, 'w') as f:
        f.writelines(lines_to_keep)
        


if __name__ == "__main__":
    _main()
