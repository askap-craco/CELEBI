# imports
import argparse, os, sys
from os.path import isfile


def get_args():
    """
    Args
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", help = "Directory of FRB output data", type = str)
    parser.add_argument("-l", help = "FRB label", type = str)
    parser.add_argument("--cfreq", help = "Central frequency", type = float, default = None)
    parser.add_argument("--bw", help = "Bandwidth", type = float, default = None)

    return parser.parse_args()






def _compile(args):
    """
    Compile all output data into a single neat txt file + yaml file
    
    """
    justlen = 30
    ofile = open(f"{args.l}_summary.txt", "w")


    # Write filepath for HTR data
    ofile.write("# HTR DATA FILEPATHS #\n")
    htr_file = os.path.join(args.d, f"htr/dynspec_fnames.txt")
    if isfile(htr_file):
        with open(htr_file, "r") as file:
            lines = file.readlines()
        for i, S in enumerate("IQUV"):
            ofile.write(f"ds{S}: ".ljust(justlen) + f"{lines[i]}")

    # General data
    ofile.write("\n# GENERAL DATA #\n")
    
    # FRB name
    ofile.write("FRB name: ".ljust(justlen) + f"{args.l}\n")
    
    # Central freq and bandwidth
    if args.cfreq is None:
        ofile.write("cfreq: ".ljust(justlen) + "NONE".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")
    else:
        ofile.write("cfreq: ".ljust(justlen) + f"{args.cfreq}".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")
    
    if args.bw is None:
        ofile.write("bw: ".ljust(justlen) + "NONE".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")
    else:
        ofile.write("bw: ".ljust(justlen) + f"{args.bw}".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")
    
    if (args.cfreq is not None) and (args.bw is not None):
        ofile.write("DM_ref_freq: ".ljust(justlen) + f"{args.cfreq-args.bw/2}".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")
    else:
        ofile.write("DM_ref_freq: ".ljust(justlen) + "NONE".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")

    # polyco file
    polyco_file = os.path.join(args.d, "binconfigs/craftfrb.polyco")
    if isfile(polyco_file):
        with open(polyco_file, "r") as file:
            lines = file.readlines()
        ofile.write("corr_DM: ".ljust(justlen) + f"{float(lines[0].split()[4])}".ljust(justlen) + "".ljust(justlen) + "[pc/cm^3]\n")
        ofile.write("corr_ref_freq: ".ljust(justlen) + f"{float(lines[1].split()[5])}".ljust(justlen) + "".ljust(justlen) + "[MHz]\n")
        ofile.write("corr_MJD: ".ljust(justlen) + f"{float(lines[0].split()[3])}".ljust(justlen) + "".ljust(justlen) + "[days]\n")
        
    # geometric delay
    geo_file = os.path.join(args.d, "binconfigs/geo_delay.txt")
    if isfile(geo_file):
        with open(geo_file, "r") as file:
            geo_delay = float(file.readline())
        ofile.write("Geocentric delay: ".ljust(justlen) + f"{geo_delay}".ljust(justlen) + "".ljust(justlen) + "[s]\n")

    # beamforming MJD
    antMJD_file = os.path.join(args.d, f"binconfigs/{args.l}_corrected_ant_MJD.txt")
    if isfile(antMJD_file):
        with open(antMJD_file) as file:
            lines = file.readlines()
    
        # get max MJD start time
        bf_MJD = []
        for line in lines:
            bf_MJD += [float(line.split(':')[1].strip())]
        min_ind = bf_MJD.index(max(bf_MJD))
        ofile.write("bform_MJD: ".ljust(justlen) + f"{bf_MJD[min_ind]}".ljust(justlen) + "".ljust(justlen) + "[days]\n")
        ofile.write("bform_ref_ant: ".ljust(justlen) + f"{lines[min_ind].split(':')[0]}\n")

    # crop MJD
    cropMJD_file = os.path.join(args.d, f"binconfigs/frb_crop_MJD.txt")
    if isfile(cropMJD_file):
        with open(cropMJD_file) as file:
            crop_MJD = file.readline()
            ofile.write("crop_MJD: ".ljust(justlen) + crop_MJD.ljust(justlen) + "".ljust(justlen) + "[days]\n")



    # position information
    ofile.write("\n# POSITION #\n")

    # FRB corrected positions
    pos_file = os.path.join(args.d, f"position/{args.l}_final_position.txt")
    if isfile(pos_file):
        with open(pos_file, "r") as file:
            lines = file.readlines()
        for i in range(23,25):
            line = lines[i].split()
            ofile.write((line[0] + ":").ljust(justlen) + line[1].ljust(justlen) + f"+/- {line[3]}\n")

        # uncertainty ellipse
        ofile.write("\nUNCERTAINTY ELLIPSE\n")
        line = lines[28].split(',')
        ofile.write("Major Axis: ".ljust(justlen) + line[0].split()[3].ljust(justlen) + "".ljust(justlen) + "[arcsec]\n")
        ofile.write("Minor Axis: ".ljust(justlen) + line[1].split()[3].ljust(justlen) + "".ljust(justlen) + "[arcsec]\n")
        ofile.write("PA: ".ljust(justlen) + line[2].split()[2].ljust(justlen) + "".ljust(justlen) + "[deg]\n")





    # Now do polcal solutions
    ofile.write("\n# POLCAL SOLUTIONS #\n")

    # get polcal name
    polname_file = os.path.join(args.d, f"polcal/polcal_name.txt")
    if isfile(polname_file):
        with open(polname_file, "r") as file:
            name = file.readline()
        ofile.write("POLCAL name: ".ljust(justlen) + name)
    
    # get polcal solutions
    polcal_file = os.path.join(args.d, f"polcal/{args.l}_polcal_solutions.txt")
    if isfile(polcal_file):
        with open(polcal_file, "r") as file:
            lines = file.readlines()
        pol_vars = ["psi", "RM", "tau", "phi", "L_scale_factor", "V_scale_factor",
                    "alpha"]
        pol_units = ["[rad]", "[rad/m^2]", "[us]", "[rad]", "", "", "[rad]"]
        for i, line in enumerate(lines[:-1]):
            pol_val = float(line.split(':')[1])
            pol_err = float(line.split(':')[2])
            ofile.write(f"{pol_vars[i]}: ".ljust(justlen) + f"{pol_val}".ljust(justlen) + f"+/- {pol_err}".ljust(justlen) + f"{pol_units[i]}\n")
        ellipt_bool = int(lines[-1].split(':')[1])
        if ellipt_bool:
            bool_str = "TRUE"
        else:
            bool_str = "FALSE"
        ofile.write("ELLIPTICITY: ".ljust(justlen) + bool_str + "\n")


    









if __name__ == "__main__":
    # main block of code

    args = get_args()

    _compile(args)

