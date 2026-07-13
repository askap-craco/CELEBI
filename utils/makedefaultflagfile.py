#	-----------------------------------------------------------------------------
#	A script to make a default flag file 
#	
#		Identified common antennas in calibrators and target
#		Flags missing antennas 
#		Flags list of antennas from the config file
#
#								AB, November 25
#	-----------------------------------------------------------------------------

import sys
import os
import argparse
import glob

#	-----------------------------------------------------------------------------
#	Function to read command line arguments

def get_args() -> argparse.Namespace:
    #Parse command line arguments
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fcaldata",
        help="Flux cal data directory",
        required=True,
    )
    parser.add_argument(
        "--pcaldata",
        default=None,
        help="Pol cal data directory",
    )
    parser.add_argument(
        "--frbdata",
        help="FRB data directory",
        required=True,
    )
    parser.add_argument(
        "--exclants",
        help="Comma separated list of bad antennas",
        required=True,
    )    
    parser.add_argument(
        "--outname",
        default="CRAFT",
        help="Base name for the output files",
    )
    parser.add_argument(
        "--maxants",
        type=int,
        default=36,
        help="Total number of antennas",
    )    

    args = parser.parse_args()
    verify_args(args, parser)
    return args

#	---------------------------------------------------------------------------------
#	Function to verify that the provided command line arguments are valid
def verify_args(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> None:

    # Check that provided paths exist
    if not os.path.exists(args.fcaldata):
        parser.error(
            f"Target directory {args.fcaldata} doesn't exist"
        )
    else:
    	print("Flux cal data path",args.fcaldata)
        
    if not os.path.exists(args.frbdata):
        parser.error(
            f"Target directory {args.frbdata} doesn't exist"
        )   
    else:
    	print("FRB data path",args.frbdata)
        
#	------------------------------------------------------------------
#	The main part of the script
    
args 		= get_args()
    
fcalants 	= sorted(glob.glob(f"{args.fcaldata}/ak*"))
frbants 	= sorted(glob.glob(f"{args.frbdata}/ak*"))

for i in range(0,len(fcalants)):
	fcalants[i] = fcalants[i][-2:]
	
for i in range(0,len(frbants)):
	frbants[i] = frbants[i][-2:]

print("Flux cal antenna count = ",len(fcalants))
print("FRB antenna count      = ",len(frbants))

commants	= sorted(list(set(fcalants) & set(frbants)))

pcalants	= []
if (args.pcaldata is not None):
	pcalants 	= sorted(glob.glob(f"{args.pcaldata}/ak*"))
	for i in range(0,len(pcalants)):
		pcalants[i] = pcalants[i][-2:]
	print("Pol cal antenna count  = ",len(pcalants))
	commants	= sorted(list(set(pcalants) & set(commants)))

#toex		= args.exclants.split(",")
#for exl in toex:
#	if (exl in commants):
#		commants.remove(exl)

# Split and strip whitespace just in case there are spaces between commas
toex = [x.strip() for x in args.exclants.split(",") if x.strip()]

for exl in toex:
    if exl in commants:
        commants.remove(exl)
    if exl in fcalants:
        fcalants.remove(exl)
    if exl in frbants:
        frbants.remove(exl)
    if exl in pcalants:
        pcalants.remove(exl)

exants		= []
for i in range(1,args.maxants+1):
	if (str(i).zfill(2) not in commants):
		exants.append(str(i).zfill(2))

with open(args.outname+"_antcount_fcal_pcal_frb_common.txt", "w") as fl:
    fl.write(str(len(fcalants)) + "\n")
    fl.write(str(len(pcalants)) + "\n")
    fl.write(str(len(frbants)) + "\n")
    fl.write(str(len(commants)) + "\n")

with open(args.outname+"_good_ants.txt", "w") as fl:
    for ant in commants:
        fl.write(str(ant) + "\n")

with open(args.outname+"_exants.txt", "w") as fl:
    for ant in exants:
        fl.write("antennas="+str(ant)+" bchan=0 echan=0 timerang=0,0,0,0,0,23,59,59 reason='EXCLUDED' / \n")
























