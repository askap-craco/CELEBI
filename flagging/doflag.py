import os
import sys
import subprocess


def run(cmd):
    """
    A drop in replacement for os.system, that will check the return code of the
    executed command and report failure/quit if the command fails

    cmd : str
       The system command to run.
    """
    print(f'{__file__}$> {cmd}')
    # Run the command and capture all the relevant details
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # forward stdout/stderr to the relelevant streams
    if result.stdout:
        print(result.stdout, end='')
    if result.stderr:
        print(result.stderr, file=sys.stderr, end='')

    # quit if the return code is not zero, with a note
    if result.returncode != 0:
        print(f"{__file__}$> ERR:FAILED:exitcode:{result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


#-----------------------------------------------------------------------------------------------
#   Flagging script
#
#   Inputs are - <input_fits> <output_fits> <badchan_file> <flagmode> <log_file> <badant_file>
#-----------------------------------------------------------------------------------------------

if(len(sys.argv)<7):
    print("Arguments are - <input_fits> <output_fits> <badchan_file> <flagmode> <log_file> <badant_file>")
    sys.exit()

# Read the ankdir from the environment or use a default for ozstar
# ankdir = os.environ.get("ANKDIR", "/fred/oz313/src/ankflag_craft/")
ankdir = os.path.split(os.path.abspath(__file__))[0]+'/'

infits		= sys.argv[1]
outfits		= sys.argv[2]
badchanfile	= sys.argv[3]
flagmode	= sys.argv[4]
logfile		= sys.argv[5]
badantfile	= sys.argv[6]	

print("doflag running with -- "+infits+" "+outfits+" "+badchanfile+" "+flagmode+" "+logfile+" "+badantfile+"\n")

print("copying goutfile")
run("cp "+ankdir+"glogout.dat .")

if(flagmode=='proper'):
    run("python3 "+ankdir+"runank.py "+infits+" temp_1.fits 1 "+badchanfile+" 1  | tee -a "+logfile)
    run("python3 "+ankdir+"runank.py temp_1.fits temp_2.fits 2 none 0  | tee -a "+logfile)
    run("python3 "+ankdir+"runank.py temp_2.fits temp_3.fits 3 none 0  | tee -a "+logfile)
    run("python3 "+ankdir+"runank.py temp_3.fits "+outfits+" 4 none 0  | tee -a "+logfile)
    run("python3 "+ankdir+"print_badant.py "+outfits+" "+badantfile+"  | tee -a "+logfile)
    run("rm -rf temp_*.fits")
else:
    run("python3 "+ankdir+"runank.py "+infits+" "+outfits+" "+badchanfile+" 1  | tee -a "+logfile)

























