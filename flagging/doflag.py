import os,sys
import subprocess

#-----------------------------------------------------------------------------------------------
#   Flagging script
#
#   Inputs are - <input_fits> <output_fits> <badchan_file> <flagmode> <log_file> <badant_file>
#-----------------------------------------------------------------------------------------------

def run(cmd):
    ret = subprocess.run(cmd, shell=True).returncode
    if ret != 0:
        print("ERR:FAILED COMMAND:{cmd}")
        sys.exit(ret)


if(len(sys.argv)<7):
    print("Arguments are - <input_fits> <output_fits> <badchan_file> <flagmode> <log_file> <badant_file>")
    sys.exit()

# Read the ankdir from the environment or use a default for ozstar
ankdir = os.environ.get("ANKDIR", "/fred/oz313/src/ankflag_craft/")

infits		= sys.argv[1]
outfits		= sys.argv[2]
badchanfile	= sys.argv[3]
flagmode	= sys.argv[4]
logfile		= sys.argv[5]
badantfile	= sys.argv[6]	

print("doflag running with -- "+infits+" "+outfits+" "+badchanfile+" "+flagmode+" "+logfile+" "+badantfile+"\n")

print("copying goutfile -- cp "+ankdir+"glogout.dat .")
run("cp "+ankdir+"glogout.dat .")

if(flagmode=='proper'):
    print("doflag - python3 "+ankdir+"runank.py "+infits+" temp_1.fits 1 "+badchanfile+" 1 >> "+logfile)
    run("python3 "+ankdir+"runank.py "+infits+" temp_1.fits 1 "+badchanfile+" 1 >> "+logfile)
    print("doflag - python3 "+ankdir+"runank.py temp_1.fits temp_2.fits 2 none 0 >> "+logfile)
    run("python3 "+ankdir+"runank.py temp_1.fits temp_2.fits 2 none 0 >> "+logfile)
    print("doflag - python3 "+ankdir+"runank.py temp_2.fits temp_3.fits 3 none 0 >> "+logfile)
    run("python3 "+ankdir+"runank.py temp_2.fits temp_3.fits 3 none 0 >> "+logfile)
    print("doflag -python3 "+ankdir+"runank.py temp_3.fits "+outfits+" 4 none 0 >> "+logfile )
    run("python3 "+ankdir+"runank.py temp_3.fits "+outfits+" 4 none 0 >> "+logfile)
    print("doflag - python3 "+ankdir+"print_badant.py "+outfits+" "+badantfile+" >> "+logfile)
    run("python3 "+ankdir+"print_badant.py "+outfits+" "+badantfile+" >> "+logfile)
    run("rm -rf temp_*.fits")
else:
    run("python3 "+ankdir+"runank.py "+infits+" "+outfits+" "+badchanfile+" 1 >> "+logfile)

























