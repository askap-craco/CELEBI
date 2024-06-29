# import sys
# import os
from glob import glob

workdir = raw_input("Enter work directory containing .fits visibility file to be flagged: ")
fitsfile = glob(workdir+"/*fits")[0]
vis = fitsfile[:-4] + "ms"#

print("fitsfile = "+fitsfile)
print("vis      = "+vis)

basename = fitsfile[len(workdir)+1:-5]
suggested_flagfile = "flagfiles/"+basename+"_flags.txt"

print("Suggested flagfile: "+suggested_flagfile)
flag_input = raw_input("Press ENTER to accept, or enter alternative flagfile: ")

flagfile = flag_input if flag_input != "" else suggested_flagfile

print("flagfile = "+flagfile)

importuvfits(fitsfile=fitsfile, vis=vis)

flagged_channels = []

print("Flagging channels")
print("For each iteration, input one of:")
print("    > A single channel number to flag")
print("    > A channel range to be flagged (inclusive)")
print("    > -1 to stop")
while True:
    plotms(vis=vis, antenna="*&", xaxis="channel")
    flags = raw_input("> ")
    flags = tuple(map(int, flags.split("~")))

    if len(flags) == 1:
        if flags[0] < 0:
            break
        else:
            flagdata(vis=vis, mode="manual", spw="0:{}".format(*flags))
            flagged_channels.append(flags)
    elif len(flags) == 2:
        flagdata(vis=vis, mode="manual", spw="0:{0}~{1}".format(*flags))
        flagged_channels.append(flags)

print("")
print("Flagging antennas")
print("Plotting amp vs chan, iterating over antennas.")
print("Check them, and input in a comma-separated list")
print("any that are bad.")
plotms(vis=vis, antenna="*&", xaxis="channel", iteraxis="antenna")
bad_ants = raw_input("> ")
if bad_ants == "":
    bad_ants = ()
else:
    bad_ants = tuple(map(int, bad_ants.split(",")))

print("")
print("Writing flag file: {}".format(flagfile))
with open(flagfile, "w") as f:
    f.write("dtimrang = 1  timeoff = 0\n\n")

    basestr = (
        "antennas=0 bchan={0} echan={1} timerang=0,0,0,0,0,23,59,59 reason='RFI' /\n"
    )
    for chans in flagged_channels:
        if len(chans) == 1:
            f.write(basestr.format(chans[0] + 1, chans[0] + 2))
        else:
            f.write(basestr.format(chans[0] + 1, chans[1] + 2))

    if len(bad_ants) > 0:
        basestr = (
            "antennas={} bchan=0 echan=0 timerang=0,0,0,0,0,23,59,59 reason='RFI' /\n"
        )
        for ant in bad_ants:
            f.write(basestr.format(ant))

print("Done")
