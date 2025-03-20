import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--farfieldfcm', help='Far-field FCM file (must already exist)', required=True)
parser.add_argument('-n', '--nearfieldfcm', help='Near-field FCM file (must not yet exist)', required=True)
parser.add_argument('-c', '--corrections', help='Text file with near field delay corrections', required=True)
args = parser.parse_args()

if not os.path.exists(args.farfieldfcm):
    parser.error("Far field fcm file {0} doesn't exist".format(args.farfieldfcm))

if os.path.exists(args.nearfieldfcm):
    parser.error("Near field fcm file {0} already exists".format(args.nearfieldfcm))

if not os.path.exists(args.corrections):
    parser.error("Near field corrections file {0} doesn't exist".format(args.corrections))

corrlines = open(args.corrections).readlines()
farfieldfcmlines = open(args.farfieldfcm).readlines()

with open(args.nearfieldfcm, "w") as output:
    for line in farfieldfcmlines:
        if "delay" in line and "ns" in line:
            splitline = line.split()
            antenna = splitline[0].split('.')[2]
            currentdelay = float(splitline[-1][:-2])
            found = False
            for corrline in corrlines:
                if corrline.split(':')[0] == antenna:
                    found = True
                    delay = currentdelay + float(corrline.split(':')[-1].strip())
            if not found:
                print("Couldn't find antenna", antenna, "in the corrections file")
                sys.exit()
            output.write("common.antenna.{0}.delay = {1:.5f}ns\n".format(antenna, delay))
        else:
            output.write(line)
