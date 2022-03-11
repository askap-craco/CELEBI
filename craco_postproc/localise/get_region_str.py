#!/usr/bin/env python3
import sys

def parse_jmfit(fname):
    vals = {}
    with open(fname, "r") as f:
        for line in f:
            line = line[:-1] # trim newline
            line = line.split(":")
            vals[line[0].strip()] = ":".join(line[1:]).strip()
    return vals
    
def create_region_str(fname):
    vals = parse_jmfit(fname)
    return f"fk5;ellipse({vals['Actual RA']}, {vals['Actual Dec']}, {float(vals['Est. RA error (mas)'])/1e3}\", {float(vals['Est. Dec error (mas)'])/1e3}\", 0)"
    
print(create_region_str(sys.argv[1]))