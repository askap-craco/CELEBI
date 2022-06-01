# import casa analysis scripts for converting component lists into a position and corresponding pixel
sys.path.append("/home/ubuntu/analysis_scripts/")  # adjust if this changes!
import analysisUtils as au

allinput = raw_input().split(",")
target = allinput[0]
nsources = int(allinput[1])
cutoff = float(allinput[2])
posfile = allinput[3]

def get_pixels_from_field(target, nsources, cutoff, posfile):
    """Write out pixels corresponding to continuum sources in a field
    image for calibration

    :param target: field image (.image or .fits format) to search in
    :type target: str
    :param nsources: maximun number of sources to get back
    :type nsources: int
    :param cutoff: limit in flux relative to strongest in field
        e.g. 0.05 (5%)
    :type imsize: float
    :param posfile: File to write pixels of identified sources to
    :type posfile: str
    """

    # open file, use casa task "findsources" to identify bright continuum sources
    ia.open(target)
    clrec = ia.findsources(nmax=nsources, point=False, cutoff=cutoff)
    pos_pix = "#ra,dec\n"
    # loop through each source, extract relevant position, convert it
    for i in range(len(clrec.keys()) - 1):
        cl1 = clrec["component" + str(i)]
        m0 = cl1["shape"]["direction"]["m0"]["value"]  # in rad
        m1 = cl1["shape"]["direction"]["m1"]["value"]  # in rad
        pos = str(au.rad2radec(m0, m1, hmsdms=True))
        radec = au.findRADec(target, pos, round=True)
        if radec != None:
            ra, dec = radec
            new_pos_pix = str(ra) + "," + str(dec) + "\n"
            pos_pix += new_pos_pix

    # write to file
    with open(posfile, "w") as f:
        f.write(pos_pix)

get_pixels_from_field(target, nsources, cutoff, posfile)