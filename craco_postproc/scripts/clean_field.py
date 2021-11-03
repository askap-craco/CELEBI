reqvars = [
    "vis",
    "imagename",
    "imsize",
    "pixsize",
    "phasecenter",
    "posfile",
    "fitsimage",
]
if not all([var in globals() for var in reqvars]):
    print("===== How to use this file =====")
    print("In an INTERACTIVE CASA session, define ALL of the following vars:")
    for var in reqvars:
        print(f"\t{var}")
    print("Then, do:")
    print('\texecfile("clean_field.py", globals())')
    exit()

tclean(
    vis=vis,
    imagename=imagename,
    imsize=[imsize],
    cell=[pixsize],
    phasecenter=phasecenter,
    pblimit=-1,  # everything below here is to be left for now
    selectdata=True,
    field="",
    spw="",
    timerange="",
    uvrange="",
    antenna="",
    scan="",
    observation="",
    intent="",
    datacolumn="corrected",
    stokes="I",
    projection="SIN",
    startmodel="",
    specmode="mfs",
    reffreq="",
    nchan=-1,
    start="",
    width="",
    outframe="LSRK",
    veltype="radio",
    restfreq=[],
    interpolation="linear",
    perchanweightdensity=True,
    gridder="widefield",
    facets=1,
    psfphasecenter="",
    chanchunks=1,
    wprojplanes=-1,
    vptable="",
    mosweight=True,
    aterm=True,
    psterm=False,
    wbawp=True,
    conjbeams=False,
    cfcache="",
    usepointing=False,
    computepastep=360.0,
    rotatepastep=360.0,
    pointingoffsetsigdev=[],
    normtype="flatnoise",
    deconvolver="multiscale",
    scales=[],
    nterms=2,
    smallscalebias=0.0,
    restoration=True,
    restoringbeam=[],
    pbcor=False,
    outlierfile="",
    weighting="natural",
    robust=0.5,
    noise="1.0Jy",
    npixels=0,
    uvtaper=[],
    niter=200,
    gain=0.1,
    threshold=0.0,
    nsigma=0.0,
    cycleniter=20,
    cyclefactor=1.0,
    minpsffraction=0.05,
    maxpsffraction=0.8,
    interactive=True,
    usemask="user",
    mask="",
    pbmask=0.0,
    sidelobethreshold=3.0,
    noisethreshold=5.0,
    lownoisethreshold=1.5,
    negativethreshold=0.0,
    smoothfactor=1.0,
    minbeamfrac=0.3,
    cutthreshold=0.01,
    growiterations=75,
    dogrowprune=True,
    minpercentchange=-1.0,
    verbose=False,
    fastnoise=True,
    restart=True,
    savemodel="modelcolumn",
    calcres=True,
    calcpsf=True,
    parallel=False,
)

sources = []

print("Cleaning finished, opening cleaned image")
viewer(f"{imagename}.image")
print(
    "Identify sources in the field, and type their pixel coordinates as\n"
    + "comma-separated integers. When all are done, input -1"
)
while True:
    source_xy = eval(input("> "))
    source_xy = tuple(map(int, source_xy.split(",")))

    if len(source_xy) == 1:
        if source_xy[0] < 0:
            break
    elif len(source_xy) == 2:
        sources.append(source_xy)

print("")
print(f"Writing source position file: {posfile}")
with open(posfile, "w") as f:
    basestr = "{0},{1}\n"
    for source in sources:
        f.write(basestr.format(source[0], source[1]))

print(f"Converting {imagename}.image to {fitsimage}")
exportfits(imagename=imagename + ".image", fitsimage=fitsimage)

print("Done")
