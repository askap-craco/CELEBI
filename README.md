# CELEBI: The CRAFT Effortless Localisation and Enhanced Burst Inspection Pipeline

[CELEBI](https://arxiv.org/abs/2301.13484) is an automated data processing pipeline for producing sub-arcsecond precision localisations and high-time resolution polarimetric measurements of fast radio bursts (FRBs) from voltages obtained with the Australian Square Kilometre Array Pathfinder (ASKAP).

CELEBI operates on the VCRAFT data format, and is designed to be run on a supercomputer. Once the dependencies have been installed, you should set up a config file based on the [template](configs/template.config) for the data you are processing, and then run [main.nf](pipelines/main.nf):
```
nextflow /path/to/CELEBI/pipelines/main.nf -c [config file]
```
It is recommended that you do this from a separate processing directory, and make use of the `-with-report` and `-w` Nextflow configuration options.

## Options
Running CELEBI without any additional flags will run everything up until completion of correlation. To perform flux calibration, imaging, and localisation, add `--calibrate`. To perform polarisation calibration and beamforming after this, add `--beamform`. Both of these flags can be provided together to run the entire pipeline in one go.

You can omit the FRB and polcal workflows from running with `--nofrb` and `--nopolcal` respectively. If the polcal is omited, the beamforming will substitute zeros in its polarisation calibration solutions.

RFI subtraction can be skipped with `--skiprfi`.

Visibility flagging can be skipped with `--noflag`. You can provide custom AIPS flag files with `--fieldflagfile`, `--polflagfile`, and `--fluxflagfile`. These can be provided alongside using automatic flagging.

## Dependencies
- [Nextflow](https://nextflow.io/)
- [AIPS](https://doi.org/10.1007/0-306-48080-8_7)
- [ParselTongue](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue)
- [CASA](https://dx.doi.org/10.1088/1538-3873/ac9642)
- [CASA Analysis Utilities](https://casaguides.nrao.edu/index.php/Analysis_Utilities)
- [DiFX](http://dx.doi.org/10.1086/658907)
- [psrvlbireduce](https://github.com/dingswin/psrvlbireduce)
- Numpy
- Scipy
- Matplotlib
- Astropy
- [Astroquery](https://dx.doi.org/10.3847/1538-3881/aafc33)
