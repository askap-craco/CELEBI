import ast
import configparser
import re
import sys
from pathlib import Path

import pandas as pd
from astropy.time import Time
from craftpy import utils
from craftpy2 import templates
from dotty_dict import dotty
from natsort import natsorted
from scipy import constants

FFT_SEC_PER_CHANNEL = 27 / 32 / constants.mega  # noqa:WPS432
NCHAN_MIN = 128


def load_fcm(fcmpath):
    """
    Read the java style properties antenna fcm file
    """

    with open(fcmpath) as fp:
        properties_str = "[fcm]\n" + fp.read()
    config = configparser.ConfigParser()
    config.read_string(properties_str)
    dot = dotty()
    for key, value in config.items("fcm", raw=True):
        dot[key] = value
    return dot


def get_delays(fpga_delay_file, fpga):
    fpga_delays = {}
    with open(fpga_delay_file) as fp:
        for line in fp:
            line = line.split("#")[0]  # Remove comments
            if line == "" or line.isspace():
                continue
            keys = line.split()
            if len(keys) != 3:
                sys.stderr.write("Cannot parse : " + line)
                continue
            ant = keys[0]
            thisFPGA = keys[1]
            delay = keys[2]
            if thisFPGA != fpga:
                continue
            fpga_delays[ant] = delay
    return fpga_delays


class DIFXSetup:
    def __init__(
        self,
        srcname,
        srcra,
        srcdec,
        startmjd,
        stopmjd,
        fcm_file,
        binconfig_file,
        eoptxt_file,
        targetants,
        integration,
        nchan=128,
        framesize=8064,
        bits=1,
        outdir="",
        fpga=None,
        fpga_delay_file=None,
        outputbasename="craftfrb",
    ):
        self.fcm_file = fcm_file
        self.binconfig_file = binconfig_file
        self.eoptxt_file = eoptxt_file
        self.fcm = load_fcm(fcm_file)
        self.targetants = targetants
        self.integration = integration
        self.nchan = max(nchan, NCHAN_MIN)
        self.framesize = framesize
        self.bits = bits
        self._outdir = outdir
        self.outputbasename = outputbasename

        self.fpga_delay = {}
        if fpga and fpga_delay_file:
            if not Path(fpga_delay_file).is_file():
                self.logger.info(
                    f"FPGA Delay file '{fpga_delay_file}' does not exist."
                )
            self.fpga_delay = get_delays(fpga_delay_file, fpga)

        self._set_obs_antennas()
        self._set_integration()

    @property
    def freq_file(self):
        outpath = Path(self._outdir) / "askapfreq.dat"
        return outpath.as_posix()

    @property
    def station_file(self):
        outpath = Path(self._outdir) / "askapstation.dat"
        return outpath.as_posix()

    @property
    def key_file(self):
        outpath = Path(self._outdir) / f"{self.outputbasename}.key"
        return outpath.as_posix()

    @property
    def v2d_file(self):
        outpath = Path(self._outdir) / f"{self.outputbasename}.v2d"
        return outpath.as_posix()

    @property
    def vex_file(self):
        outpath = Path(self._outdir) / f"{self.outputbasename}.vex"
        return outpath.as_posix()

    @property
    def calc_file(self):
        outpath = Path(self._outdir) / f"{self.outputbasename}.calc"
        return outpath.as_posix()

    @property
    def sched_logfile(self):
        outpath = Path(self._outdir) / "sched_log.txt"
        return outpath.as_posix()

    @property
    def vex2difx_logfile(self):
        outpath = Path(self._outdir) / "vex2difxlog.txt"
        return outpath.as_posix()

    @property
    def difxcalc_logfile(self):
        outpath = Path(self._outdir) / "difxcalclog.txt"
        return outpath.as_posix()

    @property
    def obs_antennas(self):
        return self._obs_antennas

    @property
    def nobs_ants(self):
        return len(self._obs_antennas)

    @property
    def twoletternames(self):
        return [ant.get("twolettername") for ant in self.obs_antennas]

    @property
    def nfft_chan(self):
        return self.nchan

    @property
    def tint(self):
        return self._tint

    @property
    def subint_nsec(self):
        return self._subint_nsec

    def run(self):
        self.write_freq_file()
        self.write_station_file()
        self.write_key_file()
        self.write_v2d_file()

        self._run_sched()
        self._run_vex2difx()
        self._run_difxcalc()

    def write_freq_file(self):
        freqout = [
            templates.sched_freq_block_template.format(
                ant_code=ant.get("antcode")
            )
            for ant in self.obs_antennas
        ]
        freqoutstr = "\n".join(freqout)
        with open(self.freq_file, "w+") as fp:
            fp.write(freqoutstr)

    def write_station_file(self):
        statout = [
            templates.sched_freq_file_template.format(
                ant_code=ant.get("antcode"),
                twolettername=ant.get("twolettername"),
                itrfpos_x=ast.literal_eval(ant.get("location").get("itrf"))[0],
                itrfpos_y=ast.literal_eval(ant.get("location").get("itrf"))[1],
                itrfpos_z=ast.literal_eval(ant.get("location").get("itrf"))[2],
            )
            for ant in self.obs_antennas
        ]
        statoutstr = "\n".join(statout)
        with open(self.station_file, "w+") as fp:
            fp.write(statoutstr)

    def write_key_file(self):
        start_time = Time(self.startmjd, format="mjd", scale="utc")
        stop_time = Time(self.stopmjd, format="mjd", scale="utc")
        if self.npol:
            dbe = "rdbe_ddc"
            bbfilt = 8
            freqoff = "-608.0, -600.0, -592.0, -584.0, -576.0, -568.0, -560.0, -552.0"
            pol = ", ".join(list("L" * 8))

        else:
            dbe = "rdbe_pfb"
            bbfilt = 32
            freqoff = (
                "-1008.0, -976.0, -944.0, -912.0, -880.0, -848.0, -816.0, -784.0, "
                * 2
            )
            pol = ", ".join(list("L" * 8) + list("R" * 8))

        keyout = templates.sched_key_file_template.format(
            station_file=self.station_file,
            freq_file=self.freq_file,
            srcname=self.srcname,
            srcra=self.srcra,
            srcdec=self.srcdec,
            dbe=dbe,
            nchan=self.nchan,
            bbfilt=bbfilt,
            freqoff=freqoff,
            pol=pol,
            startyear=start_time.ymdhms.year,
            startmonth=start_time.ymdhms.month,
            startday=start_time.ymdhms.day,
            start=start_time.datetime.time().isoformat(timespec="seconds"),
            stations=", ".join(self.twoletternames),
            duration=round((stop_time - start_time).sec),
        )
        with open(self.key_file, "w+") as fp:
            fp.write(keyout)

    def write_v2d_file(self):
        v2dout = []
        v2dout_pre = templates.difx_v2d_file_pre_template.format(
            vex_file=self.vex_file,
            startseries=0,
            antennas=", ".join(self.twoletternames),
        )
        v2dout.append(v2dout_pre)

        for ant in self.obs_antennas:
            clock_offset = (
                -float(re.sub("ns", "", ant.get("delay"))) * constants.milli
            )
            if ant in self.fpga_delay:
                clock_delays = (
                    [str(int(self.fpga_delay.get(ant)) * 6.75)]
                    * 4
                    * 2
                    * self.npol
                )
                freq_clockoffs = ", ".join(clock_delays)
            twolettername = ant.get("twolettername")

            v2dout_ant = templates.difx_v2d_file_ant_template.format(
                twolettername=twolettername,
                antname=ant.get("name"),
                clock_offset=clock_offset,
                freqClockOffs=freq_clockoffs,
                framesize=self.framesize,
                bits=self.nbits,
            )
            v2dout.append(v2dout_ant)

            if self.npol > 1:
                v2dout.append(
                    f"  datastreams={twolettername}-P0,{twolettername}-P1"
                )
                v2dout.append("}\n")
                for ipol in range(self.npol):
                    v2dout_stream = (
                        templates.difx_v2d_file_stream_template.format(
                            twolettername=twolettername,
                            ipol=ipol,
                            codif_file=self.codif_file,
                            framesize=self.framesize,
                            bits=self.nbit,
                        )
                    )
                    v2dout.append(v2dout_stream)
            else:
                v2dout.append(f"  file = {self.codif_file}")
                v2dout.append("}")

        v2dout_post = templates.difx_v2d_file_post_template.format(
            tint=self.tint,
            subint_nsec=self.subint_nsec,
            nfft_chan=self.nfft_chan,
            nchan=self.nchan,
            binconfig_file=self.binconfig_file,
            srcname=self.srcname,
        )
        v2dout.append(v2dout_post)

        with open(self.eoptxt_file) as eopf:
            eoplines = eopf.read().rstrip()

        for eopline in [_f for _f in re.split(r"\n+", eoplines) if _f]:
            if eopline.startswith("#"):
                continue
            v2dout.append(eopline)

        v2doutstr = "\n".join(v2dout)
        with open(self.v2d_file, "w+") as fp:
            fp.write(v2doutstr)

    def _set_obs_antennas(self):
        obs_antennas = {}
        for iant, ant in enumerate(
            natsorted(list(self.fcm["common.antenna"].keys()))
        ):
            if ant == "ant":
                continue
            ant_dict = self.fcm["common.antenna"][ant]
            ant_name = ant_dict.get("name")
            if ant_name not in self.targetants:
                self.logger.info(
                    f"Skipping antenna '{ant_name}' from fcm file. Was not requested."
                )
                continue
            stcode = str(iant) if iant < 10 else chr(ord("A") + iant - 10)
            obs_antennas[ant] = {
                "name": ant_name,
                "delay": ant_dict.get("delay", "0.0ns"),
                "location": ant_dict.get("location"),
                "antcode": re.sub("ak", "", ant_name),
                "twolettername": f"A{stcode}",
            }

        self._obs_antennas = obs_antennas

    def _set_integration(self):
        tfft = FFT_SEC_PER_CHANNEL * self.nfft_chan
        tint = self.integration
        if self.integration > 0.5:
            tsubint = 0.013824
        elif self.integration < 0.02:
            tint = round(tint / tfft) * tfft
            tsubint = tint
        else:
            tsubint = tint / round(tint / 0.014)

        tsubint = round(tsubint / tfft) * tfft
        self._subint_nsec = round(tsubint * constants.giga)
        self._tint = round(tint / tsubint) * tsubint

    def _run_sched(self):
        # Run Sched
        self.logger.info("Running sched")
        sched_result = utils.run_command(
            f"sched < {self.key_file}", logfile=self.sched_logfile
        )
        if sched_result.exit_code != 0:
            self.logger.info(
                f"Sched failed: Check '{self.sched_logfile}' for log output. Exiting"
            )
            sys.exit()

    def _run_vex2difx(self):
        # Run vex2difx (to generate .calc, .input, .flag)
        self.logger.info("Running vex2difx")
        sched_result = utils.run_command(
            f"vex2difx {self.v2d_file}", logfile=self.vex2difx_logfile
        )
        if sched_result.exit_code != 0:
            self.logger.info(
                f"vex2difx failed: Check '{self.vex2difx_logfile}' for log output. Exiting"
            )
            sys.exit()

    def _run_difxcalc(self):
        # Run difxcalc (to generate the interferometer model .im using built-in CALC 11 model)
        self.logger.info("Running difxcalc")
        if not Path(self.calc_file).is_file():
            self.logger.info(f"Calc file '{self.calc_file}'' does not exist.")
        sched_result = utils.run_command(
            f"difxcalc {self.calc_file}", logfile=self.difxcalc_logfile
        )
        if sched_result.exit_code != 0:
            self.logger.info(
                f"difxcalc failed: Check '{self.difxcalc_logfile}' for log output. Exiting"
            )
            sys.exit()
