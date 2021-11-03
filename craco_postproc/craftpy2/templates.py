difx_v2d_file_pre_template = """
#  Template v2d file for DiFX correlation of craftfrb

vex = {vex_file}

startSeries = {startseries:d}
minLength = 1
allowAllClockOffsets = True
tweakIntTime = True
exhaustiveAutocorrs = True

visBufferLength = 8
antennas = {antennas}
"""

difx_v2d_file_ant_template = """
ANTENNA {twolettername}
{
  name = {antname}
  clockOffset={clock_offset:.6f}
  freqClockOffs={freqClockOffs}
  clockRate=0
  clockEpoch=57000.0
  phaseCalInt=0
  toneSelection=none
  format=CODIFC/27/{framesize}/{bits}
  sampling=COMPLEX_DSB



"""

difx_v2d_file_stream_template = """
DATASTREAM {twolettername}-P{ipol:d}
{
  file = {codif_file}
  format=CODIFC/27/{framesize}/{bits}
  sampling=COMPLEX_DSB
}
"""

difx_v2d_file_post_template = """
# The nChan should never be less than 128.
# For numbers of channels < 128, set specAvg so nChan/specAvg
# gives the desired number of channels
SETUP default
{
  tInt =  {tint:.9f}
  subintNS = {subint_nsec}
  nFFTChan = {nfft_chan}
  nChan = {nchan}
  doPolar = True # Full stokes
  binConfig = {binconfig_file} # Pulsar Setup
}

# This, along with SETUP default above, should always be done\n")
RULE default
{
  setup = default
}

# SETUP place holders (commented)
# SETUP askap.set
# {
# }

# Sources (pointing centers) with recorded data but no offset pointing centers:
SOURCE {srcname} { }

"""


sched_key_file_template = """
version  = 1
expt     = 'craft'
expcode  = 'craftfrb'
obstype  = 'VLBI'

piname   = 'A.T. Deller'
address1 = 'Swinburne'
address2 = ''
address3 = 'Australia'
email    = 'adeller@astro.swin.edu.au'
phone    = '+61-3-9214-5307 (w)'
obsphone = '+61-3-9214-5307 (w)'
obsmode  = 'ASKAP 300 channel'
note1    = 'ASKAP'
! ================================================================
!       Correlator section
! ================================================================
correl   = 'Socorro'
coravg   = 2
corchan  = 16
cornant  = 3
corpol   = 'on'
corwtfn  = 'uniform'
corsrcs  = 'from .sum file only'
cortape  = 'ftp'
cornote1 = 'This is special ASKAP correlation'

! ================================================================
!       Catalogs (special askap versions)
! ================================================================
stafile  = '{station_file}'
freqfile = '{freq_file}'
overwrite
srccat /
EQUINOX = J2000
SOURCE='{srcname}' RA={srcra} DEC={srcdec} REMARKS='Beam centre for dumped voltage data' /
endcat /

setinit = askap.set /
 dbe      = '{dbe}'
 format   = 'vdif'          !  Sched doesn't understand CODIF, so lie and say VDIF.
 nchan    = {nchan}               !  Put in {nchan}x8 MHz as placeholder, overwrite later
 bbfilt   = {bbfilt}
 netside  = U
 bits     = 2
 firstlo  = 2100.0
 freqref  = 2100.0
 freqoff  = {freqoff}
 pol      = {pol}
 pcal     = 'off'
   /
endset /

year     = {startyear:d}
month    = {startmonth:d}
day      = {startday:d}
start    = {start}

stations = {stations}
setup  = askap.set
minpause = 5

source = '{srcname}'  dur = {duration:d}  gap = 0   /

"""

sched_freq_file_template = """
Name = v20cm_1  Station = ASKAP{ant_code}    Priority = 1
  rf1 = 500, 500  ifname = A,    C
  rf2 = 2000, 2000  fe  = '20cm', '20cm'
  pol =  RCP,  LCP  lo1 =  2100,  2100  syn(2) = 2.1
/

"""

sched_station_file_template = """
  STATION=ASKAP{ant_code}   STCODE={twolettername}  CONTROL=VLBA
    DBNAME = {ant_code}-ASKAP
        MOUNT=ALTAZ  AX1LIM=-90,450 AX2LIM=2.25,90 AX1RATE=83.6 AX2RATE=29.0
        AX1ACC=0.75  AX2ACC=0.25
        TSETTLE=6 TLEVSET=5 MINSETUP=5  DAR=RDBE2  NBBC=16
        DISK=MARK5C   MEDIADEF=DISK    TSCAL=CONT
        HOR_AZ =   0,  5, 10, 15, 25, 30, 40, 45, 70, 75,120,125,130,135,
                 155,160,185,190,195,220,225,235,240,245,250,255,265,270,
                 275,300,305,310,315,330,335,340,345,350,360
        HOR_EL =   2,  2,  3,  2,  2,  3,  3,  4,  4,  5,  5,  4,  4,  3,
                   3,  2,  2,  3,  4,  4,  3,  3,  4,  4,  5,  6,  6,  5,
                   6,  6,  5,  6,  5,  5,  4,  4,  3,  2,  2
        AXISOFF=  0.0
        X= {itrfpos_x}  Y= {itrfpos_y}  Z=  {itrfpos_z}
        DXDT= 0.0  DYDT=  0.0  DZDT= 0.0  EPOCH=54466
        FRAME='FROM FCM'
      /

"""
