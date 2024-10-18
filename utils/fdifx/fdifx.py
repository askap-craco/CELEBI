#python libraries
import numpy as np
import math
# import parseDiFX
import glob, struct, os, sys
import pandas as pd
from util import *
import copy

#c libraries
from ctypes import *
cDiFX = CDLL(os.path.dirname(__file__)+'/cDiFX.so') #later


#TODO - IMPROVEMENTS TO POSSIBLY BE MADE LATER - although stop worrying about this for now
#1. READ FILES IN C, THEN PASS MUTABLE BUFFERS THROUGH PYTHON API
#2. pass in complex data array instead


#initialize c functions in cDiFX
freqtype = np.ctypeslib.ndpointer(np.int32, ndim=1, flags='C')

#load data function
cDiFX.byte2vis.argtypes = (POINTER(c_char),POINTER(c_char),POINTER(c_char),freqtype,freqtype,freqtype,freqtype)
cDiFX.byte2vis.restype = c_int

#additonal function in c
cDiFX.get_baselines.argtypes = (c_int,freqtype)

#save data function
cDiFX.vis2byte.argtypes = (POINTER(c_char),c_int,POINTER(c_char),c_int,POINTER(c_char),freqtype)









#############################################
######    ADDITIONAL CLASSES     ############ 
#############################################

#baseline data
class BaselineInfo:
    """ Baseline class
    """
    num = 0
    offset = 0
    numa = 0

    def print(self):
        print("Baseline Info:")
        print("--------------")
        print("# baselines: {:d}".format(self.num))
        print("Baseline offset: {:d}".format(self.offset))


#frequency data
class FreqInfo:
    """ frequency bands class
    """
    num = 0
    numa = 0
    numActive = 0
    freq = []
    allfreq = []
    bandsActive = []
    nchan = []

    def print(self):
        print("Frequency Info:")
        print("---------------")
        print("# Bands: {:d}".format(self.num))
        print("# Active Bands: {:d}".format(self.numActive))
        print("Active Bands")
        print("------------")
        for i in range(self.numActive):
            print("{:d}".format(self.bandsActive[i]))
        print("Frequency (MHz) | # channels")
        print("----------------------------")
        for i in range(self.num):
            print("{:.2f}          | {:d}".format(self.freq[i],self.nchan[i]))
    

class CommonInfo:
    """ File information class
    """
    calcfile = ""
    threadsfile = ""
    outputfile = ""
    nbits = 0

    def print(self):
        print("calcfile: "+self.calcfile)
        print("threadsfile: "+self.threadsfile)
        print("outputfile: "+self.outputfile)


class TelescopeInfo:
    """ Telescope class
    """
    num = 0


class Header():
    """ This Header class holds all visibility header
        Information. The Choice over a pandas data frame
        is that it is slightly faster, easier for the user
        to slice and control.
    """

    def __init__(self,vissize,tnum,bands):
        self.syncword = b'\x00\xff\x00\xff'
        self.binaryformat = 1
        self.config_idx = 0
        self.source_idx = 0
        self.pulsar_idx = 0
        self.headersize = vissize

        fnum = len(bands) # number of active bands
        bnum = tnum*(tnum+1)//2 # number of unique baselines

        #### [BASELINES] ####
        baselines = np.zeros(bnum,dtype = np.int32) # create baselines
        cDiFX.get_baselines(tnum,baselines) # create antennas
        self.baselines = np.repeat(baselines,fnum*4) # 32 bit int

        #### [ANTENNAS] ####
        ant1 = baselines % 256
        ant2 = (baselines - ant1) // 256
        self.ant1 = np.repeat(ant1,fnum*4) # 32 bit int
        self.ant2 = np.repeat(ant2,fnum*4) # 32 bit int 

        #### [CORRELATION] ####
        cor = ant1 == ant2
        cor = np.repeat(cor,fnum*4)
        corval = np.empty(bnum*fnum*4,dtype = str)
        corval[:] = "C"
        corval[cor] = "A"
        self.corr = corval # string

        #### [FREQUENCY] ####
        self.freq = np.tile(bands,4*bnum).astype(dtype = np.int32) # 32 bit int

        #### [POLARISATION] ####
        self.pol = np.tile(np.repeat(["XX","XY","YX","YY"],fnum),bnum) # 16 (8 + 8) bit utf-8 encoded char

        #### [TIME] ####
        self.MJD = np.zeros(vissize,dtype = np.int32) # 32 bit int
        self.seconds = np.zeros(vissize,dtype = np.float64) # 64 bit float

        #### [U,V,W] ####
        self.position = np.zeros((vissize,3),dtype = np.float64) # 64 bit float

        #### [WEIGHT] ####
        self.weight = np.zeros(vissize,dtype = np.float64) # 64 bit float

        #### [# CHANNELS] ####
        self.channels = np.zeros(vissize,dtype = np.int32) # 32 bit int

        #### [@] ####
        self.exist = np.zeros(vissize,dtype = bool) # boolean


    #create dataframe from data
    def print(self,full=False):
        ## create DataFrame ##
        pdout = pd.DataFrame(data={"freq":self.freq,"baseline":self.baselines,
                                    "ant1":self.ant1,"ant2":self.ant2,"pol":self.pol,
                                    "MJD":self.MJD,"seconds":self.seconds,"corr":self.corr,
                                    "U (m)":self.position[:,0],"V (m)":self.position[:,1],
                                    "W (m)":self.position[:,2],"weight":self.weight,
                                    "# chan":self.channels,"@":self.exist})
        
        ## print out DataFrame ##
        if full:
            print(pdout.to_string())
        else:
            print(pdout)
        

    #copy function
    def copy(self):
        #NOTE THIS WAS CHANGED
        return copy.deepcopy(self)
                                    





#############################################
###### ADDITIONAL UTIL FUNCTIONS ############ 
#############################################

#function to quickly read inputfile data
def inputfile_quickread(inputfile):

    #returns a dict with split up tables
    with open(inputfile,'r') as f:
        filedata = f.read().split('!')#assuming ! is delimeter to each table header

        f.close()
    
    output = {} #specify empty string

    name = filedata[0]
    for i in range(1,len(filedata)):
        name = name.replace("#","") #remove #
        name = name.strip() #gets rid of white space on the sides
        currentstr = np.array(filedata[i].split("\n"))
        output[name] = currentstr[1:-1] #append to dictionary
        name = currentstr[-1]
    
    #return dictionary
    return output





#Get freq and baseline info 
#TODO - in theory this could b accelerated assuming the text file is in the same format
#for now, just search.
def get_inputfile_info(inputfile,baselinedata,frequencydata,telescopedata,commondata):
    """ This function is a bit messy and could be improved, however, it works.

        The function loades in an inputfile, splits it up into it's different tables
        and looks through each table, grabbing important information.

        I was likely overzelous making this function. I was too focussed on making it as
        fast as possible. Indeed it is quite fast, but not in the grander scheme of things.
    """
    
    #load in input file
    inputdata = inputfile_quickread(inputfile)

    #read common info
    """ Go to the Common Settings table and look for the:
        .calc filename,
        .threads filename,
        output filename 
    
    """

    calc_id = np.flatnonzero(np.char.find(inputdata["COMMON SETTINGS"],"CALC FILENAME") != -1)[0]
    thread_id = np.flatnonzero(np.char.find(inputdata["COMMON SETTINGS"],"CORE CONF FILENAME") != -1)[0]
    output_id = np.flatnonzero(np.char.find(inputdata["COMMON SETTINGS"],"OUTPUT FILENAME") != -1)[0]
    bits_id = np.flatnonzero(np.char.find(inputdata["DATASTREAM TABLE"],"QUANTISATION BITS") != -1)[0]
    commondata.calcfile = inputdata["COMMON SETTINGS"][calc_id].split(':')[1].strip()
    commondata.threadsfile = inputdata["COMMON SETTINGS"][thread_id].split(':')[1].strip()
    commondata.outputfile = inputdata["COMMON SETTINGS"][output_id].split(':')[1].strip()
    commondata.nbits = inputdata["DATASTREAM TABLE"][bits_id].split(':')[1].strip()


    #read baseline info
    """ Go to the baseline table and look for:
        num of freqency band of each baseline (find all bands being used) -> num bands,
        num of telescope,
        maximum number of possible baselines,
        get baseline offset from first baseline stream.

    
    """

    baseline_freqs_id = np.flatnonzero(np.char.find(inputdata["BASELINE TABLE"],"NUM FREQS") != -1)
    baseline_freqs = np.char.split(inputdata["BASELINE TABLE"][baseline_freqs_id],sep = ":")
    frequencydata.numActive = np.max(np.array([xline[1].strip() for xline in baseline_freqs],dtype = int))

    numtelescope = int(inputdata["TELESCOPE TABLE"][0].split(':')[1].strip())
    baselinedata.num = numtelescope * (numtelescope - 1) // 2
    baselinedata.numa = numtelescope * (numtelescope + 1) // 2
    telescopedata.num = numtelescope

    #read first baseline, get baseline offset
    basestreamA = int(inputdata["BASELINE TABLE"][1].split(':')[1].strip())
    basestreamB = int(inputdata["BASELINE TABLE"][2].split(':')[1].strip())
    baselinedata.offset = 256*(basestreamA+1) + basestreamB + 1



    #read frequency info
    """ Look through Frequency table and find:
        get num of channels in each frequency band,
        get frequency (in MHz) for each frequency band,
        get active frequency bands,
    
    """
    frequencydata.num = int(inputdata["FREQ TABLE"][0].split(':')[1].strip())
    chan_id = np.flatnonzero(np.char.find(inputdata["FREQ TABLE"],"NUM CHANNELS") != -1)
    chanstring = np.char.split(inputdata["FREQ TABLE"][chan_id],sep = ":")
    avg_id = np.flatnonzero(np.char.find(inputdata["FREQ TABLE"],"CHANS TO AVG") != -1)
    avgstring = np.char.split(inputdata["FREQ TABLE"][avg_id],sep = ":")
    chan_val = np.array([yline[1].strip() for yline in chanstring],dtype = int)
    avg_val =  np.array([yline[1].strip() for yline in avgstring],dtype = int)

    if np.any(chan_val % avg_val):
        print("Channel averaging returns non-interger number of channels, Aborting...")
        return
    frequencydata.nchan = np.asarray(chan_val // avg_val,dtype = int)

    freqs_id = np.flatnonzero(np.char.find(inputdata["FREQ TABLE"],"FREQ (MHZ)") != -1)
    freqstring = np.char.split(inputdata["FREQ TABLE"][freqs_id],sep = ":")
    frequencydata.allfreq = np.array([zline[1].strip() for zline in freqstring],dtype = float)

    #get index of active frequency bands
    recfreq_ind = np.flatnonzero(np.char.find(inputdata["DATASTREAM TABLE"],"REC FREQ INDEX") != -1)
    recfreqstring = np.char.split(inputdata["DATASTREAM TABLE"][recfreq_ind],sep = ":")
    recfreqarr = np.array([wline[1].strip() for wline in recfreqstring],dtype = int)
    frequencydata.bandsActive = np.unique(recfreqarr)

    frequencydata.freq = frequencydata.allfreq[frequencydata.bandsActive]

    if len(frequencydata.bandsActive) != frequencydata.numActive:
        print("Number of recorded bands does not match number of Active bands... Aborting.")
        return

    



#difx class
class difx:

    def __init__(self,inputfile):
        
        
        ######### data parameters ###########
        self.nvis = 0 # number of visibility's 
        self.bin = 0 # current bin
        self.pbin = 0 # previous bin
        self.nbins = 0 # number of bins
        self.npol = 4 #assumed by default
        self.nbits = 4 # assumed by default, but will be updated when loading in
        

        ############ data and header #############
        self.vis = [] # visibility data
        self.header: [] # header data (buffer)
        self.vischanlen = [] 
        self.visshape = [0,0] # DIFX shape (vis,nchan) MAX
        self.missingvis = [] #depreciated


        ############# file parameters #############
        self.currentinputfile = inputfile
        self.currentdifxfile = ""
        self.bytesize = 0
        

        ############ input file stats ############
        self.base_info = BaselineInfo() # baseline info
        self.freq_info = FreqInfo() # frequency info
        self.tele_info = TelescopeInfo() # telescope info
        self.common = CommonInfo() # common input file info 
        get_inputfile_info(inputfile,self.base_info,self.
                        freq_info,self.tele_info,self.common)
        self.difxdir = self.common.outputfile

        # set nbits
        self.nbits = self.common.nbits
        

        ######### get list of bin #############
        self.binfiles = np.array(sorted(glob.glob(self.difxdir + "/DIFX*"))) # list of bin filepaths
        self.nbins = len(self.binfiles)
        self.binarr = np.array([int(binstr[-4:]) for binstr in self.binfiles]) # list of bin numbers


        ######### visibility and header arrays ################
        self.visshape[0] = (self.base_info.num + self.tele_info.num) * 4 * self.freq_info.numActive 
        self.visshape[1] = np.max(self.freq_info.nchan[self.freq_info.bandsActive])#create MAX shape bounds

        self.vis = np.zeros(self.visshape,dtype = "complex64") # create Visibility array
        self.header = np.empty((self.visshape[0],74),dtype = np.dtype('b')) # create Header array

        #create buffer arrays
        self.visbuffer = bytes(b'\x00'*self.visshape[0]*self.visshape[1]*8)
        self.headbuffer = bytes(b'\x00'*self.visshape[0]*74)


        ######### Background varaibles used ##############
        self.bandActiveindex = np.zeros(self.freq_info.num,dtype = int) # for use in c-loading function
        for i in range(self.freq_info.numActive):
            self.bandActiveindex[self.freq_info.bandsActive[i]] = i 


        ######### initialize header data frame ###############
        self.headerTable = Header(self.visshape[0],self.tele_info.num,self.freq_info.bandsActive)

        ##NOTE - The difx file data is not complete, but under the assumption that it is near complete, then
        #ideal indexing would use boolean arrays. This is even more true for larger arrays due to boolean value
        #caching. Hence the choice of using boolean indexing is 2 fold. First, It is faster for 'near' complete 
        #data, second, This speedup will increase for larger arrays, assuming near completeness in data.

        return



    ################################################
    ############## LOAD BIN FUNCTION ###############

    def loadbin(self, bin = 1,init = False):

        ######## reset variables ########
        # number of channels in each visibility
        vischanlen = np.zeros(self.visshape[0],dtype = np.int32)

        # visibility's recorded (True or False)
        vispresent = np.ones(self.visshape[0],dtype = bool)


        ######## get DIFX file ########
        if init: #get bin file based on bin ID 
            #load last bin, is used to avoid loading in bin 0
            binID = np.where(self.binarr == self.binarr[-1])[0]
            bin = int(self.binarr[-1][-4:])


        else:
            
            binID = np.where(self.binarr == bin)[0]
            if not len(binID):
                print("DiFX file for bin {:d} not found... Aborting".format(bin))
                return


        #load in file using os.read
        self.bin = bin
        self.currentdifxpath = self.binfiles[binID[0]] # save file metadata
        self.currentdifxfile = os.path.basename(self.currentdifxpath)
        self.bytesize = os.stat(self.currentdifxpath).st_size

        ######## Load in DiFX File ########
        #NOTE: I found that this method is on average 1.5-2X faster then
        #python open.
        difx_f = os.open(self.currentdifxpath,os.O_RDONLY) # open file
        bytearr = os.read(difx_f,self.bytesize)
        os.close(difx_f) # close file

        ######## create varaibles to pass into c-loading function ########
        """ number of file bytes, 
            baseline index offset, 
            number of actice bands,
            max num of channels,
            max num of visibility's,
            num of antenna,
            nbits,
        """
        var = np.array([self.bytesize,self.base_info.offset,self.freq_info.numActive,
                        self.visshape[1],self.visshape[0],self.tele_info.num],dtype = np.int32)

        

        ######## split databuffer into vis and header through (c) ########
        """
            Function first resets header and visibility buffer arrays from any previous 
            difx file, then incrementally splits the header and vis data into a 
            predetermined formatted structure. i.e. each visibility has a predetermined 
            index, the method is robustenough to work even if visibility's are missing 
            or the number of channels do not match between visibilities.
        """
        self.nvis = cDiFX.byte2vis(bytearr,self.visbuffer,self.headbuffer,
                    np.asarray(self.freq_info.nchan,dtype = np.int32),vischanlen,
                    np.asarray(self.bandActiveindex,dtype = np.int32),var)



        ######### convert buffer arrays to numpy arrays ########
        #NOTE: the bytes array is immutable in python. Even though this can be
        #bypassed in c by manipulated the raw data (which it is), the python object 
        #cannot be changed within the Python API. Hence if we want to manipulate the 
        #visibility data or header data, we need to create a bytesarray object (which IS immutable) 
        #from the buffer data before typing it, else the immultability will pass over
        #to the data structures and writing permission will be FALSE (and uncahangable).
        self.vis = np.frombuffer(bytearray(self.visbuffer),dtype = "complex64").reshape(self.visshape)
        self.header = np.frombuffer(bytearray(self.headbuffer),dtype = np.dtype('b')).reshape(self.visshape[0],74)

        

        ######## get visibility's that exist in DiFX file ########
        missingvis = np.where(vischanlen == 0)[0] # get index of missing vis
        vispresent[missingvis] = False


        ######## Fill Header Table ########
        self.headerTable.MJD = np.frombuffer(self.header[:,12:16].reshape(self.header.shape[0]*4),dtype = np.int32)
        self.headerTable.seconds = np.frombuffer(self.header[:,16:24].reshape(self.header.shape[0]*8),dtype = np.float64)
        self.headerTable.weight = np.frombuffer(self.header[:,42:50].reshape(self.header.shape[0]*8),dtype = np.float64)
        self.headerTable.exist = vispresent
        self.headerTable.channels = vischanlen

        #### fill in UVW coords ####
        corrds = np.frombuffer(self.header[0::self.freq_info.numActive*4,50:74].reshape(24*self.base_info.numa),
                                dtype = np.float64).reshape(self.base_info.numa,3)

        # corrds = 
        self.headerTable.position = np.repeat(corrds,self.freq_info.numActive*4,axis = 0)


        ######## set flags ########
        self.pbin = bin



    ################################################
    ############## GET VIS FUNCTION ################

    def get_vis_index(self,freq = None,corr=None,baseline=None,
                        ant1=None,ant2=None,pol=None,weight=False,exist=True):
        """Get visibility's given specifications,
            This function goes through each visibility property
            and matches them with inputs.


            freq: frequency band ID
            corr: correlation [C (cross),A (auto)]
            baseline: baseline ID
            ant1: First antenna number
            ant2: Second antenna number
            pol: Polarisation [XX, XY, YX, YY]
            weight: vis header weight 
            exist: True if visibility was recorded in loaded difx file
        """
        visind = np.arange(self.visshape[0])
        transf_ind = []

        if freq is not None:
            transf_ind = self.headerTable.freq[visind] == freq
            visind = visind[transf_ind]
        
        if exist:
            transf_ind = self.headerTable.exist[visind]
            visind = visind[transf_ind]

        if weight:
            transf_ind = self.headerTable.weight[visind] > 0
            visind = visind[transf_ind]

        if corr is not None:
            transf_ind = self.headerTable.corr[visind] == corr
            visind = visind[transf_ind]

        if baseline is not None:
            transf_ind = self.headerTable.baselines[visind] == baseline
            visind = visind[transf_ind]

        if ant1 is not None:
            transf_ind = self.headerTable.ant1[visind] == ant1
            visind = visind[transf_ind]

        if ant2 is not None:
            transf_ind = self.headerTable.ant2[visind] == ant2
            visind = visind[transf_ind]

        if pol is not None:
            transf_ind = self.headerTable.pol[visind] == pol
            visind = visind[transf_ind]

        
        return visind

    
    #get freq function
    def get_freq_index(self,freq):
        """ Return boolean array vis entries in subband "freq"
        """
        return self.headerTable.freq == freq

    
    #get freq weights
    def get_freq_weights(self):
        #returns all weights in each active freq band

        freqweights = np.zeros(self.freq_info.numActive,dtype = float)
        idx = np.arange(self.visshape[0])


        for i in range(self.freq_info.numActive):
            freqweights[i] = self.headerTable.freq[i::self.freq_info.numActive][0]

        #return
        return freqweights


    #get nonzero vis
    def filter_nonzero_index(self):
        #return indexes of non zero visibility's
        return np.where(np.count_nonzero(self.vis,axis=1) > 0)[0]


    #get zero vis
    def filter_zero_index(self):
        #return indexes of zero visibility's
        return np.where(np.count_nonzero(self.vis,axis=1) == 0)[0]


    ################################################
    ############## GET CAS FUNCTION ################

    def get_cas(self,freq = None,corr = "C"):
        #get cross amplitude sum (CAS)
        #Ideal to use this function due to zero padding

        return np.sum(np.abs(np.mean(self.vis[self.get_vis_index(freq = freq,
                exist = True,corr = corr)],axis = 1)))


    
    ################################################
    ############## SAVE BIN FUNCTION ###############
    
    def savebin(self,filename):
        """ Takes a filename for the new difx file.
            The header data is reconstructed from the 
            "Header" object.

            Which 

        """

        ######## create buffers ########
        visid = self.filter_nonzero_index()
        chans = self.headerTable.channels[visid] # get number of channels with non-zero data
        num = len(visid)

        headbuffer = np.empty((num,74),dtype = np.byte) # construct header buffer

        bytelen = int(np.sum(chans)*8 + num*74)
        outbuffer = create_string_buffer(b'\x00'*bytelen,bytelen) # construct output buffer to write to file

        ######## build header buffer ########
        ## syncword ##
        headbuffer[:,0:4] = bytearray(self.headerTable.syncword)
        ## binary format ##
        headbuffer[:,4:8] = bytearray(self.headerTable.binaryformat.to_bytes(4,'little'))
        ## baseline ##
        headbuffer[:,8:12] = np.array(bytearray(self.headerTable.baselines[visid].tobytes())
                                        ,dtype = np.byte).reshape(num,4)
        ## MJD ##
        headbuffer[:,12:16] = np.array(bytearray(self.headerTable.MJD[visid].tobytes())
                                        ,dtype = np.byte).reshape(num,4)
        ## seconds ##
        headbuffer[:,16:24] = np.array(bytearray(self.headerTable.seconds[visid].tobytes())
                                        ,dtype = np.byte).reshape(num,8)
        ## config index ##
        headbuffer[:,24:28] = bytearray(self.headerTable.config_idx.to_bytes(4,'little'))
        ## source index ##
        headbuffer[:,28:32] = bytearray(self.headerTable.source_idx.to_bytes(4,'little'))
        ## freq index ##
        headbuffer[:,32:36] = np.array(bytearray(self.headerTable.freq[visid].tobytes())
                                        ,dtype = np.byte).reshape(num,4)
        ## polarisation pair ##
        headbuffer[:,36:38] = np.array(bytearray(bytes(np.char.encode(self.headerTable.pol[visid],'utf-8')))
                                        ,dtype = np.byte).reshape(num,2)               
        ## pulsar bun ##
        headbuffer[:,38:42] = bytearray(self.headerTable.pulsar_idx.to_bytes(4,'little'))
        ## Weights ##
        headbuffer[:,42:50] = np.array(bytearray(self.headerTable.weight[visid].tobytes())
                                        ,dtype = np.byte).reshape(num,8)
        ## U,V,W ##
        headbuffer[:,50:74] = np.array(bytearray(self.headerTable.position[visid].tobytes())
                                        ,dtype = np.byte).reshape(num,24)
                                    
        
        ######## build bytearray for c function ########
        outvis = self.vis[visid].tobytes()
        outheader = headbuffer.tobytes()


        ######## run c function ########
        """ Takes vis buffer and chan len array and removes zero padding,
            create c-contigous array of header and vis bytes "outbuffer".
        
        """
        cDiFX.vis2byte(outbuffer,bytelen,outvis,self.visshape[1],
                                outheader,chans)


        ######## write buffer to file ########
        with open(filename,'wb') as fileout:
            fileout.write(outbuffer)
            fileout.close()

        return


    # ################################################
    # ############ PYTHON UTIL FUNCTIONS #############

    # #print method
    def __str__(self):
        stringreturn = """
        InputFile: {:s}
        Bin: {:d}
        Active Bands:""".format(self.currentinputfile,
                                self.bin) + np.array2string(self.freq_info.freq) + """
        # vis: {:d}
        """.format(self.nvis)
        return stringreturn
        #print info of bin loaded in
    

    #copy function
    def copy(self):
        """ The class is complex, hence a deepcopy is nessesary to
            Copy all contents. 
        
        """
        #return copy of difx instance
        return copy.deepcopy(self)
    

    


    # ################################################################################################################################################
    # ############ EXTENDED UTILITIES ################################################################################################################

