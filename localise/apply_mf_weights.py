##===============================================##
##===============================================##
## Author: Tyson Dial
## Email: tdial@swin.edu.au
## Date Created: 01/07/2023
## Last Updated: 01/10/2024 
##
##
## This script takes a wmask structure containing
## time(freq-) dependent weights as well as a 
## finder and rfi mask of the correlated data, then
## weights the onpulse finder bins, subtracts the 
## rfi bins and scrunches the bins together.
##
## NOTE: The time-dependent weights are applied 
## here whilst the freq-dependent weights are
## applied to the difx header weights which are
## implemented during imaging (specifically - 
## gridding).
##
##===============================================##
##===============================================##


## Imports
import numpy as np
import argparse, os, sys, glob
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from time import time
from fdifx import fdifx,util


#wmask class
class wmask:

    def __init__(self,wmaskfilename):
        
        #load in wmask and get different weights
        with open(wmaskfilename,"rb") as file:
            self.tw = np.load(file)           #time dependent weights (formed by frequency scrunching)
            self.fw = np.load(file)           #freq dependent weights (formed by time scrunching)
            self.tmask = np.load(file)        #time mask
            self.findermask = np.load(file)   #finder mask
            self.rfimask = np.load(file)      #rfi mask

            #fieldmask: ndarray          #field mask (FOR NOW JUST TAKE ZEROTH BIN)







def create_new_difx_dirs(in_dir, out_dir):
    """
    Create new c*_f* directories to store weighted and scrunched .DIFX bin,
    this will also make copies of the .input files with proper Header infomation.


    Parameters
    ----------
    in_dir: str
        original directory of c*_f* dirs
    out_dir : str
        new directory to put new c*_f* dirs

    
    Returns
    -------
    cf_inputs: list
        list of input files to load in .difx data
    cf_outputs: list
        list of input files to save weighted/scrunched .difx bin data

    """

    # get cf inputs that are sorted properly
    cf_inputs = sorted(glob.glob(os.path.join(in_dir,"c*_f*/*D2D.input"))) 

    print(cf_inputs)
    
    cf_outputs = [""] * len(cf_inputs)     # for output directory .input files

    # create new output dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    for i, cf_input in enumerate(cf_inputs):
        # create new c*_f* dir
        cf_dir = os.path.join(out_dir,cf_input.split(os.path.sep)[-2])
        os.mkdir(cf_dir)

        # make new input file by copying contents to new directory
        cf_outputs[i] = os.path.join(cf_dir,os.path.basename(cf_input))
        with open(cf_input, 'r') as file:
            lines = file.readlines()
        for j,line in enumerate(lines):
            if "OUTPUT FILENAME" in line:
                lines[j] = "OUTPUT FILENAME:    " + os.path.abspath(cf_outputs[i])[:-5] + "difx\n"
                break
        
        # copy to new file
        with open(cf_outputs[i], 'w') as file:
            file.writelines(lines)
        
        # create new difx dir
        os.mkdir(cf_outputs[i][:-5] + "difx")
    
    return cf_inputs, cf_outputs














#args
def get_args():

    #arguments
    parser = argparse.ArgumentParser()

    #data directory
    parser.add_argument('-d', help = "Input directory for data (difx bin folder)", type = str)

    #path to weight file .npy
    parser.add_argument('-w','--weights',help = "htr weights to apply to correlated data")

    #freq channels to flag
    parser.add_argument('--flagbins',help = "flag frequency bins")

    #plot data
    parser.add_argument('-p','--plot',action="store_true",help = "Save a plot of the weighting and scrunching for checking purposes")
    args = parser.parse_args()

    # args.inputsuffix = "D2D.input"
    # args.outputsuffix = "D2DW.input"

    #get additional parameters
    #initialise data path directories
    # args.cfinputs = util.get_inputfiles(args.d,args.inputsuffix)
    # args.cfoutputs = util.create_difx_output_dir(args.d,
    #                 args.cfinputs,args.inputsuffix,args.outputsuffix)
    args.cfinputs, args.cfoutputs = create_new_difx_dirs(args.d, os.path.join(args.d, "D2DW"))

    #get frequencies of cards and index ordering
    args.freq = []
    for i,inp in enumerate(args.cfinputs):
        difx = fdifx.difx(inp)
        args.freq.append(difx.freq_info.freq[0])
        args.freq.append(difx.freq_info.freq[1])

                                            
    args.freq = np.array(args.freq) #make into numpy array
    print(args.freq)

    #get size of weight array (and data array)
    args.tbins = difx.nbins
    args.fbins = len(args.cfinputs)*2 #assume 4MHz coarse channels


    return args




















def f_flag(args,wmask):
    
    #check flag
    if args.flagbins is not None:

        #take range of frequencies
        if '~' in args.flagbins:
            edges = args.flagbins.split('~')
            flaggedbins = np.arange(int(edges[0]),int(edges[1])+1)
        
        else:
            #evaluate a stringed array
            flaggedbins = np.array(eval(args.flagbins))
    
        print("Flagged Frequency Bins: ")

        #update tmask + rfimask elements (which updates elements in global scope)
        wmask.tmask[flaggedbins] = 0
        wmask.rfimask[flaggedbins] = 0






















#main weighting function 
def weight_scrunch(args,wmask):

    #global time stamp
    t1 = time()

    ##=======================##
    ## Additional functions  ##
    ##=======================##

    def concat_header_table(tab1,tab2):
        #concatenate info from tab2 into tab 1
        #NOTE: only concatenate possible entries that are missing in tab1
        m_idx = tab1.exist == False # get missing entries

        #update 
        tab1.MJD[m_idx] = tab2.MJD[m_idx]               # MJD (in units of days)
        tab1.seconds[m_idx] = tab2.seconds[m_idx]       # seconds (to add to MJD)
        tab1.position[m_idx] = tab2.position[m_idx]     # [U,V,W] coords
        tab1.channels[m_idx] = tab2.channels[m_idx]     # fine-channel num
        tab1.exist[m_idx] = tab2.exist[m_idx]           # FLAG for if visibility is missing 

        return






    ##======================##
    ## Create f-corrections ##
    ##======================##

    #freq array is not in ascending order, but order
    #at which they appear in the card fpgas (not nessesarily the same).
    #get card fpga indices relative to band of ascending order, used for weighting.
    f_idx = util.get_inverse_sort_index(args.freq)[::-1]

    #get frequency amplitude corrections and corresponding weight corrections
    f_amp_corr = np.zeros(args.fbins)
    f_weight_corr = np.zeros(args.fbins)
    print(wmask.fw)
    print(wmask.tmask)
    for i,fweight in enumerate(wmask.fw):
        if fweight <= 0:
            f_amp_corr[i] = 0
            f_weight_corr[i] = 0
        else:
            # Amplitude correction is 1 / intensity (power domain flattening)
            f_amp_corr[i] = 1.0 / fweight # This will make the FRB spectrally flat
            # The visibility weights (representing inverse variance) must scale by 1 / (amplitude_correction)^2
            f_weight_corr[i] = 1.0 / (f_amp_corr[i] ** 2)
            # Ultimately, the product of the amplitude correction and the visibility weight leads to effectively multiplying by the original fw
    
    #NOTE: Each card fpga is made of 2 4MHz coarse channels, hence
    #each inputfile iteration will need two sets of varaibles, one 
    #for each coarse channel.

    ## ALSO create additional variables
    # cross apmlitude sum of difx data [unweighted,weighted,weighted+rfi subtracted]
    CAS = np.zeros((3,args.fbins,args.tbins-1),dtype = float)

    hw_difx = []                         #header weights from loaded difx correlations 

    hwt = np.zeros((2,args.tbins))       #time dependent weighted header weights [chan1,chan2]

    #set weights in rfi mask to -1 regardless of !=0 value (anything 0 is ignored)
    wmask.rfimask[wmask.rfimask != 0] = -1


    ##==========================##
    ## WEIGHTING AND SCRUNCHING ##
    ##==========================##

    for i,inputfile in enumerate(args.cfinputs):

        cardfpga = os.path.basename(os.path.dirname(inputfile)) #get c*_f* specifier

        print("Loading "+cardfpga+": {:d}/{:d}".format(i+1,len(args.cfinputs))) #log current cardfpga

        #create difx classes
        difx = fdifx.difx(inputfile) #current card fpga
        rfidifx = difx.copy()        #rfi
        outdifx = difx.copy()        #weighted and scrunched difx
        rfidifx_temp = difx.copy()

        #reset variables
        hwt *= 0


        #NOTE: The algorithm specified below looks at the tmask from wmask and determines
        #which time bins (.DIFX) to load in (if w_t > 0). Since two coarse channels exist
        #within each cardfpga, we load in all possible time bins for which w_t in either
        #f channel is >0.

        ## GET FINDER BINS TO LOAD IN ##
        bins_chan1 = np.where(wmask.tmask[f_idx[i*2]] != 0)[0]              #bins for chan 1
        bins_chan2 = np.where(wmask.tmask[f_idx[i*2+1]] != 0)[0]            #bins for chan 2

        finder_bins = np.unique(np.concatenate((bins_chan1,bins_chan2)))    #complete list of bins to load


        # GET RFI BINS TO LOAD IN ##
        bins_chan1 = np.where(wmask.rfimask[f_idx[i*2]] != 0)[0]
        bins_chan2 = np.where(wmask.rfimask[f_idx[i*2+1]] != 0)[0]

        rfi_bins = np.unique(np.concatenate((bins_chan1,bins_chan2)))       #complete list of rfi bins
        

        ##===============##
        ## RFI ESTIMATE  ##
        ##===============##

        rfihw = np.zeros((2,rfi_bins.size))               #rfi header weights from difx [chan1,chan2]

        for k,bin in enumerate(rfi_bins):
            difx.loadbin(bin = bin)

            #get freq header weights and save to array
            rfihw[:,k] = difx.get_freq_weights()

            #weight
            difx.vis[0::2] *= wmask.rfimask[f_idx[i*2],bin] * rfihw[0,k]
            difx.vis[1::2] *= wmask.rfimask[f_idx[i*2+1],bin] * rfihw[1,k]

            #add to rfi difx class (t-scrunching rfi data then averaging)
            rfidifx.vis += difx.vis

        
        #Take average of data (knowing that zeros will be present because only a handful
        # of rfi bins was chosen.)
        rfiN = abs(np.sum(wmask.rfimask[f_idx[i*2:i*2+2]],axis = 1)) #sum current 2 4MHz freq channels
        rfihw_tot = np.mean(rfihw,axis=1)                            #average header weights

        #NOTE: the reason for averaging the header weights and dividing through the rfi visibility's
        #is to account for the header weights of the rfi data when calculating the average rfi to subtract.

        if rfiN[0] > 0: #chan 1
            rfidifx.vis[0::2] /= rfiN[0] * rfihw_tot[0]
        
        if rfiN[1] > 0: #chan 2
            rfidifx.vis[1::2] /= rfiN[1] * rfihw_tot[1]


        ##==================##
        ## FINDER ESTIMATE  ##
        ##==================##

        for j,bin in enumerate(finder_bins):
            difx.loadbin(bin = bin)

            hw_difx = difx.get_freq_weights()                               #load current difx weights
            hwt[0,j] = hw_difx[0] * wmask.tmask[f_idx[i*2],bin]             #update new chan header weight
            hwt[1,j] = hw_difx[1] * wmask.tmask[f_idx[i*2+1],bin]

            #calculate CAS of loaded data (without weighting)
            #band active is just band present in difx data, i.e. bandsactive = [8,9]
            CAS[0,i*2,bin-1] = difx.get_cas(difx.freq_info.bandsActive[0])  
            CAS[0,i*2+1,bin-1] = difx.get_cas(difx.freq_info.bandsActive[1])

            #weight difx data
            difx.vis[0::2] *= hwt[0,j]
            difx.vis[1::2] *= hwt[1,j]

            #add to total visibility's
            outdifx.vis += difx.vis

            #calculate CAS of weighted difx data
            CAS[1,i*2,bin-1] = CAS[0,i*2,bin-1] * hwt[0,j]
            CAS[1,i*2+1,bin-1] = CAS[0,i*2+1,bin-1] * hwt[1,j]

            #calculate CAS of weighted and RFI subtracted visibility's (subtract complex vis first)
            rfidifx_temp.vis = rfidifx.vis * 1.0
            rfidifx_temp.vis[0::2] *= hwt[0,j]
            rfidifx_temp.vis[1::2] *= hwt[1,j]

            rfidifx_temp.vis += difx.vis                                #subtract current difx data
            rfidifx_temp.vis[difx.get_vis_index(corr = "A")] *= 0       #remove auto correlations

            CAS[2,i*2,bin-1] = np.sum(np.abs(np.mean(rfidifx_temp.vis[0::2],axis = 1)))
            CAS[2,i*2+1,bin-1] = np.sum(np.abs(np.mean(rfidifx_temp.vis[1::2],axis = 1)))

            #NOTE: now the header table is filled in, the purpose of thise function is to make sure
            # we consider possible missing visibility's across .DIFX files and frequencies. Assuming that the 
            # data is near complete (i.e. there is the possiblility of only a handful of visibility's missing
            # in a handful of .DIFX files) this method is fine. The function concatenates the header info
            # from multiple .DIFX files.
            concat_header_table(outdifx.headerTable,difx.headerTable)



        ##===============================##
        ## RFI subtraction and averaging
        ##===============================##

        #average header weights, then average t-scrunched visibility's
        hw_tot = np.sum(hwt,axis = 1)

        #weight rfi before subtraction
        rfidifx.vis[0::2] *= hw_tot[0]
        rfidifx.vis[1::2] *= hw_tot[1]

        #make sure visibility's exist in the data also before rfi subtraction
        rfidifx.vis[outdifx.headerTable.exist == False] *= 0

        #do rfi subtraction
        outdifx.vis += rfidifx.vis

        #average
        if hw_tot[0] > 0:
            outdifx.vis[0::2] /= hw_tot[0]
        
        if hw_tot[1] > 0:
            outdifx.vis[1::2] /= hw_tot[1]
        



        ##=========================##
        ## do frequency correction ##
        ##=========================##

        #apply amplitude correction factor to flatten spectrum
        outdifx.vis[0::2] *= f_amp_corr[f_idx[i*2]]
        outdifx.vis[1::2] *= f_amp_corr[f_idx[i*2+1]]

        #update weights in header table based on amplitude correction factor
        outdifx.headerTable.weight[0::2] = hw_tot[0] * f_weight_corr[f_idx[i*2]]
        outdifx.headerTable.weight[1::2] = hw_tot[1] * f_weight_corr[f_idx[i*2+1]]

        #remove autocorrelations
        ac_idx = outdifx.headerTable.corr == "A"
        outdifx.headerTable.weight[ac_idx] = 0
        outdifx.vis[ac_idx] *= 0




        ##==================================##
        ## Save weighted + t-scrunched data
        ##==================================##

        filename_out = args.cfoutputs[i][:-5] + "difx/" + difx.currentdifxfile[:-4] + "0000"

        outdifx.savebin(filename_out)


    
    print("Weighting and Scrunching completed...")
    print("Execution time: {:.2f} s".format(time() - t1))

    return CAS                          


























#diagnostic plot for time weighting and rfi subtraction
def diagplot_weighting(data,t_data,tsub_data):

    #create figure and axes
    fig = plt.figure("Diagnostic Plot: Time weighting and RFI subtraction",figsize = (12,10))
    ax1 = fig.add_axes([0.06,0.07,0.29,0.78])
    ax2 = fig.add_axes([0.37,0.07,0.29,0.78])
    ax3 = fig.add_axes([0.68,0.07,0.29,0.78])
    ax4 = fig.add_axes([0.06,0.90,0.91,0.05])


    ax1.imshow(data/np.max(data),aspect = 'auto',vmin = 0, vmax = 1)
    ax1.set_title("CAS",fontsize = 16)
    ax1.set_ylabel("Frequency Bin",fontsize = 16)

    ax2.imshow(t_data/np.max(t_data),aspect = 'auto',vmin = 0, vmax = 1)
    ax2.set_title("Time Weighted CAS",fontsize = 16)
    ax2.set_xlabel("Time Bin",fontsize = 16)
    ax1.get_yaxis().set_visible(False)

    ax3.imshow(tsub_data/np.max(tsub_data),aspect = 'auto',vmin = 0,vmax = 1)
    ax3.set_title("Weighted CAS [RFI subtracted]",fontsize = 16)
    ax1.get_yaxis().set_visible(False)

    #create colorbar
    colbar = np.linspace(0,1,512).reshape(1,512)
    ax4.imshow(colbar,aspect = 'auto',vmin = 0,vmax = 1,extent = [0,1,0,1])
    ax4.get_yaxis().set_visible(False)
    ax4.set_xticks([0,1.0])

    #save file to png 
    filename = os.getcwd() + "/applied_masking.png"
    print("Saving diagnostic plot of applied masking/weights to: "+filename+"...")

    plt.savefig(filename)


    AX = [ax1,ax2,ax3,ax4]

    return (fig, AX)

























#diagnostic plot to check alignment of weights and data
def diagplot_align(data,weights,suffix = ""):
    
    
    #create figure
    fig = plt.figure("Diagnostic Plot: Alignment",figsize = (12,10))
    ax1 = fig.add_axes([0.06,0.07,0.91,0.36])
    ax2 = fig.add_axes([0.06,0.48,0.91,0.36])
    ax3 = fig.add_axes([0.06,0.90,0.91,0.05])

    #Axes
    ## Data ##
    ax2.imshow(data/np.max(data),aspect = 'auto',vmin = 0,vmax = 1)
    ax2.get_xaxis().set_visible(False)
    ax2.set_ylabel("Frequency Bin",fontsize = 12)
    ax2.set_title("t-Weighted (RFI subtracted)",fontsize = 16)

    ## finder bins ##
    ax1.imshow(weights[:,1:]/np.max(weights[:,1:]),aspect = 'auto',vmin = 0,vmax = 1)
    ax1.set_xlabel("Time Bin",fontsize = 12)
    ax1.set_ylabel("Frequency Bin",fontsize = 12)
    ax1.set_title("Finder Mask",fontsize = 16)


    ## Color Bar ##
    colbar = np.linspace(0,1,512).reshape(1,512)
    ax3.imshow(colbar,aspect = 'auto',vmin = 0,vmax = 1,extent = [0,1,0,1])
    ax3.get_yaxis().set_visible(False)
    ax3.set_xticks([0,1.0])


    #save file to png 
    filename = os.getcwd() + f"/data_alignment{suffix}.png"
    print("Saving diagnostic plot of alignment of correlated and HTR data to: "+filename+"...")

    plt.savefig(filename)


    AX = [ax1,ax2,ax3]

    return (fig, AX)

























if __name__ == "__main__":
    ##run main code##

    #get arguments
    args = get_args()


    #load in weighting masks
    wmask = wmask(args.weights)
    
    # # for testing
    # wmask.fw = np.ones(wmask.fw.size)
    # wmask.tmask[wmask.tmask != 0] = 1.0


    #check if weights and difx data are same size
    if wmask.tmask.shape != (args.fbins,args.tbins):
        print("NUMBER OF T-BINS: {:d}".format(args.tbins))
        print("NUMBER OF F-BINS: {:d}".format(args.fbins))
        print(wmask.tmask.shape)
        print("weights loaded from ["+args.weights+"] do not match input time and freqeuncy bin dimensions! Aborting...")
        sys.exit()


    #apply flagging
    f_flag(args,wmask)


    #apply weighting
    CAS = weight_scrunch(args,wmask)


    #save diagnostic plots
    fsort_idx = util.get_sort_index(args.freq)[::-1] # get index to sort CAS data


    CAS = CAS[:,fsort_idx]

    #TODO: save diagnostics plots to output files, maybe include option for name?

    if args.plot:
        #diagnostic plot for weighting 
        diagplot_weighting(CAS[0],CAS[1],CAS[2])


        #diagnostic plot for alignment
        diagplot_align(CAS[2],wmask.findermask,"_rfisub_weighted")
        diagplot_align(CAS[0],wmask.findermask,"_norifsub_noweights")
        diagplot_align(CAS[0],wmask.tmask,"_norfisub_tmask")



    print("apply_mf.py Completed successfully!")


    # END OF SCRIPT...
