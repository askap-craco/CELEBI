#imports
import glob,sys,shutil,os
import numpy as np

#append / to dir if needed
def fix_dir(dir):
    if not dir[-1] == "/":
        return dir + "/"
    return dir

#construct fullfile
def fullfile(dir,filename):
    return fix_dir(dir) + filename



    

#function to get input files
def get_inputfiles(dir,suf="D2D.input"):
    cfinputs = glob.glob(fix_dir(dir) + "c*_f*/*" + suf)
    return sorted(cfinputs)

#function to create output difx files and input files
def create_difx_output_dir(idir,ifiles,isuf,osuf):
    #create output .input files
    ofiles = []
    isuf_len = len(isuf)
    for i,inputfile in enumerate(ifiles):
        ofiles.append(inputfile[:-isuf_len] + osuf)

    #copy file to new directory from old directory
    for i,outputfile in enumerate(ofiles):
        odifx = ofiles[i][:-5] + "difx"
        #check if dir exists
        if not os.path.isdir(odifx):
            os.mkdir(odifx)

        with open(ifiles[i],'r') as fin:
            lines = fin.readlines()
            for j,line in enumerate(lines):
                if "OUTPUT FILENAME" in line:
                    lines[j] = "OUTPUT FILENAME:    " + odifx + "\n"

                    #read lines to a new file
                    with open(ofiles[i],'w') as fout:
                        fout.writelines(lines)
                        fout.close()

                    break

    #create output difx files
    return ofiles

def create_difx_output_dir2(idir,ifiles,isuf,osuf,odir):
    #create output .input files
    ofiles = []
    isuf_len = len(isuf)
    for i,inputfile in enumerate(ifiles):
        ofiles.append(fullfile(odir,fullfile(os.path.basename(os.path.dirname(inputfile)),os.path.basename(inputfile)))[:-isuf_len] + osuf)

    #copy file to new directory from old directory
    for i,outputfile in enumerate(ofiles):
        odifx = ofiles[i][:-5] + "difx"
        #check if dir exists
        if not os.path.isdir(os.path.dirname(odifx)):
            os.mkdir(os.path.dirname(odifx))

        if not os.path.isdir(odifx):
            os.mkdir(odifx)

        with open(ifiles[i],'r') as fin:
            lines = fin.readlines()
            for j,line in enumerate(lines):
                if "OUTPUT FILENAME" in line:
                    lines[j] = "OUTPUT FILENAME:    " + odifx + "\n"

                    #read lines to a new file
                    with open(ofiles[i],'w') as fout:
                        fout.writelines(lines)
                        fout.close()

                    break

    #create output difx files
    return ofiles




#find sorted freq index
def get_inverse_sort_index(arr):
    #this function IS NOT A SIMPLE SORT IN ORDER, the idea of this function 
    #is that we start with a 
    idxarr = np.zeros(len(arr),dtype = int)

    #sort
    sortarr = np.sort(arr)

    #TODO - not really important, but could do this faster... probably
    for i,f in enumerate(arr):
        idxarr[i] = np.where(sortarr == f)[0][0]
    
    return idxarr


def get_sort_index(arr):

    return sorted(range(len(arr)), key = lambda k: arr[k])


