#import <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <string.h>

/*
CODE TO EXECUTE IN TERMINAL TO COMPILE FILE (move to /fred/oz002/tdial/dlib (or wherever the file is))

gcc -shared -Wl,-soname,cDiFX -o cDiFX.so -fPIC cDiFX.c

*/


//inverse factorial
int invfact(int num,int stop){
    int sum = 0;

    if (stop < 1) {
        return sum;
    }

    sum += num;

    for (int i = 1; i < stop; i++) {

        sum += num - i;

    }

    return sum;
}


int byte2vis(unsigned char* buffer,unsigned char* vis,unsigned char* header,int* freqs,int* vischan,int* zffind,int* var) {

    //initialise variables to use
    int freq_i = 0;//freq index
    int chan_i = 0;//num of channels
    int nvis = 0;//num of visibility's
    int count = 0;//buffer counter
    unsigned char* tp = vis;//temp vis pointer
    unsigned char* hp = header;//temp header pointer

    //unpack var
    int bytelen = var[0]; //length of vis buffer in bytes
    int baseoffset = var[1]; //baseline offset
    int fnum = var[2]; //number of active bands
    int tnum = var[5]; //number of telescope
    int max_chan = var[3]; //max num of channels
    int max_vis = var[4]; //max number of visibility's
    int blocknum = 0; // referes to the T1 baseline block
    int baseline_i = 0; //baseline of ith vis
    int pol_i = 0; //pol of ith vis
    int blocksize = tnum*fnum*4;
    int visINDEX = 0; //current vis index
    int visINDEXp = 0; //previous vis index

    //reset arrays
    memset(header,0,74*var[4]);
    memset(vis,0,8*var[3]*var[4]);


    //create block indexes
    int *block_n = malloc(tnum*4);
    for (int i = 0; i < tnum; i++) {
        block_n[i] = invfact(tnum,i) * fnum * 4;
    }
    int block_i = 0;


    while (count < bytelen) {

        //get freq index from header infomation
        buffer += 8;
        baseline_i = *(int*)buffer;
        buffer += 24;
        freq_i = *(int*)buffer;
        
        chan_i = freqs[freq_i]; // get number of channels in visibility
        buffer += 4;
        
        //get pol_i
        switch (*(short int*)buffer) {

            case 22616: pol_i = 0; break; //[X,X]
            case 22872: pol_i = 1; break; //[X,Y]
            case 22617: pol_i = 2; break; //[Y,X]
            case 22873: pol_i = 3; break; //[Y,Y]
        }
        buffer -= 36;
        //get visibility index
        block_i = (baseline_i - baseoffset)/256;
        

        visINDEX = block_n[block_i] + fnum*pol_i + zffind[freq_i] + fnum*4*(baseline_i 
                    - 256 * (block_i + 1) - 1 - block_i);

        //copy header data
        hp += 74*(visINDEX-visINDEXp);
        memmove(hp,buffer,74);
        buffer += 74;
        
        //set channel length to chan array
        vischan[visINDEX] = chan_i;

        //copy visibility data
        tp += 8*max_chan*(visINDEX-visINDEXp);
        memmove(tp,buffer,8*chan_i);
        
        buffer += 8*chan_i;
        count += 74 + 8*chan_i;

        visINDEXp = visINDEX;
        
        nvis++;
    }
    
    return nvis;
}


//function to write to single buffer
void vis2byte(unsigned char* buffer,int bytelen,unsigned char* vis,int max_chan,unsigned char* header,int* chans) {


    int count = 0; //while loop stopping condition
    unsigned char* tp = buffer; 
    int i = 0;
    int chan_i = 0;

    while (count < bytelen) {

        chan_i = chans[i];

        //move header memory
        memmove(tp,header,74);
        tp +=74;
        header += 74;

        //move visibility memory
        memmove(tp,vis,8*chan_i);
        tp += 8*chan_i;
        vis += 8*max_chan;

        count += 74 + 8*chan_i;


        i++;
    }


    return;


}



//Function to get Polarisation, baselines and antenna from fnum and tnum
void get_baselines(int tnum,int* baselines) {

    // loop over tnum
    int base_ind = 0;
    for (int i = 0; i < tnum; i++) {

        for (int j = i; j < tnum; j++) {

            baselines[base_ind] = (i+1)*256 + (j+1);
            base_ind++;

        }

    }

    return;

}