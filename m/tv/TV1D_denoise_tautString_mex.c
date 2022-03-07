/*
 #  File            : TV1D_denoise_tautString_mex.c, copied from condat_fast_tv.c 
 #                       (based of v 2.0 Feb '13 of condat_fast_tv.c)
 #
 #  Version History : 1.0, Oct. 11, 2016
 #
 #  Author          :   mex interface: Stephen Becker, U. Colorado
                        Laurent Condat, PhD, CNRS research fellow in France.
 #
 #  Description     : This file contains a mex wrapper
 #                    to the implementation in the C language
 #                    of algorithms described in the research paper:
 #	
 #                    L. Condat, "A Direct Algorithm for 1D Total Variation
 #                    Denoising", preprint hal-00675043, 2012.
 #
 #                    This implementation comes with no warranty: due to the
 #                    limited number of tests performed, there may remain
 #                    bugs. In case the functions would not do what they are
 #                    supposed to do, please email the author (contact info
 #                    to be found on the web).
 #
 #                    If you use this code or parts of it for any purpose,
 #                    the author asks you to cite the paper above or, in 
 #                    that event, its published version. Please email him if 
 #                    the proposed algorithms were useful for one of your 
 #                    projects, or for any comment or suggestion.
 #
 #  Usage rights    : Copyright Laurent Condat.
 #                    This file is distributed under the terms of the CeCILL
 #                    licence (compatible with the GNU GPL), which can be
 #                    found at the URL "http://www.cecill.info".
 #
 #  This software is governed by the CeCILL license under French law and
 #  abiding by the rules of distribution of free software. You can  use,
 #  modify and or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL :
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
*/


#include <math.h>
#include <time.h>			
#include <stdio.h>
#include <stdlib.h>


#if defined(__GNUC__) && !(defined(__clang__)) && defined(NEEDS_UCHAR)
#include <uchar.h>
#endif
#include "mex.h"


#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define floatType double
/* 
This function implements 1D total variation denoising by
the taut string algorithm. It was adapted from the Matlab 
code written by Lutz Duembgen, which can be found at 
http://www.imsv.unibe.ch/content/staff/personalhomepages/
duembgen/software/multiscale_densities/index_eng.html
Note: numerical blow-up if floatType=float for N>=10^6
because the algorithm is based on the running sum of the 
signal values.
*/
void TV1D_denoise_tautstring(floatType* input, floatType* output, int width, const floatType lambda) {
	width++;
	int* index_low=mxCalloc(width,sizeof(int));
	floatType* slope_low=mxCalloc(width,sizeof(floatType));
	int* index_up=mxCalloc(width,sizeof(int));
	floatType* slope_up=mxCalloc(width,sizeof(floatType));
	int* index=mxCalloc(width,sizeof(int));
	floatType* z=mxCalloc(width,sizeof(floatType));
	floatType* y_low=mxMalloc(width*sizeof(floatType));
	floatType* y_up=mxMalloc(width*sizeof(floatType));
	int s_low=0;
	int c_low=0;
	int s_up=0;
	int c_up=0;
	int c=0;
	int i=2;
	y_low[0]=y_up[0]=0;
	y_low[1]=input[0]-lambda;
	y_up[1]=input[0]+lambda;
	for (;i<width;i++) {
		y_low[i]=y_low[i-1]+input[i-1];
		y_up[i]=y_up[i-1]+input[i-1];
	}
	y_low[width-1]+=lambda;
	y_up[width-1]-=lambda;
	slope_low[0] = INFINITY;
	slope_up[0] = -INFINITY;
	z[0]=y_low[0];
	for (i=1;i<width;i++) {
		index_low[++c_low] = index_up[++c_up] = i;
		slope_low[c_low] = y_low[i]-y_low[i-1];
		while ((c_low > s_low+1) && (slope_low[MAX(s_low,c_low-1)]<=slope_low[c_low])) {
			index_low[--c_low] = i;
			if (c_low > s_low+1) 
				slope_low[c_low] = (y_low[i]-y_low[index_low[c_low-1]]) /
					(i-index_low[c_low-1]);
			else
				slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c]);
		}
		slope_up[c_up] = y_up[i]-y_up[i-1];
		while ((c_up > s_up+1) && (slope_up[MAX(c_up-1,s_up)]>=slope_up[c_up])) {
			index_up[--c_up] = i;
			if (c_up > s_up+1)
				slope_up[c_up] = (y_up[i]-y_up[index_up[c_up-1]]) /
					(i-index_up[c_up-1]);
			else
				slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c]);
		}
		while ((c_low==s_low+1) && (c_up>s_up+1) && (slope_low[c_low]>=slope_up[s_up+1])) {
			index[++c] = index_up[++s_up];
			z[c] = y_up[index[c]];
			index_low[s_low] = index[c];
			slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c]);
		}
		while ((c_up==s_up+1) && (c_low>s_low+1) && (slope_up[c_up]<=slope_low[s_low+1])) {
			index[++c] = index_low[++s_low];
			z[c] = y_low[index[c]];
			index_up[s_up] = index[c];
			slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c]);
		}
	}
	for (i=1;i<=c_low-s_low;i++) 
		z[c+i]=y_low[index[c+i]=index_low[s_low+i]];
	c = c + c_low-s_low;
	int j=0;
	float a;
	i=1;
	while (i<=c) {
		a = (z[i]-z[i-1]) / (index[i]-index[i-1]);
		while (j<index[i]) {
			output[j] = a;
			j++;
		}
		i++;
	}
	mxFree(index_low); mxFree(slope_low);
	mxFree(index_up);  mxFree(slope_up);
	mxFree(index);     mxFree(z);
	mxFree(y_low);     mxFree(y_up);
}

void printUsage() {
    mexPrintf("z = TV1D_denoise_mex_tautString(y,lambda)\n");
    mexPrintf("  returns z = argmin_{x} .5||x-y||^2 + lambda ||x||_TV\n");
    mexPrintf("  using the C code of Laurent Condat and Lutz Duembgen\n");
    mexPrintf("  and the taut string algorithm\n");
}

#define IN_X 0
#define IN_LAMBDA 1
#define OUT_X 0
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] ) {
    floatType   *input, *output, lambda;
    mwSize      m, n; 
    int         width;

    if ( nrhs != 2 ) {
        printUsage();
        mexErrMsgIdAndTxt(  "MATLAB:mexFile:invalidNumInputs",
                "Two input arguments required.");
    } else if (nlhs != 1) {
        printUsage();
        mexErrMsgIdAndTxt( "MATLAB:mexFile:maxlhs",
                "Too many output arguments, needs 1 output.");
    }
    m           = mxGetM(prhs[IN_X]);
    n           = mxGetN(prhs[IN_X]);
    width       = (int)m*(int)n;
    plhs[OUT_X] = mxCreateDoubleMatrix( m, n, mxREAL);
    output      = (floatType *)mxGetPr( plhs[OUT_X] );
    input       = (floatType *)mxGetPr( prhs[IN_X] );
    lambda      = (floatType)mxGetScalar( prhs[IN_LAMBDA] );

    TV1D_denoise_tautstring(input,output,width,lambda); 
}
