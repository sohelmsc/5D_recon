/*
  * pocs5D.cpp
 *
 *  Created on: Mar 28, 2015
 *     Author: Mafijul Bhuiyan
 */

#include <cmath>
#include <algorithm>
using ::std::sort;
#include <stdlib.h>
#include <vector>
#include <iostream>
using std::cout;
using std::ofstream;
using std::endl;
#include <fstream>
#include <string>
#include <complex.h>
#include "fftw3.h"

#ifndef PI
#define PI (3.141592653589793)
#endif

#include "recon.hpp"
#include "params.hpp"

// d ----> is the 2D data with time X (nx1*nx2*nx3*nx4)
// wd_no_pad ---> Padding zeros in different dimensions
// fmax ------> Maximum frequency
// p = 1 for FISTA
// aa2 aa3 aa4 aa5 ---> are for anti aliasing
// debias ----> factor

int compare (const void * a, const void * b);

int recon::pocs(params *params)
{

  int it,ix,iw,ntfft,nx1fft,nx2fft,nx3fft,nx4fft,nw,nk,if_low,if_high,padfactor,ix_no_pad,ix1,ix2,ix3,ix4,iter,*n,nzero;
  float perc,perci,percf,*wd1,thres,*in1,*out2,f_low,f_high,*amp,*b1,*b2,median,bias,**M,**Mshift;
  fftwf_complex *out1,*in2;
  complex float **D,**Dobs,czero,*freqslice;
  fftwf_plan p1,p2,p3,p4;
  int ix_shift,ix1_shift,ix2_shift,ix3_shift,ix4_shift;
  perci = 1.0;
  percf = 0.0;

  czero = 0.0+0.0*I;
  findSampling(params);

  if (((float)params->sum_wd/(float)(params->n2*params->n3*params->n4*params->n5)) < 0.03)
  {
    return 0;
  }

  /* copy data from input to FFT array and pad with zeros */
  ntfft = params->padFactort*params->n1;
  nx1fft = params->padFactorx*params->n2;
  nx2fft = params->padFactorx*params->n3;
  nx3fft = params->padFactorx*params->n4;
  nx4fft = params->padFactorx*params->n5;
  if(params->n2==1) nx1fft = 1;
  if(params->n3==1) nx2fft = 1;
  if(params->n4==1) nx3fft = 1;
  if(params->n5==1) nx4fft = 1;
  nw=ntfft/2+1;
  nk=nx1fft*nx2fft*nx3fft*nx4fft;

  wd1 = allocPointerFloat(nk);
  D = alloc2DComPointer(nw,nk);
  Dobs = alloc2DComPointer(nw,nk);
  M = alloc2DPointerFloat(nw,nk);
  Mshift = alloc2DPointerFloat(nw,nk);

  for (ix=0;ix<nk;ix++)
	  for (iw=0;iw<nw;iw++)
		  D[ix][iw] = czero;

  for (ix=0;ix<nk;ix++)
	  for (iw=0;iw<nw;iw++)
		  Dobs[ix][iw] = czero;

  for (ix=0;ix<nk;ix++)
	  for (iw=0;iw<nw;iw++)
		  M[ix][iw] = 0.0;

  for (ix=0;ix<nk;ix++)
	  for (iw=0;iw<nw;iw++)
		  Mshift[ix][iw] = 0.0;

  in1 = allocPointerFloat(ntfft);

  out1 = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nw);
  p1 = fftwf_plan_dft_r2c_1d(ntfft, in1, out1, FFTW_ESTIMATE);

  for (ix1=0;ix1<nx1fft;ix1++)
  {
	  for (ix2=0;ix2<nx2fft;ix2++)
	  {
		  for (ix3=0;ix3<nx3fft;ix3++)
		  {
			  for (ix4=0;ix4<nx4fft;ix4++)
			  {
				if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5)
				{
				  ix_no_pad = ix1*params->n3*params->n4*params->n5 + ix2*params->n4*params->n5 + ix3*params->n5 + ix4;
				  ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
				  for (it=0; it<params->n1; it++)
					  in1[it] = params->d[ix_no_pad][it];
				  for (it=params->n1;it<ntfft;it++)
					  in1[it] = 0.0;

				  fftwf_execute(p1);

				  for(iw=0;iw<nw;iw++)
				  {
					  D[ix][iw] = out1[iw];
				  	  Dobs[ix][iw] = out1[iw];
				  }
				}
			  }
		  }
	  }
  } // End of tx to f-w

  for (ix=0;ix<nk;ix++)
	  wd1[ix] = 0.0;

  for (ix1=0;ix1<nx1fft;ix1++)
  {
	  for (ix2=0;ix2<nx2fft;ix2++)
	  {
		  for (ix3=0;ix3<nx3fft;ix3++)
		  {
			  for (ix4=0;ix4<nx4fft;ix4++)
			  {
				  if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5)
				  {
					  ix_no_pad = ix1*params->n3*params->n4*params->n5 + ix2*params->n4*params->n5 + ix3*params->n5 + ix4;
					  ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
					  wd1[ix] = params->wd[ix_no_pad];
				  }
			  }
		  }
	  }
  }

  f_low = 0.1;
  f_high = params->fmax;

  if(f_low>0)
	  if_low = trunc(f_low*params->d1*ntfft);
  else
	  if_low = 0;

  if(f_high*params->d1*ntfft<nw)
	  if_high = trunc(f_high*params->d1*ntfft);
  else
	  if_high = 0;

  freqslice = allocComPointer(nk);

  n = allocPointerInt(4); n[0] = nx1fft; n[1] = nx2fft; n[2] = nx3fft; n[3] = nx4fft;

  p2 = fftwf_plan_dft(4, n, (fftwf_complex*)freqslice, (fftwf_complex*)freqslice, FFTW_FORWARD, FFTW_ESTIMATE);
  p3 = fftwf_plan_dft(4, n, (fftwf_complex*)freqslice, (fftwf_complex*)freqslice, FFTW_BACKWARD, FFTW_ESTIMATE);

  amp = allocPointerFloat(nk*nw);

  for (iter=0;iter<params->niter;iter++)
  {

    for (iw=if_low;iw<if_high;iw++)
    {
      for (ix1=0;ix1<nx1fft;ix1++)
      {
    	 for (ix2=0;ix2<nx2fft;ix2++)
    	 {
    		for (ix3=0;ix3<nx3fft;ix3++)
    		{
    			for (ix4=0;ix4<nx4fft;ix4++)
    			{
    				ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
    				if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5)
    				{
    					freqslice[ix] = D[ix][iw];
    				}
    				else
    				{
    					freqslice[ix] = czero;
    				}
    			}
    		}
    	 }
      }
      fftwf_execute(p2);
      for (ix=0;ix<nk;ix++)
    	  D[ix][iw] = freqslice[ix];
    }
    for (ix=0;ix<nk;ix++)
    	for (iw=0;iw<nw;iw++)
    		M[ix][iw] = cabsf(D[ix][iw]);

    // radial smoothing to dealias the spectrum
    if (params->aa2>0 || params->aa3>0 || params->aa4>0 || params->aa5>0)
    {
		// map M to Mshift
		for (iw=if_low;iw<if_high;iw++)
		{
		  for (ix1=0;ix1<nx1fft;ix1++)
		  {
			  for (ix2=0;ix2<nx2fft;ix2++)
			  {
				  for (ix3=0;ix3<nx3fft;ix3++)
				  {
					  for (ix4=0;ix4<nx4fft;ix4++)
					  {
						  if (ix1 < (int) truncf(nx1fft/2)) ix1_shift = ix1 + (int) truncf(nx1fft/2);
						  else ix1_shift = ix1 - ((int) truncf(nx1fft/2));
						  if (ix2 < (int) truncf(nx2fft/2)) ix2_shift = ix2 + (int) truncf(nx2fft/2);
						  else ix2_shift = ix2 - ((int) truncf(nx2fft/2));
						  if (ix3 < (int) truncf(nx3fft/2)) ix3_shift = ix3 + (int) truncf(nx3fft/2);
						  else ix3_shift = ix3 - ((int) truncf(nx3fft/2));
						  if (ix4 < (int) truncf(nx4fft/2)) ix4_shift = ix4 + (int) truncf(nx4fft/2);
						  else ix4_shift = ix4 - ((int) truncf(nx4fft/2));
						  ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
						  ix_shift = ix1_shift*nx2fft*nx3fft*nx4fft + ix2_shift*nx3fft*nx4fft + ix3_shift*nx4fft + ix4_shift;
						  Mshift[ix_shift][iw] = M[ix][iw];
					  }
				  }
			  }
		 }
    }

    if (params->aa2>0)
    {
      radial_filter_gathers(Mshift,0,1,nw,0,1,nx1fft,0,1,nx2fft,0,1,nx3fft,0,1,nx4fft,2,params->aa2);
    }
    if (params->aa3>0)
    {
      radial_filter_gathers(Mshift,0,1,nw,0,1,nx1fft,0,1,nx2fft,0,1,nx3fft,0,1,nx4fft,3,params->aa3);
    }
    if (params->aa4>0)
    {
      radial_filter_gathers(Mshift,0,1,nw,0,1,nx1fft,0,1,nx2fft,0,1,nx3fft,0,1,nx4fft,4,params->aa4);
    }
    if (params->aa5>0)
    {
      radial_filter_gathers(Mshift,0,1,nw,0,1,nx1fft,0,1,nx2fft,0,1,nx3fft,0,1,nx4fft,5,params->aa5);
    }

    // map Mshift to M
    for (iw=if_low;iw<if_high;iw++)
    {
      for (ix1=0;ix1<nx1fft;ix1++)
      {
    	  for (ix2=0;ix2<nx2fft;ix2++)
    	  {
    		  for (ix3=0;ix3<nx3fft;ix3++)
    		  {
    			  for (ix4=0;ix4<nx4fft;ix4++)
    			  {
    				  if (ix1 < (int) truncf(nx1fft/2)) ix1_shift = ix1 + (int) truncf(nx1fft/2);
    				  else ix1_shift = ix1 - ((int) truncf(nx1fft/2));
    				  if (ix2 < (int) truncf(nx2fft/2)) ix2_shift = ix2 + (int) truncf(nx2fft/2);
    				  else ix2_shift = ix2 - ((int) truncf(nx2fft/2));
    				  if (ix3 < (int) truncf(nx3fft/2)) ix3_shift = ix3 + (int) truncf(nx3fft/2);
    				  else ix3_shift = ix3 - ((int) truncf(nx3fft/2));
    				  if (ix4 < (int) truncf(nx4fft/2)) ix4_shift = ix4 + (int) truncf(nx4fft/2);
    				  else ix4_shift = ix4 - ((int) truncf(nx4fft/2));
    				  ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
    				  ix_shift = ix1_shift*nx2fft*nx3fft*nx4fft + ix2_shift*nx3fft*nx4fft + ix3_shift*nx4fft + ix4_shift;
    				  M[ix][iw] = Mshift[ix_shift][iw];
    			  }
    		  }
    	  }
      }
    }
  }

	// sort non-zero amplitudes
	nzero = 0;
	for (iw=0;iw<nw;iw++){
	  for (ix=0;ix<nk;ix++){
		amp[ix*nw + iw] = M[ix][iw];
		if (amp[ix*nw + iw] > 0.000001) nzero++;
	  }
	}
	qsort(amp, nk*nw, sizeof(*amp),compare);

	perc = perci + iter*((percf-perci)/(params->niter-1));
	thres = amp[(int) truncf(nk*nw - 1 - (1-perc)*nzero)];
	for (iw=if_low;iw<if_high;iw++){
	  for (ix=0;ix<nk;ix++){
		if (M[ix][iw]<thres || cabsf(D[ix][iw])<thres) D[ix][iw] = czero;
		else D[ix][iw] = D[ix][iw]*(1 - powf(thres/(cabsf(D[ix][iw]) + 0.0000001),params->p));
	  }
	}

	for(iw=if_low;iw<if_high;iw++)
	{
	  for (ix=0;ix<nk;ix++)
		  freqslice[ix] = D[ix][iw];
	  fftwf_execute(p3);
	  for (ix1=0;ix1<nx1fft;ix1++){ for (ix2=0;ix2<nx2fft;ix2++){ for (ix3=0;ix3<nx3fft;ix3++){ for (ix4=0;ix4<nx4fft;ix4++)
	  {
		ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
		if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5){
		  D[ix][iw] = freqslice[ix]/nk;
		}
		else
		{
		  D[ix][iw] = czero;
		}
	  }
	 }
	}
	}
	}

    for (iw=if_low;iw<if_high;iw++)
    	for (ix=0;ix<nk;ix++)
    		D[ix][iw] = params->PocsAlpha*Dobs[ix][iw] + (1-params->PocsAlpha*wd1[ix])*D[ix][iw];

    std::cout << "Iteration: " << iter+1 <<std::endl;

  }   /// End of iteration


  delete [] (freqslice);

  in2 = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ntfft);

  out2 = allocPointerFloat(ntfft);

  p4 = fftwf_plan_dft_c2r_1d(ntfft, (fftwf_complex*)in2, out2, FFTW_ESTIMATE);

  /// Data conversion from f-w to t-x ****************************////
  for (ix1=0;ix1<nx1fft;ix1++)
  {
	 for (ix2=0;ix2<nx2fft;ix2++)
	 {
		 for (ix3=0;ix3<nx3fft;ix3++)
		 {
			 for (ix4=0;ix4<nx4fft;ix4++)
			 {
				 if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5)
				 {
					 ix_no_pad = ix1*params->n3*params->n4*params->n5 + ix2*params->n4*params->n5 + ix3*params->n5 + ix4;
					 ix = ix1*nx2fft*nx3fft*nx4fft + ix2*nx3fft*nx4fft + ix3*nx4fft + ix4;
					 for (iw=0; iw<ntfft; iw++)
						 in2[iw] = D[ix][iw];
						 fftwf_execute(p4);
						 for(it=0;it<params->n1;it++)
							 params->d[ix_no_pad][it] = out2[it]/ntfft;
				 }
			 }
		 }
	 }
  }

  //// Starting of the debaiasing, not required for FISTA Algorithm *********************************//
  if (params->debias)
  {
    // scale RMS amplitude to median value of observed traces
    b1 = allocPointerFloat(params->nx);
    b2 = allocPointerFloat(params->nx);
    nzero = 0;
    for (ix=0;ix<params->nx;ix++){
      b1[ix] = b2[ix] = 0.0;
      if (params->wd[ix]){
        for(it=0;it<params->n1;it++) b1[ix] += params->d[ix][it]*params->d[ix][it];
        nzero++;
      }
      for(it=0;it<params->n1;it++) b2[ix] += params->d[ix][it]*params->d[ix][it];
      b1[ix] = sqrtf(b1[ix]);
      b2[ix] = sqrtf(b2[ix]);
    }
    qsort (b1, params->nx, sizeof(*b1), compare);
    median = b1[(int) truncf(params->nx - 1 - 0.5*nzero)];
    for (ix=0;ix<params->nx;ix++)
    {
      if (params->wd[ix]==0 && b2[ix] > 1e-8) bias = median/b2[ix];
      else bias = 1.0;
      if (bias < 2.0) for(it=0;it<params->n1;it++) params->d[ix][it] *= bias;
      else            for(it=0;it<params->n1;it++) params->d[ix][it] *= 2.0;
    }
    // smooth RMS amplitude laterally
    for (ix=0;ix<params->nx;ix++)
    {
      b1[ix] = b2[ix] = 0.0;
      for(it=0;it<params->n1;it++) b1[ix] += params->d[ix][it]*params->d[ix][it];
      for(it=0;it<params->n1;it++) b2[ix] += params->d[ix][it]*params->d[ix][it];
      b1[ix] = sqrtf(b1[ix]);
      b2[ix] = sqrtf(b2[ix]);
    }

    int smooth1,smooth2,smooth3,smooth4;
    if (params->n2>50) smooth1 = 5;
    else        smooth1 = 0;
    if (params->n3>50) smooth2 = 5;
    else        smooth2 = 0;
    if (params->n4>50) smooth3 = 5;
    else        smooth3 = 0;
    if (params->n5>50) smooth4 = 5;
    else        smooth4 = 0;
    smooth4d(b1,wd1,params->n2,smooth1,params->n3,smooth2,params->n4,smooth3,params->n5,smooth4,5,params->verbose);
    for (ix=0;ix<params->nx;ix++)
    {
      if (params->wd[ix]==0 && b2[ix] > 1e-8) for(it=0;it<params->n1;it++) params->d[ix][it] *= b1[ix]/b2[ix];
    }
    delete [] (b1);
    b1 = NULL;
    delete [] (b2);
    b2 = NULL;
  }

  delete [] (wd1);
  wd1 = NULL;
  delete [] (amp);
  amp = NULL;
  delete [] (in1);
  in1 = NULL;
  delete [] (out2);
  out2 = NULL;
  delete [] (out1);
  out1 = NULL;
  delete [] (in2);
  in2 = NULL;

  deallocComPointer(nk, D);
  deallocComPointer(nk, Dobs);
  deallocPointer(nk, M);
  deallocPointer(nk, Mshift);
  fftwf_destroy_plan(p1);
  fftwf_destroy_plan(p2);
  fftwf_destroy_plan(p3);
  fftwf_destroy_plan(p4);

  return 1;

}

int compare (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}

void recon::smooth4d(float *d, float *wd,
              int nx1,int smooth1,
              int nx2,int smooth2,
              int nx3,int smooth3,
              int nx4,int smooth4,
              int nrepeat,bool verbose)
{
  float *din,*d_gather1,*d_gather2,*d_gather3,*d_gather4;
  int irepeat,ix,ix1,ix2,ix3,ix4;

  din = allocPointerFloat(nx1*nx2*nx3*nx4);
  d_gather1 = allocPointerFloat(nx1);
  d_gather2 = allocPointerFloat(nx2);
  d_gather3 = allocPointerFloat(nx3);
  d_gather4 = allocPointerFloat(nx4);

  for (ix=0;ix<nx1*nx2*nx3*nx4;ix++) din[ix] =  d[ix];
  for (irepeat=0;irepeat<nrepeat;irepeat++){
  if (smooth1>1){
    for (ix2=0;ix2<nx2;ix2++){ for (ix3=0;ix3<nx3;ix3++){ for (ix4=0;ix4<nx4;ix4++){
      for (ix1=0;ix1<nx1;ix1++) d_gather1[ix1] = d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1];
      mean_filter(d_gather1,nx1,smooth1,1);
      for (ix1=0;ix1<nx1;ix1++) d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1] = d_gather1[ix1];
    }}}
  }
  if (smooth2>1){
    for (ix1=0;ix1<nx1;ix1++){ for (ix3=0;ix3<nx3;ix3++){ for (ix4=0;ix4<nx4;ix4++){
      for (ix2=0;ix2<nx2;ix2++) d_gather2[ix2] = d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1];
      mean_filter(d_gather2,nx2,smooth2,1);
      for (ix2=0;ix2<nx2;ix2++) d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1] = d_gather2[ix2];
    }}}
  }
  if (smooth3>1){
    for (ix1=0;ix1<nx1;ix1++){ for (ix2=0;ix2<nx2;ix2++){ for (ix4=0;ix4<nx4;ix4++){
      for (ix3=0;ix3<nx3;ix3++) d_gather3[ix3] = d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1];
      mean_filter(d_gather3,nx3,smooth3,1);
      for (ix3=0;ix3<nx3;ix3++) d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1] = d_gather3[ix3];
    }}}
  }
  if (smooth4>1){
    for (ix1=0;ix1<nx1;ix1++){ for (ix2=0;ix2<nx2;ix2++){ for (ix3=0;ix3<nx3;ix3++){
      for (ix4=0;ix4<nx4;ix4++) d_gather4[ix4] = d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1];
      mean_filter(d_gather4,nx4,smooth4,1);
      for (ix4=0;ix4<nx4;ix4++) d[ix4*nx3*nx2*nx1 + ix3*nx2*nx1 + ix2*nx1 + ix1] = d_gather4[ix4];
    }}}
  }
  for (ix=0;ix<nx1*nx2*nx3*nx4;ix++) d[ix] =  wd[ix]*din[ix] + (1-wd[ix])*d[ix];
  }

  delete [] (din);
  delete [] (d_gather1);
  delete [] (d_gather2);
  delete [] (d_gather3);
  delete [] (d_gather4);

  return;
}

void recon::mean_filter(float *trace, int nt, int ntw, int nrepeat)
/*< mean filter in 1d. ntw is window length and nrepeat allows to repeat the smoothing >*/
{
  int irepeat,it,itw,index1;
  float sum,nsum;
  float *tracein;
  tracein = allocPointerFloat(nt);
  /* calculate mean value of a window and assign it to the central index */
    for (irepeat=0;irepeat<nrepeat;irepeat++){
      for (it=0;it<nt;it++) tracein[it] = trace[it];
      for (it=0;it<nt;it++){
        sum  = 0.0;
        nsum = 0.0;
        for (itw=0;itw<ntw;itw++){
          index1 = it - (int) truncf(ntw/2) + itw;
          if (index1>=0 && index1<nt){
            sum += tracein[index1];
            nsum += 1.0;
          }
        }
        trace[it] = sum/nsum;
      }
    }
  delete [] (tracein);
  return;
}

void recon::radial_filter_gathers(float **d,
                           float o1,float d1,int n1,
                           float o2,float d2,int n2,
                           float o3,float d3,int n3,
                           float o4,float d4,int n4,
                           float o5,float d5,int n5,
                           int axis,int L)
{
  float **d_gather;
  int i1,i2,i3,i4,i5;
  // process gathers
  if (axis==2){
    d_gather = alloc2DPointerFloat(n1,n2);
    for (i3=0;i3<n3;i3++){ for (i4=0;i4<n4;i4++){ for (i5=0;i5<n5;i5++){
      for (i2=0;i2<n2;i2++) for (i1=0;i1<n1;i1++) d_gather[i2][i1] = d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1];
      radial_filter(d_gather,n1,n2,L);
      for (i2=0;i2<n2;i2++) for (i1=0;i1<n1;i1++) d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1] = d_gather[i2][i1];
    }}}
    deallocPointer(n2, d_gather);
  }
  else if (axis==3){
    d_gather = alloc2DPointerFloat(n1,n3);
    for (i2=0;i2<n2;i2++){ for (i4=0;i4<n4;i4++){ for (i5=0;i5<n5;i5++){
      for (i3=0;i3<n3;i3++) for (i1=0;i1<n1;i1++) d_gather[i3][i1] = d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1];
      radial_filter(d_gather,n1,n3,L);
      for (i3=0;i3<n3;i3++) for (i1=0;i1<n1;i1++) d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1] = d_gather[i3][i1];
    }}}
    deallocPointer(n3, d_gather);
  }
  else if (axis==4){
    d_gather = alloc2DPointerFloat(n1,n4);
    for (i2=0;i2<n2;i2++){ for (i3=0;i3<n3;i3++){ for (i5=0;i5<n5;i5++){
      for (i4=0;i4<n4;i4++) for (i1=0;i1<n1;i1++) d_gather[i4][i1] = d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1];
      radial_filter(d_gather,n1,n4,L);
      for (i4=0;i4<n4;i4++) for (i1=0;i1<n1;i1++) d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1] = d_gather[i4][i1];
    }}}
    deallocPointer(n4, d_gather);
  }
  else if (axis==5){
    d_gather = alloc2DPointerFloat(n1,n5);
    for (i2=0;i2<n2;i2++){ for (i3=0;i3<n3;i3++){ for (i4=0;i4<n4;i4++){
      for (i5=0;i5<n5;i5++) for (i1=0;i1<n1;i1++) d_gather[i5][i1] = d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1];
      radial_filter(d_gather,n1,n5,L);
      for (i5=0;i5<n5;i5++) for (i1=0;i1<n1;i1++) d[i5*n4*n3*n2 + i4*n3*n2 + i3*n2 + i2][i1] = d_gather[i5][i1];
    }}}
    deallocPointer(n5, d_gather);
  }

  return;
}

void recon::radial_filter(float **d, int nt, int nx, int nL)
{
  int iL,it0,ix0,it,ix;
  float dx,dt,oL,dL,L,x0,t0,x,t,theta,t1,t2,x1,x2,w1,w2,w3,w4;
  float **dout;

  dt = 1/(float) nt;
  dx = 2/(float) nx;
  dL = dt/4;//(dt < dx) ? dt : dx;
  oL = -dL*(((float) nL-1)/2);
  dout = alloc2DPointerFloat(nt,nx);
  for (ix0=0;ix0<nx;ix0++){
    for (it0=0;it0<nt;it0++){
      dout[ix0][it0] = 0.0;
      x0 = ix0*dx - 1;
      t0 = it0*dt;
      if (t0>0) theta = atanf(fabsf(x0)/t0);
      else theta = 90;
      for (iL=0;iL<nL;iL++){
        L = iL*dL + oL;
        //bilinear interpolation to get 1 value from 4 surrounding points.
        t = t0 + L*cosf(theta);
        it = (int) truncf(t/dt);
        t1 = it*dt;
        t2 = t1 + dt;
        x = x0 + signf(x0)*L*sinf(theta);
        ix = (int) truncf((x + 1)/dx);
        x1 = ix*dx - 1;
        x2 = x1 + dx;
        w1 = (t2-t)*(x2-x)/((t2-t1)*(x2-x1));
        w2 = (t-t1)*(x2-x)/((t2-t1)*(x2-x1));
        w3 = (t2-t)*(x-x1)/((t2-t1)*(x2-x1));
        w4 = (t-t1)*(x-x1)/((t2-t1)*(x2-x1));
        if (it >= 0 && it+1 < nt && ix >= 0 && ix+1 < nx && w1 >= 0 && w2 >= 0 && w3 >= 0 && w4 >= 0){
          dout[ix0][it0] += (w1*d[ix][it] + w2*d[ix][it+1] + w3*d[ix+1][it] + w4*d[ix+1][it+1])/nL;
        }
        //MARK
        //dout[ix0][it0] += d[ix][it]/nL;
        //MARK
      }
    }
  }
  for (ix=0;ix<nx;ix++) for (it=0;it<nt;it++) d[ix][it] = dout[ix][it];
  deallocPointer(nx, dout);
  return;
}

void recon::findSampling(params* params)
{
    params->nx = params->n2*params->n3*params->n4*params->n5;

    // Allocate 2D float data in the memory, n1 is the fastest dimension
    params->d = alloc2DPointerFloat(params->n1,params->nx);

    // Allocate 1D float data in the memory, n1 is the fastest dimension
    params->wd = allocPointerFloat(params->nx);

    // For finding the sampling of the data
    for (params->i2=0; params->i2<params->nx; params->i2++)
    {
      for (params->i1=0; params->i1<params->n1; params->i1++)
    	  params->d[params->i2][params->i1] = 0;
      params->wd[params->i2] = 0;
    }

    params->sum_wd = 0;
    for (params->ix=0; params->ix<params->nx; params->ix++)
    {
      params->sum = 0;
      for (params->i1=0; params->i1<params->n1; params->i1++)
      {
    	  params->sum += params->data[(params->ix*params->n1)+params->i1]*params->data[(params->ix*params->n1)+params->i1];
    	  params->d[params->ix][params->i1] = params->data[(params->ix*params->n1)+params->i1];
      }
      if (params->sum)
      {
        params->wd[params->ix] = 1;
        params->sum_wd++;
      }
    }
    delete [] params->data;
}

inline float recon::signf(float a)
/*< sign of a float >*/
{
 float b;
 if (a>0)      b = 1.0;
 else if (a<0) b =-1.0;
 else          b = 0.0;
 return b;
}

/***** Allocate 2D pointer ***************************************/
inline float ** recon::alloc2DPointerFloat(int n1, int n2)
{
	float **data;
	if (!(data = new float *[n2])) {
		std::cout << " Error out of memory " << std::endl;
		exit(1);
	} else {
		for (int i = 0; i < n2; i++)
			data[i] = new float[n1];
	}
	return data;
}

/***** Allocate 1D pointer ***************************************/
inline float * recon::allocPointerFloat(int n1)
{
	float *data;
	data = new float[n1];
	return data;
}

/***** Allocate 1D pointer ***************************************/
inline int * recon::allocPointerInt(int n1)
{
	int *data;
	data = new int[n1];
	return data;
}

inline float complex ** recon::alloc2DComPointer(int n1, int n2)
{
	complex float **data;
	data = (float complex **) fftwf_malloc(sizeof(float complex *) * n2);

	for(int i=0; i<n2; i++)
		data[i] = (float complex *) fftwf_malloc(sizeof(float complex) * n1);

	return data;
}

inline float complex * recon::allocComPointer(int n1)
{
	complex float *data;
	data = (float complex *) fftwf_malloc(sizeof(float complex) * n1);
	return data;
}

void recon::deallocPointer(int n1, float **data)
{
	int i;
	for (i = 0; i < n1; i++)
	{
		delete [] data[i];
		data[i] = NULL;
	}
	delete [] data;
	data = NULL;
}

void recon::deallocComPointer(int n1, complex float **data)
{
	int i;
	for (i = 0; i < n1; i++)
	{
		delete [] data[i];
		data[i] = NULL;
	}
	delete [] data;
	data = NULL;
}


