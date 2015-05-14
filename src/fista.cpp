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

using std::vector;

#ifndef PI
#define PI (3.141592653589793)
#endif

#include "recon.hpp"
#include "params.hpp"

int recon::fista(params *params)
{
	int it,ix,iw,ntfft,nx1fft,nx2fft,nx3fft,nx4fft,nw,nk,if_low,if_high,padfactor,ix_no_pad,ix1,ix2,ix3,ix4,iter,*n;
	float *wd1,*in1,*out2,f_low,f_high,*amp;
	complex float czero,*freqslice, **tempModel1, **tempData, **model, **tempModel2, **tempModel3;

	fftwf_complex *out1,*in2;
	fftwf_plan p1,p2,p3,p4;
	czero = 0.0+0.0*I;

	// Find sampling of 5D observed data
	std::cout<< "Step 1" <<std::endl;
	findSampling(params);


	// Calculating threshold value
	params->thresh = params->lambda/(2*params->FistaAlpha);

	// Checking whether the percentage of data is greater than 3
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

	/* Length of frequency axis after considering symmetry ******////
	nw=ntfft/2+1;

	nk=nx1fft*nx2fft*nx3fft*nx4fft;

	/* Allocating memory ***/
	wd1 = allocPointerFloat(nk);
	params->D = alloc2DComPointer(nw,nk);
	tempModel1 = alloc2DComPointer(nw,nk);
	tempModel2 = alloc2DComPointer(nw,nk);
	tempModel3 = alloc2DComPointer(nw,nk);
	model = alloc2DComPointer(nw,nk);
	tempData = alloc2DComPointer(nw,nk);

	std::cout<< "Step 2" <<std::endl;

	// Observed data in f-x Domain ******************************////
	for (ix=0;ix<nk;ix++)
		for (iw=0;iw<nw;iw++)
			params->D[ix][iw] = czero;

	//// Temp data in f-x Domain ******************************////
	for (ix=0;ix<nk;ix++)
		for (iw=0;iw<nw;iw++)
			tempData[ix][iw] = czero;

	//// Temp model1 in f-w Domain *******************************////
	for (ix=0;ix<nk;ix++)
		for (iw=0;iw<nw;iw++)
			tempModel1[ix][iw] = czero;

	//// Temp model2 in f-w Domain *******************************////
	for (ix=0;ix<nk;ix++)
		for (iw=0;iw<nw;iw++)
			tempModel2[ix][iw] = czero;

	//// Temp model3 in f-w Domain *******************************////
	for (ix=0;ix<nk;ix++)
		for (iw=0;iw<nw;iw++)
			tempModel3[ix][iw] = czero;

	/// Spatial sampling with pad factor *************************////
	for (ix=0;ix<nk;ix++)
		wd1[ix] = 0.0;

	/// Convert data from tx to fx *******************************////
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

						wd1[ix] = params->wd[ix_no_pad];

						for (it=0; it<params->n1; it++)
							in1[it] = params->d[ix_no_pad][it];

						for (it=params->n1;it<ntfft;it++)
							in1[it] = 0.0;

						fftwf_execute(p1);

						for(iw=0;iw<nw;iw++)
							params->D[ix][iw] = out1[iw];
					}
				}
			}
		}
    } // End of tx to fx

	std::cout<< "Step 3" <<std::endl;

	///// Applying bandpass filter *************************////
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

	/// To store Fourier coefficients of 4D spatial data  //////////
	freqslice = allocComPointer(nk);
	n = allocPointerInt(4);
	n[0] = nx1fft;
	n[1] = nx2fft;
	n[2] = nx3fft;
	n[3] = nx4fft;

	p2 = fftwf_plan_dft(4, n, (fftwf_complex*)freqslice, (fftwf_complex*)freqslice, FFTW_FORWARD,  FFTW_ESTIMATE);
	p3 = fftwf_plan_dft(4, n, (fftwf_complex*)freqslice, (fftwf_complex*)freqslice, FFTW_BACKWARD, FFTW_ESTIMATE);

	float tmpt, temp, t=1.0;

	///***Start FISTA iterations *************************************************************////
	for (iter=0;iter<params->niter;iter++)
	{
		for(ix=0;ix<nk;ix++)
		  	for(iw=0;iw<nw;iw++)
		  		tempModel1[ix][iw] = model[ix][iw];

		//*** Convert from model(f-k) to data(f-x) *******************************************////
		forwardOperator(params, p3, nx1fft, nx2fft, nx3fft, nx4fft, nk, if_low, if_high, freqslice, tempModel2, tempData, wd1);

		for (ix=0;ix<nk;ix++)
			for (iw=if_low;iw<if_high;iw++)
					tempData[ix][iw] = (params->D[ix][iw] - tempData[ix][iw]);

		//*** Convert from data(f-x) to model(f-k) *******************************************////
		adjointOperator(params, p2, nx1fft, nx2fft, nx3fft, nx4fft, nk,if_low,if_high,freqslice,tempData,tempModel3,wd1);

		for (ix=0;ix<nk;ix++)
			for (iw=if_low;iw<if_high;iw++)
				tempModel3[ix][iw] = tempModel2[ix][iw]+(tempModel3[ix][iw]/params->FistaAlpha);

		/// *** Applying thresholding in f-k domain
		wthresh(if_low,if_high,nk,params->thresh, tempModel3, model);

		tmpt = t;
		t = (1+sqrtf(1+4*powf(t,2)))/2;
		temp = ((tmpt-1)/t);

		for (ix=0;ix<nk;ix++)
			for (iw=if_low;iw<if_high;iw++)
				tempModel2[ix][iw] = model[ix][iw]+ temp*(model[ix][iw] - tempModel1[ix][iw]);

		std::cout << "Iteration: " << iter+1 <<std::endl;

	  }   /// End of FISTA iteration

	  /// Make the fully sampled matrix
	  for (ix=0;ix<nk;ix++)
	  	  wd1[ix] = 1.0;

	  std::cout << "Done0" << std::endl;

	  //***Convert from model(f-k) to data(f-x) *******************************************////
	  forwardOperator(params, p3, nx1fft, nx2fft, nx3fft, nx4fft, nk, if_low, if_high, freqslice, model, tempData, wd1);

	  std::cout << "Done1" << std::endl;

	  in2 = (complex float *) fftwf_malloc(sizeof(fftwf_complex) * ntfft);

	  out2 = allocPointerFloat(ntfft);

	  p4 = fftwf_plan_dft_c2r_1d(ntfft, (fftwf_complex*)in2, out2, FFTW_ESTIMATE);

	  /// Data conversion from f-x to t-x ****************************////
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
						 {
							 in2[iw] = tempData[ix][iw];
						 }
						 fftwf_execute(p4);
						 for(it=0;it<params->n1;it++)
						 {
							 params->d[ix_no_pad][it] = out2[it]/ntfft;
						 }
					 }
				 }
			 }
		 }
	  }

	  std::cout << "Done2" << std::endl;

	  delete [] freqslice;
	  freqslice = NULL;
	  delete [] (wd1);
	  wd1 = NULL;
	  delete [] (in1);
	  in1 = NULL;
	  delete [] (out2);
	  out2 = NULL;
	  delete [] (out1);
	  out1 = NULL;
	  delete [] (in2);
	  in2 = NULL;

	  std::cout << "Done3" << std::endl;

	  deallocComPointer(nk, params->D);
	  deallocComPointer(nk, tempData);
	  deallocComPointer(nk, tempModel1);
	  deallocComPointer(nk, tempModel2);
	  deallocComPointer(nk, tempModel3);
	  deallocComPointer(nk, model);
	  fftwf_destroy_plan(p1);
	  fftwf_destroy_plan(p2);
	  fftwf_destroy_plan(p3);
	  fftwf_destroy_plan(p4);

	  return 1;
}

//****** Forward Operator using 4D FFT (Y = AX) ****************************************//
void recon::forwardOperator(params* params, fftwf_plan p3, int n1, int n2, int n3, int n4, int nk, int if_low, int if_high, complex float *freqslice, complex float **Model, complex float **data, float *wd1)
{
	int ix1,ix2,ix3,ix4,ix,iw;

	for (iw=if_low;iw<if_high;iw++)
	{

	  for (ix=0;ix<nk;ix++)
		  freqslice[ix] = Model[ix][iw];

	  fftwf_execute(p3);

	  for (ix1=0;ix1<n1;ix1++)
	  {
		  for (ix2=0;ix2<n2;ix2++)
		  {
			  for (ix3=0;ix3<n3;ix3++)
			  {
				  for (ix4=0;ix4<n4;ix4++)
				  {
					  ix = ix1*n2*n3*n4 + ix2*n3*n4 + ix3*n4 + ix4;
					  if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5)
						  data[ix][iw] = wd1[ix]*(freqslice[ix]/nk);
					  else
						  data[ix][iw] = 0.0+0.0*I;
				  }
			  }
		  }
	  }
   } /// [ End of Forward modelling: Y = AX ]

}

//****** Adjoint Operator using 4D iFFT (X = A^H Y) ****************************************//
void recon::adjointOperator(params* params, fftwf_plan p2, int n1, int n2, int n3, int n4, int nk, int if_low, int if_high, complex float *freqslice, complex float **Data, complex float **Model, float *wd1)
{

	int ix1,ix2,ix3,ix4,ix,iw;

	//// Considering frequencies within bandpass filter (Hz, not Wavenumber);
	//// Convert f-x(4D) to f-w ( X = A^H Y ) *******************////////////
	for (iw=if_low;iw<if_high;iw++)
	{
	  for (ix1=0;ix1<n1;ix1++)
	  {
		 for (ix2=0;ix2<n2;ix2++)
		 {
			for (ix3=0;ix3<n3;ix3++)
			{
				for (ix4=0;ix4<n4;ix4++)
				{
					ix = ix1*n2*n3*n4 + ix2*n3*n4 + ix3*n4 + ix4;
					if (ix1 < params->n2 && ix2 < params->n3 && ix3 < params->n4 && ix4 < params->n5)
					{
						freqslice[ix] = wd1[ix]*Data[ix][iw];
					}
					else
					{
						freqslice[ix] = 0.0+0.0*I;
					}
				}
			}
		 }
	  }

	  fftwf_execute(p2);
	  for (ix=0;ix<nk;ix++)
		  Model[ix][iw] = freqslice[ix];
	}/// End of Adjoint operator; X = A^H Y
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

//*** Thresholding function for complex number *******************//
inline void recon::wthresh(int lFreq, int hFreq, int n2,  float thresh, complex float **Model, complex float **threshModel)
{
    float temp;
    complex float sgn_data;

	for(int j=0; j<n2; j++)
	{
		for(int k=lFreq; k<hFreq; k++)
		{
			temp = cabsf(Model[j][k]) - thresh;
			temp = (temp+fabs(temp))/2;
			sgn_data = (Model[j][k]/cabsf(Model[j][k]));
			threshModel[j][k] = sgn_data*temp;
		}
	}
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
