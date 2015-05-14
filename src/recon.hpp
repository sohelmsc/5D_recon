/*
  * pocs5D.hpp
 *
 *  Created on: Mar 28, 2015
 *      Author: Mafijul Bhuiyan
 */

#ifndef RECON_HPP_
#define RECON_HPP_
#include "params.hpp"
#include <complex.h>
#include <vector>
#include "fftw3.h"

using std::vector;

class recon{

public:
	int pocs(params *params);
	void deallocPointer(int n1, float **data);
	int fista(params *params);

//// Private: to make more secured ***************************************/
private:

fftwf_plan p3, p2;

void deallocComPointer(int n1, complex float **data);

inline void wthresh(int lFreq, int hFreq, int n2,  float thresh, complex float **Model, complex float **threshModel);

void adjointOperator(params* params, fftwf_plan p2, int n1, int n2, int n3, int n4, int nk, int if_low, int if_high, complex float *freqslice, complex float **Data, complex float **Model, float *wd1);

void forwardOperator(params* params, fftwf_plan p3, int n1, int n2, int n3, int n4, int nk, int if_low, int if_high, complex float *freqslice, complex float **Model, complex float **data, float *wd1);

void smooth4d(float *d, float *wd,
              int nx1,int smooth1,
              int nx2,int smooth2,
              int nx3,int smooth3,
              int nx4,int smooth4,
              int nrepeat,bool verbose);

void mean_filter(float *trace, int nt, int ntw, int nrepeat);

void radial_filter_gathers(float **d,
                           float o1,float d1,int n1,
                           float o2,float d2,int n2,
                           float o3,float d3,int n3,
                           float o4,float d4,int n4,
                           float o5,float d5,int n5,
                           int axis,int L);

void radial_filter(float **d, int nt, int nx, int nL);

float signf(float a);

float ** alloc2DPointerFloat(int n1, int n2);
inline float * allocPointerFloat(int n1);
inline int * allocPointerInt(int n1);
inline float complex * allocComPointer(int n1);
float complex ** alloc2DComPointer(int n1, int n2);
void findSampling(params* params);
};

#endif /* POCS5D_HPP_ */



