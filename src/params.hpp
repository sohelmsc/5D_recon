/*
  * params.hpp
 *
 *  Created on: Mar 28, 2015
 *      Author: Mafijul Bhuiyan
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <iostream>
#include <complex.h>
using std::cout;
using std::ostream;
using std::endl;

enum reconType{POCS, FISTA, MWNI, SVD, HOSVD}; /// Different reconstruction types

class params{

public:
	// These parameters depend on the seismic data and velocity model *******************/
	complex float **D;
	float *data;									/// 2D data
	float *wd;
	float *trace;
	float **d;
	int ix;
	int nx;
	int n1;
	int n2;
	int n3;
	int n4;
	int n5;
	int i1;
	int i2;
	int padFactort;
	int padFactorx;
	float d1;
	float d2;
	float d3;
	float d4;
	float d5;
	float o1;
	float o2;
	float o3;
	float o4;
	float o5;
	float thresh;
	float lambdaTemp[5];
	float lambda;
	float sum;
	int niter;
	float PocsAlpha;
	float FistaAlpha;
	float fmax;
	int sum_wd;
	bool verbose;
	float p;
	bool debias;
	int aa2;
	int aa3;
	int aa4;
	int aa5;
	std::string inFileName1,inFileName2, outFileName;
	reconType recType;

	void reset()
	{
		delete [] wd;
	}


};



#endif /* PARAMS_HPP_ */
