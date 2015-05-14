  /* recon.cpp
 *
 *  Created on: Mar 28, 2015
 *      Author: Mafijul Bhuiyan
 */

#include "params.hpp"
#include "recon.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char* argv[])
{
	clock_t start, finish;
	long int i,numDim=10, fsize;
	FILE *fp1, *fp2, *fp3;
	int count=0, out=0, initPatchNum, finalPatchNum;
	float *data, *dim;
	std::stringstream ss;

    /** Creating the objects of the dependence classes **********************************/
    params params;
    recon rec;

    /** Assigning values for the parameters *********************************************/
    params.inFileName1 = "/home/entropy/workspace/pocs5D/src/short_patch_98.bin";
    params.outFileName = "/home/entropy/workspace/pocs5D/src/Recon_";

    params.recType = POCS;
    params.d1 = 0.001;
    params.o1 = 0.0;
    params.d2 = 1;
    params.o2 = 0.0;
    params.d3 = 1;
    params.o3 = 0.0;
    params.d4 = 1;
    params.o4 = 0.0;
    params.d5 = 1;
    params.o5 = 0.0;
    params.niter = 220;
    params.PocsAlpha = 1.0;
    params.fmax = 180;
    if (params.fmax > 0.5/params.d1) params.fmax = 0.5/params.d1;
    params.sum_wd = 0;
    params.verbose = true;
    params.p = 1.0; //1 ==> soft thresholding, >1 ==> hard thresholding
    params.debias = false;
    params.aa2 = 0;
    params.aa3 = 0;
    params.aa4 = 0;
    params.aa5 = 0;
    params.padFactort = 2;
    params.padFactorx = 2;
    params.FistaAlpha = 1.1; //// FistaAlpha >= max(eig(H'*H)); Lipschtiz constant

    params.lambdaTemp[0] = 4.1;
    params.lambdaTemp[1] = 4.0;
    params.lambdaTemp[2] = 4.3;
    params.lambdaTemp[3] = 4.6;
    params.lambdaTemp[4] = 5.8;

    //  Read the dimension of all the patches ******************//
    //  fp1 = fopen(params.inFileName1.c_str(), "rb");
    //  dim = (float *) malloc(sizeof(float) * 960);
    //  fread(dim, sizeof(float), 960, fp1);
    //  fclose(fp1);

    initPatchNum = -1;
    finalPatchNum = 0;

    for(i=initPatchNum; i<finalPatchNum; i++)
	{
        params.inFileName2 = "/home/entropy/workspace/pocs5D/src/small_patch.bin";
        params.outFileName = "/home/entropy/workspace/pocs5D/src/small_Recon_patch_";
        params.lambda = params.lambdaTemp[i+1];

        ss << i+1;
    	std::string str = ss.str();
//    	params.inFileName2 = params.inFileName2 + str;
//    	params.inFileName2 = params.inFileName2 + std::string(".bin");
    	fp2 = fopen(params.inFileName2.c_str(), "rb");

    	params.outFileName = params.outFileName + str;
    	params.outFileName = params.outFileName + std::string(".bin");

		///* dimension of each data patch *************/
//		params.n1 = dim[(i+1)*5-1]; /// Time Axis
//		params.n2 = dim[(i+1)*5-5]; /// CMPX
//		params.n3 = dim[(i+1)*5-4]; /// CMPY
//		params.n4 = dim[(i+1)*5-3]; /// OFFSET
//		params.n5 = dim[(i+1)*5-2]; /// AZIMUTH
		///* dimension of data patch *************/

	    params.n1 = 301;
	    params.n2 = 5;
	    params.n3 = 5;
	    params.n4 = 10;
	    params.n5 = 10;

		fsize =params.n1*params.n2*params.n3*params.n4*params.n5;

		params.nx = params.n2*params.n3*params.n4*params.n5;

		params.data = (float *) malloc(sizeof(float) * fsize);

		fread(params.data,sizeof(float), fsize, fp2);

		start = clock();

		if(params.recType == POCS){out = rec.pocs(&params);} // Reconstructing data using POCS in f-x domain

		else if(params.recType == FISTA){out = rec.fista(&params);} // Reconstructing data using FISTA in f-x domain

		else if(params.recType == MWNI){/*rec.mwni(&params);*/} // Not implemented yet

		finish = clock();

		std::cout << "Time: " << (finish-start)/double(CLOCKS_PER_SEC) << " Seconds " <<std::endl;

		std::cout << "Size: " << fsize << " Seconds " <<std::endl;

		if(out!=0)
		{
			std::cout<< "Reconstruction completed at: " << i+1 << std::endl;

			params.data = (float *) malloc(sizeof(float) * fsize);

			//// Writing the file in the output file
			for (params.ix=0; params.ix<params.nx; params.ix++)
			   for (params.i1=0; params.i1<params.n1; params.i1++)
				  params.data[(params.ix*params.n1)+params.i1] = params.d[params.ix][params.i1];

			//// Write a binary file in the out
			fp3 = fopen(params.outFileName.c_str(), "wb");
			fwrite(params.data, sizeof(float), fsize, fp3);

			rec.deallocPointer(params.nx, params.d);
		}
		else
			std::cout<< "Nothing Done at: " << i+1 << std::endl;

		fclose(fp2);
		fclose(fp3);
		ss.clear();
    	ss.str(std::string());
	}

    params.reset();

  return 0;
}
