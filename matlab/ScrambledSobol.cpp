#include <iostream>
#include <cassert>

#include "MatousekAffineOwen.h"

int main(int argc, char *argv[])
{

	MatousekAffineOwenInitialize();

	mxArray *d_a = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(d_a) = (double)atoi(argv[1]);

	mxArray *skip_a = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(skip_a) = (double)atoi(argv[2]);

	mxArray *npoints_a = mxCreateDoubleMatrix(1, 1, mxREAL);
	*mxGetPr(npoints_a) = 5.0;

	mxArray *sobol_set, *sobol_points;


	/*
	mlfMatousekAffineOwen(int nargout, mxArray** x, mxArray* dim, mxArray* skip, mxArray* Npoints);
	*/

	assert(mlfMatousekAffineOwen(1, &sobol_points, d_a, skip_a, npoints_a));
	std::cout << sobol_points << std::endl;
	std::cout << mxGetM(sobol_points) <<"x" <<mxGetN(sobol_points) << std::endl;

	MatousekAffineOwenTerminate();


}
