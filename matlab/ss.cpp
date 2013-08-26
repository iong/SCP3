#include <iostream>

#include "MatousekAffineOwen.h"

int main(int argc, char *argv[])
{

	double d = atoi(argv[1]);
	mxArray *d_a = mxCreateDoubleScalar(d);
	mxArray *skip_a = mxCreateDoubleScalar((double)(1<<atoi(argv[2])));
	mxArray *npoints_a = mxCreateDoubleScalar(5.0);

	mxArray *sobol_set, *sobol_points;

	MatousekAffineOwenInitialize();
	mlfScrambledSobolset(1, &sobol_set, d_a, skip_a);
	mlfGetPoints(1, &sobol_points, sobol_set, npoints_a);
	MatousekAffineOwenTerminate();

	std::cout << mxGetM(sobol_points) <<"x"
		<<mxGetN(sobol_points) << std::endl;



}
