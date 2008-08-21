#include "jProfiler.h"

#include "tensor.h"

#include <stdio.h>
#include <conio.h>

// number of dimensions to test with (handy for making sure something works and also for profiling)
#define nD 3

/*
	functions to make the metric with
*/

double x0(vector<double,nD>& p)
{
	return p[0]+p[1]+p[2]+p[3]+1;
}

double x1(vector<double,nD>& p)
{
	return p[0]*p[0]+p[3]*p[3]-2*p[1]*p[2]-1;
}

int main()
{	
	jProfiler profiler = jProfiler(L"<global>::main");
	// create a fairly curvy metric for testing stuff with
	tensor<double,nD> tl = Eta<double,nD>();
	
	tl[0][0] = x1;
	for(int i = 1; i < nD; i++)
	{
		tl[0][i] = tl[i][0] = x0;
	}
	
	printf("Covariant metric\n");
	tl.printDebug();
	
	// tensor<> u = Gamma<double,4>(t);
	tensor<double,nD> tu = RaiseMetricIndices(tl);
	printf("Contravariant metric\n");
	tu.printDebug();

	/*
		u.printDebug();
	
		u = u.raiseIndex(0,tu);
	
		u.printDebug();
	*/
	
	// really really slow for some reason, probably all the compositions of functions, and derivatives etc..
	tensor<double,nD> R = Riemann<double,nD>(tl);
	
	R.printDebug();
	
	while(!_getch()) continue;
	
	return 0;
}