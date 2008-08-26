#include "foo\tensor.h"
#include <stdio.h>

int main()
{
	Function<4> F = 1;
	Function<4> t = Function<4>(Vector<4>::T);
	Function<4> x = Function<4>(Vector<4>::X);
	Function<4> y = Function<4>(Vector<4>::Y);
	Function<4> z = Function<4>(Vector<4>::Z);
	Function<4> G = t*t + x*x + y*y + z*z;

	Vector<4> pos1 = 1.0;
	Vector<4> pos2 = sqrt(2) * 0.5;


	printf("G = ");
	G.DebugPrint(true);
	printf("\r\n");
	printf("dG/dx = ");
	G.Derivative(x).DebugPrint(true);
	printf("\r\n");
	
	Function<4> H = 3*x*x*x;
	printf("H = ");
	H.DebugPrint(true);
	printf("\r\n");
	printf("dH/dx = ");
	H.Derivative(x).DebugPrint(true);
	printf("\r\n");
	
	return 0;
}