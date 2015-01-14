#include "tensor.h"
#include <stdio.h>

int main()
{
	Function<4> t = Function<4>(Vector<4>::T);
	Function<4> x = Function<4>(Vector<4>::X);
	Function<4> y = Function<4>(Vector<4>::Y);
	Function<4> z = Function<4>(Vector<4>::Z);
	
	TensorField<4> g = TensorField<4>(1);
	g[0] = -t*t;
	g[1] = x*x;
	g[2] = y*y;
	g[3] = z*z;
	
	printf("g_u = ");
	g.DebugPrint(true);
	printf("\r\n");
	
	printf("g_u,v = ");
	(g, NewIndex()).DebugPrint(true);
	printf("\r\n");
	
	return 0;
}