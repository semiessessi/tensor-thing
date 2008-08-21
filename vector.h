/* vector template class

has the following operators:

x =	y			- assignment
x[i]			- ith element of x
x +	y			- addition
x -	y			- subtraction
x *	y			- dot product
Abs(x)			- length (absolute value)
AbsSquared(x)	- square of length

*/

#ifndef _VECTOR_TEMPLATE_H_
#define _VECTOR_TEMPLATE_H_

#pragma warning( disable : 4068 )

#ifndef __cplusplus
#error vector types require C++ functionality
#endif

#include <math.h>

template <class c = double, int n = 4>
class vector
{
private:
	c	a[n];
public:
	vector<c,n>()
	{
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			a[i] = 0;
		}
	}

	vector<c,n>(c x)
	{
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			a[i] = x;
		}
	}

	vector<c,n>(c x, c y)
	{
		int tmp = n - n%2;

		#pragma omp parallel for
		for(int i = 0; i < tmp; i+=2)
		{
			a[i] = x;
			a[i+1] = y;
		}

		if(n!=tmp)
		{
			a[tmp] = x;
		}
	}

	vector<c,n>(c x, c y, c z)
	{
		int tmp = n - n%3;

		#pragma omp parallel for
		for(int i = 0; i < tmp; i+=3)
		{
			a[i] = x;
			a[i+1] = y;
			a[i+2] = z;
		}

		tmp = n - tmp;

		if(tmp==1)
		{
			a[n-1] = x;
		}
		else if(tmp==2)
		{
			a[n-2] = x;
			a[n-1] = y;
		}
	}

	vector<c,n>(c x, c y, c z, c w)
	{
		int tmp = n - n%4;

		#pragma omp parallel for
		for(int i = 0; i < tmp; i+=4)
		{
			a[i] = x;
			a[i+1] = y;
			a[i+2] = z;
			a[i+3] = w;
		}

		tmp = n - tmp;

		if(tmp==1)
		{
			a[n-1] = x;
		}
		else if(tmp==2)
		{
			a[n-2] = x;
			a[n-1] = y;
		}
		else if(tmp==3)
		{
			a[n-1] = x;
			a[n-2] = y;
			a[n-3] = z;
		}
	}

	vector<c,n>(c* x)
	{
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			a[i] = x[i];
		}
	}

	vector<c,n>(vector<c,n>& x)
	{
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			a[i] = x[i];
		}
	}

	c* __getDataPointer()
	{
		return a;
	}

	vector<c,n>& operator =(vector<c,n> x)
	{
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			a[i] = x[i];
		}

		return (*this);
	}

	c& operator [](int i)
	{
		return a[i];
	}

	vector<c,n> operator +(vector<c,n> x)
	{
		vector<c,n> ret;

		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			ret[i] = a[i] + x[i];
		}

		return ret;
	}

	vector<c,n> operator -(vector<c,n> x)
	{
		vector<c,n> ret;

		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			ret[i] = a[i] - x[i];
		}

		return ret;
	}

	vector<c,n> operator *(c x)
	{
		vector <c,n> ret;
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			ret[i] = a[i] * x;
		}
		
		return ret;
	}
	
	friend vector<c,n> operator *(c x, vector<c,n> v)
	{
		return v * x;
	}
	
	// dot
	c operator *(vector<c,n> x)
	{
		c ret = 0;
		
		#pragma omp parallel for
		for(int i = 0; i < n; i++)
		{
			ret += a[i] * x[i];
		}

		return ret;
	}

	friend c Abs(vector<c,n> x)
	{
		return (c)sqrt((double)(x*x));
	}

	friend c AbsSquared(vector<c,n> x)
	{
		return (x*x);
	}
};

typedef vector<long double,2> vector2_long_double;
typedef vector<double,2> vector2_double;
typedef vector<float,2> vector2_float;
typedef vector<int,2> vector2_int;
typedef vector<unsigned int,2> vector2_unsigned_int;
typedef vector<long,2> vector2_long;
typedef vector<unsigned long,2> vector2_unsigned_long;
typedef vector<long long,2> vector2_long_long;
typedef vector<unsigned long long,2> vector2_unsigned_long_long;

typedef vector<long double,3> vector3_long_double;
typedef vector<double,3> vector3_double;
typedef vector<float,3> vector3_float;
typedef vector<int,3> vector3_int;
typedef vector<unsigned int,3> vector3_unsigned_int;
typedef vector<long,3> vector3_long;
typedef vector<unsigned long,3> vector3_unsigned_long;
typedef vector<long long,3> vector3_long_long;
typedef vector<unsigned long long,3> vector3_unsigned_long_long;

typedef vector<long double,4> vector4_long_double;
typedef vector<double,4> vector4_double;
typedef vector<float,4> vector4_float;
typedef vector<int,4> vector4_int;
typedef vector<unsigned int,4> vector4_unsigned_int;
typedef vector<long,4> vector4_long;
typedef vector<unsigned long,4> vector4_unsigned_long;
typedef vector<long long,4> vector4_long_long;
typedef vector<unsigned long long,4> vector4_unsigned_long_long;

#endif