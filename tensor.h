#ifndef __cplusplus
#error TensorField and the classes it depend on require C++ functionality
#endif

#ifndef TENSOR_H
#define TENSOR_H

#include <math.h>

#define debuggy printf("line %d, %s\r\n", __LINE__, __FILE__)
#define PI		3.1415926535897932384626433832795
#define HALFPI 1.5707963267948966192313216916398

template <int d> class Vector;
template <int d> class Function;
template <int d> class TensorReference;
template <int d> class TensorField;

/*
		dummy class for comma derivatives, etc...
*/

class NewIndex
{
private:
public:
    NewIndex() { }
    ~NewIndex() { }
};

template <int d = 4>
class Vector
{
private:
	double	C[d];
public:
	Vector<d>()
	{
		for(unsigned int i = 0; i < d; ++i) C[i] = 0.0;
	}

	Vector<d>(const double& D)
	{
		for(unsigned int i = 0; i < d; ++i) C[i] = D;
	}
	
	~Vector<d>() {}

	double& operator [](unsigned int i)
	{
		return C[i];
	}

	static double T(Vector<d> Parameter) { return Parameter[0]; }
	static double X(Vector<d> Parameter) { return Parameter[1]; }
	static double Y(Vector<d> Parameter) { return Parameter[2]; }
	static double Z(Vector<d> Parameter) { return Parameter[3]; }
	static double W(Vector<d> Parameter) { return Parameter[4]; }

	static double X0(Vector<d> Parameter) { return Parameter[0]; }
	static double X1(Vector<d> Parameter) { return Parameter[1]; }
	static double X2(Vector<d> Parameter) { return Parameter[2]; }
	static double X3(Vector<d> Parameter) { return Parameter[3]; }
	static double X4(Vector<d> Parameter) { return Parameter[4]; }
	static double X5(Vector<d> Parameter) { return Parameter[5]; }
	static double X6(Vector<d> Parameter) { return Parameter[6]; }
	static double X7(Vector<d> Parameter) { return Parameter[7]; }
	static double X8(Vector<d> Parameter) { return Parameter[8]; }
	static double X9(Vector<d> Parameter) { return Parameter[9]; }
};

template <int d = 4>
class Function
{
private:
	double					C;												// if the function represents a constant this holds the value for speed
	double					(*VF)(Vector<d> Parameter);						// the function if it is a variable function
	double					(Function<d>::*F)(Vector<d> Parameter);			// the function pointer to be called

	Function<d>*			Child1;											// to allow building of trees of functions for function composition
	Function<d>*			Child2;

	int						CompositionParameter;
	double					(Function<d>::*Composition)(Vector<d> Parameter);

	// functions
	static double vfNull(Vector<d> Parameter)
	{
		return 0.0;
	}

	/*
	double fNull(Vector<d> Parameter)
	{
		return 0.0;
	}
	*/

	double fConstant(Vector<d> Parameter)
	{
		return C;
	}

	double fFunction(Vector<d> Parameter)
	{
		return VF(Parameter);
	}

	double fComposition(Vector<d> Parameter)
	{
		return (this->*Composition)(Parameter);
	}

	// compositions
	double cNull(Vector<d> Parameter)
	{
		debuggy;
		return 0.0;
	}

	double cAdd(Vector<d> Parameter)
	{
		return (*Child1)(Parameter) + (*Child2)(Parameter);
	}

	double cSub(Vector<d> Parameter)
	{
		return (*Child1)(Parameter) - (*Child2)(Parameter);
	}

	double cMul(Vector<d> Parameter)
	{
		return (*Child1)(Parameter) * (*Child2)(Parameter);
	}

	double cDiv(Vector<d> Parameter)
	{
		return (*Child1)(Parameter) / (*Child2)(Parameter);
	}

	double cPow(Vector<d> Parameter)
	{
		return powf((*Child1)(Parameter), (*Child2)(Parameter));
	}

	double cSin(Vector<d> Parameter)
	{
		return sinf((*Child1)(Parameter));
	}

	double cCos(Vector<d> Parameter)
	{
		return cosf((*Child1)(Parameter));
	}

	double cTan(Vector<d> Parameter)
	{
		return tanf((*Child1)(Parameter));
	}

	double cCot(Vector<d> Parameter)
	{
		return tanf(HALFPI - (*Child1)(Parameter));
	}

	double cSec(Vector<d> Parameter)
	{
		return 1.0/cosf((*Child1)(Parameter));
	}

	double cCsc(Vector<d> Parameter)
	{
		return 1.0/sinf((*Child1)(Parameter));
	}

public:
	Function<d>( ) : C( 0.0 ), VF( vfNull ), F( &Function<d>::fConstant ), Child1( 0 ), Child2( 0 ), Composition( &Function<d>::cNull ), CompositionParameter( 0 ) {}
	Function<d>( const double& dValue ) : C( dValue ), VF( vfNull ), F( &Function<d>::fConstant ), Child1( 0 ), Child2( 0 ), Composition( &Function<d>::cNull ), CompositionParameter( 0 ) {}
	Function<d>( double( *f )( Vector<d> Parameter ) ) : C( 0.0 ), F( &Function<d>::fFunction ), Child1( 0 ), Child2( 0 ), Composition( &Function<d>::cNull ), CompositionParameter( 0 )
	{
		VF = f;
		if(VF == Vector<d>::T) VF = Vector<d>::X0;
		if(VF == Vector<d>::X) VF = Vector<d>::X1;
		if(VF == Vector<d>::Y) VF = Vector<d>::X2;
		if(VF == Vector<d>::Z) VF = Vector<d>::X3;
		if(VF == Vector<d>::W) VF = Vector<d>::X4;
	}
	
	Function<d>(const Function<d>& f) : C(f.C), VF(f.VF), F(f.F), Composition(f.Composition), CompositionParameter(f.CompositionParameter)
	{
		if(f.Child1) Child1 = new Function<d>(*(f.Child1));
		else Child1 = 0;
		if(f.Child2) Child2 = new Function<d>(*(f.Child2));
		else Child2 = 0;
	}
	
	Function<d>( const Function<d>& F1, const Function<d>& F2, double ( Function<d>::*c )( Vector<d> Parameter ), int i = 0 ) : C( 0 ), VF( vfNull ), F( &Function<d>::fComposition ), Composition( c ), CompositionParameter( i )
	{
		Child1 = new Function<d>(F1);
		Child2 = new Function<d>(F2);
	}
	
	Function<d>(const TensorReference<d>& R)
	{
		*this = R.GetFunction();
	}

	~Function<d>()
	{
		if(Child1) delete Child1;
		if(Child2) delete Child2;
	}

	Function<d>& operator =(const Function<d>& f)
	{
		C = f.C;
		VF = f.VF;
		F = f.F;
		Child1 = new Function<d>(*(f.Child1));
		Child2 = new Function<d>(*(f.Child2));
		Composition = f.Composition;
		CompositionParameter = f.CompositionParameter;
		return *this;
	}

	double operator ()(Vector<d> Parameter)
	{
		return (this->*F)(Parameter);
	}
	
	Function<d> operator -()
	{
		return 0.0 - *this;
	}

	Function<d>& operator +=(const Function<d> &X)
	{
		return (*this = *this + X);
	}

	Function<d> operator +(const Function<d> &X) const
	{
		if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return C + X.C;
		}
		else if( ( F == &Function<d>::fFunction ) && ( X.F == &Function<d>::fFunction ) && ( VF == X.VF ) )
		{
			// X + X = 2X
			return 2*X;
		}
		else if( ( F == &Function<d>::fFunction ) && ( X.F == &Function<d>::fComposition ) )
		{
			// aX + X = (a+1)X
			if( X.Composition == &Function<d>::cMul )
			{
				if( ( X.Child1->F == &Function<d>::fConstant ) && ( X.Child2->F == &Function<d>::fFunction ) && ( X.Child2->VF == VF ) )
				{
					return (X.Child1->C + 1) * *(X.Child2);
				}
				else if( ( X.Child1->F == &Function<d>::fFunction ) && ( X.Child1->VF == VF ) && ( X.Child2->F == &Function<d>::fConstant ) )
				{
					return (X.Child2->C + 1) * *(X.Child1);
				}
			}
		}
		else if( ( X.F == &Function<d>::fFunction ) && ( F == &Function<d>::fComposition ) )
		{
			// Xa + X = (a+1)X
			if( Composition == &Function<d>::cMul )
			{
				if( ( Child1->F == &Function<d>::fConstant ) && ( Child2->F == &Function<d>::fFunction ) && ( Child2->VF == X.VF ) )
				{
					return (Child1->C + 1) * *(Child2);
				}
				else if( ( Child1->F == &Function<d>::fFunction ) && ( Child1->VF == X.VF ) && ( Child2->F == &Function<d>::fConstant ) )
				{
					return (Child2->C + 1) * *(Child1);
				}
			}
		}
		else if( ( F == &Function<d>::fComposition ) && ( X.F == &Function<d>::fComposition ) )
		{
			if( ( Composition == &Function<d>::cMul ) && ( X.Composition == &Function<d>::cMul ) )
			{
				if( ( Child1->F == &Function<d>::fFunction ) && ( X.Child1->F == &Function<d>::fFunction ) && ( Child1->VF == X.Child1->VF ) )
				{
					// X*a + X*b = (a+b)*X
					return (*(Child2) + *(X.Child2))*(*Child1);
				}
				else if( ( Child2->F == &Function<d>::fFunction ) && ( X.Child1->F == &Function<d>::fFunction ) && ( Child2->VF == X.Child1->VF ) )
				{
					// a*X + X*b = (a+b)*X
					return (*(Child1) + *(X.Child2))*(*Child2);
				}
				else if( ( Child1->F == &Function<d>::fFunction ) && ( X.Child2->F == &Function<d>::fFunction ) && ( Child1->VF == X.Child2->VF ) )
				{
					// X*a + b*X = (a+b)*X
					return (*(Child2) + *(X.Child1))*(*Child1);
				}
				else if( ( Child2->F == &Function<d>::fFunction ) && ( X.Child2->F == &Function<d>::fFunction ) && ( Child2->VF == X.Child2->VF ) )
				{
					// a*X + b*X = (a+b)*X
					return (*(Child1) + *(X.Child1))*(*Child2);
				}
			}
		}
		else if (F == &Function<d>::fConstant)
		{
			if(C == 0) return X; // X + 0 = X
			if(C < 0) return X - -C; // X + -a = X - a
			else if( X.F == &Function<d>::fComposition )
			{
				if( X.Composition == &Function<d>::cAdd )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return (X.Child1->C + C) + *(X.Child2);
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return *(X.Child1) + (X.Child2->C + C);
					}
				}
				else if( X.Composition == &Function<d>::cSub )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return (X.Child1->C + C) - *(X.Child2);
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return *(X.Child1) - (X.Child2->C - C);
					}
				}
			}
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if(X.C == 0) return *this;
			if(X.C < 0) return *this - -X.C;
			else if( F == &Function<d>::fComposition )
			{
				if( Composition == &Function<d>::cAdd )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return (Child1->C + X.C) + *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 + (Child2->C + X.C);
					}
				}
				else if( Composition == &Function<d>::cSub )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return (Child1->C + X.C) - *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 - (Child2->C - X.C);
					}
				}
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cAdd );
	}

	Function<d>& operator -=(const Function<d> &X)
	{
		return (*this = *this - X);
	}

	Function<d> operator -(const Function<d> &X) const
	{
		if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return C - X.C;
		}
		else if( ( F == &Function<d>::fFunction ) && ( X.F == &Function<d>::fFunction ) && ( VF == X.VF ) )
		{
			return 0.0;
		}
		else if( F == &Function<d>::fConstant )
		{
			if(C == 0) return X;
			else if( X.F == &Function<d>::fComposition )
			{
				if( X.Composition == &Function<d>::cAdd )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return (C - X.Child1->C) - *(X.Child2);
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return (C - *(X.Child1)) - X.Child2->C;
					}
				}
				else if( X.Composition == &Function<d>::cSub )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return (C - X.Child1->C) + *(X.Child2);
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return (C - *(X.Child1)) + X.Child2->C;
					}
				}
			}
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if(X.C == 0) return *this;
			if(X.C < 0) return *this + -X.C;
			else if( F == &Function<d>::fComposition )
			{
				if( Composition == &Function<d>::cAdd )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return (Child1->C - X.C) + *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 + (Child2->C - X.C);
					}
				}
				else if( Composition == &Function<d>::cSub )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return (Child1->C - X.C) - *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 - (Child2->C + X.C);
					}
				}
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cSub );
	}

	Function<d>& operator *=(const Function<d> &X)
	{
		return (*this = *this * X);
	}

	Function<d> operator *(const Function<d> &X) const
	{
		if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return Function<d>(C * X.C);
		}
		else if( F == &Function<d>::fConstant )
		{
			if(C == 0) return 0.0;
			else if(C == 1) return X;
			else if( X.F == &Function<d>::fComposition )
			{
				if( X.Composition == &Function<d>::cAdd )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return Function<d>( Function<d>( X.Child1->C * C ), Function<d>( Function<d>( C ), *( X.Child2 ), &Function<d>::cMul ), &Function<d>::cAdd );
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return Function<d>( Function<d>( Function<d>( C ), *( X.Child1 ), &Function<d>::cMul ), Function<d>( X.Child2->C * C ), &Function<d>::cAdd );
					}
				}
			}
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if(X.C == 0) return 0.0;
			else if(X.C == 1) return *this;
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cMul );
	}

	Function<d>& operator /=(const Function<d> &X)
	{
		return (*this = *this / X);
	}

	Function<d> operator /(const Function<d> &X) const
	{
		if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return Function<d>(C / X.C);
		}
		else if( ( F == &Function<d>::fFunction ) && ( X.F == &Function<d>::fFunction ) && ( VF == X.VF ) )
		{
			return 1.0;
		}
		else if( F == &Function<d>::fConstant )
		{
			if(C == 0) return 0.0;
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if( X.C == 0 )
			{
				// SE: we want to fuck ourselves here.
				const unsigned int uNAN = 0x7F800000;
				return *reinterpret_cast< const float* >( &uNAN ); // we're really fucked now
			}
			else if( X.C == 1 )
			{
				return *this;
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cDiv );
	}

	Function<d> Derivative(const Function& wrt) const
	{
		if( wrt.F == &Function<d>::fFunction )
		{
			return Derivative(wrt.VF);
		}
		else return 0.0;
	}
	Function<d> Derivative(double (*wrt)(Vector<d> Parameter)) const
	{
		if(wrt == Vector<d>::T) wrt = Vector<d>::X0;
		if(wrt == Vector<d>::X) wrt = Vector<d>::X1;
		if(wrt == Vector<d>::Y) wrt = Vector<d>::X2;
		if(wrt == Vector<d>::Z) wrt = Vector<d>::X3;
		if(wrt == Vector<d>::W) wrt = Vector<d>::X4;
		
		if( F == &Function<d>::fConstant )
		{
			return 0.0;
		}
		else if( F == &Function<d>::fFunction )
		{
			if(VF == wrt)
			{
				return 1.0;
			}
		}
		else if( F == &Function<d>::fComposition )
		{
			if( Composition == &Function<d>::cAdd )
			{
				return Child1->Derivative(wrt) + Child2->Derivative(wrt);
			}
			else if( Composition == &Function<d>::cSub )
			{
				return Child1->Derivative(wrt) - Child2->Derivative(wrt);
			}
			else if( Composition == &Function<d>::cMul ) // product rule
			{
				return (*Child2)*(Child1->Derivative(wrt)) + (*Child1)*(Child2->Derivative(wrt));
			}
			else if( Composition == &Function<d>::cDiv ) // quotient rule
			{
				return ((*Child2)*(Child1->Derivative(wrt)) - (*Child1)*(Child2->Derivative(wrt)))/((*Child2)*(*Child2));
			}
		}

		return 0.0;
	}

	void DebugPrint(bool txyz = false)
	{
		if( F == &Function<d>::fConstant )
		{
			printf("%g", C);
		}
		else if( F == &Function<d>::fFunction )
		{
			if((VF == Vector<d>::X0) || (VF == Vector<d>::T))
			{
				printf(txyz ? "t" : "x0");
			}
			else if((VF == Vector<d>::X1) || (VF == Vector<d>::X))
			{
				printf(txyz ? "x" : "x1");
			}
			else if((VF == Vector<d>::X2) || (VF == Vector<d>::Y))
			{
				printf(txyz ? "y" : "x2");
			}
			else if((VF == Vector<d>::X3) || (VF == Vector<d>::Z))
			{
				printf(txyz ? "z" : "x3");
			}
			else if((VF == Vector<d>::X4) || (VF == Vector<d>::W))
			{
				printf(txyz ? "w" : "x4");
			}
			else if(VF == Vector<d>::X5)
			{
				printf("x5");
			}
			else if(VF == Vector<d>::X6)
			{
				printf("x6");
			}
			else if(VF == Vector<d>::X7)
			{
				printf("x7");
			}
			else if(VF == Vector<d>::X8)
			{
				printf("x8");
			}
			else if(VF == Vector<d>::X9)
			{
				printf("x9");
			}
		}
		else if( F == &Function<d>::fComposition )
		{
			if( Composition == &Function<d>::cAdd )
			{
				//if(Child1->F == fComposition) printf("(");
				Child1->DebugPrint(txyz);
				//if(Child1->F == fComposition) printf(")");
				printf("+");
				//if(Child2->F == fComposition) printf("(");
				Child2->DebugPrint(txyz);
				//if(Child2->F == fComposition) printf(")");
			}
			else if( Composition == &Function<d>::cSub )
			{
				//if(Child1->F == fComposition) printf("(");
				Child1->DebugPrint(txyz);
				//if(Child1->F == fComposition) printf(")");
				printf("-");
				if( Child2->F == &Function<d>::fComposition ) printf( "(" );
				Child2->DebugPrint(txyz);
				if( Child2->F == &Function<d>::fComposition ) printf( ")" );
			}
			else if( Composition == &Function<d>::cMul )
			{
				if( Child1->F == &Function<d>::fConstant )
				{
					printf("%g", Child1->C);
				}
				else
				{
					if( Child1->F == &Function<d>::fComposition ) printf( "(" );
					Child1->DebugPrint(txyz);
					if( Child1->F == &Function<d>::fComposition ) printf( ")" );
					printf("*");
				}
				if( Child2->F == &Function<d>::fComposition ) printf( "(" );
				Child2->DebugPrint(txyz);
				if( Child2->F == &Function<d>::fComposition ) printf( ")" );
			}
			else if( Composition == &Function<d>::cDiv )
			{
				if( Child1->F == &Function<d>::fComposition ) printf( "(" );
				Child1->DebugPrint(txyz);
				if( Child1->F == &Function<d>::fComposition ) printf( ")" );
				printf("/");
				if( Child2->F == &Function<d>::fComposition ) printf( "(" );
				Child2->DebugPrint(txyz);
				if( Child2->F == &Function<d>::fComposition ) printf( ")" );
			}
		}
	}

	friend Function<4> operator +(double dValue, const Function<4> f)
	{
		return Function<4>( dValue ) +f;
	}

	friend Function<4> operator -( double dValue, const Function<4> f )
	{
		return Function<4>( dValue ) -f;
	}

	friend Function<4> operator *( double dValue, const Function<4> f )
	{
		return Function<4>( dValue ) * f;
	}

	friend Function<4> operator /( double dValue, const Function<4> f )
	{
		return Function<4>( dValue ) / f;
	}
};

/*
		tensor field, can also be used to represent individual tensors by using constants for all the components
*/

template <int d = 4>
class TensorReference
{
private:
	TensorField<d>*		Parent;
	unsigned int		D;
	unsigned int*		C;
public:
	TensorReference<d>() : Parent(0), D(0), C(0) {}
	
	TensorReference<d>(TensorField<d>* P, unsigned int numDims, unsigned int* components) : Parent(P), D(numDims)
	{
		C = new unsigned int[D];
		for(unsigned int i = 0; i < D; ++i) C[i] = components[i];
	}
	
	TensorReference<d>(TensorReference<d>& T) : Parent(T.Parent), D(T.D)
	{
		C = new unsigned int[D];
		for(unsigned int i = 0; i < D; ++i) C[i] = T.C[i];
	}
	
	~TensorReference<d>()
	{
		if(C) delete[] C;
	}
	
	TensorReference<d>& operator =(const TensorReference<d>& T)
	{
		Parent = T.Parent;
		D = T.D;
		C = new unsigned int[D];
		for(unsigned int i = 0; i < D; ++i) C[i] = T.C[i];
		return *this;
	}
	
	TensorReference<d>& operator =(const Function<d>& F)
	{
		debuggy;
		printf("%X %u %u", Parent, Parent->GetRank(), D);
		GetFunction() = F;
		debuggy;
		return *this;
	}
	
	TensorReference<d> operator [](unsigned int i)
	{
		if(D < Parent->GetRank())
        {
            // create array of indices
            unsigned int* p = new unsigned int[D + 1];
            
            // copy from old array and add i to the end
            for(unsigned int j = 0; j < D; ++j) p[j] = C[j];
            p[D] = i;
            
            // create return reference and delete array of indices
            TensorReference<d> ret = TensorReference<d>(Parent, D + 1, p);
            delete[] p;
            
            return ret;
        }
        // extra references should do nothing
        else
        {
            return (*this);
        }
	}
	
	Function<d>* operator ->()
	{
		return &(GetFunction());
	}
	
	Function<d>& GetFunction()
	{
		if(Parent->GetRank() == D)
		{
			debuggy;
			return Parent->FunctionFromIndices(C);
		}
		else
		{
			debuggy;
			return Parent->FirstEntry();
		}
	}
};

template <int d = 4>
class TensorField
{
protected:
    unsigned int            Rank;								// rank of tensor
	unsigned int            CovariantRank;
	unsigned int            ContravariantRank;
    Function<d>**  			Components;                     	// array of components
    unsigned int            Size;                     			// number of entries in total
    bool*                   ContravariantIndices;           	// which indices are contravariant

public:
	TensorField<d>(unsigned int rank = 2, bool* contravariantIndices = 0, double defaultValue = 0.0)
	{
		Rank = rank;
		CovariantRank = Rank;
		ContravariantRank = 0;

		if(contravariantIndices)
		{
			ContravariantIndices = new bool[Rank];
			for(unsigned int i = 0; i < Rank; ++i)
			{
				ContravariantIndices[i] = contravariantIndices[i];
				if(ContravariantIndices[i])
				{
					++ContravariantRank;
					--CovariantRank;
				}
			}
		}
		else ContravariantIndices = 0;

		Size = 1;
        for(unsigned int i = 0; i < Rank; ++i) Size *= d;

		Components = new Function<d>*[Size];
		for(unsigned int i = 0; i < Size; ++i) Components[i] = new Function<d>(defaultValue);
	}

	/*
	TensorField<d>(const double& d)
	{
		CovariantRank = 0;
		ContravariantRank = 0;
		Rank = 0;
		Size = 1;
		Components = new Function<d>*[1];
		Components[0] = new Function<d>(d);
		ContravariantIndices = 0;
	}
	*/

	TensorField<d>(TensorField& t)
	{
		Rank = t.Rank;
		CovariantRank = t.CovariantRank;
		ContravariantRank = t.ContravariantRank;
		Size = t.Size;

		ContravariantIndices = new bool[Rank];
		for(unsigned int i = 0; i < Rank; ++i) ContravariantIndices[i] = t.ContravariantIndices[i];

		Components = new Function<d>*[Size];
		for( unsigned int i = 0; i < Size; ++i )
		{
			Components[ i ] = new Function<d>( *( t.Components[ i ] ) );
		}
	}

	~TensorField<d>()
	{
		if(ContravariantIndices) delete[] ContravariantIndices;
		for(unsigned int i = 0; i < Size; ++i) delete Components[i];
		delete[] Components;
	}

	TensorField<d>& operator =(const TensorField<d>& t)
	{
		Rank = t.Rank;
		CovariantRank = t.CovariantRank;
		ContravariantRank = t.ContravariantRank;
		Size = t.Size;

		ContravariantIndices = new bool[Rank];
		for(unsigned int i = 0; i < Rank; ++i) ContravariantIndices[i] = t.ContravariantIndices[i];

		Components = new Function<d>*[Size];
		for(unsigned int i = 0; i < Size; ++i) Components[i] = new Function<d>(t.Components[i]);

		return *this;
	}

	TensorField<d> operator +(const TensorField<d>& t)
	{
		if(t.Rank != Rank) return TensorField<d>(0.0);
		if(t.ContravariantRank != ContravariantRank) return TensorField<d>(0.0);
		if(t.CovariantRank != CovariantRank) return TensorField<d>(0.0);

		TensorField<d> ret = *this;

		for(unsigned int i = 0; i < Size; ++i)
			ret.Components[i] += t.Components[i];

		return ret;
	}

	TensorField<d> operator -(const TensorField<d>& t)
	{
		if(t.Rank != Rank) return TensorField<d>(0.0);
		if(t.ContravariantRank != ContravariantRank) return TensorField<d>(0.0);
		if(t.CovariantRank != CovariantRank) return TensorField<d>(0.0);

		TensorField<d> ret = *this;

		for(unsigned int i = 0; i < Size; ++i)
			ret.Components[i] -= t.Components[i];

		return ret;
	}
	
	TensorField<d> operator ,(NewIndex _haxpile)
	{
		bool* b = new bool[Rank + 1];
		if(ContravariantIndices)
		{
			for(unsigned int i = 0; i < Rank; ++i) b[i] = ContravariantIndices[i];
		}
		else
		{
			for(unsigned int i = 0; i < Rank; ++i) b[i] = false;
		}
		b[Rank] = false;
		TensorField<d> ret = TensorField<d>(Rank + 1, b, 0.0);
		unsigned int Indices[d] = {0};
		for(unsigned int i = 0; i < d; ++i) Indices[i] = 0;
		while(Indices[0] < Rank)
		{
			printf("Comma derivative:\r\nRank, Dimension: %d, %d\r\nIndices: ", Rank, d);
			for(unsigned int k = 0; k < Rank; ++k) printf("%d ", Indices[k]);
			printf("\r\n");
			
			for(unsigned int i = 0; i < d; ++i)
			{
				double (*f)(Vector<d>);
				switch(i)
				{
					default:
					case 0:
						f = Vector<d>::X0;
						break;
					case 1:
						f = Vector<d>::X1;
						break;
					case 2:
						f = Vector<d>::X2;
						break;
					case 3:
						f = Vector<d>::X3;
						break;
					case 4:
						f = Vector<d>::X4;
						break;
					case 5:
						f = Vector<d>::X5;
						break;
					case 6:
						f = Vector<d>::X6;
						break;
					case 7:
						f = Vector<d>::X7;
						break;
					case 8:
						f = Vector<d>::X8;
						break;
					case 9:
						f = Vector<d>::X9;
						break;
				}
				debuggy;
				ret.ReferenceFromIndices(Indices, Rank)[i] = FunctionFromIndices(Indices).Derivative(f);
				debuggy;
			}
			
			++Indices[Rank-1];
			
			printf("Comma derivative:\r\nRank, Dimension: %d, %d\r\nIndices: ", Rank, d);
			for(unsigned int k = 0; k < Rank; ++k) printf("%d ", Indices[k]);
			printf("\r\n");
			
			for(unsigned int i = Rank - 1; i >= 0; --i)
			{
				if(Indices[i] >= d)
				{
					Indices[i] = 0;
					if(i > 0) ++Indices[i - 1];
				}
			}
		}
		
		return ret;
	}
	
	TensorReference<d> operator [](unsigned int i)
	{
		unsigned int* p = new unsigned int;
        *p = i;
        TensorReference<d> ret = TensorReference<d>(this, 1, p);
        delete p;
        return ret;
	}
	
	unsigned int GetRank() { return Rank; }
	
	Function<d>& FirstEntry()
	{
		return *(Components[0]);
	}
	
	TensorReference<d> ReferenceFromIndices(unsigned int* P, unsigned int D)
    {
		return TensorReference<d>(this, D, P);
    }
	
	Function<d>& FunctionFromIndices(unsigned int* P)
    {
        int j = 0;
        int k = 1;
        for(int i = 0; i < static_cast< int >( Rank ); ++i)
        {
            j += k*P[i];
            k *= d;
        }
		printf("Function from indices, index: %d\r\n", j);
        return *(Components[j]);
    }
	
	void DebugPrint(bool txyz = false)
	{
		if(Rank == 0)
		{
			(*(Components[0])).DebugPrint(txyz);
		}
		else if(Rank == 1)
		{
			printf("(");
			for(unsigned int i = 0; i < d - 1; ++i)
			{
				(*(Components[i])).DebugPrint(txyz);
				printf(", ");
			}
			(*(Components[d-1])).DebugPrint(txyz);
			printf(")");
		}
		else if(Rank == 2)
		{
			printf("[");
			for(unsigned int j = 0; j < d; ++j)
			{
				printf("(");
				for(unsigned int i = 0; i < d - 1; ++i)
				{
					(*this)[i][j]->DebugPrint(txyz);
					printf(", ");
				}
				(*this)[d-1][j]->DebugPrint(txyz);
				printf(")");
				if(j < (d - 1)) printf(",");
			}
			printf("]");
		}
	}
};

#endif