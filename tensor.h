#ifndef __cplusplus
#error TensorField and the classes it depends on require C++ functionality
#endif

#ifndef TENSOR_H
#define TENSOR_H

#include "jProfiler.h"

#include <math.h>

#define debuggy printf("line %d, %s\r\n", __LINE__, __FILE__)
#define PI        3.1415926535897932384626433832795
#define HALFPI 1.5707963267948966192313216916398

template <int d> class Vector;
template <int d> class Function;
template <int d> class TensorReference;
template <int d> class TensorField;

#define Swap(a, b) a ^= b; b ^= a; a ^= b

/*
dummy class for comma derivatives, etc...
*/

class NewIndex {};

template <int d = 4>
class Vector
{
private:
	double    C[ d ];
public:
	Vector<d>( )
	{
		for( unsigned int i = 0; i < d; ++i ) C[ i ] = 0.0;
	}

	Vector<d>( const double& D )
	{
		for( unsigned int i = 0; i < d; ++i ) C[ i ] = D;
	}

	double& operator []( unsigned int i )
	{
		return C[ i ];
	}

	static double T( Vector<d> Parameter ) { return Parameter[ 0 ]; }
	static double X( Vector<d> Parameter ) { return Parameter[ 1 ]; }
	static double Y( Vector<d> Parameter ) { return Parameter[ 2 ]; }
	static double Z( Vector<d> Parameter ) { return Parameter[ 3 ]; }
	static double W( Vector<d> Parameter ) { return Parameter[ 4 ]; }

	static double X0( Vector<d> Parameter ) { return Parameter[ 0 ]; }
	static double X1( Vector<d> Parameter ) { return Parameter[ 1 ]; }
	static double X2( Vector<d> Parameter ) { return Parameter[ 2 ]; }
	static double X3( Vector<d> Parameter ) { return Parameter[ 3 ]; }
	static double X4( Vector<d> Parameter ) { return Parameter[ 4 ]; }
	static double X5( Vector<d> Parameter ) { return Parameter[ 5 ]; }
	static double X6( Vector<d> Parameter ) { return Parameter[ 6 ]; }
	static double X7( Vector<d> Parameter ) { return Parameter[ 7 ]; }
	static double X8( Vector<d> Parameter ) { return Parameter[ 8 ]; }
	static double X9( Vector<d> Parameter ) { return Parameter[ 9 ]; }

	//template <int i = 0>
	//friend double Mass(Vector<d> Parameter) { return GetMass(i); }

	static double( *GetCoordinate( unsigned int Dimension ) )( Vector<d> )
	{
		switch( Dimension )
		{
			default:
			case 0:
				return Vector<d>::X0;
			case 1:
				return Vector<d>::X1;
			case 2:
				return Vector<d>::X2;
			case 3:
				return Vector<d>::X3;
			case 4:
				return Vector<d>::X4;
			case 5:
				return Vector<d>::X5;
			case 6:
				return Vector<d>::X6;
			case 7:
				return Vector<d>::X7;
			case 8:
				return Vector<d>::X8;
			case 9:
				return Vector<d>::X9;
		}
	}
};

template <int d = 4>
class Function
{
private:
	double                      C;                                                  // if the function represents a constant this holds the value for speed
	double( *VF )( Vector<d> Parameter );                        // the function if it is a variable function
	double                      ( Function<d>::*F )( Vector<d> Parameter );            // the function pointer to be called

	Function<d>*                Child1;                                            // to allow building of trees of functions for function composition
	Function<d>*                Child2;

	int                         CompositionParameter;
	double                      ( Function<d>::*Composition )( Vector<d> Parameter );

	// functions
	static double vfNull( Vector<d> Parameter )
	{
		return 0.0;
	}

	/*
	double fNull(Vector<d> Parameter)
	{
	return 0.0;
	}
	*/

	double fConstant( Vector<d> Parameter )
	{
		return C;
	}

	double fFunction( Vector<d> Parameter )
	{
		return VF( Parameter );
	}

	double fComposition( Vector<d> Parameter )
	{
		return ( this->*Composition )( Parameter );
	}

	// compositions
	double cNull( Vector<d> Parameter )
	{
		return 0.0;
	}

	double cAdd( Vector<d> Parameter )
	{
		return ( *Child1 )(Parameter)+( *Child2 )( Parameter );
	}

	double cSub( Vector<d> Parameter )
	{
		return ( *Child1 )(Parameter)-( *Child2 )( Parameter );
	}

	double cMul( Vector<d> Parameter )
	{
		return ( *Child1 )(Parameter)* ( *Child2 )( Parameter );
	}

	double cDiv( Vector<d> Parameter )
	{
		return ( *Child1 )( Parameter ) / ( *Child2 )( Parameter );
	}

	double cPow( Vector<d> Parameter )
	{
		return pow( ( *Child1 )( Parameter ), ( *Child2 )( Parameter ) );
	}

	double cSin( Vector<d> Parameter )
	{
		return sin( ( *Child1 )( Parameter ) );
	}

	double cCos( Vector<d> Parameter )
	{
		return cos( ( *Child1 )( Parameter ) );
	}

	double cTan( Vector<d> Parameter )
	{
		return tan( ( *Child1 )( Parameter ) );
	}

	double cCot( Vector<d> Parameter )
	{
		return tan( HALFPI - ( *Child1 )( Parameter ) );
	}

	double cSec( Vector<d> Parameter )
	{
		return 1.0 / cos( ( *Child1 )( Parameter ) );
	}

	double cCsc( Vector<d> Parameter )
	{
		return 1.0 / sin( ( *Child1 )( Parameter ) );
	}

	static double ZeroHax( ) { return 0; }

public:
	Function<d>( ) : C( 0.0 ), VF( vfNull ), F( &Function<d>::fConstant ), Child1( 0 ), Child2( 0 ), Composition( &Function<d>::cNull ), CompositionParameter( 0 ) {}
	Function<d>( const double& D ) : C( D ), VF( vfNull ), F( &Function<d>::fConstant ), Child1( 0 ), Child2( 0 ), Composition( &Function<d>::cNull ), CompositionParameter( 0 ) {}

	Function<d>( double( *f )( Vector<d> Parameter ) ) : C( 0.0 ), F( &Function<d>::fFunction ), Child1( 0 ), Child2( 0 ), Composition( &Function<d>::cNull ), CompositionParameter( 0 )
	{
		VF = f;
		if( VF == Vector<d>::T ) VF = Vector<d>::X0;
		if( VF == Vector<d>::X ) VF = Vector<d>::X1;
		if( VF == Vector<d>::Y ) VF = Vector<d>::X2;
		if( VF == Vector<d>::Z ) VF = Vector<d>::X3;
		if( VF == Vector<d>::W ) VF = Vector<d>::X4;
	}

	Function<d>( const Function<d>& f ) : C( f.C ), VF( f.VF ), F( f.F ), Composition( f.Composition ), CompositionParameter( f.CompositionParameter )
	{
		if( f.Child1 ) Child1 = new Function<d>( *( f.Child1 ) );
		else Child1 = 0;
		if( f.Child2 ) Child2 = new Function<d>( *( f.Child2 ) );
		else Child2 = 0;
	}

	Function<d>( const Function<d>& F1, const Function<d>& F2, double ( Function<d>::*c )( Vector<d> Parameter ), int i = 0 ) : C( 0 ), VF( vfNull ), F( &Function<d>::fComposition ), Composition( c ), CompositionParameter( i )
	{
		Child1 = new Function<d>( F1 );
		Child2 = new Function<d>( F2 );
	}

	Function<d>( const TensorReference<d>& R )
	{
		*this = R.GetFunction( );
	}

	~Function<d>( )
	{
		if( Child1 ) delete Child1;
		if( Child2 ) delete Child2;
	}

	bool IsConstant( )
	{
		return F == &Function<d>::fConstant;
	}

	Function<d>& operator =( const Function<d>& f )
	{
		C = f.C;
		VF = f.VF;
		F = f.F;
		if( f.Child1 ) Child1 = new Function<d>( *( f.Child1 ) );
		else Child1 = 0;
		if( f.Child2 ) Child2 = new Function<d>( *( f.Child2 ) );
		else Child2 = 0;
		Composition = f.Composition;
		CompositionParameter = f.CompositionParameter;
		return *this;
	}

	double operator ()( Vector<d> Parameter )
	{
		return ( this->*F )( Parameter );
	}

	bool operator ==( const Function& f ) const
	{
		if( ( F == &Function<d>::fFunction ) && ( f.F == &Function<d>::fFunction ) && ( VF == f.VF ) ) return true;
		else if( ( F == &Function<d>::fConstant ) && ( f.F == &Function<d>::fConstant ) && ( C == f.C ) ) return true;
		else if( ( F == &Function<d>::fComposition ) && ( f.F == &Function<d>::fComposition ) && ( Composition == f.Composition ) )
		{
			if( ( Composition == &Function<d>::cAdd ) || ( Composition == &Function<d>::cMul ) )
			{
				return ( ( *Child1 == *( f.Child1 ) ) && ( *Child2 == *( f.Child2 ) ) ) ||
					( ( *Child1 == *( f.Child2 ) ) && ( *Child2 == *( f.Child1 ) ) );
			}
			else return ( *Child1 == *( f.Child1 ) ) && ( *Child2 == *( f.Child2 ) );
		}

		else return false;
	}

	bool operator !=( const Function& f ) const
	{
		return !( *this == f );
	}

	Function<d> operator -( ) const
	{
		return -1.0*( *this );
	}

	Function<d>& operator +=( const Function<d> &X )
	{
		return ( *this = *this + X );
	}

	Function<d> operator +( const Function<d> &X ) const
	{
		if( X == 0.0 ) return *this;
		else if( *this == 0.0 ) return X;
		else if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return C + X.C;
		}
		else if( *this == X )
		{
			// X + X = 2X
			return 2 * X;
		}
		else if( X.F == &Function<d>::fComposition )
		{
			// X + aX = (a+1)X
			if( X.Composition == &Function<d>::cMul )
			{
				if( *( X.Child2 ) == *this )
				{
					return ( *( X.Child1 ) + 1 ) * *( X.Child2 );
				}
				else if( *( X.Child1 ) == *this )
				{
					return ( *( X.Child2 ) + 1 ) * *( X.Child1 );
				}
			}
		}
		else if( F == &Function<d>::fComposition )
		{
			// aX + X = (a+1)X
			if( Composition == &Function<d>::cMul )
			{
				if( *Child2 == X )
				{
					return ( *Child1 + 1 ) * *( Child2 );
				}
				else if( *Child1 == X )
				{
					return ( *Child2 + 1 ) * *( Child1 );
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
					return ( *(Child2)+*( X.Child2 ) )*( *Child1 );
				}
				else if( ( Child2->F == &Function<d>::fFunction ) && ( X.Child1->F == &Function<d>::fFunction ) && ( Child2->VF == X.Child1->VF ) )
				{
					// a*X + X*b = (a+b)*X
					return ( *(Child1)+*( X.Child2 ) )*( *Child2 );
				}
				else if( ( Child1->F == &Function<d>::fFunction ) && ( X.Child2->F == &Function<d>::fFunction ) && ( Child1->VF == X.Child2->VF ) )
				{
					// X*a + b*X = (a+b)*X
					return ( *(Child2)+*( X.Child1 ) )*( *Child1 );
				}
				else if( ( Child2->F == &Function<d>::fFunction ) && ( X.Child2->F == &Function<d>::fFunction ) && ( Child2->VF == X.Child2->VF ) )
				{
					// a*X + b*X = (a+b)*X
					return ( *(Child1)+*( X.Child1 ) )*( *Child2 );
				}
			}
		}
		else if( F == &Function<d>::fConstant )
		{
			if( C < 0 ) return X - -C; // X + -a = X - a
			else if( X.F == &Function<d>::fComposition )
			{
				if( X.Composition == &Function<d>::cAdd )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return ( X.Child1->C + C ) + *( X.Child2 );
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return *( X.Child1 ) + ( X.Child2->C + C );
					}
				}
				else if( X.Composition == &Function<d>::cSub )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return ( X.Child1->C + C ) - *( X.Child2 );
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return *( X.Child1 ) - ( X.Child2->C - C );
					}
				}
			}
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if( X.C < 0 ) return *this - -X.C;
			else if( F == &Function<d>::fComposition )
			{
				if( Composition == &Function<d>::cAdd )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return ( Child1->C + X.C ) + *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 + ( Child2->C + X.C );
					}
				}
				else if( Composition == &Function<d>::cSub )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return ( Child1->C + X.C ) - *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 - ( Child2->C - X.C );
					}
				}
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cAdd );
	}

	Function<d>& operator -=( const Function<d> &X )
	{
		return ( *this = *this - X );
	}

	Function<d> operator -( const Function<d> &X ) const
	{
		if( *this == X )
		{
			return 0.0;
		}
		else if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return C - X.C;
		}
		else if( F == &Function<d>::fConstant )
		{
			if( C == 0.0 ) return -X;
			else if( X.F == &Function<d>::fComposition )
			{
				if( X.Composition == &Function<d>::cAdd )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return ( C - X.Child1->C ) - *( X.Child2 );
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return ( C - *( X.Child1 ) ) - X.Child2->C;
					}
				}
				else if( X.Composition == &Function<d>::cSub )
				{
					if( X.Child1->F == &Function<d>::fConstant )
					{
						return ( C - X.Child1->C ) + *( X.Child2 );
					}
					else if( X.Child2->F == &Function<d>::fConstant )
					{
						return ( C - *( X.Child1 ) ) + X.Child2->C;
					}
				}
			}
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if( X.C == 0 ) return *this;
			if( X.C < 0 ) return *this + -X.C;
			else if( F == &Function<d>::fComposition )
			{
				if( Composition == &Function<d>::cAdd )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return ( Child1->C - X.C ) + *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 + ( Child2->C - X.C );
					}
				}
				else if( Composition == &Function<d>::cSub )
				{
					if( Child1->F == &Function<d>::fConstant )
					{
						return ( Child1->C - X.C ) - *Child2;
					}
					else if( Child2->F == &Function<d>::fConstant )
					{
						return *Child1 - ( Child2->C + X.C );
					}
				}
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cSub );
	}

	Function<d>& operator *=( const Function<d> &X )
	{
		return ( *this = *this * X );
	}

	Function<d> operator *( const Function<d> &X ) const
	{
		if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return Function<d>( C * X.C );
		}
		else if( F == &Function<d>::fConstant )
		{
			if( C == 0 ) return 0.0;
			else if( C == 1 ) return X;
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
			if( X.C == 0 ) return 0.0;
			else if( X.C == 1 ) return *this;
		}
		else if( F == &Function<d>::fComposition )
		{
			if( X.F == &Function<d>::fComposition )
			{
				if( ( Composition == &Function<d>::cDiv ) && ( X.Composition == &Function<d>::cDiv ) )
				{
					return ( *Child1 * *( X.Child1 ) ) / ( *Child2 * *( X.Child2 ) );
				}
			}
			else if( Composition == &Function<d>::cDiv )
			{
				return ( X * *Child1 ) / *Child2;
			}
		}
		else if( X.F == &Function<d>::fComposition )
		{
			if( X.Composition == &Function<d>::cDiv )
			{
				return ( *this * *( X.Child1 ) ) / *( X.Child2 );
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cMul );
	}

	Function<d>& operator /=( const Function<d> &X )
	{
		return ( *this = *this / X );
	}

	Function<d> operator /( const Function<d> &X ) const
	{

		if( *this == X )
		{
			return 1.0;
		}
		else if( ( F == &Function<d>::fConstant ) && ( X.F == &Function<d>::fConstant ) )
		{
			return Function<d>( C / X.C );
		}
		else if( F == &Function<d>::fConstant )
		{
			if( C == 0.0 ) return 0.0;
			else if( X.F == &Function<d>::fComposition )
			{
				if( X.Composition == &Function<d>::cDiv )
				{
					return ( C* *( X.Child2 ) ) / *( X.Child1 );
				}
			}
		}
		else if( X.F == &Function<d>::fConstant )
		{
			if( X.C == 0.0 ) return 1 / Function<d>::ZeroHax( );//1.0/0.0; // we're really fucked now
			else if( X.C == 1.0 ) return *this;
		}
		else if( F == &Function<d>::fComposition )
		{
			if( X.F == &Function<d>::fComposition )
			{
				if( ( Composition == &Function<d>::cDiv ) && ( X.Composition == &Function<d>::cDiv ) )
				{
					return ( *Child1 * *( X.Child2 ) ) / ( *Child2 * *( X.Child1 ) );
				}
			}

			if( Composition == &Function<d>::cDiv )
			{
				return *Child1 / ( *Child2 * X );
			}

			// try cancellation
			unsigned int N1 = CountMultiplicationChain( ) + 1;
			unsigned int N2 = X.CountMultiplicationChain( ) + 1;

			if( ( N1 > 1 ) || ( N2 > 1 ) )
			{
				Function<d>* C1 = new Function<d>[ N1 ];
				Function<d>* C2 = new Function<d>[ N2 ];

				for( unsigned int i = 0; i < N1; ++i ) C1[ i ] = GetFromMultiplicationChain( i );
				for( unsigned int i = 0; i < N2; ++i ) C2[ i ] = X.GetFromMultiplicationChain( i );

				for( unsigned int i = 0; i < N1; ++i )
				for( unsigned int j = 0; j < N2; ++j )
				{
					if( C1[ i ] == C2[ j ] )
					{
						Function<d> Numerator = 1.0;
						Function<d> Denominator = 1.0;
						for( unsigned int n = 0; n < N1; ++n )
						{
							if( n != i ) Numerator *= C1[ n ];
						}
						for( unsigned int m = 0; m < N2; ++m )
						{
							if( m != j ) Denominator *= C2[ m ];
						}
						delete[ ] C1;
						delete[ ] C2;
						return Numerator / Denominator;
					}
				}

				delete[ ] C1;
				delete[ ] C2;
			}
		}

		Function<d> f = *this;
		return Function<d>( f, X, &Function<d>::cDiv );
	}

	friend Function<d> Sin( const Function<d>& F )
	{
		return Function<d>( F, 0.0, &Function<d>::cSin );
	}

	friend Function<d> Cos( const Function<d>& F )
	{
		return Function<d>( F, 0.0, &Function<d>::cCos );
	}

	friend Function<d> Tan( const Function<d>& F )
	{
		return Function<d>( F, 0.0, &Function<d>::cTan );
	}

	friend Function<d> Cot( const Function<d>& F )
	{
		return Function<d>( F, 0.0, &Function<d>::cCot );
	}

	friend Function<d> Sec( const Function<d>& F )
	{
		return Function<d>( F, 0.0, &Function<d>::cSec );
	}

	friend Function<d> Csc( const Function<d>& F )
	{
		return Function<d>( F, 0.0, &Function<d>::cCsc );
	}

	Function<d> Derivative( const Function& wrt ) const
	{
		if( wrt.F == &Function<d>::fFunction )
		{
			return Derivative( wrt.VF );
		}
		else return 0.0;
	}

	Function<d> Derivative( double( *wrt )( Vector<d> Parameter ) ) const
	{
		if( F == &Function<d>::fConstant )
		{
			return 0.0;
		}
		else if( F == &Function<d>::fFunction )
		{
			if( VF == wrt )
			{
				return 1.0;
			}
		}
		else if( F == &Function<d>::fComposition )
		{
			if( Composition == &Function<d>::cAdd )
			{
				return Child1->Derivative( wrt ) + Child2->Derivative( wrt );
			}
			else if( Composition == &Function<d>::cSub )
			{
				return Child1->Derivative( wrt ) - Child2->Derivative( wrt );
			}
			else if( Composition == &Function<d>::cMul ) // product rule
			{
				return ( *Child2 )*( Child1->Derivative( wrt ) ) + ( *Child1 )*( Child2->Derivative( wrt ) );
			}
			else if( Composition == &Function<d>::cDiv ) // quotient rule
			{
				return ( ( *Child2 )*( Child1->Derivative( wrt ) ) - ( *Child1 )*( Child2->Derivative( wrt ) ) ) / ( ( *Child2 )*( *Child2 ) );
			}
			else if( Composition == &Function<d>::cSin ) // chain rule for sin
			{
				return ( Child1->Derivative( wrt ) )*Cos( *Child1 );
			}
			else if( Composition == &Function<d>::cCos ) // chain rule for cos
			{
				return -( Child1->Derivative( wrt ) )*Sin( *Child1 );
			}
			else if( Composition == &Function<d>::cTan ) // chain rule for tan
			{
				return ( Child1->Derivative( wrt ) )*Sec( *Child1 )*Sec( *Child1 );
			}
			else if( Composition == &Function<d>::cCot ) // chain rule for cot
			{
				return -( Child1->Derivative( wrt ) )*Csc( *Child1 )*Csc( *Child1 );
			}
			else if( Composition == &Function<d>::cSec ) // chain rule for sec
			{
				return ( Child1->Derivative( wrt ) )*Sec( *Child1 )*Tan( *Child1 );
			}
			else if( Composition == &Function<d>::cCsc ) // chain rule for csc
			{
				return -( Child1->Derivative( wrt ) )*Csc( *Child1 )*Cot( *Child1 );
			}
		}

		return 0.0;
	}

	unsigned int CountMultiplicationChain( ) const
	{
		if( ( F == &Function<d>::fComposition ) && ( Composition == &Function<d>::cMul ) ) return 1 + Child1->CountMultiplicationChain( ) + Child2->CountMultiplicationChain( );
		else return 0;
	}

	Function<d> GetFromMultiplicationChain( unsigned int i ) const
	{
		if( ( F == &Function<d>::fComposition ) && ( Composition == &Function<d>::cMul ) )
		{
			unsigned int n = Child1->CountMultiplicationChain( ) + 1;
			if( i < n )
			{
				return Child1->GetFromMultiplicationChain( i );
			}
			else
			{
				return Child2->GetFromMultiplicationChain( i - n );
			}
		}
		else
		{
			return *this;
		}
	}

	void DebugPrint( bool txyz = false )
	{
		if( F == &Function<d>::fConstant )
		{
			printf( "%g", C );
		}
		else if( F == &Function<d>::fFunction )
		{
			if( ( VF == Vector<d>::X0 ) || ( VF == Vector<d>::T ) )
			{
				printf( txyz ? "t" : "x0" );
			}
			else if( ( VF == Vector<d>::X1 ) || ( VF == Vector<d>::X ) )
			{
				printf( txyz ? "x" : "x1" );
			}
			else if( ( VF == Vector<d>::X2 ) || ( VF == Vector<d>::Y ) )
			{
				printf( txyz ? "y" : "x2" );
			}
			else if( ( VF == Vector<d>::X3 ) || ( VF == Vector<d>::Z ) )
			{
				printf( txyz ? "z" : "x3" );
			}
			else if( ( VF == Vector<d>::X4 ) || ( VF == Vector<d>::W ) )
			{
				printf( txyz ? "w" : "x4" );
			}
			else if( VF == Vector<d>::X5 )
			{
				printf( "x5" );
			}
			else if( VF == Vector<d>::X6 )
			{
				printf( "x6" );
			}
			else if( VF == Vector<d>::X7 )
			{
				printf( "x7" );
			}
			else if( VF == Vector<d>::X8 )
			{
				printf( "x8" );
			}
			else if( VF == Vector<d>::X9 )
			{
				printf( "x9" );
			}
		}
		else if( F == &Function<d>::fComposition )
		{
			if( Composition == &Function<d>::cAdd )
			{
				//if(Child1->F == &Function<d>::fComposition) printf("(");
				Child1->DebugPrint( txyz );
				//if(Child1->F == &Function<d>::fComposition) printf(")");
				printf( "+" );
				//if(Child2->F == &Function<d>::fComposition) printf("(");
				Child2->DebugPrint( txyz );
				//if(Child2->F == &Function<d>::fComposition) printf(")");
			}
			else if( Composition == &Function<d>::cSub )
			{
				//if(Child1->F == &Function<d>::fComposition) printf("(");
				Child1->DebugPrint( txyz );
				//if(Child1->F == &Function<d>::fComposition) printf(")");
				printf( "-" );
				if( Child2->F == &Function<d>::fComposition ) printf( "(" );
				Child2->DebugPrint( txyz );
				if( Child2->F == &Function<d>::fComposition ) printf( ")" );
			}
			else if( Composition == &Function<d>::cMul )
			{
				if( Child1->F == &Function<d>::fConstant )
				{
					if( Child1->C == -1.0 ) printf( "-" );
					else printf( "%g", Child1->C );
				}
				else
				{
					if( Child1->F == &Function<d>::fComposition ) printf( "(" );
					Child1->DebugPrint( txyz );
					if( Child1->F == &Function<d>::fComposition ) printf( ")" );
					printf( "*" );
				}
				if( Child2->F == &Function<d>::fComposition ) printf( "(" );
				Child2->DebugPrint( txyz );
				if( Child2->F == &Function<d>::fComposition ) printf( ")" );
			}
			else if( Composition == &Function<d>::cDiv )
			{
				if( Child1->F == &Function<d>::fComposition ) printf( "(" );
				Child1->DebugPrint( txyz );
				if( Child1->F == &Function<d>::fComposition ) printf( ")" );
				printf( "/" );
				if( Child2->F == &Function<d>::fComposition ) printf( "(" );
				Child2->DebugPrint( txyz );
				if( Child2->F == &Function<d>::fComposition ) printf( ")" );
			}
			else if( Composition == &Function<d>::cSin )
			{
				printf( "sin(" );
				Child1->DebugPrint( txyz );
				printf( ")" );
			}
			else if( Composition == &Function<d>::cCos )
			{
				printf( "cos(" );
				Child1->DebugPrint( txyz );
				printf( ")" );
			}
			else if( Composition == &Function<d>::cTan )
			{
				printf( "tan(" );
				Child1->DebugPrint( txyz );
				printf( ")" );
			}
			else if( Composition == &Function<d>::cCot )
			{
				printf( "cot(" );
				Child1->DebugPrint( txyz );
				printf( ")" );
			}
			else if( Composition == &Function<d>::cSec )
			{
				printf( "sec(" );
				Child1->DebugPrint( txyz );
				printf( ")" );
			}
			else if( Composition == &Function<d>::cCsc )
			{
				printf( "csc(" );
				Child1->DebugPrint( txyz );
				printf( ")" );
			}
		}
	}

	friend Function<d> operator +( double D, const Function<d>& f )
	{
		return Function<d>( D ) +f;
	}

	friend Function<d> operator -( double D, const Function<d>& f )
	{
		return Function<4>( D ) -f;
	}

	friend Function<d> operator *( double D, const Function<d>& f )
	{
		return Function<d>( D ) * f;
	}

	friend Function<d> operator /( double D, const Function<d>& f )
	{
		return Function<d>( D ) / f;
	}
};

/*
tensor field, can also be used to represent individual tensors by using constants for all the components
*/

template <int d = 4>
class TensorReference
{
private:
	TensorField<d>*        Parent;
	unsigned int        D;
	unsigned int*        C;
public:
	TensorReference<d>( ) : Parent( 0 ), D( 0 ), C( 0 ) {}

	TensorReference<d>( TensorField<d>* P, unsigned int numDims, unsigned int* components ) : Parent( P ), D( numDims )
	{
		C = new unsigned int[ D ];
		for( unsigned int i = 0; i < D; ++i ) C[ i ] = components[ i ];
	}

	TensorReference<d>( const TensorReference<d>& T ) : Parent( T.Parent ), D( T.D )
	{
		C = new unsigned int[ D ];
		for( unsigned int i = 0; i < D; ++i ) C[ i ] = T.C[ i ];
	}

	~TensorReference<d>( )
	{
		if( C ) delete[ ] C;
	}

	TensorReference<d>& operator =( const TensorReference<d>& T )
	{
		Parent = T.Parent;
		D = T.D;
		if( C ) delete[ ] C;
		C = new unsigned int[ D ];
		for( unsigned int i = 0; i < D; ++i ) C[ i ] = T.C[ i ];
		return *this;
	}

	TensorReference<d>& operator =( const Function<d>& F )
	{
		GetFunction( ) = F;
		return *this;
	}

	TensorReference<d> operator []( unsigned int i ) const
	{
		if( D < Parent->GetRank( ) )
		{
			// create array of indices
			unsigned int* p = new unsigned int[ D + 1 ];

			// copy from old array and add i to the end
			for( unsigned int j = 0; j < D; ++j ) p[ j ] = C[ j ];
			p[ D ] = i;

			// create return reference and delete array of indices
			TensorReference<d> ret = TensorReference<d>( Parent, D + 1, p );
			delete[ ] p;

			return ret;
		}
		// extra references should do nothing
		else
		{
			return TensorReference<d>( *this );
		}
	}

	Function<d>* operator ->( ) const
	{
		return &( GetFunction( ) );
	}

	bool operator ==( const Function<d>& F ) const
	{
		return GetFunction( ) == F;
	}

	bool operator !=( const Function<d>& F ) const
	{
		return GetFunction( ) != F;
	}

	Function<d> operator -( ) const
	{
		return -Function<d>( GetFunction( ) );
	}

	Function<d> operator +( const Function<d>& F ) const
	{
		return GetFunction( ) + F;
	}

	Function<d> operator -( const Function<d>& F ) const
	{
		return GetFunction( ) - F;
	}

	Function<d> operator *( const Function<d>& F ) const
	{
		return GetFunction( ) * F;
	}

	Function<d> operator /( const Function<d>& F ) const
	{
		return GetFunction( ) / F;
	}

	double operator ()( Vector<d> Parameter )
	{
		return ( GetFunction( ) )( Parameter );
	}

	Function<d>& GetFunction( ) const
	{
		if( Parent->GetRank( ) == D )
		{
			return Parent->FunctionFromIndices( C );
		}
		else
		{
			return Parent->FirstEntry( );
		}
	}

	friend Function<d> operator +( double D, const TensorReference<d>& R )
	{
		return Function<d>( D ) +R;
	}

	friend Function<d> operator -( double D, const TensorReference<d>& R )
	{
		return Function<4>( D ) -R;
	}

	friend Function<d> operator *( double D, const TensorReference<d>& R )
	{
		return Function<d>( D ) * R;
	}

	friend Function<d> operator /( double D, const TensorReference<d>& R )
	{
		return Function<d>( D ) / R;
	}
};

template <int d = 4>
class TensorField
{
protected:
	unsigned int            Rank;                                // rank of tensor
	unsigned int            CovariantRank;
	unsigned int            ContravariantRank;
	Function<d>**              Components;                         // array of components
	unsigned int            Size;                                 // number of entries in total
	bool*                   ContravariantIndices;               // which indices are contravariant

public:
	TensorField<d>( unsigned int rank = 2, bool* contravariantIndices = 0, Function<d> defaultValue = Function<d>( 0.0 ) )
	{
		Rank = rank;
		CovariantRank = Rank;
		ContravariantRank = 0;

		if( contravariantIndices )
		{
			ContravariantIndices = new bool[ Rank ];
			for( unsigned int i = 0; i < Rank; ++i )
			{
				ContravariantIndices[ i ] = contravariantIndices[ i ];
				if( ContravariantIndices[ i ] )
				{
					++ContravariantRank;
					--CovariantRank;
				}
			}
		}
		else
		{
			ContravariantIndices = new bool[ Rank ];
			for( unsigned int i = 0; i < Rank; ++i ) ContravariantIndices[ i ] = false;
		}

		Size = 1;
		for( unsigned int i = 0; i < Rank; ++i ) Size *= d;

		Components = new Function<d>*[ Size ];
		for( unsigned int i = 0; i < Size; ++i ) Components[ i ] = new Function<d>( defaultValue );
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

	TensorField<d>( const TensorField& t )
	{
		Rank = t.Rank;
		CovariantRank = t.CovariantRank;
		ContravariantRank = t.ContravariantRank;
		Size = t.Size;

		if( t.ContravariantIndices )
		{
			ContravariantIndices = new bool[ Rank ];
			for( unsigned int i = 0; i < Rank; ++i ) ContravariantIndices[ i ] = t.ContravariantIndices[ i ];
		}
		else ContravariantIndices = 0;

		Components = new Function<d>*[ Size ];
		for( unsigned int i = 0; i < Size; ++i ) Components[ i ] = new Function<d>( *( t.Components[ i ] ) );
	}

	~TensorField<d>( )
	{
		if( ContravariantIndices ) delete[ ] ContravariantIndices;
		for( unsigned int i = 0; i < Size; ++i ) delete Components[ i ];
		delete[ ] Components;
	}

	TensorField<d>& operator =( const TensorField<d>& t )
	{
		Rank = t.Rank;
		CovariantRank = t.CovariantRank;
		ContravariantRank = t.ContravariantRank;
		Size = t.Size;

		if( t.ContravariantIndices )
		{
			ContravariantIndices = new bool[ Rank ];
			for( unsigned int i = 0; i < Rank; ++i ) ContravariantIndices[ i ] = t.ContravariantIndices[ i ];
		}
		else ContravariantIndices = 0;

		Components = new Function<d>*[ Size ];
		for( unsigned int i = 0; i < Size; ++i ) Components[ i ] = new Function<d>( *( t.Components[ i ] ) );

		return *this;
	}

	bool operator ==( const TensorField<d>& T ) const
	{
		if( Rank != T.Rank ) return false;

		for( unsigned int i = 0; i < Rank; ++i )
		{
			if( ContravariantIndices[ i ] != T.ContravariantIndices[ i ] ) return false;
		}

		for( unsigned int i = 0; i < Size; ++i )
		{
			if( *( Components[ i ] ) != *( T.Components[ i ] ) ) return false;
		}

		return true;
	}

	bool operator !=( const TensorField<d>& T ) const
	{
		return !( *this == T );
	}

	TensorField<d> operator +( const TensorField<d>& t ) const
	{
		if( t.Rank != Rank ) return TensorField<d>( Rank );
		if( t.ContravariantRank != ContravariantRank ) return TensorField<d>( Rank );
		if( t.CovariantRank != CovariantRank ) return TensorField<d>( Rank );

		TensorField<d> ret = *this;

		for( unsigned int i = 0; i < Size; ++i )
			*( ret.Components[ i ] ) += *( t.Components[ i ] );

		return ret;
	}

	TensorField<d> operator -( const TensorField<d>& t ) const
	{
		if( t.Rank != Rank ) return TensorField<d>( Rank );
		if( t.ContravariantRank != ContravariantRank ) return TensorField<d>( Rank );
		if( t.CovariantRank != CovariantRank ) return TensorField<d>( Rank );

		TensorField<d> ret = *this;

		for( unsigned int i = 0; i < Size; ++i )
			*( ret.Components[ i ] ) -= *( t.Components[ i ] );

		return ret;
	}

	// tensor product
	TensorField<d> operator *( const TensorField<d>& t ) const
	{
		// create new index list
		unsigned int NewRank = Rank + t.Rank;
		bool* Indices = new bool[ NewRank ];
		unsigned int i = 0;
		for( unsigned int j = 0; j < Rank; ++j )
		{
			Indices[ i ] = ContravariantIndices[ j ];
			++i;
		}

		for( unsigned int j = 0; j < t.Rank; ++j )
		{
			Indices[ i ] = t.ContravariantIndices[ j ];
			++i;
		}

		TensorField<d> Ret = TensorField<d>( NewRank, Indices );

		delete[ ] Indices;

		// NOTE: repositioning these allows this function to be parallelised
		// out is marginally faster when single-threaded
		unsigned int* Indices1 = new unsigned int[ Rank ];
		unsigned int* Indices2 = new unsigned int[ t.Rank ];
		unsigned int* Indices3 = new unsigned int[ NewRank ];

		//#pragma omp parallel for
		for( unsigned int i = 0; i < Size; ++i )
		{
			unsigned int k = i;
			// work out indices for first tensor
			for( unsigned int j = 0; j < Rank; ++j )
			{
				Indices1[ j ] = k % d;
				k = k / d;
			}

			// iterate over second tensor
			for( unsigned int j = 0; j < t.Size; ++j )
			{
				k = j;
				// work out indices for second tensor
				for( unsigned int l = 0; l < t.Rank; ++l )
				{
					Indices2[ l ] = k % d;
					k = k / d;
				}

				// create indices for return tensor
				k = 0;
				for( unsigned int l = 0; l < Rank; ++l )
				{
					Indices3[ k ] = Indices1[ l ];
					++k;
				}

				for( unsigned int l = 0; l < t.Rank; ++l )
				{
					Indices3[ k ] = Indices2[ l ];
					++k;
				}

				// workout product of Components and store            
				Ret.FunctionFromIndices( Indices3 ) = *( Components[ i ] ) * t.FunctionFromIndices( Indices2 );
			}
		}

		delete[ ] Indices1;
		delete[ ] Indices2;
		delete[ ] Indices3;

		return Ret;
	}

	TensorField<d> operator +( const Function<d>& F ) const
	{
		TensorField<d> ret = *this;

		for( unsigned int i = 0; i < Size; ++i )
			*( ret.Components[ i ] ) = *( ret.Components[ i ] ) + F;

		return ret;
	}

	TensorField<d> operator -( const Function<d>& F ) const
	{
		TensorField<d> ret = *this;

		for( unsigned int i = 0; i < Size; ++i )
			*( ret.Components[ i ] ) = *( ret.Components[ i ] ) - F;

		return ret;
	}

	TensorField<d> operator *( const Function<d>& F ) const
	{
		TensorField<d> ret = *this;

		for( unsigned int i = 0; i < Size; ++i )
			*( ret.Components[ i ] ) = *( ret.Components[ i ] )*F;

		return ret;
	}

	TensorField<d> operator /( const Function<d>& F ) const
	{
		TensorField<d> ret = *this;

		for( unsigned int i = 0; i < Size; ++i )
			*( ret.Components[ i ] ) = *( ret.Components[ i ] ) / F;

		return ret;
	}

	friend TensorField<d> operator +( const Function<d>& F, const TensorField<d>& T )
	{
		return T + F;
	}

	friend TensorField<d> operator -( const Function<d>& F, const TensorField<d>& T )
	{
		return TensorField<d>( T.Rank, T.ContravariantIndices, F ) - T;
	}

	friend TensorField<d> operator *( const Function<d>& F, const TensorField<d>& T )
	{
		return TensorField<d>( T )*F;
	}

	TensorField<d> operator ,( NewIndex _haxpile ) const
	{
		bool* b = new bool[ Rank + 1 ];
		if( ContravariantIndices )
		{
			for( unsigned int i = 0; i < Rank; ++i ) b[ i ] = ContravariantIndices[ i ];
		}
		else
		{
			for( unsigned int i = 0; i < Rank; ++i ) b[ i ] = false;
		}
		b[ Rank ] = false;
		TensorField<d> ret = TensorField<d>( Rank + 1, b, 0.0 );
		unsigned int Indices[ d ] = { 0 };
		for( unsigned int i = 0; i < d; ++i ) Indices[ i ] = 0;
		while( Indices[ 0 ] < d )
		{
			for( unsigned int i = 0; i < d; ++i )
			{
				ret.ReferenceFromIndices( Indices, Rank )[ i ] = FunctionFromIndices( Indices ).Derivative( Vector<d>::GetCoordinate( i ) );
			}

			++Indices[ Rank - 1 ];

			for( unsigned int i = Rank - 1; i > 0; --i )
			{
				if( Indices[ i ] >= d )
				{
					Indices[ i ] = 0;
					if( i > 0 ) ++Indices[ i - 1 ];
				}
			}
		}

		return ret;
	}

	TensorReference<d> operator []( unsigned int i )
	{
		unsigned int* p = new unsigned int;
		*p = i;
		TensorReference<d> ret = TensorReference<d>( this, 1, p );
		delete p;
		return ret;
	}

	TensorField<d> operator ()( Vector<d> Parameter )
	{
		TensorField<d> Ret = TensorField<d>( *this );
		for( unsigned int i = 0; i < Size; ++i ) *( Ret.Components[ i ] ) = *( Ret.Components[ i ] )( Parameter );
		return Ret;
	}

	unsigned int GetRank( ) const { return Rank; }

	Function<d>& FirstEntry( ) const
	{
		return *( Components[ 0 ] );
	}

	TensorReference<d> ReferenceFromIndices( unsigned int* P, unsigned int D )
	{
		return TensorReference<d>( this, D, P );
	}

	Function<d>& FunctionFromIndices( unsigned int* P ) const
	{
		int j = 0;
		int k = 1;
		for( unsigned int i = 0; i < Rank; ++i )
		{
			j += k*P[ i ];
			k *= d;
		}
		return *( Components[ j ] );
	}

	TensorField<d> ChangeIndices( int* Indices )
	{
		// bad indices so return this
		if( !Indices ) return *this;

		// any repeated indices and we also return this
		for( unsigned int i = 0; i < Rank; ++i )
		{
			for( unsigned int j = 0; j < Rank; ++j )
			{
				if( i != j )
				{
					if( Indices[ i ] == Indices[ j ] ) return *this;
				}
			}
		}

		TensorField<d> Ret = TensorField<d>( Rank, ContravariantIndices );

		unsigned int* Indices1 = new unsigned int[ Rank ];
		unsigned int* Indices2 = new unsigned int[ Rank ];

		// fill out new tensor with new indices
		for( unsigned int i = 0; i < Size; ++i )
		{
			int k = i;
			// work out reference for new tensor and old
			for( unsigned int j = 0; j < Rank; ++j )
			{
				Indices1[ j ] = k % d;
				Indices2[ Indices[ j ] ] = Indices1[ j ];
				k = k / d;
			}

			// assign value to new tensor, with switched indices
			Ret.FunctionFromIndices( Indices2 ) = FunctionFromIndices( Indices1 );
		}

		delete[ ] Indices1;
		delete[ ] Indices2;

		return Ret;
	}

	TensorField<d> Contract( unsigned int Index1, unsigned int Index2, bool ForceIndexMatching = true )
	{
		if( ForceIndexMatching && ( ContravariantIndices[ Index1 ] == ContravariantIndices[ Index2 ] ) ) return *this;
		if( ( ( CovariantRank < 1 ) && ( ContravariantRank < 1 ) ) || ( Index1>Rank ) || ( Index2>Rank ) || ( Index1 == Index2 ) ) return *this;

		int NewRank = Rank - 2;

		// create new list of indices, with the contracted pair removed
		bool* Indices = new bool[ NewRank ];

		int j = 0;
		for( unsigned int i = 0; i < Rank; ++i )
		{
			if( ( i != Index1 ) && ( i != Index2 ) )
			{
				Indices[ j ] = ContravariantIndices[ i ];
				++j;
			}
		}

		TensorField<d> Ret = TensorField<d>( NewRank, Indices );

		delete[ ] Indices;

		unsigned int* Indices1 = new unsigned int[ NewRank ];
		unsigned int* Indices2 = new unsigned int[ Rank ];

		// fill out the contracted tensor
		for( unsigned int i = 0; i < Size; ++i )
		{
			// work out first index of dest and source
			int k = i;
			Indices2[ 0 ] = k % d;

			int l = 0;
			if( ( Index1 != 0 ) && ( Index2 != 0 ) )
			{
				Indices1[ l ] = k % d;
				++l;
			}

			// work out remaining indices
			for( unsigned int j = 1; j < Rank; ++j )
			{
				k = k / d;
				if( ( Index1 != j ) && ( Index2 != j ) )
				{
					Indices1[ l ] = k % d;
					++l;
				}
				Indices2[ j ] = k % d;
			}

			// sum over indices when the components are equal, otherwise do nothing
			if( Indices2[ Index1 ] == Indices2[ Index2 ] )
			{
				Ret.FunctionFromIndices( Indices1 ) = Ret.FunctionFromIndices( Indices1 ) + this->FunctionFromIndices( Indices2 );
			}
		}

		delete[ ] Indices1;
		delete[ ] Indices2;

		return Ret;
	}

	TensorField<d> RaiseIndex( unsigned int i, const TensorField<d>& Metric ) const
	{
		// if there are no indices to raise, return unmodified
		if( this->CovariantRank == 0 ) return ( *this );
		// make sure metric provided is contravariant and rank-2
		if( Metric.Rank != 2 ) return ( *this );
		if( Metric.ContravariantRank != 2 ) return ( *this );

		// make tensor product
		TensorField<d> Ret = ( *this ) * Metric;

		// contract the product to get the new tensor
		Ret = Ret.Contract( i, Ret.Rank - 2 );

		return Ret;
	}

	TensorField<d> LowerIndex( unsigned int i, const TensorField<d>& Metric ) const
	{
		// if there are no indices to lower, return unmodified
		if( this->ContrvariantRank == 0 ) return ( *this );
		// make sure metric provided is covariant and rank-2
		if( Metric.Rank != 2 ) return ( *this );
		if( Metric.CovariantRank != 2 ) return ( *this );

		// make tensor product
		TensorField<d> Ret = ( *this ) * Metric;

		// contract the product to get the new tensor
		Ret = Ret.Contract( i, Ret.Rank - 2 );

		return Ret;
	}

	static TensorField<d> Delta( )
	{
		bool Indices[ 2 ] = { true, false };
		TensorField<d> Ret = TensorField<d>( 2, Indices );
		for( unsigned int i = 0; i < d; ++i ) Ret[ i ][ i ] = 1;
		return Ret;
	}

	static TensorField<d> Eta( )
	{
		TensorField<d> Ret = TensorField<d>( 2 );
		Ret[ 0 ][ 0 ] = -1;
		for( unsigned int i = 1; i < d; ++i ) Ret[ i ][ i ] = 1;
		return Ret;
	}

	friend TensorField<d> Gamma( const TensorField<d>& Metric )
	{
		jProfiler profiler = jProfiler( L"Gamma" );

		TensorField<d> Ret = TensorField<d>( 3 );

		// do some quick checks to make sure the Metric is likely good
		// if its not then return some zeros
		if( Metric.Rank != 2 ) return Ret;
		if( Metric.CovariantRank != 2 ) return Ret;

		// get comma derivative of Metric
		TensorField<d> g0 = Metric;
		g0 = ( g0, NewIndex( ) );
		// make versions with different order of indices
		int Indices1[ ] = { 0, 2, 1 };
		int Indices2[ ] = { 1, 2, 0 };
		TensorField<d> g1 = g0.ChangeIndices( Indices1 );
		TensorField<d> g2 = g0.ChangeIndices( Indices2 );

		// compute Christoffel symbol
		return 0.5*( g0 + g1 - g2 );
	}

	friend TensorField<d> Riemann( const TensorField<d>& Metric )
	{
		jProfiler profiler = jProfiler( L"Riemann" );

		TensorField<d> Ret = TensorField<d>( 4 );

		// do some quick checks to make sure the Metric is likely good
		// if its not then return some zeros
		if( Metric.Rank != 2 ) return Ret;
		if( Metric.CovariantRank != 2 ) return Ret;

		// find raised metric
		TensorField<d> gu = InvertMetric( Metric );

		// create second comma derivative of metric
		TensorField<d> g0 = ( ( Metric, NewIndex( ) ), NewIndex( ) );

		// get Christoffel symbol of the first kind (L looks like upside down gamma)
		TensorField<d> L1 = Gamma( Metric );

		// get Christoffel symbol of the first kind by raising an index with the metric
		TensorField<d> L2 = L1.RaiseIndex( 0, gu );

		// create the four different orderings of the second derivatives of the metric
		int Indices0[ ] = { 0, 2, 1, 3 };
		int Indices1[ ] = { 1, 2, 0, 3 };
		int Indices2[ ] = { 0, 3, 1, 2 };
		int Indices3[ ] = { 1, 3, 0, 2 };

		TensorField<d> g1 = g0.ChangeIndices( Indices1 );
		TensorField<d> g2 = g0.ChangeIndices( Indices2 );
		TensorField<d> g3 = g0.ChangeIndices( Indices3 );
		g0 = g0.ChangeIndices( Indices0 );

		Ret = 0.5*( g0 - g1 - g2 + g3 );

		// create the two orderings of the products of the two different kinds of Christoffel symbols
		TensorField<d> L = L1*L2;

		// contract across 0 and 3

		L.Contract( 0, 3 );

		L1 = L.ChangeIndices( Indices0 );
		L2 = L.ChangeIndices( Indices2 );

		return Ret + L1 - L2;
	}

	/*
	this algorithm redcues the matrix to reduced row echelon form in one pass allowing the inverse to be found if it exists

	the idea is that we start with the first column swap the rows about so that the top left element is filled
	we then examine the other rows and remove multiples of each other so that they all start with zero
	*/

	friend TensorField<d> InvertMetric( const TensorField<d>& Metric )
	{
		// invert the indices, this holds true for all cases since we want g_ab -> g^ab, g^ab -> g_ab and g^a_b = g_a^b = delta^a_b
		bool Indices[ 2 ] = { !Metric.ContravariantIndices[ 0 ], !Metric.ContravariantIndices[ 1 ] };

		// create zeros to return if we have a bad metric
		TensorField<d> Zeros = TensorField<d>( 2, Indices );

		// check metric has correct rank
		if( Metric.Rank != 2 ) return Zeros;

		TensorField<d> Ret = TensorField<d>::Delta( );
		TensorField<d> T = Metric;

		// check if metric indices are equal and check g^a_b = delta^a_b
		if( T.ContravariantIndices[ 0 ] != T.ContravariantIndices[ 1 ] )
		{
			if( T != Ret ) return Zeros;

			Ret.ContravariantIndices[ 0 ] = Indices[ 0 ];
			Ret.ContravariantIndices[ 1 ] = Indices[ 1 ];

			return Ret;
		}

		// check each row and each column to see if it is all zeros and spit out a pile of zeroes if it is
		for( unsigned int i = 0; i < d; ++i )
		{
			unsigned int NumZeros = 0;
			for( unsigned int j = 0; j < d; ++j )
			{
				if( T[ i ][ j ] == 0.0 ) ++NumZeros;
			}
			if( NumZeros == d ) return Zeros;

			NumZeros = 0;
			for( unsigned int j = 0; j < d; ++j )
			{
				if( T[ j ][ i ] == 0.0 ) ++NumZeros;
			}
			if( NumZeros == d ) return Zeros;
		}

		// done with basic consistency checks

		// this array keeps track of the row swapping to save having to swap rows round in memory
		unsigned int Rows[ d ] = { 0 };
		for( unsigned int i = 0; i < d; ++i ) Rows[ i ] = i;

		// examine d x d then d x (d-1) then d x (d-2) matrices
		Function<d> F = 0.0;
		unsigned int u = 0;
		for( unsigned int c = 0; c < d; ++c )
		{
			// make pivot element non zero
			u = c;
			while( T[ Rows[ u ] ][ c ] == 0.0 )
			{
				T[ Rows[ u ] ][ c ]->DebugPrint( true );
				++u;

			}
			if( c != u )
			{
				Swap( Rows[ c ], Rows[ u ] );
			}

			// remove leading entry in each row by row subtraction with the relevant row
			for( unsigned int r = 0; r < c; ++r )
			{
				F = T[ Rows[ r ] ][ c ] / T[ Rows[ c ] ][ c ];
				for( unsigned int s = c; s < d; ++s )
				{
					T[ Rows[ r ] ][ s ] = T[ Rows[ r ] ][ s ] - F*T[ Rows[ c ] ][ s ];
					Ret[ Rows[ r ] ][ s ] = Ret[ Rows[ r ] ][ s ] - F*Ret[ Rows[ c ] ][ s ];
				}
			}
			// skip r == c
			for( unsigned int r = c + 1; r < d; ++r )
			{
				F = T[ Rows[ r ] ][ c ] / T[ Rows[ c ] ][ c ];
				for( unsigned int s = c; s < d; ++s )
				{
					T[ Rows[ r ] ][ s ] = T[ Rows[ r ] ][ s ] - F*T[ Rows[ c ] ][ s ];
					Ret[ Rows[ r ] ][ s ] = Ret[ Rows[ r ] ][ s ] - F*Ret[ Rows[ c ] ][ s ];
				}
			}

			// normalise the pivot row
			F = 1 / T[ Rows[ c ] ][ c ];
			for( unsigned int s = c; s < d; ++s )
			{
				T[ Rows[ c ] ][ s ] = F*T[ Rows[ c ] ][ s ];
				Ret[ Rows[ c ] ][ s ] = F*Ret[ Rows[ c ] ][ s ];
			}

			/*
			T.DebugPrint(true);
			printf("\r\n");
			Ret.DebugPrint(true);
			//*/
		}

		// put rows in correct places and normalise
		T = Ret;
		for( unsigned int r = 0; r < d; ++r )
		{
			for( unsigned int c = 0; c < d; ++c )
			{
				Ret[ r ][ c ] = T[ Rows[ r ] ][ c ].GetFunction( );
			}
		}

		return Ret;
	}

	void DebugPrint( bool txyz = false )
	{
		if( Rank == 0 )
		{
			( *( Components[ 0 ] ) ).DebugPrint( txyz );
		}
		else if( Rank == 1 )
		{
			printf( "(" );
			for( unsigned int i = 0; i < d - 1; ++i )
			{
				( *( Components[ i ] ) ).DebugPrint( txyz );
				printf( ", " );
			}
			( *( Components[ d - 1 ] ) ).DebugPrint( txyz );
			printf( ")" );
		}
		else if( Rank == 2 )
		{
			printf( "[\r\n" );
			for( unsigned int j = 0; j < d; ++j )
			{
				printf( "(" );
				for( unsigned int i = 0; i < d - 1; ++i )
				{
					( *this )[ i ][ j ]->DebugPrint( txyz );
					printf( ", " );
				}
				( *this )[ d - 1 ][ j ]->DebugPrint( txyz );
				printf( ")" );
				if( j < ( d - 1 ) ) printf( ",\r\n" );
			}
			printf( "\r\n]" );
		}
		else if( Rank == 3 )
		{
			printf( "{\r\n" );
			for( unsigned int k = 0; k < d; ++k )
			{
				printf( " [\r\n" );
				for( unsigned int j = 0; j < d; ++j )
				{
					printf( " (" );
					for( unsigned int i = 0; i < d - 1; ++i )
					{
						( *this )[ i ][ j ][ k ]->DebugPrint( txyz );
						printf( ", " );
					}
					( *this )[ d - 1 ][ j ][ k ]->DebugPrint( txyz );
					printf( ")" );
					if( j < ( d - 1 ) ) printf( ",\r\n" );
				}
				printf( "\r\n ]" );
				if( k < ( d - 1 ) ) printf( ",\r\n" );
			}
			printf( "\r\n}" );
		}
		else if( Rank == 4 )
		{
			printf( "<\r\n" );
			for( unsigned int l = 0; l < d; ++l )
			{
				printf( " {\r\n" );
				for( unsigned int k = 0; k < d; ++k )
				{
					printf( "  [\r\n" );
					for( unsigned int j = 0; j < d; ++j )
					{
						printf( "  (" );
						for( unsigned int i = 0; i < d - 1; ++i )
						{
							( *this )[ i ][ j ][ k ][ l ]->DebugPrint( txyz );
							printf( ", " );
						}
						( *this )[ d - 1 ][ j ][ k ][ l ]->DebugPrint( txyz );
						printf( ")" );
						if( j < ( d - 1 ) ) printf( ",\r\n" );
					}
					printf( "\r\n  ]" );
					if( k < ( d - 1 ) ) printf( ",\r\n" );
				}
				printf( "\r\n }" );
				if( l < ( d - 1 ) ) printf( ",\r\n" );
			}
			printf( "\r\n>" );
		}
	}
};

#endif