#ifndef __cplusplus
#error TensorField and the classes it depend on require C++ functionality
#endif

#ifndef TENSOR_H
#define TENSOR_H

#define debuggy printf("line %d, %s\r\n", __LINE__, __FILE__)

template <int d> class Vector;
template <int d> class Function;
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
	double	X[d];
public:
	Vector<d>()
	{
		for(unsigned int i = 0; i < d; ++i) X[i] = 0.0;
	}
	
	Vector<d>(const double& D)
	{
		for(unsigned int i = 0; i < d; ++i) X[i] = D;
	}
	
	double& operator [](unsigned int i)
	{
		return X[i];
	}
	
	/*
	template <int c>
	static double Component(Vector<d> Parameter) { return Parameter[c]; }
	*/
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
		return Composition(Parameter);
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
	
public:
	Function<d>() : C(0.0), VF(vfNull), F(fConstant), Child1(0), Child2(0), Composition(cNull), CompositionParameter(0) {}
	Function<d>(const double& d) : C(d), VF(vfNull), F(fConstant), Child1(0), Child2(0), Composition(cNull), CompositionParameter(0) {}
	Function<d>(double (*f)(Vector<d> Parameter)) : C(0.0), VF(f), F(fFunction), Child1(0), Child2(0), Composition(cNull), CompositionParameter(0) {}
	Function<d>(const Function<d>& f) : C(f.C), VF(f.VF), F(f.F), Composition(f.Composition), CompositionParameter(f.CompositionParameter)
	{
		if(f.Child1) Child1 = new Function<d>(*(f.Child1));
		else Child1 = 0;
		if(f.Child2) Child2 = new Function<d>(*(f.Child2));
		else Child2 = 0;
	}
	Function<d>(const Function<d>& F1, const Function<d>& F2, double (Function<d>::*c)(Vector<d> Parameter), int i = 0) : C(0), VF(vfNull), F(fComposition), Composition(c), CompositionParameter(i)
	{
		Child1 = new Function<d>(F1);
		Child2 = new Function<d>(F2);
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
		return F(Parameter);
	}
	
	Function<d>& operator +=(const Function<d> &X)
	{
		if((F == fConstant) && (X.F == fConstant))
		{
			C += X.C;
			return *this;
		}
		// TODO: check if we are adding to something already adding to a constant...
		else
		{
			Function<d> f = *this;
			return (*this = Function<d>(f, X, cAdd));
		}
	}
	
	Function<d> operator +(const Function<d> &X) const
	{
		if((F == fConstant) && (X.F == fConstant))
		{
			return Function<d>(C + X.C);
		}
		// TODO: check if we are adding to something already adding to a constant...
		else
		{
			Function<d> f = *this;
			return Function<d>(f, X, cAdd);
		}
	}
	
	Function<d>& operator -=(const Function<d> &X)
	{
		if((F == fConstant) && (X.F == fConstant))
		{
			C -= X.C;
			return *this;
		}
		// TODO: check if we are adding to something already adding to a constant...
		else
		{
			Function<d> f = *this;
			return (*this = Function<d>(f, X, cSub));
		}
	}
	
	Function<d> operator -(const Function<d> &X) const
	{
		if((F == fConstant) && (X.F == fConstant))
		{
			return Function<d>(C - X.C);
		}
		// TODO: check if we are adding to something already adding to a constant...
		else
		{
			Function<d> f = *this;
			return Function<d>(f, X, cSub);
		}
	}
	
	/*
	template <int wrt>
	Function<d> Derivative() const
	{
		if(F == fConstant)
		{
			return 0.0;
		}
		else if(F == fFunction)
		{
			if(VF == Vector<d>::Component<wrt>)
			{
				return 1.0;
			}
		}
		// TODO: rest once this simple bit works, i.e. composites:
		/*
		
		D [a + b]   = D[a] + D[b]			1 + x -> 0 + 1 -> 1
		D [a - b]   = D[a] - D[b]			1 - x -> 0 - 1 -> -1
		D [a * b]   = b*D[a] + a*D[b]		x * x -> x*1 + x*1 -> 2x 		4 * x -> x*0 + 4*1 -> 4
		D [a / b]   = D[a]*b - a*D[b])/(b*b)	can't be fagged... shame I can't use the old school trick of using chain rule with x^-1 and the product rule to find it as D[a * (1/b)]
		
		*/
	}
	*/
};

/*
		tensor field, can also be used to represent individual tensors by using constants for all the components
*/

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
		
		Size = 1;
        for(unsigned int i = 0; i < Rank; ++i) Size *= d;
		
		Components = new Function<d>*[Size];
		for(unsigned int i = 0; i < Size; ++i) Components[i] = new Function<d>(defaultValue);
	}
	
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
	
	TensorField<d>(TensorField& t)
	{
		Rank = t.Rank;
		CovariantRank = t.CovariantRank;
		ContravariantRank = t.ContravariantRank;
		Size = t.Size;
		
		ContravariantIndices = new bool[Rank];
		for(unsigned int i = 0; i < Rank; ++i) ContravariantIndices[i] = t.ContravariantIndices[i];
		
		Components = new Function<d>*[Size];
		for(unsigned int i = 0; i < Size; ++i) Components[i] = new Function<d>(t.Components[i]);
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
};

#endif