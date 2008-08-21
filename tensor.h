#ifndef TENSOR_H
#define TENSOR_H

#include "jProfiler.h"

#pragma warning( disable : 4311 )
#pragma warning( disable : 4313 )

#include "vector.h"

#define DIFF_STEP 0.0001

#include <stdio.h>

template <class, int> class tensor;

/*

    dummy class for comma derivatives
    probably will come in handy later for other things

*/

class newIndex
{
private:
public:
    newIndex() { }    
    ~newIndex() { }
};

/*

    tensor component class

    contains each component of the tensor as a constant or function

*/

template <class c = double, int n = 4>
class tensorComponent
{
private:
    c                        (*function)(vector<c,n>& p);      // pointer to function if needed
    c                        constant;                        // constant value (for speed)
    // to allow operators on tensors we have here two children to allow unary or binary operations between existing components
    // these are copies instead of "real" pointers, so that values don't change
    tensorComponent<c,n>*        child1;
    tensorComponent<c,n>*        child2;
    // index for calculus functionality
    int                          wrt;
        
    // function pointer for flow control optimisation
    c                        (*retf)(c val, c (*func)(vector<c,n>& p), vector<c,n> parm, tensorComponent<c,n>* a, tensorComponent<c,n>* b, c (*composition)(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm2, int w), int w);
    
    // function pointer for functions of functions
    c                        (*comp)(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w);
    
    // "functions of functions" functions (functional functions?)
    static c c_null(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        return (c)0;
    }
    
    static c c_add(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        return (*c1)(parm) + (*c2)(parm);
    }
    
    static c c_sub(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        return (*c1)(parm) - (*c2)(parm);
    }
    
    static c c_neg(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        return -(*c1)(parm);
    }
    
    static c c_mul(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        return (*c1)(parm) * (*c2)(parm);
    }
    
    static c c_div(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        return (*c1)(parm) / (*c2)(parm);
    }
    
    // derivative (differential)
    static c c_dif(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm, int w)
    {
        vector<c,n> h = vector<c,n>();
        h[w] = DIFF_STEP;
        c ih = (c)0.5 / h[w];
        
        return ((*c1)(parm + h) - (*c1)(parm - h))*ih;
    }
    
    // static helper that returns a value passed into it
    static c ret_val(c val, c (*func)(vector<c,n>& p), vector<c,n> parm, tensorComponent<c,n>* a, tensorComponent<c,n>* b, c (*composition)(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm2, int w), int w)
    {
        return val;
    }
    
    // static helper that returns a value from a function and parameter passed into it
    static c ret_func(c val, c (*func)(vector<c,n>& p), vector<c,n> parm, tensorComponent<c,n>* a, tensorComponent<c,n>* b, c (*composition)(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm2, int w), int w)
    {
        return func(parm);
    }
    
    static c ret_comp(c val, c (*func)(vector<c,n>& p), vector<c,n> parm, tensorComponent<c,n>* a, tensorComponent<c,n>* b, c (*composition)(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm2, int w), int w)
    {
        return composition(a, b, parm, w);
    }
public:
    // default constructor
    tensorComponent<c,n>(c value = (c)0)
    {
        // default to initialisation value
        retf = ret_val;
        constant = value;
        comp = c_null;
        child1 = 0;
        child2 = 0;
        function = 0;
        wrt = 0;
    }
    
    // constructor for compositions
    tensorComponent<c,n>(tensorComponent<c,n>& a, tensorComponent<c,n>& b, c (*composition)(tensorComponent<c,n>* c1, tensorComponent<c,n>* c2, vector<c,n> parm2, int wrt), int index = 0)
    {
        // set up as needed
        retf = ret_comp;
        comp = composition;
        child1 = new tensorComponent<c,n>(a);
        child2 = new tensorComponent<c,n>(b);
        constant = 0;
        function = 0;
        wrt = index;
    }
    
    // copy constructor
    tensorComponent<c,n>(tensorComponent<c,n>& copy)
    {
        function = copy.function;
        constant = copy.constant;
        retf = copy.retf;
        comp = copy.comp;
        if(copy.child1) child1 = new tensorComponent<c,n>(*(copy.child1));
        else child1 = 0;
        if(copy.child2) child2 = new tensorComponent<c,n>(*(copy.child2));
        else child2 = 0;
        wrt = copy.wrt;
    }

    // destructor
    ~tensorComponent<c,n>()
    {
        if(child1) delete child1;
        if(child2) delete child2;
    }
    
    // assignment
    tensorComponent<c,n>& operator =(tensorComponent<c,n>& copy)
    {
        function = copy.function;
        constant = copy.constant;
        retf = copy.retf;
        comp = copy.comp;
        if(child1) delete child1;
        if(child2) delete child2;
        if(copy.child1) child1 = new tensorComponent<c,n>(*(copy.child1));
        else child1 = 0;
        if(copy.child2) child2 = new tensorComponent<c,n>(*(copy.child2));
        else child2 = 0;
        wrt = copy.wrt;
        
        return (*this);
    }
    
    // assignment from a constant value
    tensorComponent<c,n>& operator =(c value)
    {
        constant = value;
        retf = ret_val;
        comp = c_null;
        child1 = 0;
        child2 = 0;
        function = 0;
        wrt = 0;
        return (*this);
    }
    
    // assignment from a function
    tensorComponent<c,n>& operator =(c (*func)(vector<c,n>& p))
    {
        function = func;
        retf = ret_func;
        constant = 0;
        comp = c_null;
        child1 = 0;
        child2 = 0;
        wrt = 0;
        return (*this);
    }
    
    // component access
    c operator ()(vector<c,n> position)
    {
        return retf(constant, function, position, child1, child2, comp, wrt);
    }
    
    tensorComponent<c,n> operator +(tensorComponent<c,n>& p)
    {
        return tensorComponent<c,n>(*this, p, c_add);
    }
    
    tensorComponent<c,n> operator +(c val)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(*this, sub, c_add);
    }
    
    friend tensorComponent<c,n> operator +(c val, tensorComponent<c,n>& t)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(sub, t, tensorComponent<c,n>::c_add);
    }
    
    tensorComponent<c,n> operator +=(tensorComponent<c,n>& p)
    {
        *this = *this + p;
        return *this;
    }
    
    tensorComponent<c,n> operator -(tensorComponent<c,n>& p)
    {
        return tensorComponent<c,n>(*this, p, c_sub);
    }
    
    tensorComponent<c,n> operator -(c val)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(*this, sub, c_sub);
    }
    
    friend tensorComponent<c,n> operator -(c val, tensorComponent<c,n>& t)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(sub, t, tensorComponent<c,n>::c_sub);
    }
    
    tensorComponent<c,n> operator -=(tensorComponent<c,n>& p)
    {
        *this = *this - p;
        return *this;
    }
    
    tensorComponent<c,n> operator -()
    {
        return tensorComponent<c,n>(*this, *this, c_neg);
    }
    
    // no *= to avoid confusion with tensorReference... besides, not a big deal anyway
    tensorComponent<c,n> operator *(tensorComponent<c,n>& p)
    {
        return tensorComponent<c,n>(*this, p, c_mul);
    }
    
    tensorComponent<c,n> operator *(c val)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(*this, sub, c_mul);
    }
    
    friend tensorComponent<c,n> operator *(c val, tensorComponent<c,n>& t)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(sub, t, tensorComponent<c,n>::c_mul);
    }
     
    tensorComponent<c,n> operator /(tensorComponent<c,n>& p)
    {
        return tensorComponent<c,n>(*this, p, c_div);
    }
    
    tensorComponent<c,n> operator /(c val)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(*this, sub, c_div);
    }
    
    friend tensorComponent<c,n> operator /(c val, tensorComponent<c,n>& t)
    {
        tensorComponent<c,n> sub = tensorComponent<c,n>(val);
        return tensorComponent<c,n>(sub, t, tensorComponent<c,n>::c_div);
    }
    
    // using , for derivative
    tensorComponent<c,n> operator ,(int i)
    {
        return tensorComponent<c,n>(*this, *this, c_dif, wrt);
    }
    
    // is this component a function?
    bool isFunction()
    {
        if(child1 && child2)
        {
            return child1->isFunction() || child2->isFunction();
        }
        else if(child1)
        {
            return child1->isFunction();
        }
        return (retf == ret_func);
    }
    
    // is this component a constant?
    bool isConstant()
    {
        if(child1 && child2)
        {
            return child1->isConstant() && child2->isConstant();
        }
        else if(child1)
        {
            return child1->isConstant();
        }
        return (retf == ret_val);
    }
    
    void printDebug()
    {
        char st1[5] = "True";
        char st2[6] = "False";
        printf("Debug info for tensor component at 0x%.8X\n", this);
        printf("  Constant: %s\n", isConstant() ? st1 : st2);    
        printf("  Function: %s\n", isFunction() ? st1 : st2);
        vector<c,n> position = vector<c,n>((c)0);
        printf("  Value with all params 0: %f\n", (*this)(position));
    }
};

/*

    tensor reference class

    allows for using [] to access tensor members without too many problems

*/

template <class c = double, int n = 4>
class tensorReference
{
private:
    tensor<c,n>*            parent;
    int                     numDimensions;
    int*                    refComponents;
    
    // helper function pointers and functions for speed
    c                       (*deref)(tensor<c,n>* p, int* di, vector<c,n> position);
    void                    (*set_c)(tensor<c,n>* p, int* di, c value);
    void                    (*set_f)(tensor<c,n>* p, int* di, c (*func)(vector<c,n>& p));
    void                    (*set_m)(tensor<c,n>* p, int* di, tensorComponent<c,n> a);
    
    static c ref_null(tensor<c,n>* p, int* di, vector<c,n> position)
    {
        return (c)0;
    }

    static void set_null_c(tensor<c,n>* p, int* di, c value)
    {
        return;
    }

    static void set_null_f(tensor<c,n>* p, int* di, c (*func)(vector<c,n>& p))
    {
        return;
    }
    
    static void set_null_m(tensor<c,n>* p, int* di, tensorComponent<c,n> a)
    {
        return;
    }
    
    static c ref_func(tensor<c,n>* p, int* di, vector<c,n> position)
    {
        return p->getValFromIndices(di)(position);
    }

    static void set_func_c(tensor<c,n>* p, int* di, c value)
    {
        p->getValFromIndices(di) = value;
    }

    static void set_func_f(tensor<c,n>* p, int* di, c (*func)(vector<c,n>& p))
    {
        p->getValFromIndices(di) = func;
    }
    
    static void set_func_m(tensor<c,n>* p, int* di, tensorComponent<c,n> a)
    {
        p->getValFromIndices(di) = a;
    }
public:
    // default constructor (annoying)
    tensorReference<c,n>()
    {
        parent = 0;
        numDimensions = 0;
        refComponents = 0;
        deref = 0;
        set_c = 0;
        set_f = 0;
        set_m = 0;
    }
    
    // constructor
    tensorReference<c,n>(tensor<c,n>* p, int ndims, int* dims)
    {
        parent = p;
        numDimensions = ndims;
        refComponents = new int[ndims];
        
        
        for(int i = 0; i < ndims; i++)
        {
            refComponents[i] = dims[i];
        }
        
        // if this is a final reference then we should call the tensor component's function to get a return value
        // for (), otherwise we should return null when () is called
        if(ndims==parent->getTotalRank())
        {
            deref = ref_func;
            set_c = set_func_c;
            set_f = set_func_f;
            set_m = set_func_m;
        }
        else
        {
            deref = ref_null;
            set_c = set_null_c;
            set_f = set_null_f;
            set_m = set_null_m;
        }
    }
    
    // copy constructor
    tensorReference<c,n>(tensorReference<c,n>& copy)
    {
        parent = copy.parent;
        numDimensions = copy.numDimensions;
        refComponents = new int[numDimensions];
        
        
        for(int i = 0; i < numDimensions; i++)
        {
            refComponents[i] = copy.refComponents[i];
        }
        
        deref = copy.deref;
        set_c = copy.set_c;
        set_f = copy.set_f;
        set_m = copy.set_m;
    }
    
    // destructor
    ~tensorReference()
    {    
        if(refComponents) delete[] refComponents;
    }
    
    // direct assignment
    tensorReference<c,n>& operator *=(tensorReference<c,n>& copy)
    {    
        if(refComponents) delete[] refComponents;
        
        parent = copy.parent;
        numDimensions = copy.numDimensions;    
        refComponents = new int[numDimensions];
        
        
        for(int i = 0; i < numDimensions; i++)
        {
            refComponents[i] = copy.refComponents[i];
        }

        deref = copy.deref;
        set_c = copy.set_c;
        set_f = copy.set_f;
        set_m = copy.set_m;
        
        return (*this);
    }
    
    // assignment - note this gets used to set the component the reference points to
    // the 'actual' assignment is above as *=
    tensorReference<c,n>& operator =(tensorReference<c,n>& copy)
    {    
        *(*this) = *copy;
        return (*this);
    }

    // set component to value
    tensorReference<c,n>& operator =(c value)
    {
        set_c(parent, refComponents, value);
        return (*this);
    }

    // set component to function
    tensorReference<c,n>& operator =(c (*func)(vector<c,n>& p))
    {
        set_f(parent, refComponents, func);
        return (*this);
    }
    
    // set component to other component
    tensorReference<c,n>& operator =(tensorComponent<c,n> a)
    {
        set_m(parent, refComponents, a);
        return (*this);
    }
    
    // overloaded [] so tensor can be used like a multi-dimensional array
    tensorReference<c,n> operator [](int i)
    {
        if(numDimensions < parent->getTotalRank())
        {
            // create array of indices
            int* p = new int[numDimensions+1];
            
            // copy from old array
            
            for(int j = 0; j < numDimensions; j++)
            {
                p[j] = refComponents[j];
            }
            
            // set last index
            p[numDimensions] = i;
            
            // create return reference and delete array of indices
            tensorReference<c,n>& ret = tensorReference<c,n>(parent,numDimensions+1,p);
            delete[] p;
            
            return ret;
        }
        // extra references should do nothing
        else
        {
            return (*this);
        }
    }
    
    // overload * and -> so this acts like a pointer
    // ** should optimise with function pointer if gets used much
    tensorComponent<c,n>& operator *()
    {
        if(numDimensions==parent->getTotalRank()) return parent->getValFromIndices(refComponents);
        else return parent->getFirstEntry();
    }
    
    tensorComponent<c,n>* operator ->()
    {
        if(numDimensions==parent->getTotalRank()) return &(parent->getValFromIndices(refComponents));
        else return &(parent->getFirstEntry());
    }
    
    // overloaded () to return tensor components
    c operator ()(vector<c,n> position)
    {
        return deref(parent, refComponents, position);
    }
    
    // arithmetic operations that should return a new tensor component
    tensorComponent<c,n> operator +(tensorReference<c,n>& r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) + (*r);
    }
    
    tensorComponent<c,n> operator +(c r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) + r;
    }
    
    tensorComponent<c,n> operator -(tensorReference<c,n>& r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) - (*r);
    }
    
    tensorComponent<c,n> operator -(c r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) - r;
    }
    
    tensorComponent<c,n> operator -()
    {
        return -(*(*this));
    }
    
    tensorComponent<c,n> operator *(tensorReference<c,n>& r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) * (*r);
    }
    
    tensorComponent<c,n> operator *(c r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) * r;
    }
    
    tensorComponent<c,n> operator /(tensorReference<c,n>& r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) / (*r);
    }
    
    tensorComponent<c,n> operator /(c r)
    {
        if(r.numDimensions!=numDimensions) return tensorComponent<c,n>();
        return (*(*this)) / r;
    }
    
    void printDebug()
    {
        printf("Debug info for tensor reference at 0x%.8X\n", this);
        printf("  Dimensions %d\n", numDimensions);    
        printf("  Components 0x%.8X\n", refComponents);
        for(int i = 0; i < numDimensions; i++) printf("  -- Component %d@0x%.8X  = %d\n", i, &(refComponents[i]), refComponents[i]);
    }
};

/*

    tensor class

    the tensor class

*/

template <class c = double, int  n = 4>
class tensor
{
private:
    int                     covariantRank;                    // number of covariant dimensions
    int                     contravariantRank;                // number of contravariant dimensions
    int                     totalRank;
    tensorComponent<c,n>**  components;                        // array of components
    int                     numEntries;                        // number of entries in total
    bool*                   isUpperIndex;                   // which dimensions are contravariant
public:
    // constructor
    tensor<c,n>(int rank = 2, bool* upIndices = 0, c defaultComponentConstant = (c)0)
    {
        totalRank = rank;
        isUpperIndex = new bool[rank];

        
        for(int i = 0; i < rank; i++)
        {
            isUpperIndex[i] = false;
        }

        covariantRank = rank;
        contravariantRank = 0; 

        if(upIndices)
        {
            
            for(int i = 0; i < rank; i++)
            {
                isUpperIndex[i] = upIndices[i];
                if(upIndices[i])
                {
                    covariantRank--;
                    contravariantRank++;
                }
            }
        }
        
        // numEntries = numDimensions ^ totalRank
        // notice that this makes a rank(0,0) tensor into a scalar
        numEntries = 1;
        for(int i = 0; i < totalRank; i++) numEntries *= n;
        
        components = new tensorComponent<c,n>*[numEntries];
        
        for(int i = 0; i < numEntries; i++)
        {
            components[i] = new tensorComponent<c,n>(defaultComponentConstant);
        }
    }
    
    // copy constructor
    tensor<c,n>(tensor<c,n>& copy)
    {
        totalRank = copy.totalRank;
        isUpperIndex = new bool[totalRank];

        
        for(int i = 0; i < totalRank; i++)
        {
            isUpperIndex[i] = copy.isUpperIndex[i];;
        }

        covariantRank = copy.covariantRank;
        contravariantRank = copy.contravariantRank; 
        
        // numEntries = numDimensions ^ totalRank
        // notice that this makes a rank(0,0) tensor into a scalar
        numEntries = copy.numEntries;
        
        components = new tensorComponent<c,n>*[numEntries];
        
        for(int i = 0; i < numEntries; i++)
        {
            components[i] = new tensorComponent<c,n>(*(copy.components[i]));
        }
    }

    ~tensor<c,n>()
    {    
        
        for(int i = 0; i < numEntries; i++)
        {
            delete components[i];
        }
        
        delete[] components;

        delete[] isUpperIndex;
    }

    // assignment
    tensor<c,n>& operator =(tensor<c,n>& copy)
    {
        // clean up memory
        
        for(int i = 0; i < numEntries; i++)
        {
            delete components[i];
        }
        
        delete[] components;

        delete[] isUpperIndex;
        
        totalRank = copy.totalRank;
        isUpperIndex = new bool[totalRank];

        
        for(int i = 0; i < totalRank; i++)
        {
            isUpperIndex[i] = copy.isUpperIndex[i];
        }

        covariantRank = copy.covariantRank;
        contravariantRank = copy.contravariantRank; 
        
        // numEntries = numDimensions ^ totalRank
        // notice that this makes a rank(0,0) tensor into a scalar
        numEntries = copy.numEntries;
        
        components = new tensorComponent<c,n>*[numEntries];
        
        for(int i = 0; i < numEntries; i++)
        {
            components[i] = new tensorComponent<c,n>(*(copy.components[i]));
        }

        return (*this);
    }
    
    // overloaded [] so tensor can be used like a multi-dimensional array
    tensorReference<c,n> operator[](int i)
    {
        int* p = new int;
        *p = i;
        tensorReference<c,n> ret; ret *= tensorReference<c,n>(this,1,p);
        delete p;
        return ret;
    }

    tensor<c,n> operator +(tensor<c,n>& t)
    {
        jProfiler profiler = jProfiler(L"tensor::operator + (tensor sum)");
        // check ranks
        if(!((t.covariantRank == covariantRank) && (t.covariantRank == covariantRank))) return *this;

        tensor<c,n> ret = tensor<c,n>(t);

        #pragma omp parallel for
        for(int i = 0; i < numEntries; i++)
        {
            *(ret.components[i]) = *(components[i]) + *(ret.components[i]);
        }

        return ret;
    }
    
    tensor<c,n> operator +=(tensor<c,n>& t)
    {
        *this = *this + t;
        return (*this);
    }

    tensor<c,n> operator -(tensor<c,n>& t)
    {
        jProfiler profiler = jProfiler(L"tensor::operator - (tensor subtraction)");
        // check ranks
        if(!((t.covariantRank == covariantRank) && (t.covariantRank == covariantRank))) return *this;

        tensor<c,n> ret = tensor<c,n>(t);

        #pragma omp parallel for
        for(int i = 0; i < numEntries; i++)
        {
            *(ret.components[i]) = *(components[i]) - *(ret.components[i]);
        }

        return ret;
    }
    
    tensor<c,n> operator -=(tensor<c,n>& t)
    {
        *this = *this - t;
        return (*this);
    }
    
    tensor<c,n> operator -()
    {
        tensor<c,n> ret = tensor<c,n>(t);

        #pragma omp parallel for
        for(int i = 0; i < numEntries; i++)
        {
            *(ret.components[i]) = -*(components[i]);
        }

        return ret;
    }

    // scalar product
    tensor<c,n> operator *(c r)
    {
        jProfiler profiler = jProfiler(L"tensor::operator * (scalar product)");
        tensor<c,n> ret = tensor<c,n>(totalRank, isUpperIndex);
        
        #pragma omp parallel for
        for(int i = 0; i < numEntries; i++)
        {
            *(ret.components[i]) = (*(components[i])) * r;
        }
        
        return ret;
    }
    
    friend tensor<c,n> operator *(c r, tensor<c,n>& t)
    {
        return t*r;
    }
    
    // tensor product
    tensor<c,n> operator *(tensor<c,n>& t)
    {
        jProfiler profiler = jProfiler(L"tensor::operator * (tensor product)");
        // create new index list
        int newRank = totalRank + t.totalRank;
        bool* ind = new bool[newRank];
        int i =0;
        for(int j = 0; j < totalRank; j++)
        {
            ind[i] = isUpperIndex[j];
            i++;
        }
        
        for(int j = 0; j < t.totalRank; j++)
        {
            ind[i] = t.isUpperIndex[j];
            i++;
        }
        
        tensor<c,n> ret = tensor<c,n>(newRank, ind);
        
        delete[] ind;
        
        // NOTE: repositioning these allows this function to be parallelised
        // out is marginally faster when single-threaded
        int* in1 = new int[totalRank];
        int* in2 = new int[t.totalRank];
        int* in3 = new int[newRank];
        
        //#pragma omp parallel for
        for(int i = 0; i < numEntries; i++)
        {
            // int* in1 = new int[totalRank];
            // int* in2 = new int[t.totalRank];
            // int* in3 = new int[newRank];
            
            int k = i;
            // work out indices for first tensor
            for(int j = 0; j < totalRank; j++)
            {
                in1[j] = k % n;
                k = k / n;
            }
            
            // iterate over second tensor
            for(int j = 0; j < t.numEntries; j++)
            {
                k = j;
                // work out indices for second tensor
                for(int l = 0; l < t.totalRank; l++)
                {
                    in2[l] = k % n;
                    k = k / n;
                }
                
                // create indices for return tensor
                k = 0;
                for(int l = 0; l < totalRank; l++)
                {
                    in3[k] = in1[l];
                    k++;
                }
                
                for(int l = 0; l < t.totalRank; l++)
                {
                    in3[k] = in2[l];
                    k++;
                }
                
                // workout product of components and store            
                ret.getValFromIndices(in3) = *(components[i]) * t.getValFromIndices(in2);
            }

            // delete[] in1;
            // delete[] in2;
            // delete[] in3;
        }
        
        delete[] in1;
        delete[] in2;
        delete[] in3;
        
        return ret;
    }
    
    // comma derivative wrt fixed index
    // this is a rank (r,s) tensor
    tensor<c,n> operator ,(int i)
    {
        jProfiler profiler = jProfiler(L"tensor::operator , (specified index)");
        tensor<c,n> ret = tensor<c,n>(totalRank, isUpperIndex);

        for(int j = 0; j < numEntries; j++)
        {
            *(ret.components[j]) = ((*(components[j])),i);
        }

        return ret;
    }
    
    // comma derivative wrt new index
    // this is a rank (r+1,s) tensor
    tensor<c,n> operator ,(newIndex& pointless_variable_for_nothing)
    {
        jProfiler profiler = jProfiler(L"tensor::operator , (new index)");
        // add new lower index
        int newRank = totalRank + 1;
        bool* ind = new bool[newRank];
        
        for(int i = 0; i < totalRank; i++) ind[i] = isUpperIndex[i];
        ind[newRank-1] = false;
        
        tensor<c,n> ret = tensor<c,n>(newRank, ind);
        
        // iterate through and fill with derivatives
        
        #pragma omp parallel for
        for(int i = 0; i < numEntries; i++)
        {
        
            int k = i;
            tensorReference<c,n> r; r *= ret[k % n];
            // work out reference for new tensor
            for(int j = 1; j < (newRank-1); j++)
            {
                k = k / n;
                r *= r[k % n];
            }
            
            // assign new row/column to new tensor
            for(int j = 0; j < n; j++)
            {
                r[j] = ((*(components[i])),j);
            }
        }
        
        return ret;
    }
    
    // switches order of indices
    tensor<c,n> changeIndices(int* indices)
    {
        jProfiler profiler = jProfiler(L"tensor::changeIndices");
        // bad indices so return this
        if(indices==0) return (*this);
        
        // any repeated indices and we also return this
        for(int i = 0; i < totalRank; i++)
        {
            for(int j = 0; j < totalRank; j++)
            {
                if(i!=j)
                {
                    if(indices[i]==indices[j]) return (*this);
                }
            }
        }
            
        tensor<c,n> ret = tensor<c,n>(totalRank, isUpperIndex);
        
        int* ind1 = new int[totalRank];
        int* ind2 = new int[totalRank];
        
        // fill out new tensor with new indices
        for(int i = 0; i < numEntries; i++)
        {
            int k = i;
            // work out reference for new tensor and old
            for(int j = 0; j < totalRank; j++)
            {
                ind1[j] = k % n;
                ind2[indices[j]] = ind1[j];
                k = k / n;
            }
            
            // assign value to new tensor, with switched indices
            ret.getValFromIndices(ind2) = *(components[i]);
        }
        
        delete[] ind1;
        delete[] ind2;
        
        return ret;
    }
    
    // raises an index using a metric
    tensor<c,n> raiseIndex(int i, tensor<c,n>& metric)
    {
        jProfiler profiler = jProfiler(L"tensor::raiseIndex");
        // if there are no indices to raise, return unmodified
        if(this->covariantRank == 0) return (*this);
        // make sure metric provided is contravariant and rank-2
        if(metric.totalRank != 2) return (*this);
        if(metric.contravariantRank != 2) return (*this);
        
        // make tensor product
        tensor<c,n> ret = (*this) * metric;
        
        // contract the product to get the new tensor
        ret = ret.Contract(i,ret.getTotalRank()-2);

        return ret;
    }
    
    tensor<c,n> lowerIndex(int i, tensor<c,n>& metric)
    {
        jProfiler profiler = jProfiler(L"tensor::lowerIndex");
        // if there are no indices to lower, return unmodified
        if(this->contravariantRank == 0) return (*this);
        // make sure metric provided is covariant and rank-2
        if(metric.totalRank != 2) return (*this);
        if(metric.covariantRank != 2) return (*this);
        
        tensor<c,n> ret = (*this) * metric;
        
        // contract the product to get the new tensor
        ret.Contract(i,ret.getTotalRank()-1);
        return ret;
    }
    
    // raises an index virtually. raising an index to make, e.g. T_uv -> T^u_v requires multiplication and contraction with metric
    // this just forces an index up or down
    tensor<c,n> forceRaiseIndex(int i)
    {
        jProfiler profiler = jProfiler(L"tensor::forceRaiseIndex");
        if(i > totalRank) return *this;
        bool* ind = new bool[totalRank];
        
        for(int j = 0; j < totalRank; j++)
        {
            ind[j] = isUpperIndex[j];
        }
        
        ind[i] = true;
        
        tensor<c,n> ret = tensor<c,n>(totalRank,ind);
        delete[] ind;
        
        for(int j = 0; j < numEntries; j++)
        {
            *(ret.components[j]) = *(components[j]);
        }
        
        return ret;
    }
    
    // does the same as above, except lowering indices
    tensor<c,n> forceLowerIndex(int i)
    {
        jProfiler profiler = jProfiler(L"tensor::forceLowerIndex");
        if(i > totalRank) return *this;
        bool* ind = new bool[totalRank];
        
        for(int j = 0; j < totalRank; j++)
        {
            ind[j] = isUpperIndex[j];
        }
        
        ind[i] = false;
        
        tensor<c,n> ret = tensor<c,n>(totalRank,ind);
        delete[] ind;
        
        for(int j = 0; j < numEntries; j++)
        {
            *(ret.components[j]) = *(components[j]);
        }
        
        return ret;
    }
        
    tensor<c,n> Contract(int index1, int index2, bool forceUpperLowerOnly = true)
    {
        jProfiler profiler = jProfiler(L"tensor::Contract");
        if(forceUpperLowerOnly && (isUpperIndex[index1] == isUpperIndex[index2])) return *this;
        if(((covariantRank < 1) && (contravariantRank < 1)) || (index1>totalRank) || (index2>totalRank) || (index1==index2)) return *this;

        int newRank = totalRank-2;

        // create new list of indices, with the contracted pair removed
        bool* ind = new bool[newRank];

        int j = 0;
        #pragma omp parallel for
        for(int i = 0; i < totalRank; i++)
        {
            if((i!=index1)&&(i!=index2))
            {
                ind[j] = isUpperIndex[i];
                j++;
            }
        }

        tensor<c,n> ret = tensor<c,n>(newRank, ind);

        delete[] ind;

        int* in = new int[newRank];
        int* in2 = new int[totalRank];
        
        // #pragma omp parallel for
        // fill out the contracted tensor
        for(int i = 0; i < numEntries; i++)
        {
            // work out first index of dest and source
            int k = i;
            in2[0] = k % n;
            
            int l = 0;
            if((index1!=0)&&(index2!=0))
            {
                in[l] = k % n;
                l++;
            }
            
            // work out remaining indices
            for(int j = 1; j < totalRank; j++)
            {
                k = k / n;
                if((index1!=j)&&(index2!=j))
                {
                    in[l] = k % n;
                    l++;
                }
                in2[j] = k % n;
            }

            // sum over indices when the components are equal, otherwise do nothing
            if(in2[index1] == in2[index2])
            {
                ret.getValFromIndices(in) = ret.getValFromIndices(in) + (*this).getValFromIndices(in2);
            }
        }

        delete[] in;
        delete[] in2;

        return ret;
    }
     
    int getCovariantRank()
    {
        return covariantRank;
    }
    
    int getContravariantRank()
    {
        return contravariantRank;
    }
    
    int getTotalRank()
    {
        return totalRank;
    }
    
    tensorComponent<c,n>& getValFromIndices(int* p)
    {
        int j = 0;
        int k = 1;
        for(int i = 0; i<totalRank; i++)
        {
            j += k*p[i];
            k *= n;
        }
        return *(components[j]);
    }

    tensorComponent<c,n>& getFirstEntry()
    {
        return *(components[0]);
    }
    
    void printDebug()
    {
        printf("Debug info for tensor at 0x%.8X\n", this);
        printf("  Total Rank %d\n", totalRank);
        printf("  Covariant Rank %d\n", covariantRank);
        printf("  Contravariant Rank %d\n", contravariantRank);
        printf("  Components 0x%.8X\n", components);
        for(int i = 0; i < numEntries; i++)
        {
            printf("  -- Component %d@0x%.8X  with all params 0 = %f\n", i, components[i], (*components[i])(vector<c,n>()));
        }
    }
};

template <class c, int n>
tensor<c,n> RaiseMetricIndices(tensor<c,n>&);
template <class c, int n>
tensor<c,n> LowerMetricIndices(tensor<c,n>&);

/*

    Kronecker delta
    
    (the identity matrix)

*/

template <class c, int n>
tensor<c,n> Delta(int rank = 2)
{
    jProfiler profiler = jProfiler(L"<global>::Delta");
    bool* ind = new bool[rank];
    
    for(int i = 1; i < rank; i++)
    {
        ind[i] = false;
    }

    ind[0] = true;
    
    tensor<c,n> ret = tensor<c,n>(rank,ind);
    
    delete[] ind;

    // fill all indices with 1 where they are equal
    for(int i = 0; i < n; i++)
    {
        tensorReference<c,n> r; r *= ret[i];
        for(int j = 1; j < rank; j++) r *= r[i];
        r = (c)1.0;
    }
    return ret;
}

/*

    Minkowski metric (eta)
    
    (metric for flat spacetime)

    ** could be easily optimised since i copy pasted delta and changed some nick nacks
    
*/

template <class c, int n>
tensor<c,n> Eta()
{
    jProfiler profiler = jProfiler(L"<global>::Eta");
    int rank = 2;
    tensor<c,n> ret = tensor<c,n>(2);

    // fill all indices with 1 where they are equal
    tensorReference<c,n> r;
    
    for(int i = 0; i < n; i++)
    {
        r *= ret[i];
        for(int j = 1; j < rank; j++) r *= r[i];
        r = (c)1.0;
    }
    
    r *= ret[0];
    for(int i = 1; i < rank; i++) r *= r[0];
    
    // -1 in the top left
    r = (c)-1.0;
    
    return ret;
}

/*

    Levi-Civita (epsilon) symbol
    
    (antisymmetric tensor)
    
*/

template <class c, int n>
tensor<c,n> Epsilon(int rank = 3)
{
    jProfiler profiler = jProfiler(L"<global>::Epsilon");
    if(rank < 3) rank = 3
    tensor<c,n> ret = tensor<c,n>(rank);
    
    // make all +ve cyclic permutations into 1
    tensorReference<c,n> r;
    for(int i = 0; i < rank; i++)
    {
            r *= ret[i];
            for(int j = 1; j < rank; j++) r *= r[(i+j)%rank];
            *r = 1;
    }    
    
    // and all negatives into -1
    for(int i = 0; i < rank; i++)
    {
            r *= ret[i];
            for(int j = (rank-1); j > 0; j--) r *= r[(i+j)%rank];
            *r = -1;
    }    
    
    return ret;
}

/*

    Christoffel (Gamma) symbol

    (rank (3,0) pseudo-tensor useful for covariant derivative)

*/

template <class c, int n>
tensor<c,n> Gamma(tensor<c,n>& metric)
{
    jProfiler profiler = jProfiler(L"<global>::Gamma");
    tensor<c,n> ret = tensor<c,n>(3);
    
    // do some quick checks to make sure the metric is likely good
    // if its not then return some zeros
    if(metric.getTotalRank() != 2) return ret;
    if(metric.getCovariantRank() != 2) return ret;
    // we could check if it is symmetric... but probably slow
    
    // get comma derivative of metric
    tensor<c,n> g0 = (metric,newIndex());
    // make versions with different order of indices
    int inds1[] = {0, 2, 1};
    int inds2[] = {1, 2, 0};
    tensor<c,n> g1 = g0.changeIndices(inds1);
    tensor<c,n> g2 = g0.changeIndices(inds2);
    
    // compute Christoffel symbol
    return 0.5*(g0 + g1 - g2);
}

/*
    Riemann (R) curvature tensor
    
    this is a rank-3,1 tensor describing the curvature of the coordinate system
*/
template <class c, int n>
tensor<c,n> Riemann(tensor<c,n>& metric)
{
    jProfiler profiler = jProfiler(L"<global>::Riemann");
    bool inds[] = {false,false,false,false};
    tensor<c,n> ret = tensor<c,n>(4,inds);
    
    // find raised metric
    tensor<c,n> gu = RaiseMetricIndices(metric);
    
    // create second comma derivative of metric
    tensor<c,n> g0 = ((metric,newIndex()),newIndex());
    
    // get Christoffel symbol of the first kind (L looks like upside down gamma)
    tensor<c,n> L1 = Gamma<c,n>(metric);
    
    // get Christoffel symbol of the first kind by raising an index with the metric
    tensor<c,n> L2 = L1.raiseIndex(0,gu);
    
    // create the four different orderings of the second derivatives of the metric
    int inds0[] = {0, 2, 1, 4};
    int inds1[] = {1, 2, 0, 4};
    int inds2[] = {0, 3, 1, 2};
    int inds3[] = {1, 3, 0, 2};
    
    tensor<c,n> g1 = g0.changeIndices(inds1);
    tensor<c,n> g2 = g0.changeIndices(inds2);
    tensor<c,n> g3 = g0.changeIndices(inds3);
    g0 = g0.changeIndices(inds0);
    
    ret = .5 * (g0 - g1 - g2 + g3);
    
    // create the two orderings of the products of the two different kinds of Christoffel symbols
    tensor<c,n> L = L1*L2;
    
    // contract across 0 and 3
    
    L.Contract(0,3);
    
    L1 = L.changeIndices(inds0);
    L2 = L.changeIndices(inds2);
    
    return ret + L1 - L2;
}

// solve equations to raise indices of a metric
template <class c, int n>
tensor<c,n> RaiseMetricIndices(tensor<c,n>& metric)
{
    jProfiler profiler = jProfiler(L"<global>::RaiseMetricIndices");
    bool con[] = {true,true};
    tensor<c,n> ret = tensor<c,n>(2, con);
    
    // do some quick checks to make sure the metric is likely good
    // if its not then return some zeros
    if(metric.getTotalRank() != 2) return ret;
    if(metric.getCovariantRank() != 2) return ret;

    /*
        g^ua g_av = delta^u_v

        this is the same as saying that g^uv is the matrix invert of g_uv
    */

    // solve using some gaussian elimination hax
    tensor<c,n> d = metric;
    ret = Delta<c,n>();
    
    // we treat delta and the metric as a big horizontal matrix nx2n
    
    // reduce metric half of matrix to an identity matrix
    // first put it into a upper triangular form then do pivoting to finish off
    // find rows beginning with zeros and swap them to the bottom
    // if we have any rows which are completely zero, we are fucked
    
    /*
    
        note, for some reason the matrix seems to come out upside down
        the reason is that the rows and columns are stored backwards to 'normal' inside the tensor
        so the list of components lists columns first, then rows
        this is so that the references into and out of the tensor are easier to maintain
        and so that the process is easier to visualise (for me) as the dimensions and rank run off to infinity
    
    */
    int bottommost = 0;
    for(int i = n; i >= 0; i--)
    {
        // check each row for i number of zeros
        for(int j = 0; j < n; j++)
        {
            int numZeros = 0;
            for(int k = 0; k < n; k++)
            {
                if(d[j][k]->isConstant() && (d[j][k](vector<c,n>()) == 0))
                {
                    numZeros++;
                }
                else break;
            }

            if(numZeros == n)
            {
                // fucked
                //printf("annoying zero -- line 1392\n");
                return tensor<c,n>(2,con);
            }

            if(numZeros == i)
            {
                // swap with bottommost row we haven't yet swapped
                bottommost++;
                for(int k = 0; k < n; k++)
                {
                    tensorComponent<c,n> buf;
                    buf = *(d[j][k]);
                    *(d[j][k]) = *(d[n-bottommost][k]);
                    *(d[n-bottommost][k]) = buf;
                    buf = *(ret[j][k]);
                    *(ret[j][k]) = *(ret[n-bottommost][k]);
                    *(ret[n-bottommost][k]) = buf;
                }
            }            
        }
    }

    /*
    
        do some row operations to simplify the 'matrix' into its lower triangular form
    
        the top left is not a zero so see if any of the other rows don't start with zero
        then go back and check rows starting with one zero to see if they need further reducing etc
        
    */
    // iterate over all rows except the last
    for(int m = 0; m < (n-1); m++)
    {
        // if we have a zero in our 'top left' we are also screwed
        if(d[m][m]->isConstant() && (d[m][m](vector<c,n>()) == 0))
        {
            // if screwed return zeros
            //printf("annoying zero -- line 1429\n");
            return tensor<c,n>(2,con);
        }
        
        // assume all rows above m are already done
        int i = (m+1);
        while(i<n)
        {
            // does row start with zero? if not then eliminate the zero
            if(!(d[i][m]->isConstant() && (d[i][m](vector<c,n>()) == 0)))
            {
                // eliminate zero using top row
                tensorComponent<c,n> a;
                a = *(d[m][m]);
                tensorComponent<c,n> b;
                b = *(d[i][m]);
                tensorComponent<c,n> r = b / a;
                
                // replace row[i] with row[i] - r*row[m]
                // this should remove the first non-zero value on row[i]
                for(int j = 0; j < n; j++)
                {
                    *(d[i][j]) = *(d[i][j]) - r * (*d[m][j]);
                    *(ret[i][j]) = *(ret[i][j]) - r * (*ret[m][j]);
                }
            }
            i++;
        }
        
        // sort matrix so that the rows with the most leading zeros are at the bottom
        // if any of the rows turn out to all be zero, then we are stuffed
        int bottommost = 0;
        for(int i = n; i >= 0; i--)
        {
            // check each row for i number of zeros
            for(int j = 0; j < n; j++)
            {
                int numZeros = 0;
                for(int k = 0; k < n; k++)
                {
                    if(d[j][k]->isConstant() && (d[j][k](vector<c,n>()) == 0))
                    {
                        numZeros++;
                    }
                    else break;
                }

                if(numZeros == n)
                {
                    // stuffed
                    //printf("annoying zero -- line 1479\n");
                    return tensor<c,n>(2,con);
                }

                if(numZeros == i)
                {
                    // swap with bottommost row we haven't yet swapped
                    bottommost++;
                    for(int k = 0; k < n; k++)
                    {
                        tensorComponent<c,n> buf;
                        buf = *(d[j][k]);
                        *(d[j][k]) = *(d[n-bottommost][k]);
                        *(d[n-bottommost][k]) = buf;
                        buf = *(ret[j][k]);
                        *(ret[j][k]) = *(ret[n-bottommost][k]);
                        *(ret[n-bottommost][k]) = buf;
                    }
                }            
            }
        }
    }
    
    /*
    // output d half for debugging
    printf("d-half of inversion matrix\n");
    d.printDebug();
    */
    
    /*
        now we need to do some 'pivot' operations (which are also row operations anyway...)
        this will reduce the d half of the 'matrix' to the identity matrix
        
        with the top row you can remove the n-1 other components by subtracting multiples
        of all the rows underneath, one after the other, then repeat for the row below, and so on...
    */
    
    // for each row except the last, lets eliminate all non-diagonal entries
    for(int i = 0; i < (n-1); i++)
    {
        // for now just test on first row
        // subtract multiples of other rows to remove remaining values
        for(int j = i + 1; j < n; j++)
        {
            // eliminate zero using top row
            tensorComponent<c,n> a;
            a = *(d[i][j]);
            tensorComponent<c,n> b;
            // note these occur on the diagonal and therfore should never be zero when this code is executed
            b = *(d[j][j]);
            tensorComponent<c,n> r = a / b;
            
            // replace row[i] with row[i] - r*row[j]
            // this should remove the first non-zero value on row[i]
            for(int k = 0; k < n; k++)
            {
                *(d[i][k]) = *(d[i][k]) - r * (*d[j][k]);
                *(ret[i][k]) = *(ret[i][k]) - r * (*ret[j][k]);
            }
        }
    }
    
    // now all that is left is to normalise each row so that it's lone value is 1
    // i.e. make the d half look exactly like the identity matrix
    
    for(int i = 0; i < n; i++)
    {
        tensorComponent<c,n> a;
        a = 1 / *(d[i][i]);
        
        // replace row[i] with a * row[i]
        for(int j = 0; j < n; j++)
        {
            *(d[i][j]) = *(d[i][j]) * a;
            *(ret[i][j]) = *(ret[i][j]) * a;
        }
    }
    
    /*
    // output d half for debugging
    printf("d-half of inversion matrix\n");
    d.printDebug();
    */
    
    return ret.forceRaiseIndex(1);
}

// the algorithm for this is the same as above, except for the order of the tensors
// i.e. we can force the indices to lowers, raise them using the above function
// then force them back down and we get the correct result
template <class c, int n>
tensor<c,n> LowerMetricIndices(tensor<c,n>& metric)
{
    jProfiler profiler = jProfiler(L"<global>::LowerMetricIndices");
    bool co[] = {false,false};
    tensor<c,n> ret = tensor<c,n>(2, co);
    
    // do some quick checks to make sure the metric is likely good
    // if its not then return some zeros
    if(metric.getTotalRank() != 2) return ret;
    if(metric.getContravariantRank() != 2) return ret;
    
    ret = metric.forceLowerIndex(0).forceLowerIndex(1);
    
    return RaiseMetricIndices(ret).forceLowerIndex(0).forceLowerIndex(1);
}
#pragma warning( default : 4311 )
#pragma warning( default : 4313 )

#endif