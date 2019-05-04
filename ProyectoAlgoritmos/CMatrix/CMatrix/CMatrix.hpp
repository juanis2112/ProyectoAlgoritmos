//  Cmatrix.hpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#ifndef CMatrix_hpp
#define CMatrix_hpp

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "../../CVector/CVector/CVector.hpp"
using namespace std;




//*************************************   FRIEND OPERATORS DEFINITIONS   **************************************

template <typename numberType> class Cmatrix;
template <typename numberType> Cmatrix<bool> operator == (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator != (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator >(const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator < (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator >= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator <= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);


template <typename numberType> Cmatrix<numberType>  operator+ (const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
template <typename numberType> Cmatrix<numberType>  operator- (const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
template <typename numberType> Cmatrix<numberType> operator * (const Cmatrix<numberType> &x, const int &y);
template <typename numberType> Cmatrix<numberType> operator / (const Cmatrix<numberType> &x, const int &y);
//template <typename numberType> Cmatrix<Cvector<numberType>>  operator^ (const Cmatrix<Cvector<numberType>> &x, const Cmatrix<Cvector<numberType>> &y);


template <typename numberType> ostream & operator<< (ostream &os, const Cmatrix<numberType> &rhs);


//*********************************************  CLASS DEFINITION   *********************************************

template<typename numberType>
class Cmatrix{
    public:
    
		//-------------------------------CONSTRUCTORS AND DESTRUCTOR----------------------------------------------
		// Empty
        Cmatrix();
        // Parametric
        Cmatrix(const Cmatrix<numberType> &x);
        // Fill
        Cmatrix(size_t size, const Cvector<numberType> &x);
        //Specialized
        Cmatrix(size_t row, size_t cod);
        // Destructor
        ~Cmatrix();

        //-----------------------------------------OPERATORS------------------------------------------------------

        // Class member operators
        Cmatrix<numberType> operator =(const Cmatrix<numberType> &rhs);
       // numberType operator [:] (size_t row, size_t col);

        // Friend operators
        friend Cmatrix<bool> operator == <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator != <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator > <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator < <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator >= <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator <= <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend ostream & operator << <> (ostream &os, const Cmatrix<numberType> &rhs);
        //friend istream & operator>>(istream &os, const vector &rhs);

       // Binary operators
        friend Cmatrix<numberType>  operator+ <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
        friend Cmatrix<numberType>  operator - <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
        friend Cmatrix<numberType> operator * <> (const Cmatrix<numberType> &x, const int &y);
        friend Cmatrix<numberType> operator / <> (const Cmatrix<numberType> &x, const int &y);
//        //friend Cmatrix<Cvector<numberType>>  operator^ <>(const Cmatrix<Cvector<numberType>> &x, const Cmatrix<Cvector<numberType>> &y);
    

	Cmatrix<int> Identity(size_t indx);

    void push (const Cvector<numberType> &value);
    void erase(size_t index);
    void insert (size_t index, const Cvector<numberType> & value);
    void clear();
    bool empty() const;
    size_t size() const;
	numberType access (size_t row, size_t cod) const;

	public:
		Cvector<numberType> * array;
		size_t capacity, count;
		void expandCapacity();
		
};

#include "CMatrix.cpp"
#endif //CMatrix_hpp

