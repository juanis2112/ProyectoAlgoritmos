//  Cmatrix.hpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#ifndef CMatrix_hpp
#define CMatrix_hpp

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <tuple>
#include "../../CVector/CVector/CVector.hpp"
using namespace std;




//*************************************   FRIEND OPERATORS DEFINITIONS   **************************************

template <typename numberType> class Cmatrix;
//template <typename numberType> numberType operator () (size_t row, size_t col);
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


template <typename numberType> ostream & operator<< (ostream &os, Cmatrix<numberType> &rhs);


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
        Cmatrix(size_t size, const Cvector<numberType> &x, bool axis = false);
        //Specialized
        Cmatrix(size_t row, size_t col, bool type);
        // Destructor
        ~Cmatrix();

        //-----------------------------------------OPERATORS------------------------------------------------------

        // Class member operators
        Cmatrix<numberType> operator =(const Cmatrix<numberType> &rhs);
       // numberType operator [:] (size_t row, size_t col);

        // Friend operators
        numberType operator () (size_t row, size_t col) const;
        numberType & operator () (size_t row, size_t col);
        //Cvector<numberType> operator () (size_t idx, bool type);

        friend Cmatrix<bool> operator == <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator != <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator > <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator < <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator >= <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend Cmatrix<bool> operator <= <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
        friend ostream & operator << <> (ostream &os, Cmatrix<numberType> &rhs);
        //friend istream & operator>>(istream &os, const vector &rhs);

       // Binary operators
        friend Cmatrix<numberType>  operator+ <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
        friend Cmatrix<numberType>  operator - <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
        friend Cmatrix<numberType> operator * <> (const Cmatrix<numberType> &x, const int &y);
        friend Cmatrix<numberType> operator / <> (const Cmatrix<numberType> &x, const int &y);
//        //friend Cmatrix<Cvector<numberType>>  operator^ <>(const Cmatrix<Cvector<numberType>> &x, const Cmatrix<Cvector<numberType>> &y);


    void push (const Cvector<numberType> &value);
    void erase(size_t index);
    void insert (size_t index, const Cvector<numberType> & value);
    void clear();
    bool empty() const;

    // Getters
    size_t numberCols() const;
    size_t numberRows() const;
	numberType access (size_t row, size_t cod) const;

    Cmatrix<numberType> transpose();
    Cmatrix<numberType> lowerTriangular();
    Cmatrix<numberType> upperTriangular();
    static Cmatrix<numberType> eye(size_t N);
    static Cmatrix<numberType> zeros(size_t rows, size_t cols);
    static Cmatrix<numberType> ones(size_t rows, size_t cols);
    static Cmatrix<numberType> random(size_t rows, size_t cols);


    void swap_r(size_t row1, size_t row2);
    void swap_c(size_t col1, size_t col2);
    Cmatrix<numberType> abs();


    tuple<Cvector<numberType>, Cmatrix<numberType>, Cmatrix<numberType>> LUP(double Tol);
    double determinant();

    ///////////////////// to double precision

    // Cmatrix<numberType> toDouble(const Cmatrix<numberType> &x);

    ///////////////////////////////
	public:
		Cvector<numberType> * array;
		size_t capacity, nRows, nCols;
		void expandCapacity();

};


// Matrix Product

template <typename numberType>
Cmatrix<numberType> product(Cmatrix<numberType> x, Cmatrix<numberType> y); // producto entre matrices

///////////////////////////////

#include "CMatrix.cpp"
#endif //CMatrix_hpp
