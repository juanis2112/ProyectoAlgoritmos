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
#include <typeinfo>
#include <list>
#include "../../CVector/CVector/CVector.hpp"
using namespace std;

//*************************************   FRIEND OPERATORS DECLARATIONS   *****************************************

template <typename numberType> class Cmatrix;
template <typename numberType> Cmatrix<bool> operator == (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator != (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator >(const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator < (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator >= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<bool> operator <= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);

template <typename numberType> Cmatrix<double>  operator+ (const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
template <typename numberType> Cmatrix<double>  operator- (const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
template <typename numberType> Cmatrix<double> operator * (const Cmatrix<numberType> &x, const numberType &y);
template <typename numberType> Cmatrix<double> operator * ( const numberType &y, const Cmatrix<numberType> &x);
template <typename numberType> Cmatrix<double> operator * (const Cmatrix<numberType> & x, const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<double> operator * (const Cmatrix<numberType> & x, const Cvector<numberType> &y);
template <typename numberType> Cmatrix<double> operator * (const Cvector<numberType> & x, const Cmatrix<numberType> &y);
template <typename numberType> Cmatrix<double> operator / (const Cmatrix<numberType> &x, const numberType &y);
template <typename numberType> Cmatrix<double> operator ^ (const Cmatrix<numberType> &x, const size_t &y);

template <typename numberType> ostream & operator<< (ostream &os, Cmatrix<numberType> &rhs);

//*********************************************  CLASS DEFINITION   *********************************************

template<typename numberType>
class Cmatrix{
public:

//-------------------------------CONSTRUCTORS AND DESTRUCTOR----------------------------------------------------
    
    // Empty Constructor
    Cmatrix();
    // Fill Constructor
    Cmatrix(size_t size, const Cvector<numberType> &x, bool axis = false);
    // Parametric Constructor
    Cmatrix(const Cmatrix<numberType> &x);
    //Specialized Constructor
    Cmatrix(size_t row, size_t col, bool type);
    // Destructor
    ~Cmatrix();

//-----------------------------------------OPERATORS------------------------------------------------------------

    // Class member operators
    Cmatrix<numberType> operator =(const Cmatrix<numberType> &rhs);
    numberType operator () (size_t row, size_t col) const;
    numberType & operator () (size_t row, size_t col);
    Cvector<numberType> operator [] (size_t idx) const;
    Cvector<numberType> & operator [](size_t idx);
    
    // Friend operators
    friend Cmatrix<bool> operator == <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
    friend Cmatrix<bool> operator != <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
    friend Cmatrix<bool> operator > <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
    friend Cmatrix<bool> operator < <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
    friend Cmatrix<bool> operator >= <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
    friend Cmatrix<bool> operator <= <> (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y);
    friend ostream & operator << <> (ostream &os, Cmatrix<numberType> &rhs);

   // Binary operators
    friend Cmatrix<double>  operator+ <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
    friend Cmatrix<double>  operator - <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> & y);
    friend Cmatrix<double> operator * <> (const Cmatrix<numberType> &x, const numberType &y);
    friend Cmatrix<double> operator * <> (const numberType &y, const Cmatrix<numberType> &x);
    friend Cmatrix<double> operator * <>(const Cmatrix<numberType> & x, const Cmatrix<numberType> &y);
    friend Cmatrix<double> operator * <>(const Cmatrix<numberType> & x, const Cvector<numberType> &y);
    friend Cmatrix<double> operator * <>(const Cvector<numberType> & x, const Cmatrix<numberType> &y);
    friend Cmatrix<double> operator / <> (const Cmatrix<numberType> &x, const numberType &y);
    friend Cmatrix<double> operator ^ <>(const Cmatrix<numberType> &x, const size_t &y);

//-------------------------------------------------CLASS METHODS------------------------------------------------------

    // Matrix Methods
    Cmatrix<double> toDouble() const;
    void push (const Cvector<numberType> &value, bool axis = 0);
    void erase(size_t index, bool axis);
    void insert (size_t index, const Cvector<numberType> & value, bool axis = 0);
    void clear();
    bool empty() const;
    size_t numberCols() const;
    size_t numberRows() const;
	numberType access (size_t row, size_t cod) const;
    
    //Easy Matrixes
    static Cmatrix<double> eye(size_t N);
    static Cmatrix<double> zeros(size_t rows, size_t cols);
    static Cmatrix<double> ones(size_t rows, size_t cols);
    static Cmatrix<double> random(size_t rows, size_t cols);
    static Cmatrix<double> diagonalize(const Cvector<numberType> &rhs);
    static Cvector<double> diagonal(Cmatrix<numberType> &m);
    static Cmatrix<double> permutationMatrix(Cvector<numberType> &v);

    // Modifiers
    void swap_r(size_t row1, size_t row2);
    void swap_c(size_t col1, size_t col2);
    void appendRows(const Cmatrix<numberType> &rhs);
    void appendCols(const Cmatrix<numberType> &rhs);
    
    // Special Matrixes
    Cmatrix<numberType> abs();
    Cmatrix<numberType> transpose();
    Cmatrix<double> inverse();
    Cmatrix<double> lowerTriangular();
    Cmatrix<double> upperTriangular();
    
    // Matrix Decompositions
    tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> LUP(double Tol);
    tuple<Cmatrix<double>, Cmatrix<double>> QR();
    
    // Matrix Properties
    Cvector<double> eigen_values(const double tol);
    double determinant();
   
// ************************************************ PRIVATE *********************************************+*
    
	private:
		size_t capacity, nRows, nCols;
        Cvector<numberType> * array;
		void expandCapacity();
        void Checkrep() const;

};

#include "CMatrix.cpp"
#endif //CMatrix_hpp
