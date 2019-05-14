//  CVector.hpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#ifndef CVector_hpp
#define CVector_hpp

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
using namespace std;

const int Initial_Capacity = 10;
//*************************************   FRIEND OPERATORS DEFINITIONS   **************************************

template <typename numberType> class Cvector;
template <typename numberType> Cvector<bool> operator == (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<bool>  operator != (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<bool>  operator > (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<bool>  operator < (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<bool>  operator >= (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<bool>  operator <= (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<numberType>  operator + (const Cvector<numberType> & x, const Cvector<numberType> & y);
template <typename numberType> Cvector<numberType>  operator - (const Cvector<numberType> & x, const Cvector<numberType> & y);
template <typename numberType> Cvector<double>  operator * (const Cvector<numberType> &x, const numberType &y);
template <typename numberType> Cvector<double>  operator * (const numberType &y, const Cvector<numberType> &x);
template <typename numberType> Cvector<double>  operator / (const Cvector<numberType> &x, const numberType &y);
template <typename numberType> Cvector<double>  operator ^ (const Cvector<numberType> &x, const int y);
template <typename numberType> ostream & operator << (ostream &os, const Cvector<numberType> &rhs);


//*********************************************  CLASS DEFINITION   *********************************************

template <typename numberType>
class Cvector{
public:

//-------------------------------CONSTRUCTORS AND DESTRUCTOR----------------------------------------------
    // Empty
    Cvector();                                       //FUNCIONA
    // Fill
    Cvector(size_t size, numberType value);          //FUNCIONA
    // Parametric
    Cvector(const Cvector &rhs);                      //FUNCIONA

    // Destructor
    ~Cvector();                                      //FUNCIONA

//-----------------------------------------OPERATORS------------------------------------------------------

    // Class member operators
    Cvector<numberType> operator =(const Cvector &rhs);   //FUNCIONA
    numberType operator [] (size_t idx) const;                 //FUNCIONA
    numberType & operator [](size_t idx);

    // Friend operators
    friend Cvector<bool> operator == <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator != <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator > <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator < <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator >= <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator <= <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend ostream & operator << <> (ostream &os, const Cvector<numberType> &rhs);
    //friend istream & operator>>(istream &os, const vector &rhs);

    // Binary operators
    friend Cvector<numberType>  operator + <>(const Cvector<numberType> & x, const Cvector<numberType> & y);
    friend Cvector<numberType>  operator - <>(const Cvector<numberType> & x, const Cvector<numberType> & y);
    friend Cvector<double>  operator * <>(const Cvector<numberType> &x, const numberType &y);
    friend Cvector<double>  operator * <>(const numberType &y, const Cvector<numberType> &x);
    friend Cvector<double>  operator / <>(const Cvector<numberType> &x, const numberType &y);
    friend Cvector<double>  operator ^ <>(const Cvector<numberType> &x, const int y);

//------------------------------------------CLASS METHODS------------------------------------------------------

    double dot(Cvector w);  // FALTA LA QUE SOLO SE PUEDA PARA R3
    Cvector<double> cross(Cvector<numberType> w);
    double norm();
    Cvector<double> normalize();
    void push(numberType value);
    void erase(size_t index);
    void insert (size_t index, numberType value);
    void clear();
    bool empty() const;
    size_t size() const;
    double angle (Cvector<numberType> &x);
    Cvector<double> proy(Cvector<numberType> &x);
    Cvector<Cvector<double>>  gram_schmidt(const Cvector<numberType> &rhs);

    ///////////////////// to double precision

    Cvector<double> toDouble()const;

    ///////////////////////////////

//---------------------------------------------PRIVATE---------------------------------------------------------


    size_t capacity = 10 , length = 0;
    numberType *array;
    void expandCapacity();
};

#include "CVector.cpp"
#endif /* CVector_hpp */
