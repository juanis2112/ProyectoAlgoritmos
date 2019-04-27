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
template <typename numberType> Cvector<numberType>  operator+ (const Cvector<numberType> & x, Cvector<numberType> & y);
template <typename numberType> Cvector<numberType>  operator- (const Cvector<numberType> & x, Cvector<numberType> & y);
template <typename numberType> Cvector<numberType>  operator* (const Cvector<numberType> &x, const Cvector<numberType> &y);
template <typename numberType> Cvector<numberType>  operator/ (const Cvector<numberType> &x, const Cvector<numberType> &y);
template <typename numberType> Cvector<numberType>  operator^ (const Cvector<numberType> &x, const Cvector<numberType> &y);


//template <typename numberType> ostream & operator<< (ostream &os, const Cvector<numberType> &rhs);


//*********************************************  CLASS DEFINITION   *********************************************

template <typename numberType>
class Cvector{
public:
    
//-------------------------------CONSTRUCTORS AND DESTRUCTOR----------------------------------------------
    // Empty
    Cvector();
    // Fill
    Cvector(size_t size, numberType value);
    // Parametric
    Cvector(const Cvector &rhs);
    
    // Destructor
    ~Cvector();
    
//-----------------------------------------OPERATORS------------------------------------------------------
    
    // Class member operators
    Cvector<numberType> operator =(const Cvector &rhs);
    numberType operator [] (size_t idx);
    
    // Friend operators
    friend Cvector<bool> operator == <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator != <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator > <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator < <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator >= <>(const Cvectopr<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator <= <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    //friend ostream & operator << <> (ostream &os, const Cvector<numberType> &rhs);
    //friend istream & operator>>(istream &os, const vector &rhs);
    
    // Binary operators
    friend Cvector<numberType>  operator+ <>(const Cvector<numberType> & x, Cvector<numberType> & y);
    friend Cvector<numberType>  operator- <>(const Cvector<numberType> & x, Cvector<numberType> & y);
    friend Cvector<numberType>  operator*<>(const Cvector<numberType> &x, const Cvector<numberType> &y);
    friend Cvector<numberType>  operator/<>(const Cvector<numberType> &x, const Cvector<numberType> &y);
    friend Cvector<numberType>  operator^ <>(const Cvector<numberType> &x, const Cvector<numberType> &y);
    
//------------------------------------------CLASS METHODS------------------------------------------------------
    
    double dot(Cvector w);
    Cvector<numberType> cross(Cvector<numberType> w);
    double norm();
    Cvector<numberType> normalize();
    void push(numberType value);
    void erase(size_t index);
    void insert (size_t index, numberType value);
    void clear();
    bool empty() const;
    size_t size() const;

//---------------------------------------------PRIVATE---------------------------------------------------------
    
private:
    size_t capacity = 10 , length = 0;
    numberType *array;
    void expandCapacity();
};

#include "CVector.cpp"
#endif /* CVector_hpp */
