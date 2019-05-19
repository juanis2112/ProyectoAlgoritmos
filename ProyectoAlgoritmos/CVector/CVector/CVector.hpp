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
#include <cassert>
#include<type_traits>
using namespace std;
const int Initial_Capacity = 10;

//*************************************   FRIEND OPERATORS DEFINITIONS   *****************************************

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
template <typename numberType> Cvector<double>  operator * (const Cvector<numberType> &x , const Cvector<numberType> &y);
template <typename numberType> Cvector<double>  operator / (const Cvector<numberType> &x, const numberType &y);
template <typename numberType> Cvector<double>  operator ^ (const Cvector<numberType> &x, const int y);
template <typename numberType> ostream & operator << (ostream &os, const Cvector<numberType> &rhs);


//*********************************************  CLASS DEFINITION   **********************************************

template <typename numberType>
class Cvector{
public:

//-------------------------------CONSTRUCTORS AND DESTRUCTOR------------------------------------------------------
    
    // Empty Constructor
    Cvector();
    // Fill Constructor
    Cvector(size_t size, numberType value);
    // Parametric Constructor
    Cvector(const Cvector &rhs);
    // Size Constructor
    Cvector(size_t Length);
    // Destructor
    ~Cvector();

//-----------------------------------------OPERATORS--------------------------------------------------------------

    // Class member operators
    Cvector<numberType> operator =(const Cvector &rhs);
    numberType operator [] (size_t idx) const;
    numberType & operator [](size_t idx);

    // Friend operators
    friend Cvector<bool> operator == <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator != <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator > <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator < <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator >= <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend Cvector<bool>  operator <= <>(const Cvector<numberType> &x , const Cvector<numberType> &y);
    friend ostream & operator << <> (ostream &os, const Cvector<numberType> &rhs);

    // Binary operators
    friend Cvector<numberType>  operator + <>(const Cvector<numberType> & x, const Cvector<numberType> & y);
    friend Cvector<numberType>  operator - <>(const Cvector<numberType> & x, const Cvector<numberType> & y);
    friend Cvector<double>  operator * <>(const Cvector<numberType> &x, const numberType &y);
    friend Cvector<double>  operator * <>(const numberType &y, const Cvector<numberType> &x);
    
    friend Cvector<double> operator * <> (const Cvector<numberType> &x, const Cvector<numberType> &y);
    friend Cvector<double>  operator / <>(const Cvector<numberType> &x, const numberType &y);
    friend Cvector<double>  operator ^ <>(const Cvector<numberType> &x, const int y);

//------------------------------------------CLASS METHODS--------------------------------------------------------
    
    // Vector Methods
    Cvector<double> toDouble()const;
    void push(numberType value);
    void erase(size_t index);
    void insert (size_t index, numberType value);
    void clear();
    bool empty() const;
    size_t size() const;
    
    // Math Methods
    double dot(Cvector<numberType> w);
    Cvector<double> cross(Cvector<numberType> w);
    double norm();
    Cvector<double> normalize();
    double angle (Cvector<numberType> &x);
    Cvector<double> proj(Cvector<numberType> &x);
    Cvector<Cvector<double>>  gram_schmidt();
    
    
//*********************************************  PRIVATE  *****************************************************

private:
    size_t capacity = Initial_Capacity ,
    length = 0;
    numberType *array;
    void expandCapacity();
    void Checkrep() const;
};

#include "CVector.cpp"
#endif /* CVector_hpp */
