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
template <typename numberType> Cvector<numberType>  operator + (const Cvector<numberType> & x, const Cvector<numberType> & y);
template <typename numberType> Cvector<numberType>  operator - (const Cvector<numberType> & x, const Cvector<numberType> & y);
template <typename numberType> Cvector<numberType>  operator * (const Cvector<numberType> &x, const double &y);
template <typename numberType> Cvector<numberType>  operator / (const Cvector<numberType> &x, const double &y);
template <typename numberType> Cvector<numberType>  operator ^ (const Cvector<numberType> &x, int y);
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
    friend Cvector<numberType>  operator * <>(const Cvector<numberType> &x, const double &y);
    friend Cvector<numberType>  operator / <>(const Cvector<numberType> &x, const double &y);
    friend Cvector<numberType>  operator ^ <>(const Cvector<numberType> &x, int y);
    
//------------------------------------------CLASS METHODS------------------------------------------------------
    
    double dot(Cvector w);  // FALTA LA QUE SOLO SE PUEDA PARA R3
    Cvector<numberType> cross(Cvector<numberType> w);
    double norm();
    Cvector<numberType> normalize();
    void push(numberType value);
    void erase(size_t index);
    void insert (size_t index, numberType value);
    void clear();
    bool empty() const;
    size_t size() const;
    float angle (Cvector<numberType> &x);
    Cvector<numberType> proy(Cvector<numberType> &x);

//---------------------------------------------PRIVATE---------------------------------------------------------
    
	
    size_t capacity = 10 , length = 0;
    numberType *array;
    void expandCapacity();
};

#include "CVector.cpp"
#endif /* CVector_hpp */
