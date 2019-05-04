//
//  CMatrix.cpp
//  CMatrix
//
//  Created by Juanita Gomez on 4/30/19.
//  Copyright © 2019 Juanita Gomez. All rights reserved.
//


#ifdef CMatrix_hpp

// -------------------------- Constructors

// Empty
template <typename numberType>
Cmatrix<numberType>::Cmatrix() {
    count = 0;
    capacity = Initial_Capacity;
    array = new Cvector<numberType> [capacity];
}

//Copy
template <typename numberType>
Cmatrix<numberType>::Cmatrix(const Cmatrix<numberType> &x){
    capacity = x.capacity + Initial_Capacity;
    count = x.count;
    array = new Cvector<numberType> [capacity];
    for (size_t i = 0; i < x.count; i++) array[i] = x.array[i];
}

//Fill
template <typename numberType>
Cmatrix<numberType>:: Cmatrix(size_t size, const Cvector<numberType> &x){
    capacity = size + Initial_Capacity;
    count = size;
    array = new Cvector<numberType> [capacity];
    for (size_t i = 0; i < size; i++) array[i] = x;
}

//Specialized
template <typename numberType>
Cmatrix<numberType>:: Cmatrix(size_t row, size_t cod){
    numberType tmp;
    Cvector<numberType> p(row,tmp);
    count = cod;
    for (size_t i = 0; i < count; i++) array[i] = p;
}


// -------------------------------------Destructor------------------------------------------

template <typename numberType>
Cmatrix<numberType>:: ~Cmatrix(){
    delete [] array;
}

//------------------------------------- OPERATORS ------------------------------------------

//-----------------------------------Non_friend operators ----------------------------------
// Operator =
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::operator=(const Cmatrix<numberType> &rhs){
    Cvector<numberType> *Oldarray = this -> array;
    this -> capacity = rhs.capacity + Initial_Capacity;
    this -> array = new Cvector<numberType> [capacity];
    for (size_t i = 0; i < rhs.count; i++) array[i] = rhs.array[i];
    this -> count = rhs.count;
    delete[] Oldarray;
    return *this;
}

//-------------------------------------Friend operators ------------------------------------------
// -------------------------- Comparisson Operators

// Operator <<
template <typename numberType>
ostream & operator<<(ostream &os, const Cmatrix<numberType> &rhs) {
    os << "[ " << endl;
    for (size_t i = 0; i < rhs.count; i++) cout << rhs.array[i] << "";
    os << "]" << endl;
    
    return os;
}

// Operator ==
template <typename numberType>
Cmatrix<bool> operator == (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = (x.array[i] == y.array[i]);
    }
    return result;
}


// Operator !=
template <typename numberType>
Cmatrix<bool> operator != (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = (x.array[i] != y.array[i]);
    }
    return result;
}

// Operator <=
template <typename numberType>
Cmatrix<bool> operator <= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = (x.array[i] <= y.array[i]);
    }
    return result;
}

// Operator >=
template <typename numberType>
Cmatrix<bool> operator >= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = (x.array[i] >= y.array[i]);
    }
    return result;
}

// Operator <
template <typename numberType>
Cmatrix<bool> operator < (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = (x.array[i] < y.array[i]);
    }
    return result;
}


// Operator >
template <typename numberType>
Cmatrix<bool> operator > (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = (x.array[i] > y.array[i]);
    }
    return result;
}


//////Binary operators

// Operator +
template <typename numberType>
Cmatrix<numberType> operator + (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    Cmatrix<numberType> result;
    if (x.count == y.count){
        for (size_t i = 0; i < x.count; i++){
            result.array[i] = x.array[i] + y.array[i];
            }
        result.count = x.count;
        }

    else if (x.count < y.count){
        for (size_t i = 0; i < x.count; i++){
            result.array[i] = x.array[i] + y.array[i];
            }

        for (size_t i = x.count; i < y.count; i++){
            result.array[i] = y.array[i];
            }
        result.count = y.count;
        }

    else if (x.count > y.count){
        for (size_t i = 0; i < y.count; i++){
            result.array[i] = x.array[i] + y.array[i];
        }

        for (size_t i = y.count; i < x.count; i++){
            result.array[i] = x.array[i];
        }
        result.count = x.count;
    }

    return result;
}


// Operator -
template <typename numberType>
Cmatrix<numberType> operator - (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    Cmatrix<numberType> result;
    if (x.count == y.count){
        for (size_t i = 0; i < x.count; i++){
            result.array[i] = x.array[i] - y.array[i];
        }
        result.count = x.count;
    }

    else if (x.count < y.count){
        for (size_t i = 0; i < x.count; i++){
            result.array[i] = x.array[i] - y.array[i];
        }

        for (size_t i = x.count; i < y.count; i++){
            result.array[i] = y.array[i] * (-1);
        }
        result.count = y.count;
    }

    else if (x.count > y.count){
        for (size_t i = 0; i < y.count; i++){
            result.array[i] = x.array[i] - y.array[i];
        }

        for (size_t i = y.count; i < x.count; i++){
            result.array[i] = x.array[i] * (-1);
        }
        result.count = x.count;
    }

    return result;
}



// Operator *
template <typename numberType>
Cmatrix<numberType> operator * (const Cmatrix<numberType> &x, const int &y){
    Cmatrix<numberType> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = numberType (x.array[i] * y);
    }
    
    return result;
}

// Operator /
template <typename numberType>
Cmatrix<numberType> operator / (const Cmatrix<numberType> &x, const int &y){
    Cmatrix<numberType> result;
    result.count = x.count;
    for (size_t i = 0; i < x.count; i++){
        result.array[i] = numberType (x.array[i] / y);
    }
    
    return result;
}

//// Operator ^
//template <typename numberType>
//Cmatrix<numberType> operator ^ (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
//    Cmatrix<numberType> result(x);
//    for (size_t i = 0; i < y.array[0]; i++){
//        result.array *= x.array;
//    }
//
//    return result;
//}


// Push
template <typename numberType>
void Cmatrix<numberType>::push(const Cvector<numberType> &value){
    if(count == capacity){
        expandCapacity();
    }
    array[count++] = value;
}

// Erase
template <typename numberType>
void Cmatrix<numberType>::erase(size_t index){
    if(index == count - 1){
        --count;
    }
    else if(index < 0 || index > count){
        cout << "This index is not in the matrix" << endl;
    }
    else{
        numberType *Oldarray = array;
        array = new numberType[capacity];
        for(size_t i = 0; i < index; i++){
            array[i] = Oldarray[i];
        }
        for(size_t i = index + 1; i < count; i++){
            array[i-1] = Oldarray[i];
        }
        delete[]Oldarray;
        --count;
    }
}

// Insert
template <typename numberType>
void Cmatrix<numberType>::insert (size_t index, const Cvector<numberType> & value){
    if(index == count) {
        push(value);
    }
    else{
        if(count == capacity) {
            expandCapacity();
        }
        Cvector<numberType> *Oldarray = array;
        array = new Cvector<numberType>[capacity];
        if(index == 0){
            array[0] = value;
            for (size_t i = 1; i <= count; i++){
                array[i] = Oldarray[i-1];
            }
        }
        else{
            for (size_t i = 0; i < index; i++) array[i] = Oldarray[i];
            array[index] = value;
            for (size_t i = index + 1; i <= count; i++) array[i] = Oldarray[i-1];
        }
        delete[]Oldarray;
        count++;
    }
}

// Clear
template <typename numberType>
void Cmatrix<numberType>::clear(){
    count = 0;
}

// Empty
template <typename numberType>
bool Cmatrix<numberType>::empty() const{
    return (count == 0);
}

// Size
template <typename numberType>
size_t Cmatrix<numberType>::size() const{
    return count;
}

// Identity
template <typename numberType>
Cmatrix<int> Cmatrix<numberType>::Identity(size_t indx){
	Cvector<int> tmp(indx,0);
	Cmatrix<int> result;
	for (size_t i = 0; i < indx; i++){
		tmp.array[i] = 1;
		result.push(tmp);
		tmp.array[i] = 0;
		}
	return result;    
}

// Access
template <typename numberType>
numberType Cmatrix<numberType>::access(size_t row, size_t cod) const{
    return this -> array[row - 1][cod - 1];
}


//--------------------------------------Expand Capacity--------------------------------------------

template <typename numberType>
void Cmatrix<numberType>::expandCapacity(){
Cvector<numberType> *Oldarray = array;
    capacity *= 2;
    array = new Cvector<numberType>[capacity];
    for(size_t i = 0; i < count; i++){
        array[i] = Oldarray[i];
    }
    delete[] Oldarray;
    
}

#endif // CMatrix_hpp
