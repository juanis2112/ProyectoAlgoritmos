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
    nRows = 0;
    nCols = 0;
    capacity = Initial_Capacity;
    array = new Cvector<numberType> [capacity];
}

//Copy
template <typename numberType>
Cmatrix<numberType>::Cmatrix(const Cmatrix<numberType> &x){
    capacity = x.capacity + Initial_Capacity;
    nRows = x.nRows;
    nCols = x.nCols;
    array = new Cvector<numberType> [capacity];
    for (size_t i = 0; i < nRows; i++) array[i] = x.array[i];
}

//Fill
template <typename numberType>
Cmatrix<numberType>:: Cmatrix(size_t size, const Cvector<numberType> &x, bool axis){

    //PREGUNTAR POR QUE TENGO QUE METER EL TMP Y TAMBIEN POR QUE TENGO QUE INICIALIZAR b!!!!!!!
    if (axis == false){
        capacity = size + Initial_Capacity;
        array = new Cvector<numberType> [capacity];
        nRows = size;
        nCols = x.size();
        for (size_t i = 0; i < nRows; i++) array[i] = x;
    }
    
    else if (axis == true){
        numberType b = 0;
        Cvector<numberType> tmp(size,b);
        capacity = x.size() + Initial_Capacity;
        array = new Cvector<numberType> [capacity];
        nRows = x.size();
        nCols = size;
        for (size_t i = 0; i < nRows; i++) array[i] = tmp;
        for (size_t i = 0; i < nRows; i++){
            for (size_t j = 0; j < nCols; j++){
                (array[i])[j] = x.array[i];
            }
        }
    }
}

//Specialized
template <typename numberType>
Cmatrix<numberType>:: Cmatrix(size_t row, size_t col, bool type){
    if (type == true) assert(row == col);
    nRows = row;
    nCols = col;
    capacity = nRows + Initial_Capacity;
    array = new Cvector<numberType> [capacity];
    Cvector<int> tmp(nCols,0);
    if (type == true){
        for (size_t i = 0; i < nRows ; i++){
            tmp.array[i] = 1;
            array[i] = tmp;
            tmp.array[i] = 0;
        }
    }
    
    else{
        for (size_t i = 0; i < nRows ; i++){
            array[i] = tmp;
            }
        }
    
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
    capacity = rhs.capacity + Initial_Capacity;
    array = new Cvector<numberType> [capacity];
    nRows = rhs.nRows;
    nCols = rhs.nCols;
    for (size_t i = 0; i < rhs.nRows; i++) array[i] = rhs.array[i];
    delete[] Oldarray;
    return *this;
}

//-------------------------------------Friend operators ------------------------------------------
// -------------------------- Comparisson Operators


// Operator ()
template <typename numberType>
numberType Cmatrix<numberType>::operator () (size_t row, size_t col) const{
    return array[row][col];
}

// Operator ()
template <typename numberType>
numberType& Cmatrix<numberType>::operator () (size_t row, size_t col){
    return array[row][col];
}



// Operator <<
template <typename numberType>
ostream & operator<<(ostream &os, const Cmatrix<numberType> &rhs) {
    os << "[ " << endl;
    for (size_t i = 0; i < rhs.nRows; i++) cout << rhs.array[i] << "";
    os << "]" << endl;
    
    return os;
}

// Operator ==
template <typename numberType>
Cmatrix<bool> operator == (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.nRows = x.nRows;
    result.nCols = x.nCols;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = (x.array[i] == y.array[i]);
    }
    return result;
}


// Operator !=
template <typename numberType>
Cmatrix<bool> operator != (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.nRows = x.nRows;
    result.nCols = x.nCols;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = (x.array[i] != y.array[i]);
    }
    return result;
}

// Operator <=
template <typename numberType>
Cmatrix<bool> operator <= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.nRows = x.nRows;
    result.nCols = x.nCols;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = (x.array[i] <= y.array[i]);
    }
    return result;
}

// Operator >=
template <typename numberType>
Cmatrix<bool> operator >= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.nRows = x.nRows;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = (x.array[i] >= y.array[i]);
    }
    return result;
}

// Operator <
template <typename numberType>
Cmatrix<bool> operator < (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.nRows = x.nRows;
    result.nCols = x.nCols;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = (x.array[i] < y.array[i]);
    }
    return result;
}


// Operator >
template <typename numberType>
Cmatrix<bool> operator > (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    result.nRows = x.nRows;
    result.nCols = x.nCols;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = (x.array[i] > y.array[i]);
    }
    return result;
}


//////Binary operators

// Operator +
template <typename numberType>
Cmatrix<numberType> operator + (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    Cmatrix<numberType> result;
    if (x.nRows == y.nRows){
        for (size_t i = 0; i < x.nRows; i++){
            result.array[i] = x.array[i] + y.array[i];
            }
        result.nRows = x.nRows;
        }

    else if (x.nRows < y.nRows){
        for (size_t i = 0; i < x.nRows; i++){
            result.array[i] = x.array[i] + y.array[i];
            }

        for (size_t i = x.nRows; i < y.nRows; i++){
            result.array[i] = y.array[i];
            }
        result.nRows = y.nRows;
        }

    else if (x.nRows > y.nRows){
        for (size_t i = 0; i < y.nRows; i++){
            result.array[i] = x.array[i] + y.array[i];
        }

        for (size_t i = y.nRows; i < x.nRows; i++){
            result.array[i] = x.array[i];
        }
        result.nRows = x.nRows;
    }

    return result;
}


// Operator -
template <typename numberType>
Cmatrix<numberType> operator - (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    Cmatrix<numberType> result;
    if (x.nRows == y.nRows){
        for (size_t i = 0; i < x.nRows; i++){
            result.array[i] = x.array[i] - y.array[i];
        }
        result.nRows = x.nRows;
    }

    else if (x.nRows < y.nRows){
        for (size_t i = 0; i < x.nRows; i++){
            result.array[i] = x.array[i] - y.array[i];
        }

        for (size_t i = x.nRows; i < y.nRows; i++){
            result.array[i] = y.array[i] * (-1);
        }
        result.nRows = y.nRows;
    }

    else if (x.nRows > y.nRows){
        for (size_t i = 0; i < y.nRows; i++){
            result.array[i] = x.array[i] - y.array[i];
        }

        for (size_t i = y.nRows; i < x.nRows; i++){
            result.array[i] = x.array[i] * (-1);
        }
        result.nRows = x.nRows;
    }

    return result;
}



// Operator *
template <typename numberType>
Cmatrix<numberType> operator * (const Cmatrix<numberType> &x, const int &y){
    Cmatrix<numberType> result;
    result.nRows = x.nRows;
    for (size_t i = 0; i < x.nRows; i++){
        result.array[i] = numberType (x.array[i] * y);
    }
    
    return result;
}

// Operator /
template <typename numberType>
Cmatrix<numberType> operator / (const Cmatrix<numberType> &x, const int &y){
    Cmatrix<numberType> result;
    result.nRows = x.nRows;
    for (size_t i = 0; i < x.nRows; i++){
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
    if(nRows == capacity){
        expandCapacity();
    }
    array[nRows++] = value;
}

// Erase
template <typename numberType>
void Cmatrix<numberType>::erase(size_t index){
    if(index == nRows - 1){
        --nRows;
    }
    else if(index < 0 || index > nRows){
        cout << "This index is not in the matrix" << endl;
    }
    else{
        numberType *Oldarray = array;
        array = new numberType[capacity];
        for(size_t i = 0; i < index; i++){
            array[i] = Oldarray[i];
        }
        for(size_t i = index + 1; i < nRows; i++){
            array[i-1] = Oldarray[i];
        }
        delete[]Oldarray;
        --nRows;
    }
}

// Insert
template <typename numberType>
void Cmatrix<numberType>::insert (size_t index, const Cvector<numberType> & value){
    if(index == nRows) {
        push(value);
    }
    else{
        if(nRows == capacity) {
            expandCapacity();
        }
        Cvector<numberType> *Oldarray = array;
        array = new Cvector<numberType>[capacity];
        if(index == 0){
            array[0] = value;
            for (size_t i = 1; i <= nRows; i++){
                array[i] = Oldarray[i-1];
            }
        }
        else{
            for (size_t i = 0; i < index; i++) array[i] = Oldarray[i];
            array[index] = value;
            for (size_t i = index + 1; i <= nRows; i++) array[i] = Oldarray[i-1];
        }
        delete[]Oldarray;
        nRows++;
    }
}

// Clear
template <typename numberType>
void Cmatrix<numberType>::clear(){
    nRows = 0;
}

// Empty
template <typename numberType>
bool Cmatrix<numberType>::empty() const{
    return (nRows == 0);
}

// Size
template <typename numberType>
size_t Cmatrix<numberType>::size() const{
    return nRows;
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
    for(size_t i = 0; i < nRows; i++){
        array[i] = Oldarray[i];
    }
    delete[] Oldarray;
    
}

#endif // CMatrix_hpp
