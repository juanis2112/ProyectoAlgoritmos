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
ostream & operator<<(ostream &os, Cmatrix<numberType> &rhs) {
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

// Operator * Matrix Multiplication
template <typename numberType>
Cmatrix<numberType> operator * (const Cmatrix<numberType> &x, Cmatrix<numberType> &y){
    Cmatrix<numberType> result;
    
    size_t xRows = x.nRows();
    size_t yRows = y.nRows();
    size_t xCols = x.nCols();
    size_t yCols = y.nCols();
    assert(xCols == yRows);
    
    for (int i = 0; i < xRows; i++)
    {
        for (int j = 0; j < yCols; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < yRows; k++)
            {
                result[i][j] += x[i][k] * y[k][j];
            }
        }
    }
    return result;
}

//// Operator * Matrix - Vector Multiplication
//
//template <typename numberType>
//Cmatrix<numberType> operator * (const Cvector<numberType> &v, Cmatrix<numberType> &x){
//    Cmatrix<numberType> result;
//
//    vRows = 1;
//    vCols = v.size();
//    xRows = x.nRows();
//    xCols = x.nCols();
//
//    assert(vCols() == xRows());
//
//    for (j = 0; j < yCols; j++)
//    {
//        result[i][j] = 0;
//        for (k = 0; k < yRows; k++)
//        {
//            result[i][j] += x[i][k] * y[k][j];
//        }
//    }
//
//    return result;
//}
//
//// Operator *
//template <typename numberType>
//Cmatrix<numberType> operator * (const Cmatrix<numberType> &x, Cvector<numberType> &y){
//    Cmatrix<numberType> result;
//    assert(x.count == y.count);
//
//    for(i = 0; i < r1; ++i)
//        for(j = 0; j < c2; ++j)
//            for(k = 0; k < c1; ++k)
//            {
//                result.array[i][j] += x.array[i][k] * y.array[k][j];
//            }
//
//}
//
//return result;
//}


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

// nCols
template <typename numberType>
size_t Cmatrix<numberType>::size() const{
    return nRows;
}

// Rows
template <typename numberType>
size_t Cmatrix<numberType>::nRows() const{
    return array[0].size();
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
numberType Cmatrix<numberType>::access(size_t row, size_t col) const{
    return this -> array[row][col];
}

template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::transpose() {
    size_t nRows = this->nCols();
    size_t nCols = this->nRows();
    
    Cmatrix<numberType> t(nRows, nCols);
    
    for (size_t i = 0; i < nCols; i++) {
        for(size_t j = 0; j < nRows; j++) {
            t.array[i][j]=this->array[j][i];
        }
    }
    
    return t;
}

// LU
template <typename numberType>
int Cmatrix<numberType>::LUP(Cmatrix<numberType> &A, double Tol) {
//Cvector<Cmatrix<double>> Cmatrix<numberType>::LUP(Cmatrix<numberType> &A, double Tol) {
    
    size_t j, k, indexMax;
    double maxPivot, absA;
    //Cmatrix<double> P;
    Cmatrix<double> L;
    Cmatrix<double> U;
    
    
    //assert(A.nRows()!=A.nCols());
    size_t N = A.nRows();
    Cvector<double> P;
    for (size_t i = 0; i <= N; i++)
        P.push(i);
    
    for (size_t i = 0; i < N; i++) {
        maxPivot = 0;
        indexMax = i;
        
        for (k = i; k < N; k++){
            absA = fabs(A.array[k][i]);
            if (absA > maxPivot) {
                maxPivot = absA;
                indexMax = k;
            }
        }
        
        //if (maxPivot < Tol) return 0; //failure, matrix is degenerate
        
        if (indexMax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[indexMax];
            P[indexMax] = j;
            
            //pivoting rows of A
            Cvector<double> PivotVector;
            PivotVector = A.array[i];
            A.array[i] = A.array[indexMax];
            A.array[indexMax] = PivotVector;
            
            //counting pivots starting from N (for determinant)
            P[N]++;
            
        }
        
        for (j = i + 1; j < N; j++) {
            A.array[j][i] /= A.array[i][i];
            
            for (k = i + 1; k < N; k++)
                A.array[j][k] -= A.array[j][i] * A.array[i][k];
        
        }
        
        
    }
    
//    Cvector<Cmatrix<double>> result;
//    result.length= 3;
//    //result.array[0] = P;
//    result.array[1] = L;
//    result.array[2] = U;
//    cout<<result;
    
    //return result ;
    cout<<"THis is P    "<<P<< endl;
    return 1;
}

// Determinant
//template <typename numberType>
//numberType Cmatrix<numberType>::LUP(Cmatrix<numberType> &A){

    
    




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
