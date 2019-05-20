//
//  CMatrix.cpp
//  CMatrix
//
//  Created by Juanita Gomez on 4/30/19.
//  Copyright © 2019 Juanita Gomez. All rights reserved.
//


#ifdef CMatrix_hpp

//-------------------------------CONSTRUCTORS AND DESTRUCTOR----------------------------------------------------

/**
  * It does not receive any parameter and creates an object of empty type Cmatrix with
  * length 0 and capacity equal to Initial Capacity.
 */

// Empty Constructor
template <typename numberType>
Cmatrix<numberType>::Cmatrix() {
    nRows = 0;
    nCols = 0;
    capacity = Initial_Capacity;
    array = new Cvector<numberType> [capacity];
    Checkrep();
}

/**
 * Receives a size_t Length and a numberType value and creates an object of type Empty
 * matrix with length length, capacity equal to Initial Capacity and each of its
 */

// Fill Constructor
template <typename numberType>
Cmatrix<numberType>:: Cmatrix(size_t size, const Cvector<numberType> &x, bool axis){
    if (axis == false){
        capacity = size + Initial_Capacity;
        array = new Cvector<numberType> [capacity];
        nRows = size;
        nCols = x.size();
        for (size_t i = 0; i < nRows; i++) array[i] = x;
    }
    
    else if (axis == true){
        numberType b = numberType{};
        Cvector<numberType> tmp(size,b);
        capacity = x.size() + Initial_Capacity;
        array = new Cvector<numberType> [capacity];
        nRows = x.size();
        nCols = size;
        for (size_t i = 0; i < nRows; i++) array[i] = tmp;
        for (size_t i = 0; i < nRows; i++){
            for (size_t j = 0; j < nCols; j++){
                (array[i])[j] = x[i];
            }
        }
    }
    Checkrep();
}

/**
 * It receives a Cmatrix rhs and creates an object of type Cmatrix by copying the length,
 * capacity and elements of rhs.
 */

// Parametric Constructor
template <typename numberType>
Cmatrix<numberType>::Cmatrix(const Cmatrix<numberType> &x){
    capacity = x.capacity + Initial_Capacity;
    nRows = x.nRows;
    nCols = x.nCols;
    array = new Cvector<numberType> [capacity];
    for (size_t i = 0; i < nRows; i++) array[i] = x.array[i];
    Checkrep();
}

/**
 * It receives a size_t row and a size_t col and a bool. If the bool is true, and the values of
 * row and col are equal, it constructs an identity matrix. If the bool is false it constructs a
 * matrix of size row x col filled with zeros.
 */

//Specialized Constructor
template <typename numberType>
Cmatrix<numberType>:: Cmatrix(size_t row, size_t col, bool type){
    if (type == true) assert(row == col);
    nRows = row;
    nCols = col;
    capacity = nRows + Initial_Capacity;
    array = new Cvector<numberType> [capacity];
    Cvector<numberType> tmp(nCols,0);
    
    if (type == true){
        for (size_t i = 0; i < nRows ; i++){
            tmp[i] = 1;
            array[i] = tmp;
            tmp[i] = 0;
        }
    }
    else{
        for (size_t i = 0; i < nRows ; i++){
            array[i] = tmp;
            }
        }
    Checkrep();
}

/**
 * 
 *
 */

// Destructor
template <typename numberType>
Cmatrix<numberType>:: ~Cmatrix(){
    nRows = 0;
    nCols = 0;
    delete [] array;
    Checkrep();
}

//------------------------------------- OPERATORS ---------------------------------------------------

//----------------------------------Class Member operators ------------------------------------------

// Operator =
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::operator=(const Cmatrix<numberType> &rhs){
    Checkrep();
    Cvector<numberType> *Oldarray = this -> array;
    capacity = rhs.capacity + Initial_Capacity;
    array = new Cvector<numberType> [capacity];
    nRows = rhs.nRows;
    nCols = rhs.nCols;
    for (size_t i = 0; i < rhs.nRows; i++) array[i] = rhs.array[i];
    delete[] Oldarray;
    Checkrep();
    return *this;
}

// Operator ()
template <typename numberType>
numberType Cmatrix<numberType>::operator () (size_t row, size_t col) const{
    Checkrep();
    return array[row][col];
}

// Operator ()
template <typename numberType>
numberType& Cmatrix<numberType>::operator () (size_t row, size_t col){
    Checkrep();
    return array[row][col];
}


// Operator []
template <typename numberType>
Cvector<numberType> Cmatrix<numberType>::operator [] (size_t idx) const {
    Checkrep();
    return array[idx];
}

// Operator []
template <typename numberType>
Cvector<numberType> & Cmatrix<numberType>::operator [](size_t idx){
    Checkrep();
    return array[idx];
}

//-------------------------------------Friend operators ------------------------------------------------

// Operator ==
template <typename numberType>
Cmatrix<bool> operator == (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert((x.nRows == y.nRows)&&(x.nCols == y.nCols));
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] == y.array[i]));
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator !=
template <typename numberType>
Cmatrix<bool> operator != (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert((x.nRows == y.nRows)&&(x.nCols == y.nCols));
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] != y.array[i]));
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator >
template <typename numberType>
Cmatrix<bool> operator > (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert((x.nRows == y.nRows)&&(x.nCols == y.nCols));
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] > y.array[i]));
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator <
template <typename numberType>
Cmatrix<bool> operator < (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert((x.nRows == y.nRows)&&(x.nCols == y.nCols));
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] < y.array[i]));
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator >=
template <typename numberType>
Cmatrix<bool> operator >= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert((x.nRows == y.nRows)&&(x.nCols == y.nCols));
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] >= y.array[i]));
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator <=
template <typename numberType>
Cmatrix<bool> operator <= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert((x.nRows == y.nRows)&&(x.nCols == y.nCols));
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] <= y.array[i]));
    }
    x.Checkrep();
    y.Checkrep();
    return result;
    
}

// Operator <<
template <typename numberType>
ostream & operator<<(ostream &os, Cmatrix<numberType> &rhs) {
    rhs.Checkrep();
    os << "[ " << endl;
    for (size_t i = 0; i < rhs.nRows; i++) cout << rhs.array[i] << "";
    os << "]" << endl;
    rhs.Checkrep();
    return os;
}

//------------------------------------- Binary operators ------------------------------------------------

// Operator +
template <typename numberType>
Cmatrix<double> operator + (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert((x.numberCols() == y.numberCols()) && (x.numberRows() == y.numberRows()));
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> copia1 = y.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] + copia1[i]);
        }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator -
template <typename numberType>
Cmatrix<double> operator - (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert((x.numberCols() == y.numberCols()) && (x.numberRows() == y.numberRows()));
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> copia1 = y.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] - copia1[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator * Matrix - Escalar Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cmatrix<numberType> &x, const numberType &y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] * double(y));
    }
    x.Checkrep();
    return result;
}

// Operator * Escalar - Matrix Multiplication
template <typename numberType>
Cmatrix<double> operator * (const numberType &y, const Cmatrix<numberType> &x){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] * double(y));
    }
    x.Checkrep();
    return result;
}

// Operator Matrix-Vector Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cmatrix<numberType> &x, const Cvector<numberType> &y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(x.numberCols() == y.size());
    Cvector<double> result(y.size());
    for (int i = 0; i < x.numberRows(); i++){
        for (int j = 0; j < x.numberCols(); j++){
            result[i] += x(i,j) * y[j];
        }
    }
    Cmatrix<double> out(1,result,1);
    x.Checkrep();
    return out;
}

// Operator Vector-Matrix Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cvector<numberType> &x, const Cmatrix<numberType> &y){
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(y.numberRows() == x.size());
    Cvector<double> result(x.size());
    for (int i = 0; i < y.numberCols(); i++){
        for (int j = 0; j < x.size(); j++){
            result[j] += y(i,j) * x[i];
        }
    }
    Cmatrix<double> out(1,result,0);
    y.Checkrep();
    return out;
}

// Operator * Matrix Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(x.numberCols() == y.numberRows());
    Cmatrix<numberType> result(x.numberRows(), y.numberCols(), false);
    for (int i = 0; i < x.numberRows(); i++){
        for (int j = 0; j < y.numberCols(); j++){
            for (int k = 0; k < y.numberRows(); k++){
                result(i,j) += x(i,k) * y(k,j);
            }
        }
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator /
template <typename numberType>
Cmatrix<double> operator / (const Cmatrix<numberType> &x, const numberType &y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(y != 0);
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] / double(y));
    }
    x.Checkrep();
    return result;
}

// Operator ^
template <typename numberType>
Cmatrix<double> operator ^ (const Cmatrix<numberType> &x, const size_t &y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cmatrix<double> Mout = x.toDouble();
    Cmatrix<double> aux = x.toDouble();
    for(size_t i = 1; i < y; i++){
        Mout = Mout * aux;
    }
    x.Checkrep();
    return Mout;
}

//-------------------------------------CLASS METHODS-------------------------------------------------

//-------------------------------------Matrix Methods------------------------------------------------

/**
 * Converts the elements in the Cmatrix to double precision.
 */

// To Double precision
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>:: toDouble() const{
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cmatrix<double> aux;
    for(size_t i = 0; i < nRows; i++){
        aux.push((this -> array[i]).toDouble());
    }
    Checkrep();
    return aux;
}

/**
 * Receives a numberType value and inserts it at the end of the matrix. It does not return anything.
 */

// Push
template <typename numberType>
void Cmatrix<numberType>::push(const Cvector<numberType> &value, bool axis){
    Checkrep();
    if (axis == false){
        if(nCols != 0) assert(value.size() == nCols);
        if(nRows == capacity) expandCapacity();
        array[nRows++] = value;
        if (nCols == 0) nCols = value.size();
    }
    else if(axis == true){
        if(nRows != 0) assert(value.size() == nRows);
        (*this) = (*this).transpose();
        if(nRows == capacity) expandCapacity();
        array[nRows++] = value;
        if (nCols == 0) nCols = value.size();
        (*this) = (*this).transpose();
    }
    Checkrep();
}

/**
 * Receives a size_t index and deletes the matrix element corresponding to that index.
 * It does not return anything.
 */

// Erase
template <typename numberType>
void Cmatrix<numberType>::erase(size_t index, bool axis){
    Checkrep();
    if (axis == false){
        assert(index <= nRows-1);
        for (size_t i = index; i < nRows-1; i++) array[i] = array[i+1];
        nRows --;
    }
    else if (axis == true){
        assert(index <= nCols-1);
        (*this) = (*this).transpose();
        for (size_t i = index; i < nRows-1; i++) array[i] = array[i+1];
        nRows --;
        (*this) = (*this).transpose();
    }
    Checkrep();
}

/**
 * Receives a size_t index and a numberType value and inserts it into the position of the matrix
 * corresponding to that index. It does not return anything.
 */

// Insert
template <typename numberType>
void Cmatrix<numberType>::insert (size_t index, const Cvector<numberType> & value, bool axis){
    Checkrep();
    if (axis == false){
        assert((index >= 0 && index <= nCols)&&(value.size() == nCols));
        if (nRows == capacity) expandCapacity();
        for (size_t i = nRows; i > index; i--) array[i] = array[i-1];
        array[index] = value;
        nRows ++;
    }
    else if(axis == true){
        assert((index >= 0 && index <= nRows)&&(value.size() == nRows));
        (*this) = (*this).transpose();
        if(nRows == capacity) expandCapacity();
        for (size_t i = nRows; i > index; i--) array[i] = array[i-1];
        array[index] = value;
        nRows ++;
        (*this) = (*this).transpose();
    }
    Checkrep();
}

/**
 * It does not receive parameters and returns the number of rows of the matrix.
 */

// Number Rows
template <typename numberType>
size_t Cmatrix<numberType>::numberCols() const{
    Checkrep();
    return nCols;
}

/**
 * It does not receive parameters and returns the number of columns of the matrix.
 */

// Number Cols
template <typename numberType>
size_t Cmatrix<numberType>::numberRows() const{
    Checkrep();
    return nRows;
}

// Clear
template <typename numberType>
void Cmatrix<numberType>::clear(){
    Checkrep();
    nRows = 0;
    nCols = 0;
    Checkrep();
}

// Empty
template <typename numberType>
bool Cmatrix<numberType>::empty() const{
    Checkrep();
    return ((nRows == 0)&&(nCols == 0));
}

// Access
template <typename numberType>
numberType Cmatrix<numberType>::access(size_t row, size_t col) const{
    Checkrep();
    return this -> array[row][col];
}

//-------------------------------------Easy Matrixes------------------------------------------------

/**
 * Creates an identity matrix of size N x N
 */

// Identity Matrix
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::eye(size_t N) {
    Cmatrix<double> t(N, N, false);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            t(i,j) = ( i == j ? 1 : 0 );
        }
    }
    return t;
}

/**
 * Creates a matrix of zeros of size Rows x Cols
 */

// Zeros Matrix
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::zeros(size_t rows, size_t cols) {
    Cmatrix<double> t(rows, cols, false);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            t(i,j) = 0;
        }
    }
    return t;
}

/**
 * Creates a matrix of ones of size Rows x Cols
 */

// Ones Matrix
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::ones(size_t rows, size_t cols) {
    Cmatrix<double> t(rows, cols, false);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            t(i,j) = 1;
        }
    }
    return t;
}

/**
 * Creates a matrix of random numbers of size Rows x Cols
 */

// Random Matrix
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::random(size_t rows, size_t cols) {
    Cmatrix<double> t(rows, cols, false);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            t(i,j) = rand();
        }
    }
    return t;
}

/**
 *
 */

// Matrix from Diagonal
template <typename numberType>
Cmatrix<double> diagonalize(const Cvector<numberType> &rhs) {
    Cvector<double> copia = rhs.toDouble();
    Cmatrix<double> result = Cmatrix<numberType>::zeros(copia.size());
    for (size_t i = 0; i < copia.size(); i++) {
        result(i,i) = copia[i];
    }
    return result;
}


/**
 * Creates a vector with the elements of the diagonal of a matrix.
 */

// Vector from Diagonal
template <typename numberType>
Cvector<double> Cmatrix<numberType>::diagonal(Cmatrix<numberType> &m) {
    m.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(m.numberRows() == m.numberCols());
    Cmatrix<double> copia = m.toDouble();
    Cvector<double> v;
    for (size_t i = 0; i < copia.numberRows(); i++) {
        v.push(copia(i,i));
    }
    m.Checkrep();
    return v;
}

/**
 * Creates a permutation matrix from a vector where each element of the vector indicates a "1" in
 * the corresponding column index.
 */

// Permutation Matrix
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::permutationMatrix(Cvector<numberType> &v) {
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<int> v_int;
    for(size_t i = 0; i < v.size(); i++){
        assert((v[i] <= v.size()));
        assert((v[i] >= 0));
        v_int.push(int(v[i]));
    }
    Cmatrix<double> p;
    Cvector<double> copia_vector = v_int.toDouble();
    p = Cmatrix<numberType>::zeros(copia_vector.size()-1, copia_vector.size()-1);
    for (size_t i = 0; i < copia_vector.size(); i++) {
        p(i,copia_vector[i]) = double(1);
    }
    return p;
}

//-------------------------------------Modifiers------------------------------------------------

/**
 * Swaps rows row1 and row2.
 */


// Swap_r
template <typename numberType>
void Cmatrix<numberType>::swap_r(size_t row1, size_t row2){
    Checkrep();
    Cvector<numberType> tmp = array[row1];
    array[row1] = array[row2];
    array[row2] = tmp;
    Checkrep();
}

/**
 * Swaps rows col1 and col2.
 */

// Swap_c
template <typename numberType>
void Cmatrix<numberType>::swap_c(size_t col1, size_t col2){
    Checkrep();
    Cvector<numberType> tmp1 = (this -> transpose())[col1];
    Cvector<numberType> tmp2 = (this -> transpose())[col2];
    for (size_t i = 0; i < nRows; i++){
        (array[i])[col1] = tmp2[i];
        (array[i])[col2] = tmp1[i];
    }
    Checkrep();
}

//-------------------------------------Special Matrixes------------------------------------------------

/**
 * Returns a matrix with the absolute value elementwise of "this".
 */

// Abs
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::abs(){
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cmatrix<numberType> result(nRows, nCols, false);
    for (size_t i = 0; i < nRows; i++){
        for (size_t j = 0; j < nCols; j++){
            (result[i])[j] = fabs((array[i])[j]);
        }
    }
    Checkrep();
    return result;
}

/**
 * Returns the transpose of "this".
 */

// Transpose
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::transpose() {
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    size_t r = this->numberCols();
    size_t c = this->numberRows();
    Cmatrix<numberType> t(r, c, false);
    for (size_t i = 0; i < nCols; i++) {
        for(size_t j = 0; j < nRows; j++) {
            t(i,j) = (*this)(j,i);
        }
    }
    Checkrep();
    return t;
}

/**
 * Returns the inverse of "this".
 */

// Inverse
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::inverse(){
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    assert(this -> determinant() != 0);
    Cmatrix<double> A((*this).toDouble());
    Cmatrix<double> Ainv, Linv, Uinv, Pinv;
    tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> lup = A.LUP(double(0.0000001));
    Cvector<double> P = get<0>(lup);
    Cmatrix<double> L = get<1>(lup);
    Cmatrix<double> U = get<2>(lup);
    size_t N = nCols;
    Linv = Cmatrix<double>::zeros(N,N);
    Uinv = Cmatrix<double>::zeros(N,N);
    Ainv = Cmatrix<double>::zeros(N,N);
    Pinv = Cmatrix<double>::permutationMatrix(P);
    // Inverse of L
    for(int k=0; k<N; k++){
        Linv(k,k) = 1/L(k,k);
        for(int i=k+1; i<N; i++){
            Linv(i,k) = (-L[i].dot((Linv.transpose())[k]))/L(i,i);
        }
    }
    // Inverse of U
    for(int k=N-1; k>=0; --k){
        Uinv(k,k) = 1/U(k,k);
        for(long i=k-1; i>=0; i--){
            Uinv(i,k) = (-U[i].dot((Uinv.transpose())[k]))/U(i,i);
        }
    }
    Ainv = Uinv*Linv*Pinv;
    Checkrep();
    return Ainv;
}

/**
 * Returns a lowerTriangular matrix copying the corresponding elements of "this".
 */

// Lower Triangular
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::lowerTriangular() {
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    size_t r = numberRows();
    size_t c = numberCols();
    Cmatrix<double> t(r, c, false);
    for (size_t i = 0; i < numberRows(); i++) {
        for (size_t j = 0; j < numberCols(); j++) {
            t(i,j) = ( i >= j ? this->array[i][j] : 0 );
        }
    }
    Checkrep();
    return t;
}

/**
 * Returns a upperTriangular matrix copying the corresponding elements of "this".
 */

// Upper Triangular
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::upperTriangular() {
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    size_t r=numberRows();
    size_t c=numberCols();
    Cmatrix<double> t(r, c, false);
    for (size_t i = 0; i < numberRows(); i++) {
        for (size_t j = 0; j < numberCols() ; j++) {
            t(i,j) = ( i <= j ? this->array[i][j] : 0 );
        }
    }
    Checkrep();
    return t;
}

//-------------------------------------Matrix Decompositions-----------------------------------------

/**
 * Computes the LUP decomposition of the matrix "this". Returns a tuple composed by the
 * matrix L, the matrix U and a vector P that corresponds to the permutation matrix of
 * the decomposition.
 */

// LUP
template <typename numberType>
tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> Cmatrix<numberType>::LUP(double Tol) {
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    Cmatrix<double> A((*this).toDouble());
    size_t j, k, indexMax;
    double maxPivot, absA;
    Cmatrix<double> L;
    Cmatrix<double> U;
    size_t N = A.numberRows();
    Cvector<double> P;
    for (size_t i = 0; i <= N; i++)
        P.push(i);

    for (size_t i = 0; i < N; i++) {
        maxPivot = 0;
        indexMax = i;

        for (k = i; k < N; k++){
            absA = fabs(A[k][i]);
            if (absA > maxPivot) {
                maxPivot = absA;
                indexMax = k;
            }
        }
        assert(maxPivot > Tol);
        if (indexMax != i) {
            //Finding Pivots
            j = P[i];
            P[i] = P[indexMax];
            P[indexMax] = j;

            // Permuting rows
            Cvector<double> PivotVector;
            PivotVector = A[i];
            A[i] = A[indexMax];
            A[indexMax] = PivotVector;

            P[N]++;
        }
        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    U = A.upperTriangular();
    L = A-A.upperTriangular()+Cmatrix<double>::eye(N);
    Checkrep();
    return make_tuple(P, L, U);
}

/**
 * Computes the QR decomposition of the matrix "this". Returns a tuple composed by the
 * matrix Q, and the matrix R.
 */

// QR Decomposition
template <typename numberType>
tuple<Cmatrix<double>, Cmatrix<double>> Cmatrix<numberType>::QR() {
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    tuple<Cvector<double>, Cmatrix<double>> QR;
    Cmatrix<double> A((*this).toDouble());
    Cvector<Cvector<double>> cols;
    Cmatrix<double> copia_Atrans = (A.transpose()).toDouble();
    for(size_t i = 0; i<nCols; i++)
        cols.push(copia_Atrans[i]);
    Cvector<Cvector<double>> Qcols;
    Qcols = cols.gram_schmidt();
    // Q
    Cmatrix<double> Q;
    for(size_t i=0; i<nCols; i++){
        Q.push(Qcols[i]);
    }
    Cmatrix<double> copia_Qtrans = (Q.transpose()).toDouble();
    Q = copia_Qtrans;
    
    // R
    Cmatrix<double> R(nRows, nCols, false);
    for (size_t i = 0; i < numberRows(); i++) {
        for (size_t j = 0; j < numberCols(); j++) {
            R(i,j) = ( i <= j ?  Qcols[i].dot(cols[j]): 0 );
        }
    }
    Checkrep();
    return make_tuple(Q, R);
}

//-------------------------------------Matrix Properties------------------------------------------------

/**
 * Returns a double corresponding to the determinant of the matrix "this".
 */

// Determinant
template <typename numberType>
double Cmatrix<numberType>::determinant(){
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    Cmatrix<double> A((*this).toDouble());
    tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> lup = A.LUP(double(0.0000000001));
    Cvector<double> P = get<0>(lup);
    Cmatrix<double> L = get<1>(lup);
    Cmatrix<double> U = get<2>(lup);

    double detU = 1;
    double exp = P[P.size()-1]-(P.size()-1);
    double detP = pow(-1,exp);
    for (size_t i = 0; i<U.numberCols(); i++) detU*=U(i,i);
    double detA = detU*detP;
    Checkrep();
    return detA;
}

/**
 * Returns a vector containing the eigen values of the matriz by iterating the QR method.
 */

// EigenValues
template <typename numberType>
Cvector<double> Cmatrix<numberType>::eigen_values(const double tol){
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(nRows == nCols);
    Cmatrix<double> A((*this).toDouble());
    double det = A.determinant();
    cout << "Det :" << det << endl;
    double i= 1;
    while(i>tol){
        tuple<Cmatrix<double>, Cmatrix<double>> qr = A.QR();
        Cmatrix<double> Q(get<0>(qr));
        Cmatrix<double> R(get<1>(qr));
        A = R*Q;
        i = fabs(A(1,0));
        
    }
    Checkrep();
    return Cmatrix<double>::diagonal(A);
}

// ************************************************ PRIVATE *********************************************+*

/**
 * Expands the number of rows of the Cmatrix by doubling it.
 */

//Expand Capacity
template <typename numberType>
void Cmatrix<numberType>::expandCapacity(){
    Checkrep();
    Cvector<numberType> *Oldarray = array;
    capacity *= 2;
    array = new Cvector<numberType>[capacity];
    for(size_t i = 0; i < nRows; i++){
        array[i] = Oldarray[i];
    }
    Checkrep();
    delete[] Oldarray;
}

/**
 * Check the invariant representation of the class.
 */

// Check Representation Invariant
template <typename numberType>
void Cmatrix<numberType>::Checkrep() const{
    assert(nRows >= 0 && nCols >= 0);
    assert(nRows <= capacity);
}




#endif // CMatrix_hpp
