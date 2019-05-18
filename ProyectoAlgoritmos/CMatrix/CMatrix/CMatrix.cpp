//
//  CMatrix.cpp
//  CMatrix
//
//  Created by Juanita Gomez on 4/30/19.
//  Copyright © 2019 Juanita Gomez. All rights reserved.
//


#ifdef CMatrix_hpp

//-------------------------------CONSTRUCTORS AND DESTRUCTOR----------------------------------------------------

// Empty Constructor
template <typename numberType>
Cmatrix<numberType>::Cmatrix() {
    nRows = 0;
    nCols = 0;
    capacity = Initial_Capacity;
    array = new Cvector<numberType> [capacity];
    Checkrep();
}

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
    return array[row][col];
}

// Operator ()
template <typename numberType>
numberType& Cmatrix<numberType>::operator () (size_t row, size_t col){
    return array[row][col];
}


// Operator []
template <typename numberType>
Cvector<numberType> Cmatrix<numberType>::operator [] (size_t idx) const {
    return array[idx];
}

// Operator []
template <typename numberType>
Cvector<numberType> & Cmatrix<numberType>::operator [](size_t idx){
    return array[idx];
}

//-------------------------------------Friend operators ------------------------------------------------

// Operator ==
template <typename numberType>
Cmatrix<bool> operator == (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] == y.array[i]));
    }
    return result;
}


// Operator !=
template <typename numberType>
Cmatrix<bool> operator != (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] != y.array[i]));
    }
    return result;
}
// Operator >
template <typename numberType>
Cmatrix<bool> operator > (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] > y.array[i]));
    }
    return result;
}

// Operator <
template <typename numberType>
Cmatrix<bool> operator < (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] < y.array[i]));
    }
    return result;
}

// Operator >=
template <typename numberType>
Cmatrix<bool> operator >= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] >= y.array[i]));
    }
    return result;
}

// Operator <=
template <typename numberType>
Cmatrix<bool> operator <= (const Cmatrix<numberType> &x , const Cmatrix<numberType> &y){
    Cmatrix<bool> result;
    for (size_t i = 0; i < x.nRows; i++){
        result.push((x.array[i] <= y.array[i]));
    }
    return result;
}

// Operator <<
template <typename numberType>
ostream & operator<<(ostream &os, Cmatrix<numberType> &rhs) {
    os << "[ " << endl;
    for (size_t i = 0; i < rhs.nRows; i++) cout << rhs.array[i] << "";
    os << "]" << endl;
    return os;
}

//------------------------------------- Binary operators ------------------------------------------------

// Operator +
template <typename numberType>
Cmatrix<double> operator + (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    assert((x.numberCols() == y.numberCols()) && (x.numberRows() == y.numberRows()));
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> copia1 = y.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] + copia1[i]);
        }
    return result;
}

// Operator -
template <typename numberType>
Cmatrix<double> operator - (const Cmatrix<numberType> &x, const Cmatrix<numberType> &y){
    assert((x.numberCols() == y.numberCols()) && (x.numberRows() == y.numberRows()));
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> copia1 = y.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] - copia1[i]);
    }
    return result;
}

// Operator * Matrix - Escalar Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cmatrix<numberType> &x, const numberType &y){
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] * double(y));
    }
    return result;
}

// Operator * Escalar - Matrix Multiplication
template <typename numberType>
Cmatrix<double> operator * (const numberType &y, const Cmatrix<numberType> &x){
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] * double(y));
    }
    return result;
}

// Operator Matrix-Vector Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cmatrix<numberType> &x, Cvector<numberType> &y){
    assert(x.numberCols() == y.size());
    Cvector<double> result(y.size());
    for (int i = 0; i < x.numberRows(); i++){
        for (int j = 0; j < x.numberCols(); j++){
            result[i] += x(i,j) * y[j];
        }
    }
    Cmatrix<double> out(1,result,1);
    return out;
}

// Operator Vector-Matrix Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cvector<numberType> &x, Cmatrix<numberType> &y){
    assert(y.numberRows() == x.size());
    Cvector<double> result(x.size());
    for (int i = 0; i < y.numberCols(); i++){
        for (int j = 0; j < x.size(); j++){
            result[j] += y(i,j) * x[i];
        }
    }
    Cmatrix<double> out(1,result,0);
    return out;
}

// Operator * Matrix Multiplication
template <typename numberType>
Cmatrix<double> operator * (const Cmatrix<numberType> &x, Cmatrix<numberType> &y){
    Cmatrix<numberType> result(x.numberRows(), y.numberCols(), false);
    assert(x.numberCols() == y.numberRows());
    for (int i = 0; i < x.numberRows(); i++){
        for (int j = 0; j < y.numberCols(); j++){
            for (int k = 0; k < y.numberRows(); k++){
                result(i,j) += x(i,k) * y(k,j);
            }
        }
    }
    return result;
}

// Operator /
template <typename numberType>
Cmatrix<double> operator / (const Cmatrix<numberType> &x, const numberType &y){
    assert(y != 0);
    Cmatrix<double> copia = x.toDouble();
    Cmatrix<double> result;
    for (size_t i = 0; i < copia.numberRows(); i++){
        result.push(copia[i] / double(y));
    }
    return result;
}

// Operator ^
template <typename numberType>
Cmatrix<double> operator ^ (const Cmatrix<numberType> &x, const size_t &y){
    Cmatrix<double> Mout = x.toDouble();
    Cmatrix<double> aux = x.toDouble();
    for(size_t i = 1; i < y; i++){
        Mout = Mout * aux;
    }
    return Mout;
}

//-------------------------------------CLASS METHODS-------------------------------------------------

//-------------------------------------Matrix Methods------------------------------------------------

// To Double precision
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>:: toDouble() const{
    Cmatrix<double> aux;
    for(size_t i = 0; i < nRows; i++){
        aux.push((this -> array[i]).toDouble());
    }
    return aux;
}

// Push
template <typename numberType>
void Cmatrix<numberType>::push(const Cvector<numberType> &value, bool axis){
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
}

// Erase
template <typename numberType>
void Cmatrix<numberType>::erase(size_t index, bool axis){
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
}

// Insert
template <typename numberType>
void Cmatrix<numberType>::insert (size_t index, const Cvector<numberType> & value, bool axis){
    if (axis == false){
        assert((index >= 0 && index <= nRows)&&(value.size() == nRows));
        if (nRows == capacity) expandCapacity();
        for (size_t i = nRows; i > index; i--) array[i] = array[i-1];
        array[index] = value;
        nRows ++;
    }
    else if(axis == true){
        assert((index >= 0 && index <= nCols)&&(value.size() == nCols));
        (*this) = (*this).transpose();
        if(nRows == capacity) expandCapacity();
        for (size_t i = nRows; i > index; i--) array[i] = array[i-1];
        array[index] = value;
        nRows ++;
        (*this) = (*this).transpose();
    }
}

// Number Rows
template <typename numberType>
size_t Cmatrix<numberType>::numberCols() const{
    return nCols;
}

// Number Cols
template <typename numberType>
size_t Cmatrix<numberType>::numberRows() const{
    return nRows;
}

// Clear
template <typename numberType>
void Cmatrix<numberType>::clear(){
    nRows = 0;
    nCols = 0;
}

// Empty
template <typename numberType>
bool Cmatrix<numberType>::empty() const{
    return ((nRows == 0)&&(nCols == 0));
}

// Access
template <typename numberType>
numberType Cmatrix<numberType>::access(size_t row, size_t col) const{
    return this -> array[row][col];
}

//-------------------------------------Easy Matrixes------------------------------------------------

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

// Matrix from Diagonal

// PENDIENTE -- Diagonalize

// Vector from Diagonal
template <typename numberType>
Cvector<double> Cmatrix<numberType>::diagonal(Cmatrix<numberType> m) {
    assert(m.numberRows() == m.numberCols());
    Cmatrix<double> copia = m.toDouble();
    Cvector<double> v;
    for (size_t i = 0; i < copia.numberRows(); i++) {
        v.push(copia(i,i));
    }
    return v;
}

// Permutation Matrix
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::permutationMatrix(Cvector<numberType> v) {
    Cmatrix<double> p;
    Cvector<double> copia_vector = v.toDouble();
    p = Cmatrix<numberType>::zeros(copia_vector.size()-1, copia_vector.size()-1);
    for (int i = 0; i < copia_vector.size(); i++) {
        p[i][copia_vector[i]] = 1;
    }
    return p;
}

//-------------------------------------Modifiers------------------------------------------------

// Swap_r
template <typename numberType>
void Cmatrix<numberType>::swap_r(size_t row1, size_t row2){
    Cvector<numberType> tmp = array[row1];
    array[row1] = array[row2];
    array[row2] = tmp;
}

// Swap_c
template <typename numberType>
void Cmatrix<numberType>::swap_c(size_t col1, size_t col2){
    Cvector<numberType> tmp1 = (this -> transpose())[col1];
    Cvector<numberType> tmp2 = (this -> transpose())[col2];
    for (size_t i = 0; i < nRows; i++){
        (array[i])[col1] = tmp2[i];
        (array[i])[col2] = tmp1[i];
    }
}

// Append Rows

// Append Cols

//-------------------------------------Special Matrixes------------------------------------------------

// Abs
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::abs(){
    Cmatrix<numberType> result(nRows, nCols, false);
    for (size_t i = 0; i < nRows; i++){
        for (size_t j = 0; j < nCols; j++){
            (result[i])[j] = fabs((array[i])[j]);
        }
    }
    return result;
}

// Transpose
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::transpose() {
    size_t r = this->numberCols();
    size_t c = this->numberRows();
    Cmatrix<numberType> t(r, c, false);
    for (size_t i = 0; i < nCols; i++) {
        for(size_t j = 0; j < nRows; j++) {
            t(i,j) = (*this)(j,i);
        }
    }
    return t;
}

// Inverse
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::inverse(){
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
    Pinv = permutationMatrix(P);
    // Inverse of L
    for(int k=0; k<N; k++){
        Linv(k,k) = 1/L(k,k);
        for(int i=k+1; i<N; i++){
            Linv(i,k) = (-L[i].dot((Linv.transpose())[k]))/L(i,i);
        }
    }
    // Inverse of U
    for(long k=N-1; k>=0; --k){
        size_t tmp = static_cast<unsigned>(k);
        Uinv(tmp,tmp) = 1/U(tmp,tmp);
        for(long i=tmp-1; i>=0; i--){
            Uinv(i,k) = (-U[i].dot((Uinv.transpose())[k]))/U(i,i);
        }
    }
    Ainv = Uinv*Linv*Pinv;
    return Ainv;
}

// Lower Triangular
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::lowerTriangular() {
    size_t r = numberRows();
    size_t c = numberCols();
    Cmatrix<double> t(r, c, false);
    for (size_t i = 0; i < numberRows(); i++) {
        for (size_t j = 0; j < numberCols(); j++) {
            t(i,j) = ( i >= j ? this->array[i][j] : 0 );
        }
    }
    return t;
}

// Upper Triangular
template <typename numberType>
Cmatrix<double> Cmatrix<numberType>::upperTriangular() {
    size_t r=numberRows();
    size_t c=numberCols();
    Cmatrix<double> t(r, c, false);
    for (size_t i = 0; i < numberRows(); i++) {
        for (size_t j = 0; j < numberCols() ; j++) {
            t(i,j) = ( i <= j ? this->array[i][j] : 0 );
        }
    }
    return t;
}

//-------------------------------------Matrix Decompositions-----------------------------------------

// LUP
template <typename numberType>
tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> Cmatrix<numberType>::LUP(double Tol) {
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
            absA = fabs(A.array[k][i]);
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
            PivotVector = A.array[i];
            A.array[i] = A.array[indexMax];
            A.array[indexMax] = PivotVector;

            P[N]++;
        }
        for (j = i + 1; j < N; j++) {
            A.array[j][i] /= A.array[i][i];

            for (k = i + 1; k < N; k++)
                A.array[j][k] -= A.array[j][i] * A.array[i][k];
        }
    }
    U = A.upperTriangular();
    L = A-A.upperTriangular()+Cmatrix<double>::eye(N);
    return make_tuple(P, L, U);
}

// QR Decomposition
template <typename numberType>
tuple<Cmatrix<double>, Cmatrix<double>> Cmatrix<numberType>::QR() {
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
    return make_tuple(Q, R);
}

// Singular Value Decomposition

//template <typename numberType>
//tuple<Cmatrix<numberType>, Cmatrix<numberType>> Cmatrix<numberType>::SVD() {
//    Cmatrix<numberType> A(*this);
//    Cmatrix<numberType> P;
//    int min;
//    if(A.numberCols>A.numberRows) {
//        P = A*A.transpose();
//        min = A.numberCols;
//    else {
//        P= A.transpose()*A;
//        min = A.numberRows;
//        }
//    tuple<Cmatrix<numberType>, Cmatrix<numberType>> qr = A.QR();
//    Cmatrix<double> Q(get<0>(qr));
//    Cmatrix<double> R(get<1>(qr));
//    A = R*Q;
//    i = fabs(A(1,0));
//}

//-------------------------------------Matrix Properties------------------------------------------------

// Determinant
template <typename numberType>
double Cmatrix<numberType>::determinant(){
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
    return detA;

}

// EigenValues
template <typename numberType>
Cvector<double> Cmatrix<numberType>::eigen_values(const double tol){
    Cmatrix<double> A((*this).toDouble());
    double i= 1;
    while(i>tol){
        tuple<Cmatrix<double>, Cmatrix<double>> qr = A.QR();
        Cmatrix<double> Q(get<0>(qr));
        Cmatrix<double> R(get<1>(qr));
        A = R*Q;
        i = fabs(A(1,0));
    }
    return Cmatrix<numberType>::diagonal(A);
}

// ************************************************ PRIVATE *********************************************+*

//Expand Capacity
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

// Check Representation Invariant
template <typename numberType>
void Cmatrix<numberType>::Checkrep() const{
    assert(nRows >= 0 && nCols >= 0);
    assert(nRows <= capacity);
}




#endif // CMatrix_hpp
