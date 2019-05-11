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
        numberType b = numberType{};
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
    Cvector<numberType> tmp(nCols,0);
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

//// Operator ()
//template <typename numberType>
//Cvector<numberType> Cmatrix<numberType>::operator () (size_t idx, bool type){
//    if (type == false){
//        return array[idx];
//    }
//    Cvector<numberType> tmp;
//    for (size_t i = 0; i < nRows; i++){
//        tmp.push(array[i][idx]);
//    }
//    return tmp;
//}

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
            result.array[i] = y.array[i] * (double(-1));
        }
        result.nRows = y.nRows;
    }

    else if (x.nRows > y.nRows){
        for (size_t i = 0; i < y.nRows; i++){
            result.array[i] = x.array[i] - y.array[i];
        }

        for (size_t i = y.nRows; i < x.nRows; i++){
            result.array[i] = x.array[i] * (double(-1));
        }
        result.nRows = x.nRows;
    }

    return result;
}



// Operator * Matrix - Escalar Multiplication
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
    Cmatrix<numberType> result(x.numberRows(), y.numberCols(), false);
    size_t xRows = x.numberRows();
    size_t yRows = y.numberRows();
    size_t xCols = x.numberCols();
    size_t yCols = y.numberCols();
    assert(xCols == yRows);

    for (int i = 0; i < xRows; i++)
    {
        for (int j = 0; j < yCols; j++)
        {
            for (int k = 0; k < yRows; k++)
            {
                result(i,j) += x(i,k) * y(k,j);
            }
        }
    }
    return result;
}

//// Operator * Matrix - Vector Multiplication

//template <typename numberType>
//Cvector<numberType> operator * (const Cvector<numberType> &v, Cmatrix<numberType> &x){
//
//    Cvector<numberType> result(x.numberRows(),0);
//    size_t vRows = 1;
//    size_t vCols = v.size();
//    size_t xRows = x.numberRows();
//    size_t xCols = x.numberCols();
//
//    assert(vCols == xRows);
//
//    for (size_t j = 0; j < vCols; j++)
//    {
//        for (size_t k = 0; k < xRows; k++)
//        {
//            result[k] += v[k] * x(k,j);
//        }
//    }
//
//    return result;
//}

//// Operator *  Vector - Matrix  Multiplication
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
//Cmatrix<numberType> operator ^ (const Cmatrix<numberType> &x, const numberType &y){
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
    nCols = value.size();
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
size_t Cmatrix<numberType>::numberRows() const{
    return nRows;
}

// nRows
template <typename numberType>
size_t Cmatrix<numberType>::numberCols() const{
    return nCols;
}

// Access
template <typename numberType>
numberType Cmatrix<numberType>::access(size_t row, size_t col) const{
    return this -> array[row][col];
}

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
    Cvector<numberType> tmp1;
    Cvector<numberType> tmp2;
    tmp1 = (this -> transpose()).array[col1];
    tmp2 = (this -> transpose()).array[col2];
    for (size_t i = 0; i < nRows; i++){
        (array[i])[col1] = tmp2[i];
        (array[i])[col2] = tmp1[i];
    }

}

// Abs
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::abs(){
    Cmatrix<numberType> result(nRows, nCols, false);
    for (size_t i = 0; i < nRows; i++){
        for (size_t j = 0; j < nCols; j++){
            (result.array[i])[j] = fabs((array[i])[j]);
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
            t.array[i][j]=this->array[j][i];
        }
    }
    return t;
}

// Upper Triangular
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::upperTriangular() {
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

// Lower Triangular
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::lowerTriangular() {
    size_t r=numberRows();
    size_t c=numberCols();
    Cmatrix<numberType> t(r, c, false);
    for (size_t i = 0; i < numberRows(); i++) {
        for (size_t j = 0; j < numberCols(); j++) {
            t(i,j) = ( i >= j ? this->array[i][j] : 0 );
        }
    }
    return t;
}

// Identity
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::eye(size_t N) {
    Cmatrix<numberType> t(N, N, false);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            t(i,j) = ( i == j ? 1 : 0 );
        }
    }
    return t;
}

// Zeros
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::zeros(size_t rows, size_t cols) {
    Cmatrix<numberType> t(rows, cols, false);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            t(i,j) = 0;
        }
    }
    return t;
}

// Ones Matrix
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::ones(size_t rows, size_t cols) {
    Cmatrix<numberType> t(rows, cols, false);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            t(i,j) = 1;
        }
    }
    return t;
}

// Random
template <typename numberType>
Cmatrix<numberType> Cmatrix<numberType>::random(size_t rows, size_t cols) {
    Cmatrix<numberType> t(rows, cols, false);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            t(i,j) = rand();
        }
    }
    return t;
}

// LUP
template <typename numberType>
tuple<Cvector<numberType>, Cmatrix<numberType>, Cmatrix<numberType>> Cmatrix<numberType>::LUP(double Tol) {

    Cmatrix<numberType> A(*this);
    size_t j, k, indexMax;
    double maxPivot, absA;
    Cmatrix<double> L;
    Cmatrix<double> U;

    //assert(A.numberRows()!=A.numberCols());
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

        //if (maxPivot < Tol) return 0; //failure, matrix is degenerate

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

// Determinant
template <typename numberType>
double Cmatrix<numberType>::determinant(){
    Cmatrix<numberType> A(*this);
    tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> lup = A.LUP(0.0001);
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

//--------------------------------------To double precision---------------------------------------

template <typename numberType>
Cmatrix<double> Cmatrix<numberType>:: toDouble(){
    Cmatrix<double> aux;
    for(unsigned i = 0; i < this -> nRows; i++){
        aux.push((this -> array[i]).toDouble());
    }
    *this = aux;
    return *this;
}

//--------------------------------------Matrix Product--------------------------------------------

template <typename numberType>
Cmatrix<numberType> product(Cmatrix<numberType> x, Cmatrix<numberType> y){

    unsigned fX = x.nRows;
    unsigned cX = x.nCols;

    unsigned fY = y.nRows;
    unsigned cY = y.nCols;

    Cmatrix<numberType> Mout;

    cout << "Las matrices a multiplicar son: " << endl;
    cout << "Matriz 1: " << fX << " x " << cX << endl;
    cout << endl;
    cout << x << endl;

    cout << "Matriz 2: " << fY << " x " << cY << endl;
    cout << endl;
    cout << y << endl;

    if(cX != fY){
        cout << "No se puede realizar el producto entre matrices. " << endl;
    }
    else{

        Cmatrix<numberType> aux = y.transpose();

        Cmatrix<numberType> twoMatr;

        for(unsigned i = 0; i < x.array -> size(); i++){
            twoMatr.push(x.array[i]);
        }

        for(unsigned i = 0; i < y.array -> size(); i++){
            twoMatr.push(aux.array[i]);
        }
        // cout << twoMatr << endl;

        // twoMatr.transpose();

        // cout << twoMatr << endl;

        int n = 0;
        // cout << Mout.array -> size() << endl;
        // cout << x.array -> size() << endl;
        // cout << y.array -> size() << endl;
        numberType d;
        while(Mout.nRows < x.nRows){
            Cvector<numberType> out;
            for(unsigned j = 0; j < aux.nRows; j++){
                numberType cont = 0;
                for(unsigned k = 0; k < twoMatr.array[0].size(); k++){
                    d = twoMatr.array[n][k] * twoMatr.array[j + (x.array -> size())][k];
                    cont += d;
                }
                out.push(cont);

                // cout << cont << endl;
                // cout << out << endl;
            }
            n++;
            Mout.push(out);
            // cout << Mout.array -> size() << endl;
            // cout << Mout << endl;
        }

        // Mout = Mout.transpose();

        cout << "El resultado es una matriz: " << fX << " x " << cY << endl;
        cout << endl;
    }

    return Mout;
}


/////////////////////////////////////////////////////////////////

#endif // CMatrix_hpp
