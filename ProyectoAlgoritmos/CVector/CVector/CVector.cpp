//  CVector.cpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#ifdef CVector_hpp

//------------------------------------- CONSTRUCTOR -----------------------------------------

// Empty Constructor
template <typename numberType>
Cvector<numberType>::Cvector(){
    capacity = Initial_Capacity;
    array = new numberType[capacity];
    length = 0;
    Checkrep();
    
}

// Parametric Constructor
template <typename numberType>
Cvector<numberType>::Cvector(const Cvector &rhs){
    capacity = rhs.capacity + Initial_Capacity;
    array = new numberType[capacity];
    for (size_t i = 0; i < rhs.length; i++){
        array[i] = rhs.array[i];
    }
    length = rhs.length;
    Checkrep();
}

// Capacity Constructor
template <typename numberType>
Cvector<numberType>::Cvector(size_t Length){
    capacity = Length + Initial_Capacity;
    array = new numberType[capacity];
    numberType tmp = numberType{};
    for (size_t i = 0; i < Length; i ++) array[i] = tmp;
    length = Length;
    Checkrep();
}

// Fill Constructor
template <typename numberType>
Cvector<numberType>::Cvector(size_t size, numberType value){
    capacity = Initial_Capacity;
    array = new numberType [capacity];
    length = size;
    for (size_t i = 0; i < length; i++) array[i] = value;
    Checkrep();
    
}

//------------------------------------- DESTRUCTOR -----------------------------------------

// Destructor
template <typename numberType>
Cvector<numberType>::~Cvector(){
    length = 0;
    delete[] array;
    Checkrep();
    
}

//------------------------------------- OPERATORS ------------------------------------------

//-----------------------------------Non_friend operators ----------------------------------
// Operator =
template <typename numberType>
Cvector<numberType> Cvector<numberType>::operator=(const Cvector &rhs){
    Checkrep();
    rhs.Checkrep();
    numberType *Oldarray = this -> array;
    capacity = rhs.capacity + Initial_Capacity;
    array = new numberType[capacity];
    for (size_t i = 0; i < rhs.length; i++){
        array[i] = rhs.array[i];
    }
    this -> length = rhs.length;
    delete[] Oldarray;
    Checkrep();
    rhs.Checkrep();
    return *this;
}

// Operator []
template <typename numberType>
numberType Cvector<numberType>::operator [] (size_t idx) const {
    Checkrep();
    return array[idx];
}

// Operator []
template <typename numberType>
numberType & Cvector<numberType>::operator [](size_t idx){
    Checkrep();
    return array[idx];
}
//-------------------------------------Friend operators ------------------------------------------
//// Comparisson Operators


// Operator <<
template <typename numberType>
ostream & operator<<(ostream &os, const Cvector<numberType> &rhs) {
    rhs.Checkrep();
    os << "[ ";
    for (size_t i = 0; i < rhs.length; i++) cout << rhs.array[i] << " ";
    os << "]" << endl;
    rhs.Checkrep();
    return os;
}

// Operator ==
template <typename numberType>
Cvector<bool> operator == (const Cvector<numberType> &x , const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<bool> result(x.size());
    for (size_t i = 0; i < x.length; i++){
        result[i] = (x.array[i] == y.array[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}
// Operator !=
template <typename numberType>
Cvector<bool> operator != (const Cvector<numberType> &x , const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<bool> result(x.size());
    for (size_t i = 0; i < x.length; i++){
        result[i] = (x.array[i] != y.array[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator <=
template <typename numberType>
Cvector<bool> operator <= (const Cvector<numberType> &x , const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<bool> result(x.size());
    for (size_t i = 0; i < x.length; i++){
        result[i] = (x.array[i] <= y.array[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
    
}

// Operator >=
template <typename numberType>
Cvector<bool> operator >= (const Cvector<numberType> &x , const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<bool> result(x.size());
    for (size_t i = 0; i < x.length; i++){
        result[i] = (x.array[i] >= y.array[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator <
template <typename numberType>
Cvector<bool> operator < (const Cvector<numberType> &x , const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<bool> result(x.size());
    for (size_t i = 0; i < x.length; i++){
        result[i] = (x.array[i] < y.array[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}


// Operator >
template <typename numberType>
Cvector<bool> operator > (const Cvector<numberType> &x , const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<bool> result(x.size());
    for (size_t i = 0; i < x.length; i++){
        result[i] = (x.array[i] > y.array[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

////Binary operators

// Operator +
template <typename numberType>
Cvector<numberType> operator + (const Cvector<numberType> &x, const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<numberType> result(x.size());
    for (size_t i = 0; i < x.length; i++)
        result[i] = x.array[i] + y.array[i];
    return result;
    x.Checkrep();
    y.Checkrep();
}

// Operator -
template <typename numberType>
Cvector<numberType> operator - (const Cvector<numberType> &x, const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<numberType> result(x.size());
    for (size_t i = 0; i < x.length; i++)
        result[i] = x.array[i] - y.array[i];
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator *
template <typename numberType>
Cvector<double> operator * (const Cvector<numberType> &x, const numberType &y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result;
    for (size_t i = 0; i < copia.length; i++){
        result.push(copia.array[i] * double (y));
    }
    x.Checkrep();
    return result;
}

// Operator *
template <typename numberType>
Cvector<double> operator * (const Cvector<numberType> &x, const Cvector<numberType> &y){
    x.Checkrep();
    y.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(x.length == y.length);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> copia1;
    copia1 = y.toDouble();
    Cvector<double> result(copia.size());
    for (size_t i = 0; i < copia.size(); i++){
        result[i] = (copia[i] * copia1[i]);
    }
    x.Checkrep();
    y.Checkrep();
    return result;
}

// Operator *
template <typename numberType>
Cvector<double> operator * (const numberType &y, const Cvector<numberType> &x){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result;
    for (size_t i = 0; i < copia.length; i++){
        result.push(copia.array[i] * double(y));
    }
    x.Checkrep();
    return result;
}

// Operator /
template <typename numberType>
Cvector<double> operator / (const Cvector<numberType> &x, const numberType &y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(y != 0);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result;
    for (size_t i = 0; i < copia.length; i++){
        result.push(copia.array[i] / y);
    }
    x.Checkrep();
    return result;
}

// Operator ^
template <typename numberType>
Cvector<double> operator ^ (const Cvector<numberType> &x, const int y){
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result(x.size());
    for (size_t i = 0; i < x.size(); i++){
        result[i] = pow(copia[i],double(y));
    }
    x.Checkrep();
    return result;
}

//---------------------------------------Class Methods------------------------------------------

// Dot Product
template <typename numberType>
double Cvector<numberType>::dot(Cvector<numberType> arr){
    Checkrep();
    arr.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(this -> length == arr.length);
    Cvector<double> copia;
    copia = arr.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    size_t size1 = copia1.size();
    size_t size2 = copia.size();
    if(size1 != size2){
        return -1;
    }

    double out = 0;
    for(size_t i = 0; i < size1; i++){
        out += copia1[i] * copia[i];
    }
    Checkrep();
    arr.Checkrep();
    return out;
}

// Cross Product
template <typename numberType>
Cvector<double> Cvector<numberType>::cross(Cvector<numberType> arr){
    Checkrep();
    arr.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    assert(this -> length == 3 && arr.length == 3);
    Cvector<double> copia;
    copia = arr.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    Cvector<double> out;
    out.push((copia1[1] * copia[2]) - (copia[1] * copia1[2]));
    double x = ((copia1[0] * copia[2]) - (copia[0] * copia1[2]));
    out.push(-1 * x);
    out.push((copia1[0] * copia[1]) - (copia[0] * copia1[1]));
    Checkrep();
    arr.Checkrep();
    return out;
}

// Norm
template <typename numberType>
double Cvector<numberType>::norm(){
    Checkrep();
    Cvector<double> copia1;
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    copia1 = this -> toDouble();
    double num = 0;
    for(size_t i = 0; i < copia1.size(); i++){
        num += pow(copia1[i], double(2));
    }
    double out = double(sqrt(num));
    Checkrep();
    return out;
}

// Normalize
template <typename numberType>
Cvector<double> Cvector<numberType>::normalize(){
    Checkrep();
    Cvector<double> copia1;
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    copia1 = this -> toDouble();
    double norm1 = copia1.norm();
    Cvector<double> out;
    for(size_t i = 0; i < copia1.size(); i++){
        out.push(copia1[i] / norm1);
    }
    Checkrep();
    return out;
}

// Push
template <typename numberType>
void Cvector<numberType>::push(numberType value){
    Checkrep();
    if(length == capacity){
        expandCapacity();
    }
    array[length++] = value;
    Checkrep();
}

// Erase
template <typename numberType>
void Cvector<numberType>::erase(size_t index){
    Checkrep();
    for (size_t i = index; i < length-1; i++) array[i] = array[i+1];
    length --;
    Checkrep();
}

// Insert
template <typename numberType>
void Cvector<numberType>::insert (size_t index, numberType value){
    Checkrep();
    if (length == capacity) expandCapacity();
    for (size_t i = length; i > index; i--) array[i] = array[i-1];
    array[index] = value;
    length++;
    Checkrep();
    
}

// Clear
template <typename numberType>
void Cvector<numberType>::clear(){
    Checkrep();
    length = 0;
    Checkrep();
    
}

// Empty
template <typename numberType>
bool Cvector<numberType>::empty() const{
    Checkrep();
    return (length == 0);
}

// Size
template <typename numberType>
size_t Cvector<numberType>::size() const{
    Checkrep();
    return length;
}


// Angle
template <typename numberType>
double Cvector<numberType>::angle(Cvector<numberType> &x){
    Checkrep();
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    double result;
    double dot =  copia1.dot(copia);
    double norms = copia1.norm() * copia.norm();
    result = acos(dot / norms);
    Checkrep();
    x.Checkrep();
    return result;
}

// Projection
template <typename numberType>
Cvector<double> Cvector<numberType>::proy(Cvector<numberType> &x){
    Checkrep();
    x.Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    Cvector<double> result;
    double T = (copia1.dot(copia) / pow(copia1.norm(),double(2)));
    for (size_t i = 0; i < copia1.size(); i++){
        result.push(copia1[i] * T);
    }
    Checkrep();
    x.Checkrep();
    return result;
}


// Gram_schmidt
template <typename numberType>
Cvector<Cvector<double>> Cvector<numberType>::gram_schmidt(){
    Checkrep();
    Cvector<Cvector<double>> copia;
    Cvector<numberType> rhs(*this);
    for (size_t i = 0; i < rhs.length; i++){
        copia.push(rhs.array[i].toDouble());
    }
    Cvector<Cvector<double>> result;
    result.push(copia[0].normalize());
    Cvector<double> tmp;
    for (size_t i = 1; i < copia.size(); i++){
        tmp = copia[i];
        for (size_t j = 0; j < result.size(); j++){
            tmp = tmp - (result[j] * (copia[i].dot(result[j])));
      }
        result.push(tmp.normalize());
    }
    Checkrep();
    return result;
}
//--------------------------------------Expand Capacity--------------------------------------------

template <typename numberType>
void Cvector<numberType>::expandCapacity(){
    Checkrep();
    numberType *Oldarray = array;
    capacity *= 2;
    array = new numberType[capacity];
    for(size_t i = 0; i < length; i++){
        array[i] = Oldarray[i];
    }
    delete[] Oldarray;
    Checkrep();
}


// --------------------------------------- To double precision -------------------------------------

template <typename numberType>
Cvector<double> Cvector<numberType>::toDouble() const{
    Checkrep();
    assert(!(is_same<numberType, char>::value));
    assert(is_arithmetic<numberType>::value);
    Cvector<double> aux;
    for(size_t i = 0; i < this -> size(); i++){
        aux.push(double(array[i]));
    }
    
    return aux;
}

template <typename numberType>
void Cvector<numberType>::Checkrep() const{
    assert(length >= 0);
    assert(length <= capacity);
}


///////////////////////////////////////////////////////////////

#endif //CVector_hpp
