//  CVector.cpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#ifdef CVector_hpp

//------------------------------------- CONSTRUCTOR -----------------------------------------

// Empty
template <typename numberType>
Cvector<numberType>::Cvector(){
    capacity = Initial_Capacity;
    array = new numberType[capacity];
    length = 0;
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
}

// Fill Constructor
template <typename numberType>
Cvector<numberType>::Cvector(size_t size, numberType value){
    capacity = Initial_Capacity;
    array = new numberType[capacity];
    length = size;
    for (size_t i = 0; i < size; i++) array[i] = value;
}

//------------------------------------- DESTRUCTOR -----------------------------------------

// Destructor
template <typename numberType>
Cvector<numberType>::~Cvector(){
    delete[] array;
}

//------------------------------------- OPERATORS ------------------------------------------

//-----------------------------------Non_friend operators ----------------------------------
// Operator =
template <typename numberType>
Cvector<numberType>  Cvector<numberType>::operator=(const Cvector &rhs){
    numberType *Oldarray = this -> array;
    this -> array = new numberType[rhs.capacity + Initial_Capacity];
    for (int i = 0; i < rhs.length; i++){
        array[i] = rhs.array[i];
    }
    this -> length = rhs.length;
    this -> capacity = rhs.capacity + Initial_Capacity;
    delete[] Oldarray;
    return *this;
}

// Operator []
template <typename numberType>
numberType Cvector<numberType>::operator [] (size_t idx){
    
    return array[idx];
}

//-------------------------------------Friend operators ------------------------------------------
//// Comparisson Operators

// Operator ==
template <typename numberType>
Cvector<bool> operator == (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    if(x.length == y.length){
        for (size_t i = 0; i < x.length; i++){
            result.array[i] = (x.array[i] == y.array[i]);
        }
    }
    return result;
}

// Operator !=
template <typename numberType>
Cvector<bool> operator != (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] != y.array[i]);
    }
    return result;
}

// Operator <=
template <typename numberType>
Cvector<bool> operator <= (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] <= y.array[i]);
    }
    return result;
}

// Operator >=
template <typename numberType>
Cvector<bool> operator >= (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] >= y.array[i]);
    }
    return result;
}

// Operator <
template <typename numberType>
Cvector<bool> operator < (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] < y.array[i]);
    }
    return result;
}

// Operator >
template <typename numberType>
Cvector<bool> operator > (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] > y.array[i]);
    }
    return result;
}

////Binary operators

// Operator +
template <typename numberType>
Cvector<numberType> operator + (const Cvector<numberType> &x, const Cvector<numberType> &y){
    Cvector<numberType> result;
    size_t len = max(x.length, y.length);
    for (size_t i = 0; i < len; i++){
        result.array[i] = x.array[i] + y.array[i];
    }
    
    return result;
}

// Operator -
template <typename numberType>
Cvector<numberType> operator - (const Cvector<numberType> &x, const Cvector<numberType> &y){
    Cvector<numberType> result;
    size_t len = max(x.length, y.length);
    for (size_t i = 0; i < len; i++){
        result.array[i] = x.array[i] - y.array[i];
    }
    
    return result;
}

// Operator *
template <typename numberType>
Cvector<numberType> operator * (const Cvector<numberType> &x, const Cvector<numberType> &y){
    Cvector<numberType> result;
    size_t len = max(x.length, y.length);
    for (size_t i = 0; i < len; i++){
        result.array[i] = x.array[i] * y.array[0];
    }
    
    return result;
}

// Operator /
template <typename numberType>
Cvector<numberType> operator / (const Cvector<numberType> &x, const Cvector<numberType> &y){
    Cvector<numberType> result;
    
    // if y= 0 ERROR
    size_t len = max(x.length, y.length);
    for (size_t i = 0; i < len; i++){
        result.array[i] = x.array[i] / y.array[0];
    }
    
    return result;
}

// Operator ^
template <typename numberType>
Cvector<numberType> operator ^ (const Cvector<numberType> &x, const Cvector<numberType> &y){
    Cvector<numberType> result(x);
    for (size_t i = 0; i < y.array[0]; i++){
        result.array *= x.array;
    }
    
    return result;
}

//---------------------------------------Class Methods------------------------------------------


// Dot Product
template <typename numberType>
double Cvector<numberType>::dot(Cvector<numberType> arr){
    size_t size1 = this -> length;
    size_t size2 = arr.length;

    if(size1 != size2){
        return -1;
    }

    double out = 0;
    for(int i = 0; i < size1; i++){
        out += this -> array[i] * arr.array[i];
    }

    return out;
}

// Cross Product
template <typename numberType>
Cvector<numberType> Cvector<numberType>::cross(Cvector<numberType> arr){
    Cvector<numberType> out;

    out.push((this -> array[1] * arr.array[2]) - (arr.array[1] * this -> array[2]));
    double x = (this -> array[0] * arr.array[2]) - (arr.array[0] * this -> array[2]);
    out.push(-1 * x);
    out.push((this -> array[0] * arr.array[1]) - (arr.array[0] * this -> array[1]));
    return out;
}

// Norm
template <typename numberType>
double Cvector<numberType>::norm(){
    int num = 0;
    for(size_t i = 0; i < this -> length; i++){
        num += pow(this -> array[i], 2);
    }
    double out = sqrt(num);
    return out;
}

// Normalize
template <typename numberType>
Cvector<numberType> Cvector<numberType>::normalize(){
    int norm1 = this -> norm();
    cout << norm1 << endl;
    Cvector<numberType> out;
    for(size_t i = 0; i < this -> length; i++){
        out.push(this -> array[i] / norm1);
    }
    return out;
}

// Push
template <typename numberType>
void Cvector<numberType>::push(numberType value){
    if(length == capacity){
        expandCapacity();
    }
    array[length++] = value;
}

// Erase
template <typename numberType>
void Cvector<numberType>::erase(size_t index){
    if(index == length - 1){
        --length;
    }
    else if(index < 0 || index > length){
        cout << "This index is not in the vector" << endl;
    }
    else{
        numberType *Oldarray = array;
        array = new numberType[capacity];
        for(size_t i = 0; i < index; i++){
            array[i] = Oldarray[i];
        }
        for(size_t i = index + 1; i < length; i++){
            array[i-1] = Oldarray[i];
        }
        delete[]Oldarray;
        --length;
    }
}

// Insert
template <typename numberType>
void Cvector<numberType>::insert (size_t index, numberType value){
    if(index == length) {
        push(value);
    }
    else{
        if(length == capacity) {
            expandCapacity();
        }
        numberType *Oldarray = array;
        array = new numberType[capacity];
        if(index == 0){
            array[0] = value;
            for (size_t i = 1; i <= length; i++){
                array[i] = Oldarray[i-1];
            }
        }
        else{
            for (size_t i = 0; i < index; i++) array[i] = Oldarray[i];
            array[index] = value;
            for (size_t i = index + 1; i <= length; i++) array[i] = Oldarray[i-1];
        }
        delete[]Oldarray;
        length++;
    }
}

// Clear
template <typename numberType>
void Cvector<numberType>::clear(){
    numberType *Oldarray = array;
    array = new numberType[capacity];
    length = 0;
    delete[]Oldarray;
}

// Empty
template <typename numberType>
bool Cvector<numberType>::empty() const{
    if(length == 0){
        return true;
    }
    else{
        return false;
    }
}

// Size
template <typename numberType>
size_t Cvector<numberType>::size() const{
    return length;
}

//--------------------------------------Expand Capacity--------------------------------------------

template <typename numberType>
void Cvector<numberType>::expandCapacity(){
    numberType *Oldarray = array;
    capacity *= 2;
    array = new numberType[capacity];
    for(size_t i = 0; i < length; i++){
        array[i] = Oldarray[i];
    }
    delete[] Oldarray;
    
}

#endif //CVector_hpp

