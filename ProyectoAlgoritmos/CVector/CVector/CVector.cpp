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
    for (size_t i = 0; i < rhs.length; i++){
        array[i] = rhs.array[i];
    }
    this -> length = rhs.length;
    this -> capacity = rhs.capacity + Initial_Capacity;
    delete[] Oldarray;
    return *this;
}

// Operator []
template <typename numberType>
numberType Cvector<numberType>::operator [] (size_t idx) const {
    return array[idx];
}

// Operator []
template <typename numberType>
numberType & Cvector<numberType>::operator [](size_t idx){
    return array[idx];
}

//// Operator []=
//template <typename numberType>
//void Cvector<numberType>::operator []= (size_t idx, numberType value){
//    array[idx] = value;
//}

//-------------------------------------Friend operators ------------------------------------------
//// Comparisson Operators

// Operator <<
template <typename numberType>
ostream & operator<<(ostream &os, const Cvector<numberType> &rhs) {
    os << "[ ";
    for (size_t i = 0; i < rhs.length; i++) cout << rhs.array[i] << " ";
    os << "]" << endl;
    
    return os;
}

// Operator ==
template <typename numberType>
Cvector<bool> operator == (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    result.length = x.length;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] == y.array[i]);
    }
    return result;
}
// Operator !=
template <typename numberType>
Cvector<bool> operator != (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    result.length = x.length;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] != y.array[i]);
    }
    return result;
}

// Operator <=
template <typename numberType>
Cvector<bool> operator <= (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    result.length = x.length;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] <= y.array[i]);
    }
    return result;
}

// Operator >=
template <typename numberType>
Cvector<bool> operator >= (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    result.length = x.length;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] >= y.array[i]);
    }
    return result;
}

// Operator <
template <typename numberType>
Cvector<bool> operator < (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    result.length = x.length;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = (x.array[i] < y.array[i]);
    }
    return result;
}

// Operator >
template <typename numberType>
Cvector<bool> operator > (const Cvector<numberType> &x , const Cvector<numberType> &y){
    Cvector<bool> result;
    result.length = x.length;
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
    if (x.length == y.length){
        for (size_t i = 0; i < x.length; i++){
            result.array[i] = x.array[i] + y.array[i];
            }
        result.length = x.length;
        }
    
    else if (x.length < y.length){
        for (size_t i = 0; i < x.length; i++){
            result.array[i] = x.array[i] + y.array[i];
            }
        
        for (size_t i = x.length; i < y.length; i++){
            result.array[i] = y.array[i];
            }
        result.length = y.length;
        }
    
    else if (x.length > y.length){
        for (size_t i = 0; i < y.length; i++){
            result.array[i] = x.array[i] + y.array[i];
        }
        
        for (size_t i = y.length; i < x.length; i++){
            result.array[i] = x.array[i];
        }
        result.length = x.length;
    }
    
    return result;
}

// Operator -
template <typename numberType>
Cvector<numberType> operator - (const Cvector<numberType> &x, const Cvector<numberType> &y){
    Cvector<numberType> result;
    if (x.length == y.length){
        for (size_t i = 0; i < x.length; i++){
            result.array[i] = x.array[i] - y.array[i];
        }
        result.length = x.length;
    }
    
    
    else if (x.length < y.length){
        for (size_t i = 0; i < x.length; i++){
            result.array[i] = x.array[i] - y.array[i];
        }
        
        for (size_t i = x.length; i < y.length; i++){
            result.array[i] = 0 - y.array[i];
        }
        result.length = y.length;
    }
    
    else if (x.length > y.length){
        for (size_t i = 0; i < y.length; i++){
            result.array[i] = x.array[i] - y.array[i];
        }
        
        for (size_t i = y.length; i < x.length; i++){
            result.array[i] = 0 - x.array[i];
        }
        result.length = x.length;
    }
    
    return result;
}

// Operator *
template <typename numberType>
Cvector<numberType> operator * (const Cvector<numberType> &x, const double &y){
    Cvector<numberType> result;
    result.length = x.length;
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = numberType (x.array[i] * y);
    }
    
    return result;
}

// Operator /
template <typename numberType>
Cvector<numberType> operator / (const Cvector<numberType> &x, const double &y){
    Cvector<numberType> result;
	result.length = x.length;
    // if y= 0 ERROR
    for (size_t i = 0; i < x.length; i++){
        result.array[i] = x.array[i] / y;
    }
    
    return result;
}

// Operator ^
template <typename numberType>
Cvector<numberType> operator ^ (const Cvector<numberType> &x,  int y){
    Cvector<numberType> result(x);
	result.length = x.length;
    for (int i = 0; i < y; i++){
        result =  (result * y);
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
    for(size_t i = 0; i < size1; i++){
        out += this -> array[i] * arr.array[i];
    }

    return out;
}

// Cross Product
template <typename numberType>
Cvector<numberType> Cvector<numberType>::cross(Cvector<numberType> arr){
    Cvector<numberType> out;
    out.length = arr.length;
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
    for (size_t i = index; i < length-1; i++) array[i] = array[i+1];
    length --;
}

// Insert
template <typename numberType>
void Cvector<numberType>::insert (size_t index, numberType value){
    if (length == capacity) expandCapacity();
    for (size_t i = length; i > index; i--) array[i] = array[i-1];
    array[index] = value;
    length++;
}

// Clear
template <typename numberType>
void Cvector<numberType>::clear(){
    length = 0;
}

// Empty
template <typename numberType>
bool Cvector<numberType>::empty() const{
    return (length == 0);
}

// Size
template <typename numberType>
size_t Cvector<numberType>::size() const{
    return length;
}


// Angle
template <typename numberType>
float Cvector<numberType>::angle(Cvector<numberType> &x){
    float result;
    numberType dot = this -> dot(x);
    result = acos((float(dot) / (float(this -> norm()) * float(x.norm()))));
    return result;
}

// Projection
template <typename numberType>
Cvector<numberType> Cvector<numberType>::proy(Cvector<numberType> &x){
    Cvector<numberType> result;
    result.length = x.length;
    double T = ((double(this -> dot(x))) / (double(pow(this -> norm(),2))));
    for (size_t i = 0; i < this -> length; i++){
        result.array[i] = numberType((this -> array[i]) * T);
    }
    return result;
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

