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
    array = new numberType [capacity];
    length = size;
    for (size_t i = 0; i < length; i++) array[i] = value;
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
    capacity = rhs.capacity + Initial_Capacity;
    array = new numberType[capacity];
    for (size_t i = 0; i < rhs.length; i++){
        array[i] = rhs.array[i];
    }
    this -> length = rhs.length;
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
Cvector<double> operator * (const Cvector<numberType> &x, const numberType &y){
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result;
    for (size_t i = 0; i < copia.length; i++){
        result.push(copia.array[i] * double (y));
    }

    return result;
}

// Operator *
template <typename numberType>
Cvector<double> operator * (const numberType &y, const Cvector<numberType> &x){
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result;
    for (size_t i = 0; i < copia.length; i++){
        result.push(copia.array[i] * double(y));
    }
    
    return result;
}

// Operator /
template <typename numberType>
Cvector<double> operator / (const Cvector<numberType> &x, const numberType &y){
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result;
    // if y= 0 ERROR
    for (size_t i = 0; i < copia.length; i++){
        result.push(copia.array[i] / double(y));
    }

    return result;
}

// Operator ^
template <typename numberType>
Cvector<double> operator ^ (const Cvector<numberType> &x, const int y){
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> result(copia);
	result.length = copia.length;
    for (size_t i = 0; i < size_t(y); i++){
        result = (result * double(y));
    }

    return result;
}

//---------------------------------------Class Methods------------------------------------------

// Dot Product
template <typename numberType>
double Cvector<numberType>::dot(Cvector<numberType> arr){
    Cvector<double> copia;
    copia = arr.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    size_t size1 = copia1.length;
    size_t size2 = copia.length;
    if(size1 != size2){
        return -1;
    }

    double out = 0;
    for(size_t i = 0; i < size1; i++){
        out += copia1.array[i] * copia.array[i];
    }
    return out;
}

// Cross Product
template <typename numberType>
Cvector<double> Cvector<numberType>::cross(Cvector<numberType> arr){
    Cvector<double> copia;
    copia = arr.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    Cvector<double> out;
    out.length = copia.length;
    out.push((copia1.array[1] * copia.array[2]) - (copia.array[1] * copia1.array[2]));
    double x = ((copia1.array[0] * copia.array[2]) - (copia.array[0] * copia1.array[2]));
    out.push(-1 * x);
    out.push((copia1.array[0] * copia.array[1]) - (copia.array[0] * copia1.array[1]));
    return out;
}

// Norm
template <typename numberType>
double Cvector<numberType>::norm(){
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    double num = 0;
    for(size_t i = 0; i < copia1.length; i++){
        num += pow(copia1.array[i], double(2));
    }
    double out = double(sqrt(num));
    return out;
}

// Normalize
template <typename numberType>
Cvector<double> Cvector<numberType>::normalize(){
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    double norm1 = copia1.norm();
    Cvector<double> out;
    for(size_t i = 0; i < copia1.length; i++){
        out.push(copia1.array[i] / norm1);
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
double Cvector<numberType>::angle(Cvector<numberType> &x){
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    double result;
    double dot =  copia1.dot(copia);
    result = acos(dot / (copia1.norm() * copia.norm()));
    return result;
}

// Projection
template <typename numberType>
Cvector<double> Cvector<numberType>::proy(Cvector<numberType> &x){
    Cvector<double> copia;
    copia = x.toDouble();
    Cvector<double> copia1;
    copia1 = this -> toDouble();
    Cvector<double> result;
    double T = (copia1.dot(copia) / pow(copia1.norm(),double(2)));
    for (size_t i = 0; i < copia1.length; i++){
        result.push(copia1.array[i] * T);
    }
    return result;
}


// Gram_schmidt
template <typename numberType>
Cvector<Cvector<double>> Cvector<numberType>::gram_schmidt(const Cvector<numberType> &rhs){
  Cvector<Cvector<double>> copia;
    for (size_t i = 0; i < rhs.length; i++){
        copia.push(rhs.array[i].toDouble());
    }
  Cvector<Cvector<double>> result;
  result.push(copia[0].normalize());
  Cvector<double> tmp;
  for (size_t i = 1; i < copia.length; i++){
    tmp = copia[i];
    for (size_t j = 0; j < result.length; j++){
        tmp = tmp - (result[j] * (copia[i].dot(result[j])));
      }
    result.push(tmp.normalize());
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


// ---------------------------------------To double precision -------------------------------------

template <typename numberType>
Cvector<double> Cvector<numberType>:: toDouble() const{
    Cvector<double> aux;
    for(size_t i = 0; i < this -> size(); i++){
        aux.push(double(array[i]));
    }
    return aux;
}

///////////////////////////////////////////////////////////////

#endif //CVector_hpp
