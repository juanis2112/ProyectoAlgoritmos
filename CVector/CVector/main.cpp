// main.cpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#include "CVector.hpp"
using namespace std;
int main(){
    Cvector<int> p(3,100.0);
    Cvector<int> q(p);
	for (int i = 0; i < 10 ; i++) q.insert(0,i);
    Cvector<int> r;
    r = p^3;
    cout << p << endl;
	cout << q << endl;
    cout << r << endl;
  
    return 0;
}
