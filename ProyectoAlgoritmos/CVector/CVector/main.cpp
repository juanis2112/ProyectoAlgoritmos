// main.cpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#include "CVector.hpp"
using namespace std;
int main(){
    Cvector<double> p(3,100.0);
    Cvector<double> q(p);
	//for (int i = 0; i < 10 ; i++) q.insert(0,i);
    Cvector<double> r;
    r = p ^ 3.0;
    cout << p << endl;
	cout << q << endl;
    cout << r << endl;
    Cvector<double> p1;
    Cvector<double> p2;
    for (int i = 0; i < 5; i++){
        p1.insert(0,double(i));
        p2.insert(0,double(i+10));
    }


    cout << p1.angle(p2) << endl;
    cout << p1.proy(p2) << endl;
    cout << "p1: " << p1 << endl;
    cout << "p2 :" << p2 << endl;
    p2.erase(2);
    cout << "p2 :" << p2 << endl;


    Cvector<int> h;
    h.push(3);
    h.push(5);
    h.push(1);
    Cvector<Cvector<int>> temp;
    temp.push(h);
    h.clear();
    h.push(2);
    h.push(1);
    h.push(0);
    temp.push(h);
    h.clear();
    h.push(0);
    h.push(1);
    h.push(0);;
    temp.push(h);
    cout << temp << endl;
    Cvector<Cvector<double>> temp1;
    temp1 = temp.gram_schmidt(temp);
    cout << temp1 << endl;

cout << "--------------------------------------------------------------------" << endl;
    Cvector<int> p10;
    p10.array[0] = 1;
    p10.array[1] = 10;
    p10.length = 2;


    Cvector<int> p20;
    p20.array[0] = 5;
    p20.array[1] = 10;
    p20.length = 2;

    cout << p10 -p20 << endl;


    return 0;
}
