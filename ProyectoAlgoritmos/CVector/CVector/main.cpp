// main.cpp

// Santiago Lopez
// Oscar Velasco
// Juanita Gomez

#include "CVector.hpp"
using namespace std;
int main(){
//    Cvector<double> p(3,100.0);
//    Cvector<double> q(p);
//    //for (int i = 0; i < 10 ; i++) q.insert(0,i);
//    Cvector<double> r;
//    r = p ^ 3.0;
//    cout << p << endl;
//      cout << q << endl;
//    cout << r << endl;
//    Cvector<double> p1;
//    Cvector<double> p2;
//    for (int i = 0; i < 5; i++){
//        p1.insert(0,double(i));
//        p2.insert(0,double(i+10));
//    }
//
//
//    cout << p1.angle(p2) << endl;
//    cout << p1.proy(p2) << endl;
//    cout << "p1: " << p1 << endl;
//    cout << "p2 :" << p2 << endl;
//    p2.erase(2);
//    cout << "p2 :" << p2 << endl;
//
//
//    Cvector<int> h;
//    h.push(3);
//    h.push(5);
//    h.push(1);
//    Cvector<Cvector<int>> temp;
//    temp.push(h);
//    h.clear();
//    h.push(2);
//    h.push(1);
//    h.push(0);
//    temp.push(h);
//    h.clear();
//    h.push(0);
//    h.push(1);
//    h.push(0);;
//    temp.push(h);
//    cout << temp << endl;
//    Cvector<Cvector<double>> temp1;
//    temp1 = temp.gram_schmidt(temp);
//    cout << temp1 << endl;
//
//cout << "--------------------------------------------------------------------" << endl;
//    Cvector<int> p10;
//    p10.array[0] = 1;
//    p10.array[1] = 10;
//    p10.length = 2;
//
//
//    Cvector<int> p20;
//    p20.array[0] = 5;
//    p20.array[1] = 10;
//    p20.length = 2;
//
//    cout << p10 -p20 << endl;
//
//    Cvector<int> t;
//    t.push(10);
//    t.push(6);
//    t.push(4);
//    t.push(11);
//    cout << t << endl;
//    cout << t.toDouble() << endl;
//
//    Cvector<float> t2;
//    t2.push(15.4);
//    t2.push(6.8);
//    t2.push(10.5);
//    cout << t2 << endl;
//    cout << t2.toDouble() << endl;

    
    
    
    Cvector<int> p;
    p.push(1);
    p.push(5);
    p.push(3);
    cout << "p : " << p << endl;
    Cvector<int> q;
    q.push(7);
    q.push(0);
    q.push(2);
    cout << "q : " << q << endl;
    Cvector<int> m;
    m.push(4);
    m.push(1);
    m.push(9);
    cout << "m : " << m << endl;
    
//    Cvector<Cvector<int>> tmp;
//    tmp.push(p);
//    tmp.push(q);
//    tmp.push(m);
//    cout << "tmp : " << tmp << endl;
//    cout << "---------------------------------------------------------------" << endl;
//    cout << "Gram_schmidt: " << tmp.gram_schmidt(tmp) << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "Angle p_q : " << q.angle(p) << endl;
    return 0;
}
