//
//  main.cpp
//  CMatrix
//
//  Created by Juanita Gomez on 4/30/19.
//  Copyright © 2019 Juanita Gomez. All rights reserved.
//

#include <iostream>
#include "CMatrix.hpp"

using namespace std;

int main() {
//    Cvector<int> tmp;
//    for (int i = 1; i < 5; i++) tmp.push(i*10);
//    Cmatrix<int> x(3,tmp);
//    cout << x << endl;
//
//    Cvector<int> tmp2;
//    for (int i = 1; i < 4; i++) tmp2.push(i);
//    Cmatrix<int> y(4,tmp2);
//    cout << y << endl;
//    Cmatrix<int> z;


    Cvector<double> v;
    v.push(1);
    v.push(4);
    v.push(3);
    Cmatrix<double> r (3,v);
    r.array[1][1] = -1;
    r.array[0][1] = -3;
    r.array[2][2] = -2;
    cout<<r<<endl;
    tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> c = r.LUP(0.001);
    Cmatrix<double> L;
    L = get<1>(c);
    Cmatrix<double> U;
    U = get<2>(c);
    cout<<"P:"<<get<0>(c)<<endl;
    cout<<"L:"<<get<1>(c)<<endl;
    cout<<"U:"<<get<2>(c)<<endl;
    
    Cmatrix<double> rinv;
    rinv = r.inverse();
    cout<<"This is the inverse of A"<<rinv;
    
    double det = r.determinant();
    cout<<det;

//    Cmatrix<double> m;
     Cmatrix<double> n;
//    Cmatrix<double> o;
//    Cmatrix<double> ra;
//    cout<<"Esto es un vector"<<r[1]<<endl;
//    Cmatrix<double> m = r.transpose();
//    cout<<"tranpose:"<<m<<endl;
//    Cmatrix<double> u = r.lowerTriangular();
//    cout<<"lower"<<u<<endl;
//    Cmatrix<double> s = r.upperTriangular();
//    cout<<"lower"<<s<<endl;
//
//    m = Cmatrix<double>::eye(8);
//    cout<<m<<endl;
//
     n = Cmatrix<double>::zeros(8,3);
//    cout<<n<<endl;
//
//    o = Cmatrix<double>::ones(8,3);
//    cout<<o<<endl;
//
//    ra = Cmatrix<double>::random(8,3);
//    cout<<ra<<endl;

    cout << endl << endl;
    // cout << r << endl;

    cout << "//////////////////////////////////////////////////" << endl;
    cout << endl;

    Cmatrix<double> nr = r.transpose();
    // cout << nr << endl;

    Cmatrix<double> a = product(r, nr);
    cout << a << endl;

    //////////////////////////////////////////////////////////////

    cout << "//////////////////////////////////////////////////" << endl;
    cout << endl;

    Cmatrix<double> t;
    t.push(v);
    // cout << t << endl;

    Cmatrix<double> b = product(r, t);
    cout << b << endl;

    cout << "///////////////////////////////////////////////////" << endl;
    cout << endl;

    Cmatrix<double> d = product(t, r);
    cout << d << endl;

    Cmatrix<double> z = d.toDouble();
    cout << z << endl;

	  return 0;
}
