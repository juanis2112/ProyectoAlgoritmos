//  main.cpp
//  CMatrix
//
//  Created by Juanita Gomez on 4/30/19.
//  Copyright © 2019 Juanita Gomez. All rights reserved.
//

#include <iostream>
#include<vector>
#include "CMatrix.hpp"

using namespace std;

int main() {
    //    // Prueba de los métodos de la clase ///////////////////////////////////////
    //
    //    cout.precision(5);
    //    cout << endl;
    //    cout << "-------------------------------------------" << endl;
    //    cout << endl;
    //
    ////    Cmatrix<int> p;
    ////    cout << "P : " << p << endl;
    ////    cout << "p.numberRows : " << p.numberRows() << endl;
    ////    cout << "p.numberCols: " << p.numberCols() << endl;
    ////    cout << "-------------------------------------------" << endl;
    ////    Cvector<int> q(100,8);
    ////    cout << "Estoy metiendo : " <<  q << endl;
    ////    p.push(q);
    ////    Cvector<int> r(100,3);
    ////    cout << "Estoy metiendo : " <<  r << endl;
    ////    p.push(r);
    ////    Cvector<int> t(100,5);
    ////    cout << "Estoy metiendo : " <<  t << endl;
    ////    p.push(t);
    ////    cout << "P Despues de hacer push : " << p << endl;
    ////    cout << "p.numberRows : " << p.numberRows() << endl;
    ////    cout << "p.numberCols: " << p.numberCols() << endl;
    ////
    ////    Cmatrix<int> k;
    ////    p.push(10);
    ////    Cmatrix<int> y = p+k;
    ////    cout << "P + K" << y << endl;
    //    Cmatrix<int> p;
    //    Cvector<int> r;
    //    r.push(9);
    //    r.push(5);
    //    r.push(7);
    //    cout << r << endl;
    //    p.push(r);
    //    p.push(r);
    //    cout << p << endl;
    //    cout << "p.numberRows : " << p.numberRows() << endl;
    //    cout << "p.numberCols : " << p.numberCols() << endl;
    //    p.swap_c(0,2);
    //    cout << p << endl;
    //    cout << "p.numberRows : " << p.numberRows() << endl;
    //    cout << "p.numberCols : " << p.numberCols() << endl;
    ////
    ////
    ////
    ////
    //
    //
    //
    //// Empty
    //// Push
    //// Transpose
    //// Ostream
    //// numberRows
    //// numberCols
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    
    
    
        // 1) Halle la matriz traspuesta de: m1 = {(5, 3, 6), (4, 10, 7), (1, 4, 2)}
//        vector<int> a1 = {5, 3, 6};
//        vector<int> b1 = {4, 10, 7};
//        vector<int> c1 = {1, 4, 2};
//        Cmatrix<int> m1;
//        Cvector<int> d1;
//        Cvector<int> e1;
//        Cvector<int> f1;
//        for(unsigned i = 0; i < a1.size(); i++){
//            d1.push(a1[i]);
//            e1.push(b1[i]);
//            f1.push(c1[i]);
//        }
//        m1.push(d1);
//        m1.push(e1);
//        m1.push(f1);
//
//        cout << "m1: " << endl;
//        cout << m1 << endl;
//        cout << "m1.numberRows : " << m1.numberRows() << endl;
//        cout << "m1.numberCols : " << m1.numberCols() << endl;
//        cout << endl;
//        cout << "Transpose of m1: " << endl;
//        Cmatrix<int> m1t = m1.transpose();
//        cout << m1t << endl;
//        cout << "m1t.numberRows : " << m1t.numberRows() << endl;
//        cout << "m1t.numberCols : " << m1t.numberCols() << endl;
//
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
//
//        // 2) Transforme la matriz m1 a triangular inferior
//        cout << "m1 to lower triangular: " << endl;
//        Cmatrix<double> m1lt = m1.lowerTriangular();
//        cout << m1lt << endl;
//        cout << "m1lt.numberRows : " << m1lt.numberRows() << endl;
//        cout << "m1lt.numberCols : " << m1lt.numberCols() << endl;
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
//
//        // 3) Transforme la matriz m1 a triangular superior
//        cout << "m1 to upper triangular: " << endl;
//
//        Cmatrix<double> m1ut = m1.upperTriangular();
//        cout << m1ut << endl;
//        cout << "m1ut.numberRows : " << m1ut.numberRows() << endl;
//        cout << "m1ut.numberCols : " << m1ut.numberCols() << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // 4) Cree la matriz identidad de orden 10
//        Cmatrix<double> m2;
//        cout << "m2: " << endl;
//        cout << m2 << endl;
//        cout << "m2.numberRows : " << m2.numberRows() << endl;
//        cout << "m2.numberCols : " << m2.numberCols() << endl;
//        cout << endl;
//        cout << "m2 to Identity: " << endl;
//
//        Cmatrix<double> m2e = m2.eye(10);
//        cout << m2e << endl;
//        cout << "m2e.numberRows : " << m2e.numberRows() << endl;
//        cout << "m2e.numberCols : " << m2e.numberCols() << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
//
//
//        // 5) Cree la matriz nula con 4 filas y 5 columnas
//
//        cout << "m2: " << endl;
//        cout << m2 << endl;
//        cout << "m2 to Zeros: " << endl;
//
//        Cmatrix<double> m2z = m2.zeros(4, 5);
//        cout << m2z << endl;
//        cout << "m2z.numberRows : " << m2z.numberRows() << endl;
//        cout << "m2z.numberCols : " << m2z.numberCols() << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // 6) Cree una matriz con 5 filas y 4 columnas, cuyos elementos sólo sean 1s
//        cout << "m2: " << endl;
//        cout << m2 << endl;
//        cout << "m2 to Ones: " << endl;
//
//        Cmatrix<double> m2o = m2.ones(5, 4);
//        cout << m2o << endl;
//        cout << "m2o.numberRows : " << m2o.numberRows() << endl;
//        cout << "m2o.numberCols : " << m2o.numberCols() << endl;
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
//
//        // 7) Cree una matriz con 3 filas y 5 columnas, cuyos valores sean random
//        cout << "m2: " << endl;
//        cout << m2 << endl;
//        cout << "m2 to random: " << endl;
//
//        Cmatrix<double> m2r = m2.random(3, 5);
//        cout << m2r << endl;
//        cout << "m2r.numberRows : " << m2r.numberRows() << endl;
//        cout << "m2r.numberCols : " << m2r.numberCols() << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // 8) Haga los siguientes swaps en la matriz m3 = {(9, 10, 5, 30), (5, 2, 6, 4), (7, 45, 1, 8)}
//        // a) fila 1 - fila 3
//        // b) columna 2 - columna 3
//        vector<int> a2 = {9, -10, 5, -30};
//        vector<int> b2 = {5, 2, -6, 4};
//        vector<int> c2 = {-7, 45, 1, -8};
//        Cmatrix<int> m3;
//        Cvector<int> d2;
//        Cvector<int> e2;
//        Cvector<int> f2;
//        for(unsigned i = 0; i < a2.size(); i++){
//            d2.push(a2[i]);
//            e2.push(b2[i]);
//            f2.push(c2[i]);
//        }
//        m3.push(d2);
//        m3.push(e2);
//        m3.push(f2);
//
//        cout << "m3: " << endl;
//        cout << m3 << endl;
//        cout << endl;
//
//        // a)
//        m3.swap_r(0, 2);
//        cout << "Changed Row 0 with Row 2 in m3: " << endl;
//        cout << m3 << endl;
//
//        // b)
//        m3.swap_c(1, 2);
//        cout << "Changed Column 1 with Column 3 in m3: " << endl;
//        cout << m3 << endl;
//
//        cout << endl;
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // 9) Apliquele valor absoluto a todos los elementos de la matriz m3
//        cout << "m3: " << endl;
//        cout << m3 << endl;
//        cout << endl;
//        Cmatrix<int> m3abs = m3.abs();
//        cout << "m3 to abs: " << endl;
//        cout << m3abs << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // 10) Halle el determinante de la matriz r
//        Cvector<double> v;
//        v.push(1);
//        v.push(4);
//        v.push(3);
//        Cmatrix<double> r (3,v);
//        r(1,1) = -1;
//        r(0,1) = -3;
//        r(2,2) = -2;
//        cout << "r: " << endl;
//        cout << r << endl;
//        tuple<Cvector<double>, Cmatrix<double>, Cmatrix<double>> c = r.LUP(0.001);
//        // cout << "P:" << get<0>(c) << endl;
//        // cout << "L:" << get<1>(c) << endl;
//        // cout << "U:" << get<2>(c) << endl;
//        double det = r.determinant();
//        cout << "The determinant of r is: " << det << endl;
//
//        cout << endl;
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // 11) Imprima:
//        // a) Numero de filas de m3
//        // b) Numero de columnas de m3
//        // c) El elemento en la fila 2, columna 3
//
//        cout << "m3: " << endl;
//        cout << m3 << endl;
//        cout << endl;
//
//        // a)
//        cout << "The Rows number of m3 is: " << m3.numberRows() << endl;
//        cout << endl;
//
//        // b)
//        cout << "The Columns number of m3 is: " << m3.numberCols() << endl;
//        cout << endl;
//
//        // c)
//        cout << "The element in row 2 and column 3 is: " << m3.access(2, 3) << endl;
//        cout << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
////        // Prueba de los operadores ////////////////////////////////////////////////
//
//    vector<int> a3 = {9, -10, 5};
//    vector<int> b3 = {5, 2, -6};
//    vector<int> c3 = {-7, 45, 1};
//    Cmatrix<int> m4;
//    Cvector<int> d3;
//    Cvector<int> e3;
//    Cvector<int> f3;
//    for(unsigned i = 0; i < a3.size(); i++){
//        d3.push(a3[i]);
//        e3.push(b3[i]);
//        f3.push(c3[i]);
//    }
//    cout << "Lo que voy a meter con push a m4 : "<< d3 << endl;
//    m4.push(d3);
//    cout << "Lo que voy a meter con push a m4 : "<< e3 << endl;
//    m4.push(e3);
//    cout <<"Lo que voy a meter con push a m4 : "<< f3 << endl;
//    m4.push(f3);
//
//    cout << "m4 : " << m4 << endl;
//    cout << "m4.numberRows : " << m4.numberRows() << endl;
//    cout << "m4.numberCols : " << m4.numberCols() << endl;
//
//    cout << "----------------------------------------------" << endl;
//    cout << endl;
////
//    vector<int> a4 = {9, -10, 5};
//    vector<int> b4 = {5, 2, -6};
//    vector<int> c4 = {-7, 45, 1};
//    Cmatrix<int> m5;
//    Cvector<int> d4;
//    Cvector<int> e4;
//    Cvector<int> f4;
//    for(unsigned i = 0; i < a4.size(); i++){
//        d4.push(a4[i]);
//        e4.push(b4[i]);
//        f4.push(c4[i]);
//    }
//
//    cout << "Lo que voy a meter con push a m5 : "<< d4 << endl;
//    m5.push(d4);
//    cout << "Lo que voy a meter con push a m5 : "<< e4 << endl;
//    m5.push(e4);
//    cout << "Lo que voy a meter con push a m5 : "<< f4 << endl;
//    m5.push(f4);
//
//    cout << "m5 : " << m5 << endl;
//    cout << "m5.numberRows : " << m5.numberRows() << endl;
//    cout << "m5.numberCols : " << m5.numberCols() << endl;
//
//    cout << "----------------------------------------------" << endl;
//    cout << endl;
//
////    // Operator == and Operator !=
//    cout << "m4 == m5: " << endl;
//    Cmatrix<bool> m4igum5 = (m4 == m5);
//    cout << m4igum5 << endl;
//
//    cout << "m4 != m5: " << endl;
//    Cmatrix<bool> m4dism5 = (m4 != m5);
//    cout << m4dism5 << endl;
////
//    m4.swap_r(0, 2);
//    m5.swap_c(0, 2);
//    cout << "----------------------------------------------" << endl;
//    cout << "DESPUES DE HACER SWAP" << endl;
//    cout << endl;
//
//    cout << "m4 : " << m4 << endl;
//    cout << "m5 : " << m5 << endl;
//
//    cout << "-------------------------------------------" << endl;
//    cout << endl;
//
////    // Operator <= and Operator >=
////
//    cout << "m4 <= m5: " << endl;
//    Cmatrix<bool> m4meim5 = (m4 <= m5);
//    cout << m4meim5 << endl;
//
//    cout << "m4 >= m5: " << endl;
//    Cmatrix<bool> m4maim5 = (m4 >= m5);
//    cout << m4maim5 << endl;
//
//    cout << "DESPUES DE HACER SWAP A M5" << endl;
//    m5.swap_c(0, 2);
//    cout << "m5 : " << m5 << endl;
//    cout << "-------------------------------------------" << endl;
//    cout << endl;
//
//    // Operator < and Operator >
//
//    cout << "m4: " << endl;
//    cout << m4 << endl;
//    cout << "m4.numberRows: " << m4.numberRows() << endl;
//    cout << "m4.numberCols: " << m4.numberCols() << endl;
//    cout << endl;
//
//    cout << "m5: " << endl;
//    cout << m5 << endl;
//    cout << "m5.numberRows: " << m5.numberRows() << endl;
//    cout << "m5.numberCols: " << m5.numberCols() << endl;
//    cout << endl;
//
//    cout << "m4 < m5: " << endl;
//    Cmatrix<bool> m4mem5 = (m4 < m5);
//    cout << m4meim5 << endl;
//
//    cout << "m4 > m5: " << endl;
//    Cmatrix<bool> m4mam5 = (m4 > m5);
//    cout << m4mam5 << endl;
//
//    cout << "-------------------------------------------" << endl;
//    cout << endl;
//////
////    // Operator + and Operator -
//    cout << "m4: " << endl;
//    cout << m4 << endl;
//    cout << "m4.numberRows: " << m4.numberRows() << endl;
//    cout << "m4.numberCols: " << m4.numberCols() << endl;
//    cout << endl;
//
//    cout << "m5: " << endl;
//    cout << m5 << endl;
//    cout << "m5.numberRows: " << m5.numberRows() << endl;
//    cout << "m5.numberCols: " << m5.numberCols() << endl;
//    cout << endl;
//
//    Cmatrix<double> m4masm5 = m4 + m5;
//    cout << "m4 + m5: " << endl;
//    cout << m4masm5 << endl;
//
////    //
//    Cmatrix<double> m4menm5 = m4 - m5;
//    cout << "m4 - m5: " << endl;
//    cout << m4menm5 << endl;
//
//    cout << "-------------------------------------------" << endl;
//    cout << endl;
//
////    //
////        // Operator * scalar and Operator / scalar
////
//         cout << "m4: " << endl;
//         cout << m4 << endl;
//         cout << "m4 * 3: " << endl;
//         Cmatrix<double> m4x3 = m4 * 3;
//         cout << m4x3 << endl;
//
//         cout << "m4 / 2: " << endl;
//         Cmatrix<double> m4d2 = m4 / 2;
//         cout << m4d2 << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//////        // Operator Cmatrix ^ n
//////
//        cout << "m4: " << endl;
//        cout << m4 << endl;
//        Cmatrix<double> m4elev4 = m4 ^ size_t(4);
//        cout << "m4 ^ 4: " << endl;
//        cout <<  m4elev4 << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
////
//        // Operator Cmatrix * CVector, Operator Cvector * Cmatrix and Operator Cmatrix
//
//        Cmatrix<double> o = r;
//        cout << "o: " << endl;
//        cout << o << endl;
//        cout << "r: " << endl;
//        cout << r << endl;
//        Cvector<double> to;
//        to.push(5);
//        to.push(2);
//        to.push(3);
//        cout << "to: " << to << endl;
//        cout << endl;
//
//        Cmatrix<double> oXto = o * to;
//        cout << "o * to: " << endl;
//        cout << oXto << endl;
//
//        Cmatrix<double> toXo = to * o;
//        cout << "to * o: " << endl;
//        cout << toXo << endl;
//
//        Cmatrix<double> rXo = r * o;
//        cout << "r * o: " << endl;
//        cout << rXo << endl;
//
//        cout << "-------------------------------------------" << endl;
//        cout << endl;
//
//        // Cmatrix<double> rele2 = r ^ 2;
//        // cout << rele2 << endl;
//
//         //Prueba de metodos basicos de Cmatrix
//
//        vector<int> a5 = {9, -10, 5, -30, 6};
//        vector<int> b5 = {5, 2, -6, 4, 0};
//        vector<int> c5 = {-7, 45, 1, -8, 4};
//        vector<int> d5 = {5, 12, 3, 19, 22};
//        Cmatrix<int> m6;
//        Cvector<int> e5;
//        Cvector<int> f5;
//        Cvector<int> g5;
//        Cvector<int> h5;
//        for(unsigned i = 0; i < a5.size(); i++){
//            e5.push(a5[i]);
//            f5.push(b5[i]);
//            g5.push(c5[i]);
//            h5.push(d5[i]);
//        }
//        m6.push(e5);
//        m6.push(f5);
//        m6.push(g5);
//        m6.push(h5);
//
//        cout << "m6: " << endl;
//        cout << m6 << endl;
//        cout << endl;
//
//         // erase()
//         cout << "m6.erase(3)" << endl;
//         m6.erase(3);
//         cout << "m6: " << endl;
//         cout << m6 << endl;
//         cout << endl;
//
//        // insert()
//        cout << "m6.insert(3, [ 9 -10 5 -30 6 ])" << endl;
//        m6.insert(3, e5);
//        cout << "m6: " << endl;
//        cout << m6 << endl;
//        cout << endl;
//
//        // clear()
//        cout << "m6.clear()" << endl;
//        m6.clear();
//        cout << "m6: " << endl;
//        cout << m6 << endl;
//        cout << endl;
//
//        // empty()
//        cout << boolalpha << "m6 is empty?: " << m6.empty() << endl;
//        cout << endl;
//
//        m6.push(h5);
//
//        // empty()
//        cout << boolalpha << "m6 is empty?: " << m6.empty() << endl;
//
//        cout << endl;
//        cout << "-------------------------------------------" << endl;
////
//
//
////    Cmatrix<double> santi = Cmatrix<double>::ones(3,3);
////    Cmatrix<double> juanis = Cmatrix<double>::ones(3,3);
////    Cvector<double> imbecil;
////    imbecil.push(2);
////    imbecil.push(1);
////    imbecil.push(-1);
////    santi(0,2) = 9;
////    juanis(1,2) = 2;
////    Cmatrix<double> mult = santi*juanis;
////    Cvector<double> mult1 = santi*imbecil;
////    Cvector<double> mult2 = imbecil*juanis;
////    Cmatrix<double> mult3 = juanis^2;
////    cout<<santi<<"X"<<juanis<<"="<<mult;
////    cout<<santi<<"X"<<imbecil<<"="<<mult1;
////    cout<<imbecil<<"X"<<juanis<<"="<<mult2;
////    cout<<juanis<<"X"<<juanis<<"="<<mult3;
    
//    Cvector<int> f1;
//    f1.push(3);
//    f1.push(3);
//    f1.push(4);
//    Cvector<int> f2;
//    f2.push(4);
//    f2.push(4);
//    f2.push(5);
//    Cvector<int> f3;
//    f3.push(3);
//    f3.push(3);
//    f3.push(4);
//    Cmatrix<int> temp;
//    temp.push(f1);
//    temp.push(f2);
//    temp.push(f3);
//    cout << "temp : " << temp << endl;
//
//    Cvector<int> f0(3,1);
//    temp.insert(1,f0,true);
//    cout << "temp_Insertado : " << temp << endl;
    
    Cvector<int> p;
    p.push(1);
    p.push(2);
    p.push(3);
    Cvector<int> q;
    q.push(4);
    q.push(5);
    q.push(6);
    Cmatrix<int> temp;
    temp.push(p,1);
    cout << temp << endl;
    return 0;
    
    
}

