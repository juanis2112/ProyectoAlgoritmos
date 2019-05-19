# include<iostream>
# include<vector>
using namespace std;

# include "CVector.hpp"

int main(){

    cout.precision(5);
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;
    cout<< "1) Halle el producto punto entre los siguientes vectores: " <<endl;
    cout<< "  a) v1 = (2, -9, 5, 7), v2 = (5, 9, -4, 6)," <<endl;  // Res = 329
    cout<< "  b) v3 = (15, 7, 12, 10, 3), v4 = (19, 24, 1, -13, 8)"<<endl; // Res = 359
    cout <<endl;
    cout<< " a) "<<endl;
    cout << endl;
    vector<int> a1 = {2, -9, 5, 7};
    vector<int> b1 = {5, 9, -4, 60};
    Cvector<int> v1;
    Cvector<int> v2;
    for(unsigned i = 0; i < a1.size(); i++){
        v1.push(a1[i]);
        v2.push(b1[i]);
    }
    cout << "  v1: " << v1 << endl;
    cout << "  v2: " << v2 << endl;
    cout << "  v1 dot v2 = " << v1.dot(v2) << endl;
    cout << "  v2 dot v1 = " << v2.dot(v1) << endl;
    cout << endl;
    cout<< " b) "<<endl;
    vector<int> c1 = {15, 7, 12, 10, 3};
    vector<int> d1 = {19, 24, 1, -13, 8};
    Cvector<int> v3;
    Cvector<int> v4;
    for(unsigned i = 0; i < c1.size(); i++){
        v3.push(c1[i]);
        v4.push(d1[i]);
    }
    cout << endl;
    cout << "  v3: " << v3 << endl;
    cout << "  v4: " << v4 << endl;
    cout << "  v3 dot v4 = " << v3.dot(v4) << endl;
    cout << "  v4 dot v3 = " << v4.dot(v3) << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< "2) Halle el producto cruz entre los siguientes vectores: "<<endl;
    cout<< " a) v5 = (2, 4, -5), v6 = (-3, -2, 1)"<<endl; // Res = (-6, 13, 8)
    cout<< " b) v7 = (1, -3, 4), v8 = (-2, 1, 1)"<<endl; // Res = (-7, -9, -5)
    cout <<endl;
    cout<< " a) "<<endl;
    cout << endl;
    vector<int> a2 = {2, 4, -5};
    vector<int> b2 = {-3, -2, 1};
    Cvector<int> v5;
    Cvector<int> v6;
    for(unsigned i = 0; i < a2.size(); i++){
        v5.push(a2[i]);
        v6.push(b2[i]);
    }
    cout << "  v5: " << v5 << endl;
    cout << "  v6: " << v6 << endl;
    cout << "  v5 x v6 = " << v5.cross(v6) << endl;
    cout << endl;
    cout<< " b) "<<endl;
    vector<int> c2 = {1, -3, 4};
    vector<int> d2 = {-2, 1, 1};
    Cvector<int> v7;
    Cvector<int> v8;
    for(unsigned i = 0; i < c2.size(); i++){
        v7.push(c2[i]);
        v8.push(d2[i]);
    }
    cout << endl;
    cout << "  v7: " << v7 << endl;
    cout << "  v8: " << v8 << endl;
    cout << "  v7 x v8 = " << v7.cross(v8) << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< "3) Halle la norma de los siguientes vectores:" <<endl;
    cout<< " a) v9 = (-5, -2, 4, 7, 5) " <<endl; // Res = sqrt(119) = 10,9087
    cout<< " b) v10 = (0, 1, 0, 0, 0, 0, 1, 1, 6) " <<endl; // Res = sqrt(39) = 6,2449
    cout <<endl;
    cout<< " a) "<<endl;
    cout << endl;
    vector<int> a3 = {-5, -2, 4, 7, 5};
    Cvector<int> v9;
    for(unsigned i = 0; i < a3.size(); i++){
        v9.push(a3[i]);
    }
    cout << "  v9: " << v9 << endl;
    cout << "  ||v9|| = " << v9.norm() << endl;
    cout << endl;
    cout<< " b) "<<endl;
    vector<int> b3 = {0, 1, 0, 0, 0, 0, 1, 1, 6};
    Cvector<int> v10;
    for(unsigned i = 0; i < b3.size(); i++){
        v10.push(b3[i]);
    }
    cout << endl;
    cout << "  v10: " << v10 << endl;
    cout << "  ||v10|| = " << v10.norm() << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< "4) Normalice los siguientes vectores: " <<endl;
    cout<< " a) v11 = (-4, 3), Res = (-0.8, 0.6) " <<endl;
    cout<< " b) v12 = (8, -8), Res = (0.7071, -0.7071) " <<endl;
    cout << endl;
    cout<< " a) "<<endl;
    cout << endl;
    vector<int> a4 = {-4, 3};
    Cvector<int> v11;
    for(unsigned i = 0; i < a4.size(); i++){
        v11.push(a4[i]);
    }
    cout << "  v11: " << v11 << endl;
    cout << "  Normalize v11: " << v11.normalize() << endl;
    cout << endl;
    cout<< " b) "<<endl;
    vector<int> b4 = {8, -8};
    Cvector<int> v12;
    for(unsigned i = 0; i < b4.size(); i++){
        v12.push(b4[i]);
    }
    cout << endl;
    cout << "  v12: " << v12 << endl;
    cout << "  Normalize of v12: " << v12.normalize() << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< "5) Encuentre el angulo en radianes entre los siguientes vectores:" <<endl;
    cout<< " a) v13 = (2, 4), v14 = (-2, 3), Res = 1,0516 rad" <<endl;
    cout<< " b) v15 = (-1, 3, 4), v16 = (5, -2, 7), Res = 1,183  rad" <<endl;
    cout << endl;
    cout<< " a) "<<endl;
    cout << endl;
    vector<int> a5 = {2, 4};
    vector<int> b5 = {-2, 3};
    Cvector<int> v13;
    Cvector<int> v14;
    for(unsigned i = 0; i < a5.size(); i++){
        v13.push(a5[i]);
        v14.push(b5[i]);
    }
    cout << "  v13: " << v13 << endl;
    cout << "  v14: " << v14 << endl;
    cout << "  Angle between v13 and v14: " << v13.angle(v14) << endl;
    cout << endl;
    cout<< " b) "<<endl;
    vector<int> c5 = {-1, 3, 4};
    vector<int> d5 = {5, -2, 7};
    Cvector<int> v15;
    Cvector<int> v16;
    for(unsigned i = 0; i < c5.size(); i++){
        v15.push(c5[i]);
        v16.push(d5[i]);
    }
    cout << endl;
    cout << "  v15: " << v15 << endl;
    cout << "  v16: " << v16 << endl;
    cout << "  Angle between v15 and v16: " << v15.angle(v16) << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< "6) Encuentre la proyección del segundo vector sobre el primero" <<endl;
    cout<< " a) v17 = (5, 6), v18 = (3, 4), Res = (195/61, 234/61) = (3.1967, 3.8360)" <<endl;
    cout<< " b) v19 = (3, 4), v20 = (5, 6), Res = (117/25, 156/25) = (4.68, 6.24)" <<endl;
    cout << endl;
    cout<< " a) "<<endl;
    cout << endl;
    vector<int> a6 = {5, 6};
    vector<int> b6 = {3, 4};
    Cvector<int> v17;
    Cvector<int> v18;
    for(unsigned i = 0; i < a6.size(); i++){
        v17.push(a6[i]);
        v18.push(b6[i]);
    }
    cout << "  v17: " << v17 << endl;
    cout << "  v18: " << v18 << endl;
    cout << "  Proyection of v18 on v17: " << v17.proj(v18) << endl;
    cout << endl;
    cout<< " b) "<<endl;
    cout << "  Proyection of v17 on v18: " << v18.proj(v17) << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< "7) Ortonormalice los siguientes conjuntos de vectores usando el proceso de Gram-Schmidt " <<endl;
    cout<< " a) B1 = {(1, 0, 1), (0, 0, 1), (-1, 1, 0)}, Res = {(sqrt(2)/2, 0, sqrt(2)/2), (-sqrt(2)/2, 0, sqrt(2)/2), (0, 1, 0)}" <<endl;
    cout<< " b) B2 = {(0, 1, -1), (1, 0, 0), (1, 0, 1)}, Res = {(0, sqrt(2)/2, -sqrt(2)/2), (1, 0, 0), (0, 0, 0)}" <<endl;
    cout << endl;
    cout<< " a) "<<endl;
    cout << endl;
    Cvector<Cvector<int>> B1;
    vector<int> a7 = {1, 0, 1};
    vector<int> b7 = {0, 0, 1};
    vector<int> c7 = {-1, 1, 0};
    Cvector<int> v19;
    Cvector<int> v20;
    Cvector<int> v21;
    for(unsigned i = 0; i < a7.size(); i++){
        v19.push(a7[i]);
        v20.push(b7[i]);
        v21.push(c7[i]);
    }
    B1.push(v19);
    B1.push(v20);
    B1.push(v21);

    cout << "  The set B1 is: " << B1 << endl;
    cout << "  The set B1 ortonormalizated is: " << B1.gram_schmidt() << endl;
    cout << endl;
    cout<< " b) "<<endl;
    Cvector<Cvector<int>> B2;
    vector<int> d7 = {0, 1, -1};
    vector<int> e7 = {1, 0, 0};
    vector<int> f7 = {1, 0, 1};
    Cvector<int> v22;
    Cvector<int> v23;
    Cvector<int> v24;
    for(unsigned i = 0; i < d7.size(); i++){
        v22.push(d7[i]);
        v23.push(e7[i]);
        v24.push(f7[i]);
    }
    B2.push(v22);
    B2.push(v23);
    B2.push(v24);

    cout << endl;
    cout << "  The set B2 is: " << B2 << endl;
    cout << "  The set B2 ortonormalizated is: " << B2.gram_schmidt() << endl;

    cout << endl;
    cout << "========================================================================" << endl;
    cout << endl;
    cout<< "OPERACIONES ENTRE VECTORES"<<endl;
    cout << endl;
    cout<< "Operator == and Operator != "<<endl;
    cout << endl;
    vector<double> aux1 = {0.5, 62, 10.3, 8};
    Cvector<double> v25;
    Cvector<double> v26;
    for(unsigned i = 0; i < aux1.size(); i++){
        v25.push(aux1[i]);
        v26.push(aux1[i]);
    }
    cout << "  v25: " << v25 << endl;
    cout << "  v26: " << v26 << endl;
    cout << "  v25 == v26: " << (v25 == v26) << endl;
    cout << "  v25 != v26: " << (v25 != v26) << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< " Operator <= and Operator >=, Operator < and Operator > "<<endl;
    cout << endl;
    vector<double> aux2 = {1.0, 5.7, 9.6, 46, 8.4, 6.1};
    vector<double> aux3 = {7.6, 5.7, 6.5, 10.9, 25, 6.74};
    Cvector<double> v27;
    Cvector<double> v28;
    for(unsigned i = 0; i < aux2.size(); i++){
        v27.push(aux2[i]);
        v28.push(aux3[i]);
    }
    cout << "  v27: " << v27 << endl;
    cout << "  v28: " << v28 << endl;
    cout << "  v27 <= v28: " << (v27 <= v28) << endl;
    cout << "  v27 >= v28: " << (v27 >= v28) << endl;
    cout << "  v27 < v28: " << (v27 < v28) << endl;
    cout << "  v27 > v28: " << (v27 > v28) << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< " Operator / "<<endl;
    cout << endl;
    vector<double> aux4 = {2, 4, 6, 8, 10, 12};
    Cvector<double> v29;
    for(unsigned i = 0; i < aux4.size(); i++){
        v29.push(aux4[i]);
    }
    cout << "  v29: " << v29 << endl;
    cout << "  v29 / 2: " << v29 / 2.0 << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<< " Operator ^ "<<endl;
    cout << endl;
    vector<int> aux5 = {2, 2, 2, 2, 2};
    Cvector<int> v30;
    for(unsigned i = 0; i < aux5.size(); i++){
        v30.push(aux5[i]);
    }
    cout << "  v30: " << v30 << endl;
    cout << "  v30 ^ 3: " << (v30 ^ 3) << endl;
    cout << "  v30 ^ 4: " << (v30 ^ 4) << endl;
    cout << endl;
    cout << "========================================================================" << endl;
    cout << endl;
    cout<<"MÉTODOS DE LA CLASE CVECTOR"<< endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<<" erase() "<< endl;
    cout << endl;
    vector<int> aux6 = {1, 8, 10, 45, 6, 4};
    Cvector<int> v31;
    for(unsigned i = 0; i < aux6.size(); i++){
        v31.push(aux6[i]);
    }
    cout << "  v31: " << v31 << endl;
    cout << "  v31.erase(3)" << endl;
    v31.erase(3);
    cout << "  v31: " << v31 << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;

    cout<<" insert() "<< endl;
    cout << endl;
    cout << "  v31.insert(3, 7)" << endl;
    v31.insert(3, 7);
    cout << "  v31: " << v31 << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;
    cout<<" clear() "<< endl;
    cout << endl;
    cout << "  v31.clear()" << endl;
    v31.clear();
    cout << "  v31: " << v31 << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;
    cout<<" empty() "<< endl;
    cout << endl;
    cout << boolalpha << "  v31 is empty?: " << v31.empty() << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;
    cout<<" push() "<< endl;
    cout << endl;
    cout << "  v31.push()" << endl;
    for(unsigned i = 0; i <= 10; i+=2){
        v31.push(i);
    }
    cout << "  v31: " << v31 << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl;
    cout<<" empty() "<< endl;
    cout << endl;
    cout << boolalpha << "  v31 is empty?: " << v31.empty() << endl;

    cout << endl;
    cout << "------------------------------------------------------------------------" << endl;

    return 0;


}
