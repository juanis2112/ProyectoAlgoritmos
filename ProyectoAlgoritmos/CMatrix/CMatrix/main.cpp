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
    Cvector<double> v;
    v.push(1);
    v.push(4);
    v.push(3);
    Cmatrix<double> r (3,v);
    r.array[1][1] = -1;
    r.array[0][1] = -3;
    r.array[2][2] = -2;
    cout<<r<<endl;
    r.LUP(r,0.001);
    cout<<r<<endl;
    
    

//    Cvector<Cmatrix<double>> result;
//    result.array[0] = r;
//    result.array[1] = r;
//    result.array[2] = r;
//    result.length = 3;
    
    
    
    //cout<<r;
	return 0;
}
