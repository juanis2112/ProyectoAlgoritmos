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
    Cvector<int> tmp;
    for (int i = 1; i < 5; i++) tmp.push(i*10);
    Cmatrix<int> r(3,tmp);
    cout << r << endl;
    Cmatrix<int> t(3,4,false);
    cout << t << endl;
    
    
	//cout << r.access(7,5) << endl; 
	return 0;
}
