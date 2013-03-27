/*
 *  TestSimplexIterator.cpp
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 2/25/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <cstdlib>
#include <iostream>

#include "SimplexIterator.h"

using namespace std;

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    int M;
    
    SimplexIterator<2>  s2(N);
    SimplexIterator<3>  s3(N);
    
    M = s3.end();
    for (int i=M; i>0; i--) {
	    s3 = i-1;
	    cout << s3 << endl;
    }
    
}
