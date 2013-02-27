/*
 *  TestSimplexIterator.cpp
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 2/25/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <iostream>

#include "SimplexIterator.h"

using namespace std;

int main(int argc, char *argv[])
{
    SimplexIterator<3>  s3(12);
    
    for (int i=0; i < 20; i++) {
        cout << s3.index[0] << " " << s3.index[1] << " " << s3.index[2] << endl;
        s3++;
    }
    
    cout << ">> " << s3.index[0] << " " << s3.index[1] << " " << s3.index[2] << endl;
    s3 = 19;
    for (int i=0; i < s3.end()-18; i++) {
        cout << i << " : " << s3.index[0] << " " << s3.index[1] << " " << s3.index[2] << endl;
        s3++;
    }
}