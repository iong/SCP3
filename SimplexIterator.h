/** @file SimplexIterator.h
 *  @author igeorges@uci.edu
 *  @date 2/25/13
 *  @brief Iterator for an N-dimensional simplex
 *
 *  ~~~~~~~~~~~~~~~~~~~~{.c}
 *  for (i=0; i<M; i++) {
 *      for (j=i+1; j<M; j++) {
 *          for (k=j+1; k<M; k++) {
 *              ...
 *  ~~~~~~~~~~~~~~~~~~~~
 */


#ifndef SCP_Double_Excitations_SimplexIterator_h
#define SCP_Double_Excitations_SimplexIterator_h

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

template <int dim>
class SimplexIterator
{
    vector<int> V[dim];
    int N, last_location;

public:
    int index[dim];

    SimplexIterator(int N_);
    SimplexIterator<dim>& operator=(int location);
    SimplexIterator operator++(int);
    SimplexIterator operator--(int);

    const int& end() {
        return V[dim-1][N-1];
    }
};

template<int dim>
ostream& operator<<(ostream& os, SimplexIterator<dim>& it)
{
    os << "(";
    for (int i=0; i<dim-1; i++) {
        os<<it.index[i] << ", ";
    }
    os << it.index[dim-1] << ")";
    
    return os;
}

template<int dim>
SimplexIterator<dim>::SimplexIterator(int N_) : N(N_)
{
    int i, j;
    
    last_location = 0;
    
    for (i=0; i<dim; i++) {
        V[i].resize(N);
        V[i][i] = 1;
    }
    
    for (i=0; i<N; i++) {
        V[0][i] = i+1;
    }
    
    for (i=1; i<dim; i++) {
        for (j=i+1; j<N; j++) {
            V[i][j] = V[i][j-1] + V[i-1][j-1];
        }
    }
    
    for (i=0; i<dim; i++) index[i] = i;
}


template<int dim>
SimplexIterator<dim>& SimplexIterator<dim>::operator=(int location)
{
    int d, j;
    if (location == last_location) {
        return *this;
    }
    else if (location == last_location + 1) {
        (*this)++;
        
        return *this;
    }
    else if (location == last_location - 1) {
        (*this)--;
        
        return *this;
    }
 
    
    last_location = location;
    
    fill(index, index+dim, 0);
    
    for (j=0, d=dim-1; d>0; d--) {
        for (; j <= N - d && location >= V[d-1][(N - 1) - (j + 1)]; j++) {
            location -= V[d-1][ (N - 1) - (j + 1) ];
        }
        index[dim - 1 - d] = j++;
    }
    index[dim - 1] = j + location;
    
    return *this;
}


template<int dim>
SimplexIterator<dim> SimplexIterator<dim>::operator++(int)
{
    int d;
    
    index[dim-1]++;
    for (d=dim-1; d>0 && index[d]>=N + d - (dim-1); d--) {
        index[d] = 0;
        index[d-1]++;
    }
    if (index[0] > N-dim) {
        cerr << "Simplex incremented out of bounds.\n";
        exit(EXIT_FAILURE);
    }
    for (d=d+1; d<dim; d++) {
        index[d] = index[d-1]+1;
    }
    
    last_location++;
    
    return *this;
}


template<int dim>
SimplexIterator<dim> SimplexIterator<dim>::operator--(int)
{
    int d;
    
    index[dim-1]--;
    for (d=dim-1; d>0 && index[d]<=index[d-1]; d--) {
        index[d] = N + d - (dim - 1) - 1;
        index[d-1]--;
    }
    if (index[0] < 0) {
        cerr << "Simplex incremented out of bounds.\n";
        exit(EXIT_FAILURE);
    }

    last_location--;
    
    return *this;
}



#endif
