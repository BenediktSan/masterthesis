#include <iostream>
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <complex>
#include<string>
#include<sys/types.h>
#include<sys/stat.h>



using cplx = std::complex<double>;
using namespace std;
using Eigen::MatrixXd;
using Observable = Eigen::MatrixXcd;
using uint = unsigned int;


// Forward Declarations
Observable compute_multi_tensor_product( std::vector<Observable> local_operators );
Observable compute_tensor_product( const Observable& A, const Observable& B );

// function : compute a tensor product of N operators
Observable compute_multi_tensor_product( std::vector<Observable> local_operators )
{
  if( local_operators.size() == 1 ) // only one operator given -> return it
  {
    return local_operators[0];
  }
  else // vector of operators given -> compute tensor product and return it
  {
    local_operators[0] = compute_tensor_product( local_operators[0], local_operators[1] ); // replace first operator by the tensor product
    local_operators.erase( local_operators.begin()+1 ); // erase second operator
    return compute_multi_tensor_product( local_operators );
  }
}

// function : compute a tensor product of 2 operators
Observable compute_tensor_product( const Observable& A, const Observable& B )
{
  uint dimA = A.rows();
  uint dimB = B.rows();
  Observable C;
  C.resize( dimA*dimB, dimA*dimB );
  for( uint rowA=0; rowA<dimA; rowA++ ){
    for( uint colA=0; colA<dimA; colA++ ){
      for( uint rowB=0; rowB<dimB; rowB++ ){
        for( uint colB=0; colB<dimB; colB++ ){
          C( dimB*rowA + rowB, dimB*colA + colB ) = A( rowA, colA ) * B( rowB, colB ); 
        }
      }
    }
  }
  return C;
}

