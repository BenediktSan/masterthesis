#include <iostream>
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <complex>
#include<string>
#include<sys/types.h>
#include<sys/stat.h>

#include "../header/functions.h"



using cplx = std::complex<double>;
using namespace std;
using Observable = Eigen::MatrixXcd;
using uint = unsigned int;
using std::vector;



//klasse um blockmatrizen zu diagonalisieren

class Block_Matrix{
private:

    vector<Observable> block_vector{};
    vector<uint> size_of_blocks{};
    
    uint matrix_size = 0;


public:

    Block_Matrix(vector<Observable> block_vector_init){
        block_vector = block_vector_init;

        for( int i= 0; i < block_vector.size(); i++ ){
            size_of_blocks.emplace_back(block_vector[i].rows());
            matrix_size = matrix_size + block_vector[i].rows();
        }

    }

    Observable diagonalization_matrix(){
        Observable unitary_transformation = Observable::Zero(matrix_size, matrix_size);

        int place_of_next_matrix =0;

        for(int i = 0; i < block_vector.size(); i++){

            Eigen::ComplexEigenSolver<Observable> Solver(block_vector[i]);
            Observable eigenvectors = Solver.eigenvectors();
            unitary_transformation.block(place_of_next_matrix, place_of_next_matrix, size_of_blocks[i], size_of_blocks[i]) = eigenvectors;
            place_of_next_matrix = place_of_next_matrix + size_of_blocks[i];

        }
        
        return unitary_transformation;
    }






};

