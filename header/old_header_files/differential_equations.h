#include <iostream>
#include <fstream>
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <complex>
#include<string>
#include<sys/types.h>
#include<sys/stat.h>

#include "../../header/matrices_and_const.h"





//################ Differential equations for semi classical approach


// input are component of spin vector 0 = x, 1 = y ...; index in vector array
// Bfield in matrices and const
Eigen::Vector3cd EoM_Zeeman_component(uint index, vector<Eigen::Vector3cd> &spin_vec){

    double B_field_strength = 1.;
    Eigen::Vector3cd B_field = B_field_direction.normalized() * B_field_strength;

    Eigen::Vector3cd dS;
    for(uint component = 0; component < 3; component ++){
        if(component == 0){
            dS(component) = B_field(1) * spin_vec[index](2) - B_field(2) * spin_vec[index](1);
        }
        else if(component == 1){
            dS(component) = B_field(2) * spin_vec[index](0) - B_field(0) * spin_vec[index](2);
        }
        else if(component == 2){
            dS(component) = B_field(0) * spin_vec[index](1) - B_field(1) * spin_vec[index](0);
        }
    }
    //dS = spin_vec[index].cross(B_field);
    return dS;
}

void EoM_Zeeman(vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, double t){

    for( uint i = 0; i < spin_vec.size(); i++){
        dxdt[i] = EoM_Zeeman_component(i, spin_vec);

    }
}



