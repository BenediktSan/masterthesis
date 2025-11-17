#include <iostream>
#include <fstream>
#pragma once
#define _USE_MATH_DEFINES
#include <boost/numeric/odeint.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <cmath>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <complex>
#include<string>
#include<sys/types.h>
#include<sys/stat.h>
#include<random>

#include "../../header/creating_spin_system/semi_classical_spins.h"




using namespace boost::numeric::odeint;
using boost::math::quadrature::trapezoidal;




// class to create 2D SC spin systems to calculate pair correlations for
// this class either uses the autocorrelation for mean field calculation

//in iterate function different differential equations can be solved
//but they have to be plugged in manually
// furthermore are all imoportant constants like seeds or the B field golbally defined in matrices_and_const.h

class spin_system_2D : public spin_system {
protected:

public:

//constructor for 2D spin system
    spin_system_2D(int amount_per_row_init, uint step_amount_init, double step_width_init, string filename_init, bool print_text_init, lattice_database& database_init)
     : spin_system(amount_per_row_init, step_amount_init, step_width_init, filename_init, print_text_init, database_init){

        dimension = 2;
        // delete alle spins to create new ones for a 2D system
        spin_vec.clear();
        //create_spins_2d();

        //special case for 3 particle system
        if(amount_per_row_init == 1){
            //cout << "\nCreating 3 particle 2D system\n" << endl;
            create_spins_3_particle();
        }
        else{
            create_spins_2d();
        }

        //check if correct database

        if(database.return_database_dimension() != dimension){

            cout << "#### ERROR: DIMENSIONS IN DATABASE AND 2D SPIN SYSTEM DO NOT MATCH ####\n" << "Database dimension is = " << database.return_database_dimension() << " Spin system dimension is = " << dimension << endl;

        }

        //reassigning values
        system_size = amount_per_row * amount_per_row; // 2D system size is amount_per_row^2

        if(amount_per_row == 1){
            system_size = 3;
        }
        
        if(filename.find("2D") == std::string::npos){
            filename = "2D_system_" + filename;
        }

    }

    //function to create the 2D spin system
    //void create_spins() override;

    void create_spins_2d(){

        Eigen::Vector3cd spin;

        for(int x = 0; x < amount_per_row; x++){
            for(int y = 0; y < amount_per_row; y++){
                //test with second method to draw spin vector
                spin_vec.emplace_back(create_Spin_Vector2( ));
                //spin_vec.emplace_back(create_Spin_Vector_non_normalized()); //create_Spin_Vector2( ));
            
            }
        }

        if(spin_vec.size() != coupling_const.size()){
            cout << "#### ERROR: SPIN VECTOR SIZE DOES NOT MATCH COUPLING CONSTANT SIZE ####\n" << "Spin vector size = " << spin_vec.size() << " Coupling constant size = " << coupling_const.size() << endl;
            cout << "#### 2D SPIN SYSTEM NOT CREATED PROPERLY ####" << endl;
        }
    }

    //function to create a 2D spin system only made up of 3 spins in a "L" shape
    void create_spins_3_particle(){

        Eigen::Vector3cd spin1 = create_Spin_Vector2();
        Eigen::Vector3cd spin2 = create_Spin_Vector2();
        Eigen::Vector3cd spin3 = create_Spin_Vector2();

        spin_vec.emplace_back(spin1);
        spin_vec.emplace_back(spin2);
        spin_vec.emplace_back(spin3);

        if(spin_vec.size() != coupling_const.size()){
            cout << "#### ERROR: SPIN VECTOR SIZE DOES NOT MATCH COUPLING CONSTANT SIZE ####\n" << "Spin vector size = " << spin_vec.size() << " Coupling constant size = " << coupling_const.size() << endl;
            cout << "#### 2D SPIN SYSTEM NOT CREATED PROPERLY ####" << endl;
        }
    }

};