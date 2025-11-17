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

//other dependancies
#include "../../header/matrices_and_const.h"
#include "../../header/lattice_database/lattice_database.h"


#include "../../header/old_header_files/differential_equations.h"

#include "../../header/TimeMeasure.h"



using namespace boost::numeric::odeint;
using boost::math::quadrature::trapezoidal;




// class to create a simple cubic lattice lattice 

//in iterate function different differential equations can be solved
//but they have to be plugged in manually
// furthermore are all imoportant constants like seeds or the B field golbally defined in matrices_and_const.h




class spin_system{

protected:



//basic constants regarding the lattice

uint dimension = 3; // default value 3

int amount_per_row;
int system_size;

uint step_amount;
double step_width;
uint current_step = 0;

//name for file
string filename;

// bool if data should be saved as a .txt
// standard value false due to flooding the terminal with text for averaging over system
bool print_text = false;


//add two new bools to determine if there is supposed to be a cut off in the calculation of pair correlations
//and if the xxz heisenberg hamiltonian has to be used instead of the isotropic one
//################IMPORTANT LINE IS ONLY EDITED IN THIS FILE#####################
//ONLY FOR 2D SYSTEM"
bool use_pair_correlation_cut_off_2D = false;
bool use_xxz_hamiltonian_2D = false;


vector<Eigen::Vector3cd> spin_vec{}; // vector of the semi classical Spins

//Vector all semi classical spin configurations get saved into, to calculate pair correlatoions in the end
//needed for every time step
vector <vector<Eigen::Vector3cd>>all_spin_vec_configurations{};



//"Matrices" filled with distances and connecting vectors. Is supposed to save runtime
// arrays are called in a way, that the first index is for the starting lattice site and the second for the site
// for which the connection is to be evaluated

//lattice database to save computation time
//is provided as a reference

lattice_database& database;

//these are filled with values from the lattice database
//data object is vector of vectors, where the first index is the starting lattice site and the second the target lattice site
//and the EIgen::Vector then is the connecting vector normalized
vector< vector<Eigen::Vector3cd> > connecting_vector_norm; //normalized connecting vector
vector< vector<double> > distances;


//coupling constant for dipolar hamiltonian
//also provided from database
vector< vector<double> > coupling_const;

//coupling constant for heisenberg hamiltonian
vector< vector<double> > coupling_const_heisenberg;


//create dxdt_Zero to make solving of EOMs more efficient

vector<Eigen::Vector3cd> dxdt_Zero{};



//Measured quantities

vector<Eigen::Vector3cd> magnetization{}; // vector of magnetization
vector<Eigen::Vector3cd> magnetization_CS{}; // vector of magnetization for central spin located at (0,0,0)
vector<double> magnetization_correlation{}; // Complex number of magnetization calculated from magnetization. Would like this to be a double

//Quantities related to pair correlations

vector<vector<double>> pair_correlations{};

//this butt ugly variable consists of a vector of all index pairs which correspond to a specific vector in the lattice. This vector then is in a vector of vecotrs corresponding to the same correlation grouped together
//and this vecotr is then included in another vector, which accomodates these vectors for all groups.
vector< vector< vector< pair<uint, uint> > > > all_correlation_groups_lattice_pairs{};


TimeMeasure calc;

public:


spin_system(int amount_per_row_init, uint step_amount_init, double step_width_init, string filename_init, bool print_text_init, lattice_database& database_init) : database(database_init){

    dimension = 3; //default value 3
    amount_per_row = amount_per_row_init;
    step_amount = step_amount_init;
    step_width = step_width_init;
    system_size = amount_per_row * amount_per_row * amount_per_row;

    filename = filename_init;
    print_text = print_text_init;

    create_spins();



    //get values from database
    database = database_init;


    connecting_vector_norm = database.return_normalized_connections();
    distances = database.return_distances();
    coupling_const = database.return_dipolar_coupling_const();
    coupling_const_heisenberg = database.return_heisenberg_coupling_const();
    dxdt_Zero = database.return_dxdt_Zero();
    all_correlation_groups_lattice_pairs = database.return_all_lattice_pairs();





    if(print_text == true){
        cout << "\nA cubic lattice with " << system_size << " lattice sites has been created"<< endl;
        cout << "Iteration is done with "<< step_amount << " steps and a step width of " << step_width << endl;
        cout << "Calculations will end at T = " << step_amount * step_width  << "µs "<< endl; 
        cout << "\"Larmor\" frequency approx " << 2. * mu_theo<< " or " << gamma_i * B_field_direction.norm() << "10^7 rad / s"<< endl;
        cout << "Magnetic field is in direction " << B_field_direction.real().transpose() << endl ;
        cout << "Amount of vector correlation groups found is " << all_correlation_groups_lattice_pairs.size() << endl<< endl;
    }


};




// function that creates spins for each lattice site. Each component is normal distributed.

void create_spins();

//Randomly drawing a spcin vector by drawing each component independently following a normal distribution

Eigen::Vector3cd create_Spin_Vector_non_normalized();

//alternative way to draw the spin vectors randomly. This time by utilizing a vector in spherical coordinates

Eigen::Vector3cd create_Spin_Vector2();







//////// Functions important for the calculations
 
//calc other quantities

//function which calculates magnetization or correlations when called and appends value to a vector

void calc_magnetization();

void calc_magnetization_central_spin();



//calculates magnetization for all time steps at the same time
//These are the functions actually used

void calc_all_magnetization_correlation();

//function that sums over pairs, to geht vector that is used to calculate pair correlations
//should be called after every ime step

vector< double > calc_pair_correlation_for_index_pairs(const vector< pair< uint, uint> >&);
    



void calc_all_pair_correlations();

void calc_auto_correlation();

void calc_auto_correlation1();


//################ Differential equations for semi classical approach



// Bfield in matrices and const



Eigen::Vector3cd EoM_Dipolar_component(uint, vector<Eigen::Vector3cd>& );



void EoM_Dipolar(vector<Eigen::Vector3cd>&, vector<Eigen::Vector3cd>& , const double );



//Wenn diese Rechnng funktioniert skaliert sie nur noch mit (n^2 + 3n)/2 statt n^2 ( sogar noch weniger wegend keiner selbstww, bin aber gerade)
// this equation saves time by using loops in a more effective was, but ist loses time due to calculations for each lattice site being done at the same moment

void EoM_Dipolar_new(vector<Eigen::Vector3cd>& , vector<Eigen::Vector3cd>&, const double);


// Equations used for a Heisenebrg system. Mostly used for 2D systems

void EoM_Heisenberg(vector<Eigen::Vector3cd>&, vector<Eigen::Vector3cd>& , const double );

//Equations used for a Heisenebrg XXZ system. Mostly used for 2D systems

void EoM_Heisenberg_XXZ(vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt , const double t );


/////////// Functions for the iteration




void iterate_over_time();


void iteration_step(function<void(vector<Eigen::Vector3cd>&, vector<Eigen::Vector3cd>&, double)>,const double );

/*

//implement iteration using a stepper with step width adjusrment

void iterate_over_time_self_adjusting();


//new rk4 stepper with step adjustment
void iteration_step_self_adjusting(double time_step);

*/

//Special case for Zeeman iteration
//only allowed if no other calculations are done in this class instance

void do_Zeeman_iteration();





/////// FUnctions which return private data
////// used to gain access to average over multiple systems 

vector<double> return_correlation(){
    if(magnetization_correlation.size() == 0){
        cout << "#### ERROR REGARDING CORRELATION: NOTHING TO RETURN ####\n" << "Correlation vector size = " << magnetization_correlation.size()<<endl;
    }
    return magnetization_correlation;
}


//also calls calculation in this function to get the possibility in multi_system_averaging to NOT calc pair correlations
vector <vector<double> > calc_and_return_pair_correlations(){

    calc_all_pair_correlations();
    //calc_auto_correlation();

    if(pair_correlations.size() == 0){
        cout << "#### ERROR REGARDING PAIR CORRELATION: NOTHING TO RETURN ####\n" << "Pair correlation vector size = " << pair_correlations.size()<<endl;
    }


    return pair_correlations;
}

/////// Functions to save Data


void save_magnetization(){

    if(magnetization[step_amount - 1](0).imag() > pow(10., -10.)){
    cout << "#### ERROR REGARDING MAGNETIZATION: SEEMS TO BE IMAGINARY ####\n" << "" << magnetization[step_amount - 1](0)<<endl;
    }


    std::ofstream myfile;

    myfile.open ("build/daten/" + filename + "Magnetization.txt");
    myfile << "M_x \t M_y \t M_z \t time\n";


    for(int i = 0; i < magnetization.size(); i++){
        myfile <<magnetization[i](0).real()<< "\t"<< magnetization[i](1).real()<< "\t"<<  magnetization[i](2).real()<< "\t"<<  i * step_width <<endl;
    }; 

    myfile.close();
    
}

void save_magnetization_central_spin(){


    std::ofstream myfile;

    myfile.open ("build/daten/" + filename + "Magnetization_CS.txt");
    myfile << "S_x \t S_y \t S_z \t time\n";

    for(int i = 0; i < magnetization.size(); i++){
        myfile <<magnetization_CS[i](0).real()<< "\t"<< magnetization_CS[i](1).real()<< "\t"<<  magnetization_CS[i](2).real()<< "\t"<<  i * step_width <<endl;
    }; 

    myfile.close();
    
}

void save_correlation(){


    std::ofstream myfile;

    myfile.open ("build/daten/" + filename + "Correlation.txt");
    myfile << "G_x time\n";

    for(int i = 0; i < magnetization_correlation.size(); i++){
        myfile <<magnetization_correlation[i]<< "\t"<<  i * step_width <<endl;
    }; 

    myfile.close();
    
}








};


