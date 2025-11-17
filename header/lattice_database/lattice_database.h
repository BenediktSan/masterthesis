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

//other dependancies
#include "../../header/matrices_and_const.h"


//#include "../../header/creating_spin_system/lattice_creation.cpp"


using namespace boost::numeric::odeint;
using boost::math::quadrature::trapezoidal;




// class to create a simple cubic lattice lattice 

//one instance should be given to the spin system class. That way not every spin system instance has to calculate 
//all distances and coupling constants. Spin system only has to 
// furthermore are all imoportant constants like seeds or the B field golbally defined in matrices_and_const.h




class lattice_database{

protected:

// is defined upon use of create database fucntion
uint dimension;

//basic constants regarding the lattice

double lattice_const; 
int amount_per_row;
int system_size;


//constant for a cutoff distance for the dipolar interaction
// this is the maximum distance for which the dipolar interaction is calculated
//if this is set to 0, all distances are calculated
double coupling_cut_off;



vector<Eigen::Vector3cd> R_vec{};   // vecotr of the lattice sites


//"Matrices" filled with distances and connecting vectors. Is supposed to save runtime
// arrays are called in a way, that the first index is for the starting lattice site and the second for the site
// for which the connection is to be evaluated

//vector< vector<double> > distance;
vector< vector<Eigen::Vector3cd> > connecting_vector;
vector< vector<Eigen::Vector3cd> > connecting_vector_norm; //normalized connecting vector
vector< vector<double> > distances;


//coupling constant for dipolar hamiltonian

vector< vector<double> > coupling_const;

vector<double> coupling_const_sorted;

//coupling constant for heisenberg hamiltonian

vector< vector<double> > coupling_const_heisenberg;

//create dxdt_Zero to make solving of EOMs more efficient

vector<Eigen::Vector3cd> dxdt_Zero{};


//######PAir correlations and corresponding vectors used for their calculation


//This list includes the different correlation vectors grouped together if contribute to the same correlation type
// after being calculatet in the class constructer they should be sorted by lenth of vectors in these groups (Verey grou has vector of the same length by construction)
vector< vector<Eigen::Vector3cd > > list_of_correlation_vector_groups{};

//This includes the multiplicity corresponding to these groups aka the size of these groups
vector< uint > multiplicity_of_correlation_groups{};


//this butt ugly variable consists of a vector of all index pairs which correspond to a specific vector in the lattice. This vector then is in a vector of vecotrs corresponding to the same correlation grouped together
//and this vecotr is then included in another vector, which accomodates these vectors for all groups.
vector< vector< vector< pair<uint, uint> > > > all_correlation_groups_lattice_pairs{};




public:


lattice_database(double lattice_const_init, int amount_per_row_init, double coupling_cut_off_init){

    //cout << "\nDatabase constructor has been called" << endl;
    lattice_const = lattice_const_init;
    amount_per_row = amount_per_row_init;
    //sytsem_size is initialised in initialize_... functions

    coupling_cut_off = coupling_cut_off_init;



};

//function that initializes system after creation of this class
// this is done, so that the create_lattice() function can be virtual and overwritten in the 2D lattice class

void initialize_database_3D(){

    dimension = 3;
    system_size = amount_per_row * amount_per_row * amount_per_row;


    create_lattice_3D();


    //fill the distance arrays

    create_distance_arrays();

    create_dipolar_coupling_const();

    create_heisenberg_coupling_const();



    //find alll correlations and group them by their type
    //this calculates the groups and multiplicity_of_groups
    //afterwards they are sorted according to the lenth of vectors in these groups

    find_correlation_groups();
    sort_correlation_groups();


    //find all lattice pairs for the pair corrleation calculations later on
    calc_all_lattice_pairs();


}

void initialize_database_2D(){

    dimension = 2;
    system_size = amount_per_row * amount_per_row ;

    if(system_size == 1){
        //special case for three particle system
        cout << "\n############### ERROR: WRONG LATTICE CREATION FUNCTION HAS BEEN CALLED ###########" << endl;
        cout << "############### CALLING 3 PARTICLE LATTICE CREATION FUNCTION WOULD HAVE BEEN REQUIRED ###########\n" << endl;
        dimension = 2;
        amount_per_row = 1;
        system_size = 3 ;
        
        create_lattice_3_particle();


        return;
    }

    create_lattice_2D();

    //fill the distance arrays

    create_distance_arrays();

    create_dipolar_coupling_const();
    create_heisenberg_coupling_const();


    //find alll correlations and group them by their type
    //this calculates the groups and multiplicity_of_groups
    //afterwards they are sorted according to the lenth of vectors in these groups

    find_correlation_groups();
    sort_correlation_groups();


    //find all lattice pairs for the pair corrleation calculations later on
    calc_all_lattice_pairs();


}

void initialize_database_3_particle(){


    dimension = 2;
    amount_per_row = 1;
    system_size = 3 ;

    

    create_lattice_3_particle();

    //fill the distance arrays

    create_distance_arrays();

    create_dipolar_coupling_const();
    create_heisenberg_coupling_const();


    //find alll correlations and group them by their type
    //this calculates the groups and multiplicity_of_groups
    //afterwards they are sorted according to the lenth of vectors in these groups

    find_correlation_groups();
    sort_correlation_groups();


    //find all lattice pairs for the pair corrleation calculations later on
    calc_all_lattice_pairs();


}

//function that gets distance, but includes periodic boundaries


//calc other quantities
//function returns connecting vector from i to j

Eigen::Vector3cd calc_connecting_vec(uint i, uint j){       
    if(i == j){
        cout << "\n\n\tERROR VEKTOR ZU SICH SELBST BERECHNET\n\n" << endl;
    } 
    double cube_length = amount_per_row * lattice_const;

    double dist_x =real( R_vec[j](0) - R_vec[i](0) );
    double dist_y =real( R_vec[j](1) - R_vec[i](1) );
    double dist_z =real( R_vec[j](2) - R_vec[i](2) );


    double tol = 1e-8;    





    if( abs(dist_x)  > cube_length * 0.5 + tol){
        if(dist_x > 0){
            dist_x = dist_x - cube_length ;
        }
        else{
            dist_x = dist_x + cube_length ;
        } 
    }
    if(abs(dist_y) > cube_length * 0.5 + tol){
        if(dist_y > 0){
            dist_y = dist_y - cube_length ;
        }
        else{
            dist_y = dist_y + cube_length ;
        }
    }
    if(abs(dist_z) > cube_length * 0.5 + tol){
        if(dist_z > 0){
            dist_z = dist_z - cube_length ;
        }
        else{
            dist_z = dist_z + cube_length ;
        }
    }

    //now consider a special case where due to periodic boundary conditions
    //different vectors map to the same lattice site
    // this only happens for even lattcie sizes
    // and leads to an underestimation of the multiplicity if not considered
    // also when calculating lattice pairs connecting lattice sites are not found,
    //because (0,0,-2) is not considered for a 4x4x4 lattice
    // but (0,0,2) is

    if(amount_per_row % 2 == 0){
        if(abs(abs(dist_x) - cube_length * 0.5) < tol && dist_x < 0){
            dist_x = -1. * dist_x;
            //cout << "Corrected dist_x for periodic BCs for vector " << dist_x/lattice_const << " " << dist_y/lattice_const << " " << dist_z/lattice_const << endl; 
        }
        if(abs(abs(dist_y) - cube_length * 0.5) < tol && dist_y < 0){
            dist_y = -1. * dist_y;
            //cout << "Corrected dist_y for periodic BCs for vector " << dist_x/lattice_const << " " << dist_y/lattice_const << " " << dist_z/lattice_const << endl; 
        }
        if(abs(abs(dist_z) - cube_length * 0.5) < tol && dist_z < 0){
            dist_z = -1. * dist_z;
            //cout << "Corrected dist_z for periodic BCs for vector " << dist_x/lattice_const << " " << dist_y/lattice_const << " " << dist_z/lattice_const << endl; 
        }
    }

    Eigen::Vector3cd vec;
    vec << dist_x, dist_y, dist_z;

    double r = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

    if( r > amount_per_row * lattice_const * sqrt(3)){
        cout << "\n\n\tERROR LENGTH OF VECTOR LONGER THAT MAXIMUM VALUE ALLOWED BY BOUNDARY CONDITIONS\n\n" << r << endl;
    }
    if (sqrt(dist_x * dist_x) > amount_per_row * lattice_const * sqrt(2) || sqrt(dist_y * dist_y) > amount_per_row * lattice_const * sqrt(2) || sqrt(dist_z * dist_z) > amount_per_row * lattice_const * sqrt(2)){
        cout << "\tERROR ONE COMPONENT OF VECTOR LONGER THAT MAXIMUM VALUE ALLOWED BY BOUNDARY CONDITIONS" << endl;
        cout << "\tdist_x " << dist_x / lattice_const << " dist_y " << dist_y / lattice_const << " dist_z " << dist_z / lattice_const << endl;
        cout << "\tAllowed maximum " << amount_per_row * sqrt(2) << endl <<endl;
    }

    return vec;
}


//////// Functions important for creation of my system
void create_lattice_3D();

void create_lattice_2D();

//function to create a 2D lattice with only 3 particles in an "L" shape
void create_lattice_3_particle();

//function that creates arrays which are filled with the distances and connecting vectors between lattice sites
// arrays are called in a way, that the first index is for the starting lattice site and the second for the site
// for which the connection is to be evaluated

void create_distance_arrays();


//create a vector with all required coupling constants for the dipolar interaction
//all required constants are in matrices_and_....
//NOTE: coupling constant vanishes if angle 54° between B-field and connecting vector is. 

void create_dipolar_coupling_const();


//create the coupling constant for a heisenberg hamiltonian

void create_heisenberg_coupling_const();


//functions to return calculated distanxes and coupling constants.
//this is being used in semi_classical_soins.h

uint return_database_dimension(){
    return dimension;
}

vector< vector<Eigen::Vector3cd> > return_normalized_connections(){
    return connecting_vector_norm;
}

vector< vector<double> > return_distances(){
    return distances;
}

vector< vector<double> > return_dipolar_coupling_const(){
    return coupling_const;
}

vector< vector<double> > return_heisenberg_coupling_const(){
    return coupling_const_heisenberg;
}

vector<Eigen::Vector3cd> return_dxdt_Zero(){
    return dxdt_Zero;
}


vector<vector<vector<pair<uint, uint>>>> return_all_lattice_pairs(){
    return all_correlation_groups_lattice_pairs;
}

int return_amount_of_correlation_groups(){
    return  multiplicity_of_correlation_groups.size();
}
//################Functions to calculate paircorrelations



//this function finds all possible pair correlations in my lattice or in my cutoff range
//returns a vecor filled with exemplary vectors which stand for a certain correlation e.g. a neraest neighbour correlation in x-direction

void find_correlation_groups();









//used to check which vectors satisfy my conditions
//also function that uses indices of functions that satisfy that to calculate correlations
//returned vector includes all index pairs for that certain wanted correlation

vector< pair<uint, uint> > find_suiting_lattice_pairs(Eigen::Vector3cd exemplary_vector);

//functions that calculates for every vector in evry correlation type group all suiting lattice pairs ans saves
//This is later returned and saved to lattice_system

void calc_all_lattice_pairs();

//function that sorts the correlations according to their connecting vector length
//also sorts the multiplicitys the same way
//mainly written by chatgpt

void sort_correlation_groups();


//function that prints correlation groups and their multiplicity to the terminal

void print_correlation_groups();






};
