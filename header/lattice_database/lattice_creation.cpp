
#include <iostream>
#include "../../header/lattice_database/lattice_database.h"




//////// Functions important for creation of my system
//this function creates the lattice
// *The first amount_per_row^2 indices in this "array" are in the z-y plane
void lattice_database::create_lattice_3D(){

    if( amount_per_row > 11){
        cout << "###### WARNING: CREATING A 3D LATTICE WITH MORE THAN 11 SITES PER ROW ######\n" << endl;
    }

    Eigen::Vector3cd R_lattice_site;

    for(int x = 0; x < amount_per_row; x++){
        for(int y = 0; y < amount_per_row; y++){
            for(int z = 0; z < amount_per_row; z++){
                R_lattice_site(0) = x * lattice_const;
                R_lattice_site(1) = y * lattice_const;
                R_lattice_site(2) = z * lattice_const;
                //cout << R_lattice_site.transpose() << endl;
                R_vec.emplace_back(R_lattice_site);
            }
        }
    }
}

//Funtion to create 2D lattice
void lattice_database::create_lattice_2D() {

        //cout << "\nCreating 2D lattice" << endl;

        Eigen::Vector3cd R_lattice_site;

        for(int x = 0; x < amount_per_row; x++){
            for(int y = 0; y < amount_per_row; y++){
                R_lattice_site(0) = x * lattice_const;
                R_lattice_site(1) = y * lattice_const;
                R_lattice_site(2) = 0;
                //cout << R_lattice_site.transpose() << endl;
                R_vec.emplace_back(R_lattice_site);
            }
        }
    }

//function to create a 2D spin system only made up of 3 spins in a "L" shape
void lattice_database::create_lattice_3_particle() {


        //cout << "\nCreating 3 particle lattice" << endl;

        Eigen::Vector3cd R_lattice_site; 

        // Create the three particles in an "L" shape
        R_lattice_site(0) = 0;
        R_lattice_site(1) = 0;
        R_lattice_site(2) = 0;
        R_vec.emplace_back(R_lattice_site);

        R_lattice_site(0) = lattice_const;
        R_lattice_site(1) = 0;
        R_lattice_site(2) = 0;
        R_vec.emplace_back(R_lattice_site);

        R_lattice_site(0) = lattice_const;
        R_lattice_site(1) = lattice_const;
        R_lattice_site(2) = 0;
        R_vec.emplace_back(R_lattice_site);

    }



// for which the connection is to be evaluated

void lattice_database::create_distance_arrays() {

    //fill the arrays witch 0 values to be able to Index in the next step
    
    vector<double> zero_scalar;
    vector<Eigen::Vector3cd> zero_vector;

    for(int i = 0; i < system_size; i++){
        //fill vector filled with 0s to be placed multiple times in distances etc.
        zero_scalar.emplace_back(0.);
        zero_vector.emplace_back(Eigen::Vector3cd::Zero());


    }

    //assigning this vector here to the experimental element dxdt_zero
    //its used to make solving the dipole EOMs more efficient

    dxdt_Zero = zero_vector;


    for(int i = 0; i < system_size; i++){

        //fill distances, connecting vectors and normalized connecting vectors with  0
        distances.emplace_back(zero_scalar);
        connecting_vector.emplace_back(zero_vector);
        connecting_vector_norm = connecting_vector;

    }

    //special case for three particle system
    if(amount_per_row == 1){
        distances[0][0] = 0;
        distances[1][1] = 0;
        distances[2][2] = 0;

        distances[0][1] = lattice_const;
        distances[1][0] = lattice_const;
        distances[1][2] = lattice_const;
        distances[2][1] = lattice_const;
        distances[0][2] = sqrt(2.) * lattice_const;
        distances[2][0] = sqrt(2.) * lattice_const;

        Eigen::Vector3cd vec_00;
        vec_00 << 0, 0, 0;
        connecting_vector[0][0] = vec_00;
        connecting_vector[1][1] = vec_00;
        connecting_vector[2][2] = vec_00;

        Eigen::Vector3cd vec_01;
        vec_01 << lattice_const, 0, 0;
        connecting_vector[0][1] = vec_01;
        connecting_vector[1][0] = -1 * vec_01;

        Eigen::Vector3cd vec_12;
        vec_12 << 0, lattice_const, 0;
        connecting_vector[1][2] = vec_12;
        connecting_vector[2][1] = -1 * vec_12;

        Eigen::Vector3cd vec_02;
        vec_02 << lattice_const, lattice_const, 0;
        connecting_vector[0][2] = vec_02;
        connecting_vector[2][0] = -1 * vec_02;

        return;
    }




    // fill the arrays with prober values

    for(int i = 0; i < system_size; i++){
        //goal of indexing to do every calculation only once
        //use this ineffictive code to catch all pair correlations later on
        for(int j= 0; j < system_size; j++){
        
        if(i == j){
            continue;
        }

        Eigen::Vector3cd vector = calc_connecting_vec(i,j);
        connecting_vector[i][j] = vector;
        //connecting_vector[j][i] = -1 * vector;        

        double distance_value = vector.norm();
        distances[i][j] = distance_value;
        //distances[j][i] = distance_value;

        //also add normalized vectors
        Eigen::Vector3cd norm_vector =vector.normalized();
        connecting_vector_norm[i][j] = norm_vector;
        //connecting_vector_norm[j][i] = -1 * norm_vector;

        }


    }

    //cout << "All required distances have been calculated\n" << endl;
}


//create a vector with all required coupling constants for the dipolar interaction
//all required constants are in matrices_and_....

//NOTE: coupling constant vanishes if angle 54° between B-field and connecting vector is. 
//This happens an easily achievable angle in a cubic lattice for in ddirection (0,0,1)
void lattice_database::create_dipolar_coupling_const(){
//fill the arrays witch 0 values to be able to Index in the next step
    
    vector<double> zero_scalar;

    for(int i = 0; i < system_size; i++){
        //fill vector filled with 0s to be placed multiple times in distances etc.
        zero_scalar.emplace_back(0.);


    }

    for(int i = 0; i < system_size; i++){

        //fill dcoupling_const with  0
        coupling_const.emplace_back(zero_scalar);

    }


    // fill the arrays with prober values

    for(int i = 0; i < system_size; i++){

        //uint that counts how many coulings i consider
        int counted_sites = 0;

        //goal of indexing to do every calculation only once
        for(int j= i + 1; j < system_size; j++){



        //check if distance is smaller than cutoff
        // summand added to account for computaion errors from distance calculation
        //this does make a difference and is therefore added
        if( coupling_cut_off > 0 && distances[i][j] - 1 * pow(10,-14)> coupling_cut_off * lattice_const){
            //cout << "cutoff distance reached" << endl;
            continue;   
        }
        counted_sites += 1;
        //calculate coupling constant only if distance is not 0

        //assign connecting vector and its length
    
        Eigen::Vector3cd vec_r_norm = connecting_vector_norm[i][j];
        double r =  distances[i][j];
        Eigen::Vector3cd B_field_direction_norm = B_field_direction.normalized();

        //calculate specific couping constant
        //used constants are defined in matrices and const.h 
        //hbar^2 = 1
        //factor 1/50 only for debugging purposes holds no physical value
        double d = ( 1. - 3. * pow( vec_r_norm.dot(B_field_direction_norm).real(), 2.))/2. * (mu_0) / (4. *M_PI) * (gamma_i * gamma_i) / pow(r,3);
        //d = ( 1. - 3. * pow( vec_r_norm.dot(B_field_direction_norm).real(), 2.))/2. * 1. / (pow(r/lattice_const,3) );

        coupling_const[i][j] = d;
        coupling_const[j][i] = d;


        

        }

        if(i == 0 && coupling_cut_off != 0){
            //cout << "Amount of coupling constants for site 0 is " << counted_sites << endl;
        }
    }        


    //cout << "nearest neighbour coupling const is " << coupling_const[0][1] << " ratio to furthest is " << coupling_const[0][int(amount_per_row / 2.)] / coupling_const[0][1] << endl;

}

//create the coupling constant for a heisenberg hamiltonian
void lattice_database::create_heisenberg_coupling_const(){
    //fill the arrays witch 0 values to be able to Index in the next step
    vector<double> zero_scalar;
    for(int i = 0; i < system_size; i++){
        //fill vector filled with 0s to be placed multiple times in distances etc.
        zero_scalar.emplace_back(0.);
    }



    for(int i = 0; i < system_size; i++){

        //fill dcoupling_const with  0
        coupling_const_heisenberg.emplace_back(zero_scalar);

    }

    // J = 1 means "FM" order
    double coupling_const_value = 1.0; // Heisenberg coupling constant value, can be adjusted

    //how may neares t neighbour coupling constants are considered
    //default value for nearest neighbour coupling constant is 1
    double coupling_cut_off = 1;

    
    for(int i = 0; i < system_size; i++){

        for(int j= i + 1; j < system_size; j++){
            
            if(i == j){
                continue;
            }
            // summand added to account for computaion errors from distance calculation
            //this does make a difference and is therefore added
            if( coupling_cut_off > 0 && distances[i][j] - 1 * pow(10,-14) > coupling_cut_off * lattice_const){
                 //cout << "cutoff distance reached" << endl;
                 continue;   
             }

         
            coupling_const_heisenberg[i][j] = coupling_const_value;
            coupling_const_heisenberg[j][i] = coupling_const_value;
        }

        
    }
}