#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include<mpi.h>
#include<random>


#include <complex.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "../header/functions.h"
#include "../header/matrices_and_const.h"

#include "../header/creating_spin_system/semi_classical_spins.h"
#include "../header/creating_spin_system/iteration.cpp"
#include "../header/creating_spin_system/EoMs.cpp"
#include "../header/creating_spin_system/calc_observables.cpp"
#include "../header/creating_spin_system/spin_creation.cpp"


#include "../header/lattice_database/lattice_database.h"
#include "../header/lattice_database/lattice_creation.cpp"
#include "../header/lattice_database/pair_correlations.cpp"

#include "../header/multi_system_averaging.h"
#include "../header/old_header_files/differential_equations.h"
#include "../header/random_generator_mpi.h"






//#include "../header/TimeMeasure.h"





using cplx = std::complex<double>;
using namespace std;
using Observable = Eigen::MatrixXcd;

//cd "d:\Dateien\Uni\Semester VI\Bachelorarbeit\Bachelorbums\testing_grounds\" ; if ($?) { g++  -std=c++17 -pedantic -Wall -Wextra -I ../header/Eigen -o testing_grounds *.cpp } ; if ($?) { .\testing_grounds }
//cd "/home/sander/work/bachelor/Bachelorbums/testing_grounds/" && g++ testing_grounds.cpp -I /net/el7/eigen/3.4.0/include/ -o testing_grounds && "/home/sander/work/bachelor/Bachelorbums/testing_grounds/"testing_ground
//cd "/data/sander/bachelor/Bachelorbums/testing_grounds/" && g++ testing_grounds.cpp -o testing_grounds && "/data/sander/bachelor/Bachelorbums/testing_grounds/"testing_grounds

//-I /net/el7/eigen/3.4.0/include/
// ssh-add -l -E sha256


//cd "/home/sander/data/Masterarbeit/testing_grounds/" && g++ testing_grounds.cpp -I /net/el7/eigen/3.4.0/include/ -o testing_grounds && "/home/sander/data/Masterarbeit/testing_grounds/"testing_grounds


int main(int argc, char* argv[]){

    // ====== Initialize MPI ======
    int world_size, my_rank; // world_size is the number of cores, my_rank is the number of "this" core 
    MPI_Init( nullptr, nullptr );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );



    if(my_rank == 0){
        build_folder();
    }

    if(my_rank == 0){
	cout << "\n\nAmount of cores used " << world_size << endl;
    }





    if(my_rank == 0){

        cout << "\n-------------------------- Create Comparison Data using MPI --------------------------------\n"<<endl;
    }

    // ====== Set up basic constants ======

    //dimensiomn of the system
    //default value 3D
    uint dimension = 3; // 3D system

    
    //initalizing constants for lattice
    vector<int> row_lengths = {6};
    vector<int> row_lengths_2D = {5}; //default value if 2D system is called at runtime

    uint steps = 850;//564;//450;
    double dt = 2.5 * pow(10,-11);

    
    //a as lattice constant defined in matrices_and_const
    double a = lattice_const;

    //amount of systems used for averaging
    // total system amount is system_amount times core amount 
    // 70 for 42000 175 for 105000 
    // 40 for 24000 and 8 for 4800
    uint system_amount = 70; 
    //bool to determine if pair correlations are calculated. 
    //This has runtime reasons
    bool calc_pair_corr = false;
    
    string dimension_string = "";


    //change step size adjustments for different magnetic field directions

    if(B_field_direction == Eigen::Vector3cd(0,1,1)){
        steps = 1090;
    }
    else if(B_field_direction == Eigen::Vector3cd(1,1,1)){
        steps = 1550;
    }


    //coupling cutoff for lattice database
    //int coupling_cut_off = 0; //0 means all distances are calculated
    double coupling_cut_off = 0; 
    
    string cut_off_name = "";
    if(coupling_cut_off != 0){
        cut_off_name = "cut_off_" + to_string(coupling_cut_off) +"_";
    }
 

    vector<string> filenames;






    // is data supposed to be saved ?
    bool print_text = false;


    //#### make it possible to pass amount of lattice sitees per row during runtime
    //I hate to use goto, but i think its fairly simple and readable so who gives a f

    if(argc > 2){
        if(string(argv[1]) == "2D"){
            if(row_lengths_2D.size() > 1 && my_rank == 0){
                cout << "Error only one value for 2D row length allowed" << endl;
                cout << "Using default value " << row_lengths_2D[0] << endl;
                goto jmp;
            }
            row_lengths_2D[0] = stoi(argv[2]);
            if(my_rank == 0){
                cout << "\nNew row length for 2D system is " << row_lengths_2D[0] << endl;
            }
        }
        else{
            if(row_lengths.size() > 1 && my_rank == 0){
                cout << "Error only one value for 3D row length allowed" << endl;
                cout << "Using default value " << row_lengths[0] << endl;
                goto jmp;
            }
            row_lengths[0] = stoi(argv[2]);
            if(my_rank == 0){
                cout << "\nNew row length for 3D system is " << row_lengths[0] << endl;
            }
        }
    }
    jmp:

    //########### Redfining certain values if 2D is passed as an argument at runtime to main()
    string comparer = "";
    if(argc > 1){
        comparer = string(argv[1]);
    }
    
    if( comparer == "2D"){
        if(my_rank == 0){
            cout << "Considering 2D system" << endl;
        }
        dimension_string = "2D_system_";

        dimension = 2;
        row_lengths = row_lengths_2D;

        //halving the step amount from 7500
        steps = 3000;        
        dt = 0.005;

        //check if other b_fields than (0,0,1) are used
        if(B_field_direction != Eigen::Vector3cd(0,0,1) & my_rank == 0){
            cout << "\n\n######### ERROR B-FIELD DIRECTION FOR 2D IS " << B_field_direction.real().transpose() << " ###########" <<  endl;
            B_field_direction << 0,0,1;
            cout << "\nThis is changed to " << B_field_direction.real().transpose() << endl;
        }

        //turn on calc of pair correlation
        //if 2D system is considered FID is basicallly trivial and the pair correlations are of interest
        calc_pair_corr = true;
    }

    string B_field_name = "[" + to_string(int(B_field_direction(0).real())) + to_string(int(B_field_direction(1).real())) + to_string(int(B_field_direction(2).real())) + "]" ;


    if(my_rank ==0){
        cout << "\n-------------------------- Calculate averaged correlation with MPI for a " << dimension << "D system --------------------------------\n"<<endl;
    }

    if(my_rank ==0){
        cout << "\nPAIR CORRELATIONS ARE EVALUATED " << boolalpha<< calc_pair_corr  << endl;
        cout << "\nCUTOFF PARAMETER FOR COUPLING IS " << coupling_cut_off<< endl;
        cout << "\nDIMENSION IS " << dimension << endl << endl;
    }
    
    for(int i = 0; i < row_lengths.size(); i++){
        //reset generator for identical initial conditions
        filenames.emplace_back(dimension_string + to_string(world_size * system_amount) + "sys_mpi_no_norm_" + cut_off_name + to_string(row_lengths[i]) + "P_comparison_" + B_field_name);
        //cout << filenames[i] << endl;

    }

    TimeMeasure calc;


    for(int i = 0; i < row_lengths.size(); i++){
        //reseed for identical initial conditions
        // ====== create random number generator on each core ======
        //seed is located in matrices_and_const.h
        std::string seed_str = to_string(seed); // "random"
        auto generated_seed = generate_seed( seed_str, my_rank );
        std::mt19937 engine{ throw_seed( generated_seed, my_rank ) };

        //still use generator gen in matrices_and_const
        //code supplied by Timo
        gen = engine;

        //create data by averaging over different initial conditions
        multi_system_averaging average(dimension, a, steps, dt, system_amount, row_lengths[i], filenames[i], my_rank, coupling_cut_off);
        average.create_averaged_data(calc_pair_corr);

        //evaluate the data and merge different cores by utilizing ALLREDUCE
        //functions are found in functions.h
        evaluate_correlations_mpi(average, filenames[i], dt, system_amount, world_size, my_rank);
        if( my_rank == 0 ){
            //calc.measure("All Calculations regarding normal Correlations");
        }
        evaluate_pair_correlations_mpi(average, filenames[i], dt, system_amount, world_size, my_rank, calc_pair_corr);
        if( my_rank == 0 ){
            //calc.measure("All PAIR Correlation Calculations");
        }
        evaluate_interim_correlation_results_mpi(average, filenames[i], dt, system_amount, world_size, my_rank);
        if( my_rank == 0 ){
            //calc.measure("All interim result Calculations");
        }
        evaluate_interim_auto_and_pair_correlation_results_mpi(average, filenames[i], dt, system_amount, world_size, my_rank);
        if( my_rank == 0 ){
            //calc.measure("All interim result Calculations");
        }
        evaluate_variance_at_certain_time(average, filenames[i], dt, system_amount, world_size, my_rank);
        if( my_rank == 0 ){
            //calc.measure("All Variance Calculations");
        }
    
   
        // ====== inform user about result ======
        if( my_rank == 0 ){
            //std::cout << "\ncomputation finished\n" <<endl; 
            calc.measure(to_string(row_lengths[i]) + "P comparison calculation done");
        }
    

    }


    //Calculate for Zeeman Hamiltonian

    /*
    print_text = true; 

    spin_system zeeman(a, 6, steps, dt,"Zeeman_", print_text);
    zeeman.do_Zeeman_iteration();

    calc.measure("Finished Calculations on Zeeman system");


    */

    // ====== Finalize MPI ======
    MPI_Finalize(); 

    calc.stop();

    return 0;

    
}
