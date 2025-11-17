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
#include<random>
#include<mpi.h>


#include "../header/creating_spin_system/semi_classical_spins.h"
#include "../header/lattice_database/lattice_database.h"
#include "../header/functions.h"


#include "../header/matrices_and_const.h"

#include "../header/2D_spin_system/2D_spin_system.h"

//#include " ../header/functions.h"

using namespace boost::numeric::odeint;
using boost::math::quadrature::trapezoidal;





// class to average over different initial conditions
// the calss spin_system ist used for all the calculations





class multi_system_averaging{

private:

//dimension of the system
//default is 3D 
uint dimension = 3;

//basic constants regarding the system and the averaging

double lattice_const; 
int amount_per_row;
int system_size;

uint step_amount;
double step_width;
uint current_step = 0;

//amounts of systems used for averaging

uint system_amount;

//name for file
string filename;

// bool used to initiate spin lattice. if true the data should of this system ist to be saved as a .txt
// generic value should be false
bool print_text = false;

//rank of current core
//used fr printing etc
int my_rank;

//database to be passed as reference to every spin_system to safe computation time
//by supplying all distances and coupling beforehand

lattice_database database;


//Vectors to save different quantities like correlation of the system from every inital conditon
//not sure if these are even used. Ill keep them around for a bit

vector< vector<double> > storage_correlation{};
vector< vector<Eigen::Vector3cd> > storage_magnetization{};

//vectors which contain all the averaged values

vector< double > averaged_correlation{};
//used for estimator of MC Error
vector< double > averaged_correlation_squared{};

//vector which contains all the different pair correlations for the different pair correlation groups

vector < vector <double> > averaged_pair_correlations{};

vector< Eigen::Vector3cd > averaged_magnetization{};

//object that saves interim correlation results
// used for one lengthy calculation with many initial condiotions to check how well the results converge
// e.g for 10000 i.c. also saves 1000,2000,...

vector< vector< double > > interim_correlation_results{};
//Inforrmation how many systems were used for a certain interim result
vector<int> interim_results_system_amount{};


// also save interim auto and pair correlations for certain correlaitons
// FID and pair correlation results share interim_resluts_system_amount

//correlation types that are considerd
//e.g 0 is for the auto correlation, 1 for c= 1 ....
vector<int> interim_pair_correlations_types = {0,1,2,3};

vector< vector< vector< double > > > interim_auto_and_pair_correlation_results{interim_pair_correlations_types.size()};






//object that for a certain iteration step saves all data data points of a all iniitial conditions for this point
//is later used to calculate variance of this point
//as point the extrema are choosen. Therefore different b-flieds need different points

//first value is the index corresponding to a certain time step, second value are the correlation values, third the squared sum of the correlation values
//index are asigned at inizialization
//corresponding indices also
tuple< vector<int>, vector< double >, vector<double> > expectation_values_at_certain_time{};


public:


multi_system_averaging(uint dimension_init, double lattice_const_init,uint step_amount_init, double step_width_init, uint system_amount_init, uint amount_per_row_init, string filename_init, int my_rank_init, double coupling_cut_off) : database(lattice_const_init, amount_per_row_init, coupling_cut_off){



    
    lattice_const = lattice_const_init;
    amount_per_row = amount_per_row_init;
    step_amount = step_amount_init;
    step_width = step_width_init;
    system_size = amount_per_row * amount_per_row * amount_per_row;

    system_amount =system_amount_init;

    filename = filename_init;

    my_rank = my_rank_init;

    if(my_rank == 0 && (dimension_init != 2 && dimension_init != 3)){
        cout << "#### ERROR: DIMENSION NOT SUPPORTED ####\n" << "Dimension = " << dimension_init << endl;
    }

    dimension = dimension_init;

    if(dimension == 2){
        system_size = amount_per_row * amount_per_row;
        if(amount_per_row == 1){
            system_size = 3;
        }
    }


    //initialize the database
    if(dimension == 2){
        if(amount_per_row == 1){
            database.initialize_database_3_particle();
        }
        else{
        database.initialize_database_2D();
        }
    }
    else{
        database.initialize_database_3D();
    }
    // when useing the extrema we nee d different points to evaluate them for different magnetic fields
    // this function assigns the correct time steps and their corrsponding indices
    find_time_steps_for_expectation_values();
    vector<double> zero_vec{};

    for(int i = 0; i < get<0>(expectation_values_at_certain_time).size(); i++){
        zero_vec.emplace_back(0.);
    }
    get<1>(expectation_values_at_certain_time) = zero_vec;
    get<2>(expectation_values_at_certain_time) = zero_vec;

    
    // already filling the averaged correlation vector with 0. to save time in later steps
    // this is done later for the pair correlation, because we do not know how many pair correlation types exist
    int T_max = int( T_max_prefactor * step_amount);


    for( int i = 0; i < T_max; i++){

        averaged_correlation.emplace_back(0.);
        averaged_correlation_squared.emplace_back(0.);
   
    }  


    //check if i want to save to more pair correlations as interim results than exist

    if(interim_pair_correlations_types.size() > database.return_amount_of_correlation_groups()){
        cout << "#### ERROR: MORE INTERIM PAIR CORRELATIONS THAN EXIST ####\n" << "Interim pair correlations = " << interim_pair_correlations_types.size() << " > " << "Existing pair correlations = " << database.return_amount_of_correlation_groups() << endl;
        cout << "Reducing amount of interim pair correlations to existing ones" << endl;
        interim_pair_correlations_types.resize(database.return_amount_of_correlation_groups());
        interim_auto_and_pair_correlation_results.resize(interim_pair_correlations_types.size());
        cout << "New amount of interim pair correlations = " << interim_pair_correlations_types.size() << endl << endl;
    }

    if(my_rank == 0){
        cout <<"\e[1mSTARTING:\e[0m \t" << filename << endl;
        cout << "\nA " << dimension <<"D cubic lattice with " << system_size << " lattice sites has been created"<< endl;
        cout << "Lattice constant is " << lattice_const << " (when processing the data I change this to 2.72325) " << endl;
        cout << "Iteration is done with " << system_amount << " initial conditions per core and "<< step_amount << " steps per initiial condition and a step width of " << step_width << endl;
        //cout << "Calculations will end at T = " << step_amount * step_width  << "µs T_max used for the correlation is therefore T_max = " << T_max_prefactor * step_amount * step_width << "µS" << endl; 
        //cout << "\"Larmor\" frequency approx " << 2. * mu_theo << " or " << gamma_i * B_field_direction.norm() << "10^7 rad / s"<< endl;
        cout << "Magnetic field is in direction [" << B_field_direction.real().transpose() << "]" << endl;
    }

};

// when useing the extrema we nee d different points to evaluate them for different magnetic fields
// this function assigns the correct time steps and their corrsponding indices

void find_time_steps_for_expectation_values(){

    Eigen::Vector3cd B_field = B_field_direction;

    vector<int> indices{};


    //indices for the extrema are choose due to prior calculations

    if(B_field == Eigen::Vector3cd(0,0,1)){
        indices = {0,120,211,297,366};

    }
    else if(B_field == Eigen::Vector3cd(0,1,1)){
        indices = {0,197,323,460};

    }
    else if(B_field == Eigen::Vector3cd(1,1, 1)){
        indices = {0,304,523,744};
  
    }
    else{
        cout << "#### ERROR: NO KNOWN MAGNETIC FIELD DIRECTION USED FOR EXPECTATION VALUES ####\n" << "B_field = " << B_field.real().transpose() << endl;
    }

    get<0>(expectation_values_at_certain_time) = indices;


}


//function that creates all the different systems and saves their data
// and if I say data i mean just the correlation

void create_averaged_data(bool calc_pair_corr){

    // Im iterationg over all the different systems. There Im already adding up all the different correlations to save time
    // IM also already normalising these to the sys1tem size. Might lead to errors and difficult to understand code



    for( int i = 0; i < system_amount; i++){

        if(i == system_amount / 2. - 1. / 2. * (system_amount % 2) && my_rank == 0 ){
            cout << "\n\tHALFWAY DONE WITH ITERATION FOR AVERAGING\n" << endl;
        }   

        unique_ptr<spin_system> lattice;

        //makue__uniwue_object is a custiom version of std::make_unique, which is not included in the gcc version on lido
        //therefore ineede a custom version, whis does not try to overwrite the std function
        if(dimension == 2){
            lattice = make_unique_object<spin_system_2D>(amount_per_row, step_amount, step_width, filename, false, database    );        }
        else{
            lattice = make_unique_object<spin_system>(amount_per_row, step_amount, step_width, filename, false, database    );        

        }
        lattice -> iterate_over_time();

        vector<double> single_system_correlation = lattice -> return_correlation();

        //calc pair correlations with the possibilitie to not do this for runtime reasons
        if(calc_pair_corr == true){
            vector< vector< double > > single_system_pair_correlation = lattice -> calc_and_return_pair_correlations();

            //filling the averaged pair correlation vector with 0. to make indexing possible
            //Done here, beacuase now we know the amount of correlation groups
            if(i == 0){
                int T_max = int( T_max_prefactor * step_amount);
                for( int j = 0; j < single_system_pair_correlation.size(); j++){
                    vector<double> zero_vec{};
                    for( int k = 0; k < T_max; k++){

                        zero_vec.emplace_back(0.);
                    }              
                    averaged_pair_correlations.emplace_back(zero_vec);
                }
            }

            //also add up all pair correlations
            //normalizing is also done in main due to MPI
            for(int j = 0; j < single_system_pair_correlation.size(); j++){
                for(int k = 0; k < single_system_pair_correlation[j].size(); k++){
                    averaged_pair_correlations[j][k] += single_system_pair_correlation[j][k];
                }
            }

        }



        if(averaged_correlation.size() != single_system_correlation.size() ){
            cout << "#### ERROR REGARDING AVERAGED CORRELATION: WRONG SIZE OF VECTORS ####\n" << "Correlation vector size = " << single_system_correlation.size() << " != " << "averaged correlation size = " << averaged_correlation.size() << endl;
            cout << "occured for the " << i << " system" <<endl;
        }


        //save quantities of interest || right now only done for correlations and not for pair correlations
        storage_correlation.emplace_back(single_system_correlation);


        //calc correlation values at certain time steps

        for(int j = 0; j < get<0>(expectation_values_at_certain_time).size(); j++){
            get<1>(expectation_values_at_certain_time)[j] += single_system_correlation[get<0>(expectation_values_at_certain_time)[j]];
            get<2>(expectation_values_at_certain_time)[j] += single_system_correlation[get<0>(expectation_values_at_certain_time)[j]] * single_system_correlation[get<0>(expectation_values_at_certain_time)[j]];
        }



        //Add all normalized correlations to average_correlation
        //normalizing is done in main due to MPI
        for(int j = 0; j < single_system_correlation.size(); j++){
            averaged_correlation[j] +=  single_system_correlation[j];
            //used for estimaton of MC Error
            averaged_correlation_squared[j] += single_system_correlation[j] * single_system_correlation[j];

        }
        calc_interim_correlation_results(i);


        if(calc_pair_corr == true){
            calc_interim_pair_correlation_results(i, calc_pair_corr);
        }


        //for(int j = 0; j < averaged_correlation.size(); j++){
        //    averaged_correlation[j] += 1. / system_amount * single_system_correlation[j];
        //}


    }




}

//function that calculates interim correlation results
// used for one lengthy calculation with many initial condiotions to check how well the results converge
// e.g for 10000 i.c. also saves 1000,2000,...
void calc_interim_correlation_results(int j){
    // no interim results for small amount of initial conditions treated
    // consider that this system amount is the number treated on one core
    if(system_amount < 5){
        return;
    }
    if((j + 1) % 5 == 0 && j != 0 && j != system_amount - 1){
        interim_correlation_results.emplace_back(averaged_correlation);
        interim_results_system_amount.emplace_back(j);
    }
}

void calc_interim_pair_correlation_results(int j, bool calc_pair_corr){
    // no interim results for small amount of initial conditions treated
    // consider that this system amount is the number treated on one core
    if(system_amount < 25){
        return;
    }

    if(calc_pair_corr == false){
        return;
    }



    if((j + 1) % 5 == 0 && j != 0 && j != system_amount - 1){
        if(interim_results_system_amount.size() == 0){
            cout << "#### ERROR: INTERIM RESULTS SYSTEM AMOUNT IS EMPTY ####\n"<< "This should have been calles before" << endl;
        }
        for( int i = 0; i < interim_pair_correlations_types.size(); i++){
            interim_auto_and_pair_correlation_results[i].emplace_back(averaged_pair_correlations[interim_pair_correlations_types[i]]);

        }
    }
    
}


//functions that returns correlation. Used for allreduce() in main.cpp

vector<double> return_averaged_correlation(){
    if(averaged_correlation.empty() == true){
        cout << "#####ERROR TRYING TO RETURN EMPTY CORRELATION#####" << endl;
    }
    
    return averaged_correlation;
}

//is used for Estimator of MC Error
vector<double> return_averaged_correlation_squared(){
    if(averaged_correlation_squared.empty() == true){
        cout << "#####ERROR TRYING TO RETURN EMPTY SQUARED CORRELATION#####" << endl;
    }
    
    return averaged_correlation_squared;
}

//function that returns interim correlation results
// used for one lengthy calculation with many initial condiotions to check how well the results converge
// e.g for 10000 i.c. also saves 1000,2000,...

pair< vector<int>, vector< vector< double > > > return_interim_correlation_results(){
    if(interim_correlation_results.size() != interim_results_system_amount.size()){
        cout << "#####ERROR INTERIM RESULTS AND ITS SYSTEM AMOUNT HAVE DIFFERENT DIMENSIONS#####" << endl;
        cout << "Interim correlation results size = " << interim_correlation_results.size() << " != " << "interim results system amount size = " << interim_results_system_amount.size() << endl;
    }
    pair< vector<int>, vector< vector< double > > > pair(interim_results_system_amount, interim_correlation_results);

    return pair;

}

// First element in tuple is vector with the different correlation types checked, second the system amounts used for the interim results
// third element is a vector of vectors containingt the correlations


tuple< vector<int>, vector<int>, vector< vector< vector< double > > > > return_interim_auto_and_pair_correlation_results(){
    if(interim_auto_and_pair_correlation_results[0].size() != interim_results_system_amount.size()){
        cout << "#####ERROR INTERIM RESULTS AND ITS SYSTEM AMOUNT HAVE DIFFERENT DIMENSIONS#####" << endl;
        cout << "Interim auto and pair correlation results size = " << interim_auto_and_pair_correlation_results[0].size() << " != " << "interim results system amount size = " << interim_results_system_amount.size() << endl;
    }
    if(interim_auto_and_pair_correlation_results.size() != interim_pair_correlations_types.size()){
        cout << "#####ERROR AMOUNT OF INTERIM RESULT GROUPS AND TYPES SUPPOSED TO CHECK HAVE DIFFERENT DIMENSIONS#####" << endl;
        cout << "Interim auto and pair correlation results size = " << interim_auto_and_pair_correlation_results.size() << " != " << "interim pair correlation types size = " << interim_pair_correlations_types.size() << endl;
    }
    return tuple< vector<int>, vector<int>, vector< vector< vector< double > > > > {interim_pair_correlations_types, interim_results_system_amount, interim_auto_and_pair_correlation_results};
}


//funtction that returns all the different averaged pair correlation types
vector <vector<double> > return_averaged_pair_correlations(){
    if(averaged_pair_correlations.size() == 0){
        cout << "#### ERROR REGARDING AVERAGED PAIR CORRELATION: NOTHING TO RETURN ####\n" << "Averaged pair correlation vector size = " << averaged_pair_correlations.size()<<endl;
    }
    return averaged_pair_correlations;
}

//function that returns the correlation values at a certain time object 

tuple< vector<int>, vector< double >, vector<double> > return_correlation_values_at_certain_time(){
    if(get<0>(expectation_values_at_certain_time).size() == 0){
        cout << "#### ERROR REGARDING CORRELATION VALUES AT CERTAIN TIME: NOTHING TO RETURN ####\n" << "Averaged pair correlation vector size = " << get<0>(expectation_values_at_certain_time).size()<<endl;
    }
    return expectation_values_at_certain_time;
}
//functions that save the averaged data
//useless dif used in MPI calculation




void save_averaged_correlation(){
    if(my_rank != 0){
        cout << "#####ERROR WRONG CORE TRYS TO WRITE DATA INTO .TXT#####" <<endl;
        return;
    }
    cout << "Data is being written into " << filename << endl; 	

    std::ofstream myfile;

    myfile.open ("build/daten/" + filename + "_averaged_Correlation.txt");
    myfile << "G_x time\n";

    for(int i = 0; i < averaged_correlation.size(); i++){
        myfile <<averaged_correlation[i]<< "\t"<<  i * step_width <<endl;
    }; 

    myfile.close();
    
}




};
