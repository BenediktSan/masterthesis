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

#include "../header/matrices_and_const.h"

//implementing thios here so i can use this on the clustr
//implemented before including multi syste averaging, otherwise this throws errors


template<typename T, typename... Args>
std::unique_ptr<T> make_unique_object(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

#include "../header/multi_system_averaging.h"





using cplx = std::complex<double>;
using namespace std;
using Observable = Eigen::MatrixXcd;
using uint = unsigned int;




//function that saves correlation after MPI calculations

void save_correlation_mpi(string filename, vector<double> averaged_correlation, vector<double> error, double step_width){
    
    std::ofstream myfile;

    //save 2D calculations to a different folder
    string postfix = "";

    if(filename.find("2D") != std::string::npos){
        postfix = "_2D";
    }

    myfile.open ("build/daten" + postfix +"/" + filename + "_averaged_Correlation.txt");
    if (myfile.is_open() == false){
    cout << "unable to open file " << filename <<endl;
    }
    myfile << "C_x err time\n";

    for(int i = 0; i < averaged_correlation.size(); i++){
        myfile << averaged_correlation[i]<< "\t"<< error[i]<< "\t"<< i * step_width <<endl;
    }; 


    myfile.close();  
    
}


//function that saves pair correlation after MPI calculations

void save_pair_correlations_mpi(string filename, vector<vector<double>> pair_correlations, double step_width){
    
    //save 2D calculations to a different folder
    string postfix = "";

    if(filename.find("2D") != std::string::npos){
        postfix = "_2D";
    }


    std::ofstream myfile;

    myfile.open ("build/daten" + postfix + "/" + filename + "_averaged_pair_Correlations.txt");
    if (myfile.is_open() == false){
    cout << "unable to open file " << filename <<endl;
    }
    string header = "";
    for( int i = 0; i < pair_correlations.size(); i++){
        header += "C_" + to_string(i) + " ";
    }
    header += "time";
    myfile << header +"\n";

    for(int i = 0; i < pair_correlations[0].size(); i++){
        
        for(int j = 0; j < pair_correlations.size(); j++){
            myfile << pair_correlations[j][i] << "\t";
        }

        myfile << i * step_width <<endl;
    }; 


    myfile.close();  
    
}

void save_interim_correlation_results_mpi(string filename, vector < int > interim_results_system_amount, vector< vector< double > > interim_correlation_results, double step_width){

    //save 2D calculations to a different folder
    string postfix = "";

    if(filename.find("2D") != std::string::npos){
        postfix = "_2D";
    }


    std::ofstream myfile;

    myfile.open ("build/daten" + postfix + "/" + filename + "_averaged_interim_Correlation_results.txt");
    if (myfile.is_open() == false){
        cout << "unable to open file " << filename <<endl;
    }
    string header = "system_amount ";
    for( int i = 0; i < interim_results_system_amount.size(); i++){
        header += to_string(interim_results_system_amount[i]) + "sys ";
    }
    header += "time";
    myfile << header +"\n";

    for(int i = 0; i < interim_correlation_results[0].size(); i++){

        if(i >= interim_results_system_amount.size()){
            myfile <<0 << "\t";
        }
        else{
            myfile << interim_results_system_amount[i] << "\t";
        }
        
        for(int j = 0; j < interim_correlation_results.size(); j++){
            myfile << interim_correlation_results[j][i] << "\t";
        }

        myfile << i * step_width <<endl;
    }; 


    myfile.close();  


}

void save_interim_auto_and_pair_correlation_results_mpi(string filename,vector< int > interim_pair_correlations_types, vector < int > interim_results_system_amount, vector< vector< vector< double > > > interim_auto_and_pair_correlation_results, double step_width){

    //save 2D calculations to a different folder
    string postfix = "";

    if(filename.find("2D") != std::string::npos){
        postfix = "_2D";
    }

    if(interim_auto_and_pair_correlation_results.size() == 0){
        cout << "#### ERROR: NO INTERIM AUTO AND PAIR CORRELATION RESULTS TO SAVE ####\n" << "Nothing to save" << endl;
        return;
    }


    std::ofstream myfile;

    for( int corr_type_index = 0; corr_type_index < interim_pair_correlations_types.size(); corr_type_index++){
        int corr_type = interim_pair_correlations_types[corr_type_index];

        if(interim_auto_and_pair_correlation_results[corr_type_index].size() == 0){
            cout << "#### ERROR: NO INTERIM AUTO AND PAIR CORRELATION RESULTS TO SAVE FOR CORRELATION TYPE C = " << corr_type << " ####\n" << "Nothing to save" << endl;
            continue;
        }

        myfile.open ("build/daten" + postfix + "/" + filename + "_averaged_interim_pair_Correlation_results_type_" + to_string(corr_type) + ".txt");


        if (myfile.is_open() == false){
            cout << "unable to open file " << filename <<endl;
        }
        vector< vector< double > > interim_pair_correlation_results = interim_auto_and_pair_correlation_results[corr_type_index];

        myfile << "Correlation type C = " << corr_type << endl;

        string header = "system_amount ";
        for( int i = 0; i < interim_results_system_amount.size(); i++){
            header += to_string(interim_results_system_amount[i]) + "sys ";
        }
        header += "time";
        myfile << header +"\n";

        for(int i = 0; i < interim_pair_correlation_results[0].size(); i++){

            if(i >= interim_results_system_amount.size()){
                myfile <<0 << "\t";
            }
            else{
                myfile << interim_results_system_amount[i] << "\t";
            }

            for(int j = 0; j < interim_pair_correlation_results.size(); j++){
                myfile << interim_pair_correlation_results[j][i] << "\t";
            }

            myfile << i * step_width <<endl;
        }; 


        myfile.close();  

    }
}


void save_variance_for_certain_step(string filename,vector<int> evaluation_indices, vector<double> exp_val, vector<double> exp_val_squared, double step_width){
    
    //save 2D calculations to a different folder
    string postfix = "";

    if(filename.find("2D") != std::string::npos){
        postfix = "_2D";
    }


    std::ofstream myfile;

    myfile.open ("build/daten" + postfix + "/" + filename + "_variance_for_step.txt");
    if (myfile.is_open() == false){
        cout << "unable to open file " << filename <<endl;
    }
    myfile << "variance, exp_val, exp_val_squared, time\n";

    for(int i = 0; i < exp_val.size(); i++){
        double variance = exp_val_squared[i] - exp_val[i] * exp_val[i];
        myfile << variance << "\t"<< exp_val[i] << "\t" << exp_val_squared[i] << "\t" << evaluation_indices[i] * step_width <<endl;
    }
    myfile.close();  

}

//Function that evaluates the calculated values for the Correlation using the magnetization

void evaluate_correlations_mpi(multi_system_averaging average, string filename, double dt, uint system_amount, int world_size, int my_rank){

    vector<double> correlation = average.return_averaged_correlation(); 

    //used to estimate Error
    vector<double> correlation_squared = average.return_averaged_correlation_squared(); 
    vector<double> estimated_error{}; 



    // ====== add up vectors on all cores ======
    //start with simple correlation calculated from magnetization
    if(world_size < 15){
        cout <<"Core Nr "<< my_rank << " waiting for other cores"<< endl;
    }    
    MPI_Barrier(MPI_COMM_WORLD);

    vector<double> buffer( correlation.size() );
    MPI_Allreduce( correlation.data(), buffer.data(), correlation.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    correlation = buffer;

    MPI_Allreduce( correlation_squared.data(), buffer.data(), correlation_squared.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    correlation_squared = buffer;


    //normalize correlation to system size
    for(int j = 0; j < correlation.size(); j++){
        correlation[j] = 1. /(world_size * system_amount) * correlation[j];

        correlation_squared[j]= 1. /pow(world_size * system_amount, 1.) * correlation_squared[j];
        double sigma_squared = 1. / (system_amount * world_size) * ( correlation_squared[j] - (correlation[j] * correlation[j]));
        //not entirely sure about sqrt
        estimated_error.emplace_back( sqrt(sigma_squared));
    }

 
    


    // ====== inform user about result ======
    if( my_rank == 0 ){
        std::cout << "evaluating correlation finished" <<endl; 
        save_correlation_mpi(filename, correlation, estimated_error, dt);
    }




}

//Function that evaluates the calculated values for pair correlations and saves them

void evaluate_pair_correlations_mpi(multi_system_averaging average, string filename, double dt, uint system_amount, int world_size, int my_rank, bool calc_pair_corr){

    //return if no pair correlations were calculated
    if(calc_pair_corr == false){
        return;
    }


    //note: pair correlations are not normalized
    vector <vector<double> > pair_correlations = average.return_averaged_pair_correlations();



    // ====== add up vectors on all cores ======
    if(world_size < 15){
        cout <<"Core Nr "<< my_rank << " waiting for other cores"<< endl;
    }        
    MPI_Barrier(MPI_COMM_WORLD);

    //calculate pair correlations

    vector<double> buffer( pair_correlations[0].size() );

    for( int i = 0; i < pair_correlations.size(); i++){
        MPI_Allreduce( pair_correlations[i].data(), buffer.data(), pair_correlations[i].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        pair_correlations[i] = buffer;
    }

    

    //normalize pair correlation to sytsem size
    for(int i = 0; i < pair_correlations.size(); i++){
        for(int j = 0; j < pair_correlations[i].size(); j++){

            pair_correlations[i][j] = 1. /(world_size * system_amount) * pair_correlations[i][j];

        }
    }


    // ====== inform user about result ======
    if( my_rank == 0 ){
        std::cout << "evaluating pair correlation finished" <<endl; 
        save_pair_correlations_mpi(filename, pair_correlations, dt);
    }

}




//oevaluate interim correlation results
// used for one lengthy calculation with many initial condiotions to check how well the results converge
// e.g for 10000 i.c. also saves 1000,2000,...
void evaluate_interim_correlation_results_mpi(multi_system_averaging average, string filename, double dt, uint system_amount, int world_size, int my_rank){
    
    pair< vector<int>, vector< vector< double > > > pair = average.return_interim_correlation_results(); 


    vector< vector< double > > interim_correlation_results = pair.second;
    //Inforrmation how many systems were used for a certain interim result
    vector<int> interim_results_system_amount = pair.first;

    if(interim_correlation_results.size() == 0){
        if (my_rank == 0){
            cout <<"No interim results to evaluate" <<endl;
        }
        return;
    }


    // ====== add up vectors on all cores ======
    //start with simple correlation calculated from magnetization
    if(world_size < 15){
        cout <<"Core Nr "<< my_rank << " waiting for other cores"<< endl;
    }    
    MPI_Barrier(MPI_COMM_WORLD);



    for(int i = 0; i< interim_correlation_results.size(); i++){
        vector<double> buffer( interim_correlation_results[i].size() );
        MPI_Allreduce( interim_correlation_results[i].data(), buffer.data(), interim_correlation_results[i].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        interim_correlation_results[i] = buffer;

    }


    vector<int> interim_results_system_amount_buffer = interim_results_system_amount;

    //normalize correlation to system size
    for(int i = 0; i< interim_correlation_results.size(); i++){
        int current_system_amount = interim_results_system_amount[i] * world_size;
        interim_results_system_amount_buffer[i] = current_system_amount;
        for(int j = 0; j < interim_correlation_results[i].size(); j++){
            interim_correlation_results[i][j] = 1. /(current_system_amount) * interim_correlation_results[i][j];
        }
    }
 
    


    // ====== inform user about result ======
    if( my_rank == 0 ){
        std::cout << "evaluating correlation finished" <<endl; 
        save_interim_correlation_results_mpi(filename, interim_results_system_amount_buffer, interim_correlation_results, dt);
    }

}

//oevaluate interim pair correlation results
// used for one lengthy calculation with many initial condiotions to check how well the results converge
// e.g for 10000 i.c. also saves 1000,2000,...
void evaluate_interim_auto_and_pair_correlation_results_mpi(multi_system_averaging average, string filename, double dt, uint system_amount, int world_size, int my_rank){
    
    tuple< vector<int>, vector<int>, vector <vector< vector< double > > > > data_tuple = average.return_interim_auto_and_pair_correlation_results(); 

    
    vector< vector< vector< double > > > all_interim_auto_and_pair_correlation_results = get<2>(data_tuple);
    //Inforrmation how many systems were used for a certain interim result
    vector<int> interim_results_system_amount = get<1>(data_tuple);
    //pair correlation types that were considered
    vector<int> interim_pair_correlations_types = get<0>(data_tuple);



    if(all_interim_auto_and_pair_correlation_results.size() == 0){
        if (my_rank == 0){
            cout <<"No interim results to evaluate" <<endl;
        }
        return;
    }

    //calculate correct system amount fo all cores used
    
    for(int i = 0; i< interim_results_system_amount.size(); i++){
        int current_system_amount = interim_results_system_amount[i] * world_size;
        interim_results_system_amount[i] = current_system_amount;
    }



    // ====== add up vectors on all cores ======
    //start with simple correlation calculated from magnetization
    if(world_size < 15){
        cout <<"Core Nr "<< my_rank << " waiting for other cores"<< endl;
    }    
    MPI_Barrier(MPI_COMM_WORLD);

    for(int corr_type_index = 0; corr_type_index < interim_pair_correlations_types.size(); corr_type_index++){
        int corr_type = interim_pair_correlations_types[corr_type_index];
        vector< vector< double > > interim_auto_and_pair_correlation_results = all_interim_auto_and_pair_correlation_results[corr_type_index];

        if(interim_auto_and_pair_correlation_results.size() == 0){
            if (my_rank == 0){
                cout <<"No interim results to evaluate for correlation type C = " << corr_type << endl;
            }
            continue;
        }


        for(int i = 0; i< interim_auto_and_pair_correlation_results.size(); i++){
            vector<double> buffer( interim_auto_and_pair_correlation_results[i].size() );
            MPI_Allreduce( interim_auto_and_pair_correlation_results[i].data(), buffer.data(), interim_auto_and_pair_correlation_results[i].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            interim_auto_and_pair_correlation_results[i] = buffer;

        }




        //normalize correlation to system size
        for(int i = 0; i< interim_auto_and_pair_correlation_results.size(); i++){
            int current_system_amount = interim_results_system_amount[i] ;
            for(int j = 0; j < interim_auto_and_pair_correlation_results[i].size(); j++){
                interim_auto_and_pair_correlation_results[i][j] = 1. /(current_system_amount) * interim_auto_and_pair_correlation_results[i][j];
            }
        }

        all_interim_auto_and_pair_correlation_results[corr_type_index] = interim_auto_and_pair_correlation_results;
    
    }


    // ====== inform user about result ======
    if( my_rank == 0 ){
        std::cout << "evaluating correlation finished" <<endl; 
        save_interim_auto_and_pair_correlation_results_mpi(filename, interim_pair_correlations_types, interim_results_system_amount, all_interim_auto_and_pair_correlation_results, dt);
    }
    
}



//this function evaluates the expecation values of my correlation at certain time steps. this is used to calculate the variance

void evaluate_variance_at_certain_time(multi_system_averaging average, string filename, double dt, uint system_amount, int world_size, int my_rank){

    //get the correlation values at certain time steps
    tuple< vector<int>, vector< double >, vector<double> > expectation_values_at_certain_time = average.return_correlation_values_at_certain_time(); 

    vector< int > indices = get<0>(expectation_values_at_certain_time);
    vector< double > exp_val = get<1>(expectation_values_at_certain_time);
    vector< double > exp_val_squared = get<2>(expectation_values_at_certain_time);


    // ====== add up vectors on all cores ======
    //start with simple correlation calculated from magnetization
    if(world_size < 15){
        cout <<"Core Nr "<< my_rank << " waiting for other cores"<< endl;
    }    
    MPI_Barrier(MPI_COMM_WORLD);

    


    vector<double> buffer( exp_val.size() );

    MPI_Allreduce( exp_val.data(), buffer.data(), exp_val.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    exp_val = buffer;

    MPI_Allreduce( exp_val_squared.data(), buffer.data(), exp_val_squared.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    exp_val_squared = buffer;


    //normalize correlation to system size

    for(int j = 0; j < exp_val.size(); j++){
        exp_val[j] = 1. /(world_size * system_amount) * exp_val[j];

        exp_val_squared[j]= 1. /(world_size * system_amount) * exp_val_squared[j];
    }

    // ====== inform user about result ======
    if( my_rank == 0 ){
        std::cout << "evaluating variance finished" <<endl; 
        save_variance_for_certain_step(filename, indices, exp_val, exp_val_squared, dt);
    }

}



//this is a function that transposes a std::vector<std::vector<double>>, effectively inverting its indices
//used in calculation of pair correlations

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty()) return {}; // Handle empty input

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j]; // Swap indices
        }
    }

    return transposed;
}



//functions that only prints something to the terminal if its root (0th core)

void root_core_cout(int my_rank){

    
}


//function that builds folder


void build_folder(){


    // ,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH
    std::string foldername1 = "build";
    char const *c = foldername1.data();
    if( mkdir( c,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH  ) == -1 ) // creating the folder failed ( because it already exists? )
    {
        std::cout << "error_build\n";
    }

    std::string foldername2 = "build/daten";
    char const *d = foldername2.data();
    if( mkdir( d,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1 ) // creating the folder failed ( because it already exists? )
    {
        std::cout << "error_build/daten\n";
    }

    std::string foldername3 = "build/daten_2D";
    char const *e = foldername3.data();
    if( mkdir( e,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1 ) // creating the folder failed ( because it already exists? )
    {
        std::cout << "error_build/daten_2D\n";
    } 

    return;
}

void vec_to_txt(string file_name, vector<cplx> vec){

    std::ofstream myfile;
    myfile.open ("build/daten/" + file_name + ".txt");
    
    myfile << "\n\t Eigenvalues:\n";

    for(int i = 0; i < vec.size(); i++){
        myfile << vec[i] <<endl;
    }; 
    

    myfile.close();
};



uint throw_seed( const std::string& seed ){
    if( seed.find("random") != std::string::npos ) // random seed
    {
        return static_cast<uint>( time(nullptr) );
    }  
    else // preset seed
    {
        return std::stoi( seed );
    }  
}



Eigen::VectorXcd mat_to_vec(Observable mat){


    Eigen::VectorXcd vec(mat.size());
    vec << (mat.transpose()).reshaped();

    return vec;
}

Observable vec_to_mat(Eigen::VectorXcd vec){

    uint len = sqrt(vec.size());
    Observable mat(len, len);

    mat << vec.reshaped(len,len).transpose();

    return mat;
}

void check_mat_conversion(){
    
    Observable A = Eigen::MatrixXd::Random(3,3);
    cout << "\n\n" << A << endl;
    cout << "\n\n"<<endl; 
    cout <<mat_to_vec(A) <<endl;;
    cout << "\n\n" << endl;
    cout << vec_to_mat(mat_to_vec(A))<<endl;

}

//Vektor mit EInheitsmatrizen für Erweiterung der Hilberträume. Spinmatrizen werden dabei über ersetzen eingefügt
//Reihenfolge: Centralspin (&Trion), Bathspins,

vector<Observable> build_identity_vec_S1(uint bath_size){
    vector<Observable> ident{};
    ident.emplace_back(SIGMA_0);

    for( uint i = 0; i < bath_size; i++){
        ident.emplace_back(SIGMA_0_S1);
    }
    return ident;
}

//dies ist der algemeinere build identity
vector<Observable> build_identity_vec(uint bath_size, float I_bath){

    Observable SIGMA_0_bath = Observable::Identity(2. * I_bath +1, 2. * I_bath +1);

    vector<Observable> ident{};
    ident.emplace_back(SIGMA_0);

    for( uint i = 0; i < bath_size; i++){
        ident.emplace_back(SIGMA_0_bath);
    }
    return ident;
}

Observable commutator(Observable A, Observable B){
    return (A * B) - (B * A);
}

vector<Observable> build_spin_vec(float I_bath){
    vector<Observable> S_vec{};
    auto S_plus = build_local_Splus(I_bath);
    auto S_minus = S_plus.adjoint();

    const Observable S_X = 1./2. * ( S_plus + S_minus);
    S_vec.emplace_back(S_X);

    const Observable S_Y = 1./(2. * cplx_i) * ( S_plus - S_minus);
    S_vec.emplace_back(S_Y);

    const Observable S_Z = 1./2. * commutator(S_plus,S_minus);
    S_vec.emplace_back(S_Z);

    return S_vec;
}


double trace_exp_val(Observable A, Observable B){
    cplx trace = (A * B).trace();
    if(imag(trace) > pow(10., -10.)){
        cout <<"\nERROR\tERWARTUNGSWERT KOMPLEX\tTERROR\n" << trace << endl;
    }
    return real(trace);
}

//Ganz ganz viele funktionen um vektoren zu printen
void print_vec(vector<Observable> vec){
    int size = vec.size();
    for(int i = 0; i<size; i++){
        cout << vec[i] << endl << endl;
    }
    cout << endl;
}

void print_vec(vector<Eigen::VectorXcd> vec){
    int size = vec.size();
    for(int i = 0; i<size; i++){
        cout << vec[i] << endl << endl;
    }
    cout << endl;
}

void print_vec(vector<double> vec){
    int size = vec.size();
    for(int i = 0; i<size; i++){
        cout << vec[i] << endl;
    }
    cout << endl;
}

void print_vec(vector<cplx> vec){
    int size = vec.size();
    for(int i = 0; i<size; i++){
        cout << vec[i] << endl;
    }
    cout << endl;
}


Observable build_local_Splus(const float j){
   uint H_dim = 2*j +1;         //initialisator
   float m = j;               //counter
   Observable S_plus = Observable::Zero(H_dim,H_dim);
   for(uint i= 0; i<H_dim-1; i++){
        m = m -1 ; 
        S_plus(i,i+1) = sqrt(j*(j+1) -m*(m+1));
   };
   return S_plus;
}

