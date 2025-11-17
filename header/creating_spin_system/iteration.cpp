
#include <iostream>
#include "../../header/creating_spin_system/semi_classical_spins.h"
//#include "../../header/creating_spin_system/calc_quantities.cpp"

// ✅ Define how Boost should handle Eigen::Vector3cd

//it is important to note that also a way to handle a vector of an Eigen Vector is handled needs to be implemented
namespace boost { namespace numeric { namespace odeint {
    template<>
    struct vector_space_norm_inf<Eigen::Vector3cd> {
        static double apply(const Eigen::Vector3cd &x) {
            return x.template lpNorm<Eigen::Infinity>();  // Max norm
        }
    };
}}}



/////////// Functions for the iterationF




void spin_system::iterate_over_time(){

    if(database.return_database_dimension() != dimension){

        cout << "#### ERROR: DIMENSIONS IN DATABASE AND SPIN SYSTEM DO NOT MATCH ####\n" << "Database dimension is = " << database.return_database_dimension() << " Spin system dimension is = " << dimension << endl;

    }


    for(uint i = 0; i < step_amount; i++){

        if(i == step_amount / 2. - 2. * (step_amount % 2) && print_text == true){
            cout << "\n\t\tHALFWAY DONE WITH SINGLE SYSTEM ITERATION\n" << endl;
        }

        //default value for EoM is the diploar EoM
        //if dim = 2 we change to Heisenberg hamiltonian
        function<void(std::vector<Eigen::Vector3cd>&, std::vector<Eigen::Vector3cd>&, double)> wrapped_EoM;

        wrapped_EoM = [this]( vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, double t) {
            this -> EoM_Dipolar(spin_vec, dxdt, t);
        };

        if(dimension == 2){
            if(use_xxz_hamiltonian_2D == true){
                wrapped_EoM = [this]( vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, double t) {
                this -> EoM_Heisenberg_XXZ(spin_vec, dxdt, t);
                };
            }
            else{
                wrapped_EoM = [this]( vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, double t) {
                    this -> EoM_Heisenberg(spin_vec, dxdt, t);
                };
            }
        }

        //using namespace std::placeholders; // for _1, _2, _3

        //auto bound_func = bind(&spin_system::EoM_Dipolar, this, placeholders::_1, placeholders::_2, placeholders::_3);

        iteration_step(wrapped_EoM, i * step_width);
        calc_magnetization();
        //cout << "magnetization = " <<  magnetization.back().transpose() <<endl;
        calc_magnetization_central_spin();
        //calc_magnetization_correlation();
        for(int j = 0; j < spin_vec.size(); j++){
            if(abs(spin_vec[j].norm() - sqrt(3.) / 2.) > 5 * pow(10., -8.)){
                cout << "#### ERROR REGARDING SPIN LENGTH AT STEP "<< j <<" ####\n" << "check spin length " << spin_vec[i].norm() *pow(10,7) << "e-7" << " != " << sqrt(3.) / 2. * pow(10,7) << "e-7" << "?\n"<< endl;
            }
        }
    }
    for(int i = 0; i < spin_vec.size(); i++){
        if(abs(spin_vec[i].norm() - sqrt(3.) / 2.) > 1 * pow(10., -8.)){
            cout << "#### ERROR REGARDING SPIN LENGTH ####\n" << "check spin length " << spin_vec[i].norm() *pow(10,7) << "e-7" << " != " << sqrt(3.) / 2. * pow(10,7) << "e-7" << "?\n"<< endl;
        }
    }
   
    // calls calculation for pair calc in return func to get the possibility in multi_system_averaging to NOT calc pair correlations
    calc_all_magnetization_correlation();


    ////function to check if data is supposed to be written to a .txt
    //if(print_text ==true){
    //    save_magnetization();
    //    save_magnetization_central_spin();
    //    save_correlation();
    //}


}



void spin_system::iteration_step(function<void(vector<Eigen::Vector3cd>&, vector<Eigen::Vector3cd>&, double)> EoM, const double time_step){

    

    // Wrap the differential equation with a lambda to pass the additional argument
    //Here the additional argument are the distance vectors and the distance
    //auto wrapped_EoM = [this]( vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, double t) {
    //    this -> EoM(spin_vec, dxdt, t);
    //};


    //could be moved outside of function to save time
    runge_kutta4< vector<Eigen::Vector3cd> > rk4;
    rk4.do_step(EoM, spin_vec, time_step, step_width);
    //cout << spin_vec[0].transpose() << endl;
    //modified_midpoint< vector<Eigen::Vector3cd> > step;
    //step.do_step(wrapped_EoM, spin_vec, time_step, step_width);
 

    //assign new system state to vector for pair correlation calculations later
    all_spin_vec_configurations.emplace_back(spin_vec);
    


}

/*

//implement iteration using a stepper with step width adjusrment

void spin_system::iterate_over_time_self_adjusting(){

    double t = 0;
    double t_end = step_amount * step_width;
    while(t < t_end){
        double t_next = t + step_width;
        if (t_next > t_end) t_next = t_end;

        if(t == (step_amount * step_width) / 2. - 2. * (step_amount % 2) && print_text == true){
            cout << "\n\t\tHALFWAY DONE WITH SINGLE SYSTEM ITERATION\n" << endl;
        }

        iteration_step_self_adjusting( t);
        calc_magnetization();
        calc_magnetization_central_spin();
        //calc_magnetization_correlation();

    }

    if(abs(spin_vec[0].norm() - sqrt(3.) / 4.) > 5 * pow(10., -8.)){
    cout << "#### ERROR REGARDING SPIN LENGTH ####\n" << "check spin length " << spin_vec[0].norm() *1000000 << "e-6" << " != " << sqrt(3.) / 4. * 1000000 << "e-6" << "?\n"<< endl;
    }
    calc_all_magnetization_correlation();


    //function to check if data is supposed to be written to a .txt
    if(print_text ==true){
        save_magnetization();
        save_magnetization_central_spin();
        save_correlation();
    }


}



//new rk4 stepper with step adjustment
void spin_system::iteration_step_self_adjusting(double time_step){


    // Wrap the differential equation with a lambda to pass the additional argument
    //Here the additional argument are the distance vectors and the distance
    auto wrapped_EoM = [this]( vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, double t) {
        this -> EoM_Dipolar(spin_vec, dxdt, t);
    };


    //could be moved outside of function to save time
    typedef runge_kutta_cash_karp54< vector<Eigen::Vector3cd> > rk4;
    auto controlled_stepper = make_controlled(1E-6, 1E-6, rk4());
    controlled_stepper.try_step(wrapped_EoM, spin_vec, time_step, step_width);



}


*/

//Special case for Zeeman iteration
//only allowed if no other calculations are done in this class instance

void spin_system::do_Zeeman_iteration(){

    if(magnetization.size() != 0){
        cout << "\nERROR: THIS CLASS INSTANCE HAS BEEN USED FOR OTHER CALCULATIONS AND NOT JUST FOR ZEEMAN\n" << endl;
        return;
    }

    runge_kutta4< vector<Eigen::Vector3cd> > rk4;

    for( int time_step = 0; time_step < step_amount; time_step ++){
        rk4.do_step(EoM_Zeeman, spin_vec, time_step, step_width);
        calc_magnetization();
        calc_magnetization_central_spin();
    }


    if(abs(spin_vec[0].norm() - sqrt(3.) / 4.) > pow(10., -10.)){
        cout << "#### ERROR REGARDING SPIN LENGTH IN ZEEMAN####\n" << "check spin length " << spin_vec[0].norm() *1000 << "e-3" << " != " << sqrt(3.) / 4. * 1000 << "e-3" << "?\n"<< endl;
    }
    calc_all_magnetization_correlation();




}



