
#include <iostream>
#include "../../header/creating_spin_system/semi_classical_spins.h"

//################ Differential equations for semi classical approach



// Bfield in matrices and const
//coupling const is calculated in latice_database


Eigen::Vector3cd spin_system::EoM_Dipolar_component(uint index, vector<Eigen::Vector3cd> &spin_vec){

    Eigen::Vector3cd dS(0,0,0);

        
    for(uint j = 0; j < spin_vec.size(); j++){

        double d =  coupling_const[index][j];

        //special cas to account for cutoff to save computational time
        //cut off is realized by setteing the dipolar coupling constant to 0
        if(d == 0.){
            continue;
        }


        if(index == j){
            continue;
        }
        //here no factor 1/2 due to douple sum in Hamilton
        dS(0) += d * ( - spin_vec[j](1) * spin_vec[index](2) - 2. * spin_vec[j](2) * spin_vec[index](1) );
        dS(1) += d * ( 2. * spin_vec[j](2) * spin_vec[index](0) + spin_vec[j](0) * spin_vec[index](2) );
        dS(2) += d * ( - spin_vec[j](0) * spin_vec[index](1) + spin_vec[j](1) * spin_vec[index](0) );

        //cout << "dS components  " << dS(0) << " " << dS(1) << " " << dS(2) <<endl;
    }
    
    //dS = spin_vec[index].cross(B_field_direction);
    return dS;
}



void spin_system::EoM_Dipolar(vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, const double t){
    
    //cout << "Check which size dxdt has. Do i have to create some new stuff?\t"<< dxdt.size() << endl;

    for( uint i = 0; i < spin_vec.size(); i++){
        dxdt[i] = EoM_Dipolar_component(i, spin_vec);

    }
}



//Wenn diese Rechnng funktioniert skaliert sie nur noch mit (n^2 + 3n)/2 statt n^2 ( sogar noch weniger wegend keiner selbstww, bin aber gerade)
// this equation saves time by using loops in a more effective was, but ist loses time due to calculations for each lattice site being done at the same moment

void spin_system::EoM_Dipolar_new(vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt, const double t){

    Eigen::Vector3cd dS(0,0,0);


    for( uint index = 0; index < spin_vec.size(); index++){
        //dxdt[index] = EoM_Dipolar_component(index, spin_vec);

        

        //here only j smaller than index due to d_ij = 0
        for(uint j = 0; j < index; j++){

            double d = coupling_const[index][j];


            if(index == j){
                continue;
            }
            
            //trying to improve this by using alternative way of writting these EoM

            Eigen::Vector3cd cross_prod = spin_vec[index].cross(spin_vec[j]);
            Eigen::Vector3cd additional_vector_1(- 3. * spin_vec[index](1) * spin_vec[j](2), + 3. * spin_vec[index](0) * spin_vec[j](2), 0.);
            Eigen::Vector3cd additional_vector_2(- 3. * spin_vec[j](1) * spin_vec[index](2), + 3. * spin_vec[j](0) * spin_vec[index](2), 0.);



            dxdt[index] += d/2. * (cross_prod + additional_vector_1);

            dxdt[j] += d/2. * (-cross_prod + additional_vector_2);


            //cout << "dS components  " << dS(0) << " " << dS(1) << " " << dS(2) <<endl;
        }
    }
}

// Equations used for a Heisenebrg system. Mostly used for 2D systems

    void spin_system::EoM_Heisenberg(vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt , const double t ){

        for(int i = 0; i < system_size; i++){
            

            // more efficient way to sum over all interactions
            for(int j =0; j < i ; j++){
                //coupling constant fr this lattice pair
                double J = coupling_const_heisenberg[i][j];
                if(coupling_const_heisenberg[i][j] != coupling_const_heisenberg[j][i]){
                    cout << "############## ERROR: COUPLING CONSTANTS FOR HEISENBERGMODEL DIFFER ############\n" << coupling_const_heisenberg[i][j] << " "  << coupling_const_heisenberg[i][j] << endl;
                }

                if(i==j  || abs(J) <  pow(10,-10) ){
                    continue;
                }
                dxdt[i] +=  - J *  spin_vec[j].cross(spin_vec[i]) ;
                dxdt[j] +=    J * spin_vec[j].cross(spin_vec[i]);


                //these EoMs have a factor 2 anisotropy in the z component
                
                //dxdt[i][0] += -J * ( spin_vec[j][1] * spin_vec[i][2] - 2. * spin_vec[j][2] * spin_vec[i][1] ) ;
                //dxdt[i][1] += -J * ( 2. * spin_vec[j][2] * spin_vec[i][0] - spin_vec[j][0] * spin_vec[i][2] ) ;
                //dxdt[i][2] += -J * ( spin_vec[j][0] * spin_vec[i][1] - spin_vec[j][1] * spin_vec[i][0] ) ;

            }
        //dxdt[i] +=  spin_vec[i].cross(B_field_direction);
        //cout << dxdt[i] << endl;

        }            
    //cout << dxdt[0] <<endl;
    //cout << "#########################################" << endl;

    //cout << dxdt[0] << endl;

    }
    
    void spin_system::EoM_Heisenberg_XXZ(vector<Eigen::Vector3cd> &spin_vec, vector<Eigen::Vector3cd> &dxdt , const double t ){

        for( uint index = 0; index < spin_vec.size(); index++){

            Eigen::Vector3cd dS(0,0,0);


            for(uint j = 0; j < spin_vec.size(); j++){

                double J =  coupling_const_heisenberg[index][j];

                //special cas to account for cutoff to save computational time
                //cut off is realized by setteing the dipolar coupling constant to 0


                if(index==j  || abs(J) <  pow(10,-10) ){
                    continue;
                }

                //here no factor 1/2 due to douple sum in Hamilton
                dS(0) += J * ( - spin_vec[j](1) * spin_vec[index](2) - 2. * spin_vec[j](2) * spin_vec[index](1) );
                dS(1) += J * ( 2. * spin_vec[j](2) * spin_vec[index](0) + spin_vec[j](0) * spin_vec[index](2) );
                dS(2) += J * ( - spin_vec[j](0) * spin_vec[index](1) + spin_vec[j](1) * spin_vec[index](0) );

                //cout << "dS components  " << dS(0) << " " << dS(1) << " " << dS(2) <<endl;
            }

            dxdt[index] = dS;
        }

        

    }

    