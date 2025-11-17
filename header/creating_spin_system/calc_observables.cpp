
#include <iostream>
#include "../../header/creating_spin_system/semi_classical_spins.h"







//////// Functions important for the calculations
 

//function which calculates magnetization or correlations when called and appends value to a vector

void spin_system::calc_magnetization(){
    
    Eigen::Vector3cd M;
    M << 0, 0, 0;
    double particle_amount = system_size;

    for(int i = 0; i < particle_amount; i++){
        M += spin_vec[i];
    }
    //maybe i caclulate the magnetization wrong here. Maybe i need a different prefactor
    //mybe use factor sqrt(N)
    M = M * 1. / sqrt(particle_amount);
    magnetization.emplace_back(M);
}

void spin_system::calc_magnetization_central_spin(){
    
    Eigen::Vector3cd M_CS;


    // save spin at (ß0,0,0) as central spin
    M_CS = spin_vec[0];
    magnetization_CS.emplace_back(M_CS);
}





//calculates magnetization for all time steps at the same time
//this is actually used

void spin_system::calc_all_magnetization_correlation(){

    if(magnetization.size() != step_amount){
        cout << "\nERROR IN CALCULATION OF CORRELATION\nmagnetization() has wrong size " << magnetization.size() << " while step amount is " << step_amount << endl << endl;
    }

    
    //correlation variable
    double G = 0;
    //could lead to errors that i use abs() here
    double G_0 = 0;//abs(magnetization[0].dot(magnetization[0]));

    //this is only for
    //double G1 = 0;
    //double G2 = 0;
    double G_test = 0;

    int T_max = int( T_max_prefactor * step_amount);

    //t is the time dependencie for the correlation function
    // This is iteated until T_max is reached, beacuse this is also the maximum time step, which can be evaluated via time averaging
    for( int t = 0; t < T_max; t++){
        //correlation variable
        G = 0;

        //this is the integration over tau
        for(int tau = 0; tau< T_max; tau++){

            G += real(magnetization[tau](0) * magnetization[t + tau](0));
            //G1 += real(magnetization[tau](1) * magnetization[t + tau](1));
            //G2 += real(magnetization[tau](2) * magnetization[t + tau](2));

            if(tau == 0 && t == 0){
                G_test = G;
                //cout << "\033[1;36m" <<"FID at t=0 w/o integration= " << G_test<<  "\033[0m" << endl;
            }
            
            
        }  
        
        if(t == 0){
            //cout << "\033[1;36m"  << "integrated magnetization FID = " << G / T_max <<  "\033[0m" << endl;
            //cout << "\033[1;36m"  << "integrated magnetization FID = " << G1 / T_max <<  "\033[0m" << endl;
            //cout << "\033[1;36m"  << "integrated magnetization FID = " << G2 / T_max <<  "\033[0m" << endl;
            //cout << "\033[1;36m"  << "Sum FID = " << G / T_max + G1 / T_max + G2 / T_max <<  "\033[0m" << endl;
            //cout << "\033[1;36m"  << "div FID = " << G_test / (G / T_max) <<  "\033[0m" << endl;

        }
        

        G = 1./ T_max * G ;    
        //case to calculate G_0
        if(t == 0){
            G_0 = G;
        }          
        //#####this normalization is a crucila stepp which need investigation
        //G = G / G_0;
        magnetization_correlation.emplace_back(G);

    }
    
}

//function that calculates my pairr correltaion for a vector of index pairs
vector< double > spin_system::calc_pair_correlation_for_index_pairs(const vector< pair< uint, uint> > &all_index_pairs){

    //can have errors for the last group, for uneven amounts of lattice sites per row

    //all configurations of spins for that correlation group
    vector< double > spin_correlations{};
    
    //calc.measure("start calculationng PAIR correlation for certain times");
    double T_max = int( T_max_prefactor * step_amount);


    //t is the time dependencie for the correlation function
    // This is iteated until T_max is reached, beacuse this is also the maximum time step, which can be evaluated via time averaging
    for(int t = 0; t < T_max; t++){
        //correlation
        double G = 0;
        //double G1 = 0;
        //double G2 = 0;
        double S_x_squared = 0;
        double S_y_squared = 0;
        double S_z_squared = 0;
        double S_i_squared_sum =0;

        //this is the integration over tau 
        for(int tau = 0; tau < T_max; tau ++ ){

            S_x_squared = 0;
            S_y_squared = 0;
            S_z_squared = 0;
            S_i_squared_sum = 0;


            //summation over pairs 
            for (int i = 0; i < all_index_pairs.size(); i ++ ){
                //Code for debbuging autocrorrealtion and checking if they yield physical results
                //calculates <S_i(tau)^2> to check if the sum over i is 0.75
                
                if(t == 0 && all_index_pairs[0].first == all_index_pairs[0].second){
                    //auto correlation
                    S_x_squared += real(all_spin_vec_configurations[tau][all_index_pairs[i].first][0] * all_spin_vec_configurations[t + tau][all_index_pairs[i].second][0]);
                    S_y_squared += real(all_spin_vec_configurations[tau][all_index_pairs[i].first][1] * all_spin_vec_configurations[t + tau][all_index_pairs[i].second][1]);
                    S_z_squared += real(all_spin_vec_configurations[tau][all_index_pairs[i].first][2] * all_spin_vec_configurations[t + tau][all_index_pairs[i].second][2]);
                }
                


                //calculation only for x component
                G += real(all_spin_vec_configurations[tau][all_index_pairs[i].first][0] * all_spin_vec_configurations[t + tau][all_index_pairs[i].second][0]);
                //G1 += real(all_spin_vec_configurations[tau][all_index_pairs[i].first][1] * all_spin_vec_configurations[t + tau][all_index_pairs[i].second][1]);
                //G2 += real(all_spin_vec_configurations[tau][all_index_pairs[i].first][2] * all_spin_vec_configurations[t + tau][all_index_pairs[i].second][2]);


            }
            if(t == 0 && all_index_pairs[0].first == all_index_pairs[0].second){

                S_i_squared_sum = S_x_squared/system_size + S_y_squared/system_size + S_z_squared/system_size;

                if(abs(S_i_squared_sum - 0.75) > 1 * pow(10., -7.) ){

                    cout << "#### ERROR REGARDING SAME TIME AUTO CORRELATIOBN ####\n" << "<S_x^2> + <S_y^2> + <S_z^2> = " << S_i_squared_sum *pow(10,7) << "e-7" << " != " << 0.75 * pow(10,7) << "e-7" << "?\n"<< endl;
                }
            }

        // Code for debbuging if auto correlations yoield physical results
        /*
        if((tau == 0 || tau == int(T_max) -1 ) && t == 0  && all_index_pairs[0].first == all_index_pairs[0].second){

            double size = all_index_pairs.size();

            double print_G;
            double print_G1;
            double print_G2;

            if(tau == 0){
                print_G = G / size;
                print_G1 = G1 / size;
                print_G2 = G2 / size;
            }
            else{
                print_G  = G /(size * int(T_max) );
                print_G1 = G1 /(size * int(T_max) );
                print_G2 = G2 / (size * int(T_max) );
                      
            }
            G_auto /= size;
            G_auto1 /= size;
            G_auto2 /= size;

            
            //cout << "index pair amount" << all_index_pairs.size() << endl;
            cout << "<S(tau)^2> at tau=" << tau << " = " << G_auto << endl;
            cout << "<S(tau)^2> at tau=" << tau << "  for y = " << G_auto1 << endl;
            cout << "<S(tau)^2> at tau=" << tau << "  for z = " << G_auto2 << endl;
            cout << "sum = " << G_auto + G_auto1 + G_auto2 << endl;
            //cout << "semiintegrated auto correlation = " << print_G  << endl;
            //cout << "semiintegrated auto correlation = " << print_G1  << endl;
            //cout << "semiintegrated auto correlation = " << print_G2 << endl;
            //cout << "sum integrated = " << print_G + print_G1 + print_G2 << endl;
            //cout << "div = " <<(print_G) / G_auto << endl;
        }
        */

        }
    /*
    if(t == 0 && all_index_pairs[0].first == all_index_pairs[0].second){
        //cout << "Auto correlation at t=0 = " << G << endl;     

        double size = all_index_pairs.size();
        G/= (size *     (T_max) );
        G1 /= (size *    (T_max) );
        G2 /= (size * (T_max) );

        cout << "integrated auto correlation = " << G  << endl;
        cout << "integrated auto correlation = " << G1  << endl;
        cout << "integrated auto correlation = " << G2 << endl;



    }
    */
    //normalozation to latice size is done in the overarching function     
    G = 1./ T_max * G;


   
    spin_correlations.emplace_back( G);
    }
    //calc.measure("finish calculationng PAIR correlation for certain times");

    return spin_correlations;
}

//Function that uses lattice pairs to calculate all kind of different pair correlation types
//each calculated entry to the vector is one pair correlation type
//NOTE: Due to calc_pair_correlation_for_pairs this beaves very different to calc_all_magnetization_correlation
//      This is due to a weird way of summing over all pairs
void spin_system::calc_all_pair_correlations(){
    if(all_spin_vec_configurations.size() != step_amount){
        cout << "\nERROR IN CALCULATION OF PAIRCORRELATIONS\nall_spin_vec_configurations() has wrong size " << all_spin_vec_configurations.size() << " while step amount is " << step_amount << endl << endl;
    }
    
    //calc.measure("start calculationng PAIR correlation");

    //can have errors for the last group, for uneven amounts of lattice sites per row
    double multiplicity = 0;
    //Note: this function gets "inverted" at the end of the function. More information there
    vector< vector<double> > pair_correlation_buffer{};
    


    //could lead to errors that i use abs() here
    //normalization constant
    double G_0 = 0;
    int T_max = int( T_max_prefactor * step_amount);


    double mult =0;
    double all_pair_amounts = 0;
    double pair_at_zero = 0;

    //summation over vectors in groups

    //intorduce cut_off fpr debuging aand runtime purposes
    //############################ IMPORTANT LINE ################
    double pair_correlation_cut_off = 0;

    //set as a default that if a 2D system is considered only the first 4 groups are calculated
    if(dimension == 2 && pair_correlation_cut_off == 0 && use_pair_correlation_cut_off_2D == true){
        pair_correlation_cut_off = 4;
    }


    if(pair_correlation_cut_off == 0){
        pair_correlation_cut_off = all_correlation_groups_lattice_pairs.size();
    }
    
    if(pair_correlation_cut_off > all_correlation_groups_lattice_pairs.size()){
        pair_correlation_cut_off = all_correlation_groups_lattice_pairs.size();
    }

    for( int group_index = 0; group_index < pair_correlation_cut_off; group_index++){

        vector< double > G_pair{};
        multiplicity =     all_correlation_groups_lattice_pairs[group_index].size();
        mult += multiplicity;
        //vector with all pairs corresponding to that group
        vector< pair< uint, uint> > all_index_pairs{};
        uint pair_amount = 0;

        for(int vec_index = 0; vec_index < multiplicity; vec_index++){
            all_index_pairs.insert(all_index_pairs.end(), all_correlation_groups_lattice_pairs[group_index][vec_index].begin(), all_correlation_groups_lattice_pairs[group_index][vec_index].end());
        }

        pair_amount = all_index_pairs.size();
        all_pair_amounts += pair_amount;

        G_pair = calc_pair_correlation_for_index_pairs(all_index_pairs);
        


        //calculate normalization constant over this loop 
        G_0 += G_pair[0];

        //Normalize using amount of pairs. ONLY USE ONE NORMALIZATION WAY!
        for(int i = 0; i < G_pair.size(); i++){
            G_pair[i] =  multiplicity* G_pair[i] / ( pair_amount); 

        }

        //this has just debugging purposes
        if(group_index == 0){
            //cout << "Integrated Auto correlation at t=0 = " << G_pair[0] << endl;
        }
        else{
            //cout << multiplicity / pair_amount << endl;
            pair_at_zero += G_pair[0];
        }

        pair_correlation_buffer.emplace_back(G_pair);


    }

    //cout << "Sum of PAIR correlation at t=0 = " << pair_at_zero << endl;
    //cout << "Sum of all pairs = " << all_pair_amounts<< endl;
    //cout << "Sum of multiplicities of all groups: " << mult << endl;
    
    //normalize pair correlations to ther sum being 1 att t=0
    //alternatively normalize them to the amount of pairs used earlier

    for(int i = 0; i < pair_correlation_buffer.size(); i++){
        for(int j = 0; j < pair_correlation_buffer[i].size(); j++){
            //pair_correlation_buffer[i][j] = 1./ G_0 * pair_correlation_buffer[i][j];
        }
    }

    pair_correlations = pair_correlation_buffer;
}



//function that calsculates auto correlation

void spin_system::calc_auto_correlation(){
    if(magnetization.size() != step_amount){
        cout << "\nERROR IN CALCULATION OF AUTO CORRELATION\nmagnetization() has wrong size " << magnetization.size() << " while step amount is " << step_amount << endl << endl;
    }
    
    //correlation variable
    //calculates the autocorrelation for all lattice sites independently
    //in the end we average over all lattice sites
    vector<vector < double>> G_lattice_site{};
    vector<double> single_site_G{};

    int T_max = int( T_max_prefactor * step_amount);

    for(int i = 0; i<T_max; i++){
	    single_site_G.emplace_back(0);
    }
	
    for( int i =0; i < system_size; i++){
        G_lattice_site.emplace_back(single_site_G);
    }
    //if(abs(spin_vec[i].norm() - sqrt(3.) / 4.) > 1 * pow(10., -7.)){
    //    cout << "#### ERROR REGARDING SPIN LENGTH ####\n" << "check spin length " << spin_vec[i].norm() *pow(10,7) << "e-7" << " != " << sqrt(3.) / 4. * pow(10,7) << "e-7" << "?\n"<< endl;
    //}
    
    //k index is the index for each lattice site
    //i index is the time step t
    //j index is the time step tau
    for( int k = 0; k < system_size; k ++){
    	for( int i = 0; i < T_max; i++){
            if(abs(all_spin_vec_configurations[i][k].norm() - sqrt(3.) / 4.) > 5 * pow(10., -8.)){
                cout << "#### ERROR REGARDING SPIN LENGTH FOUND IN AUTO CORR CALC ####\n" << "check spin length " << spin_vec[i].norm() *pow(10,7) << "e-7" << " != " << sqrt(3.) / 4. * pow(10,7) << "e-7" << "?\n"<< endl;
            }
            //correlation variable
            for(int j = 0; j< T_max; j++){	
            	G_lattice_site[k][i] += real(all_spin_vec_configurations[j][k][0] * all_spin_vec_configurations[i + j][k][0]);
	        }
        }
    }
    
    //normalizing
    for( int k = 0; k < system_size; k ++){
        for( int i = 0; i < T_max; i++){
            G_lattice_site[k][i] = G_lattice_site[k][i] / T_max;
        }
    }
    pair_correlations = G_lattice_site;
}

