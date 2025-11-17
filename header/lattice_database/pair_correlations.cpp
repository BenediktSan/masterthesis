#include <iostream>
#include "../../header/lattice_database/lattice_database.h"






//################Functions to calculate paircorrelations



//this function finds all possible pair correlations in my lattice or in my cutoff range
//returns a vecor filled with exemplary vectors which stand for a certain correlation e.g. a neraest neighbour correlation in x-direction

//Goal is to group all vectors that belong to the same pair correlation type together. 
//the amount of vectors grouped should match their multiplicity, which is taken from the paper
//grouping is done by taking the length of the vector and angle between +- the magnetic field vector and this vector
//for one of these grouped vectors the suiting lattice pairs are computed, which is then used to calculate the correaltion
//return value are these groups

void lattice_database::find_correlation_groups(){

    //Vector in which these aforementioned groups are saved
    vector < vector < Eigen::Vector3cd> > list_of_groups{};

    //Vector in which the correspondimng multiplicity is saved
    vector < uint > multiplicity_of_groups{};

    //Vector for Eigen::Vectors who form a group

    vector <Eigen::Vector3cd> pair_correlation_vectors{};



    //copying all connecting vectors for the lattice site at the origin to this vector so i can delete after assigning them to a group.
    //this way I cant double assign

    vector<Eigen::Vector3cd> connecting_vector_copy = connecting_vector[0];

    //also doing this for the distances to get a coherent index relation between a vector and a distance
    vector<double> distances_copy = distances[0];

    
    while(connecting_vector_copy.empty() == false){
        
        Eigen::Vector3cd vec = connecting_vector_copy[0];
        double reference_angle_to_Bfield = acos(real(vec.dot(B_field_direction)) / real(B_field_direction.dot(B_field_direction) * vec.dot(vec)));
        double other_reference_angle_to_Bfield = acos(real(vec.dot(-B_field_direction)) / real(B_field_direction.dot(B_field_direction) * vec.dot(vec)));

        uint multiplicity = 1;
        vector<Eigen::Vector3cd> correlation_group{vec};
        //cout << vec.transpose()* 1./lattice_const <<endl;
        //cout << "Angle is " << reference_angle_to_Bfield << " or " << other_reference_angle_to_Bfield <<endl;
        //cout << "\tadded " << reference_angle_to_Bfield + other_reference_angle_to_Bfield << endl;
        
        for(int i = 1; i < connecting_vector_copy.size();i++){
            
            double angle_to_Bfield = acos(real(connecting_vector_copy[i].dot(B_field_direction) )/ real(B_field_direction.dot(B_field_direction) * connecting_vector_copy[i].dot(connecting_vector_copy[i])));

            //if( i < 4){
            //    cout << "\nchecking for vector " << vec.transpose()* 1./lattice_const << " against " << connecting_vector_copy[i].transpose() * 1./lattice_const <<endl;
            //    cout << "distances[i] = " << distances_copy[i] << "\t" << "condition = " <<  (abs(angle_to_Bfield -reference_angle_to_Bfield) < pow(10,-14)) << " " <<   (abs(angle_to_Bfield - other_reference_angle_to_Bfield) < pow(10,-15)) << endl;
            //    cout << "Angle is " << std::setprecision (15) << reference_angle_to_Bfield << " or " << other_reference_angle_to_Bfield << " compared to " << angle_to_Bfield << endl;
            //    cout << "diff " << other_reference_angle_to_Bfield -angle_to_Bfield <<endl;
            //    cout << "Index " << i <<endl;
            //}


            if( abs(vec.norm() - distances_copy[i]) < pow(10,-14) && ( (abs(angle_to_Bfield - reference_angle_to_Bfield) < pow(10,-14))|| (abs(angle_to_Bfield - other_reference_angle_to_Bfield) < pow(10,-14)))){


                //cout << "\nchecking for vector " << vec.transpose()* 1./lattice_const << " against " << connecting_vector_copy[i].transpose() * 1./lattice_const <<endl;
                //cout << "distances[i] = " << distances_copy[i] << "\t" << "condition = " << (vec.norm() == distances_copy[i] && ( angle_to_Bfield == reference_angle_to_Bfield || angle_to_Bfield == other_reference_angle_to_Bfield)) ;
                //cout << " Index " << i <<endl;

                //adding vector to group
                correlation_group.emplace_back(connecting_vector_copy[i]);
                multiplicity += 1;

                //deleting from copy vectors
                connecting_vector_copy.erase(connecting_vector_copy.begin() + i);
                distances_copy.erase(distances_copy.begin() + i);

                //reseting i so we still check every index, even thoughwe deleted an entry
                i -= 1;




                //now consider a special case where due to periodic boundary conditions
                //different vectors map to the same lattice site
                // this only happens for even lattcie sizes
                // and leads to an underestimation of the multiplicity if not considered
                // also when calculating lattice pairs connecting lattice sites are not found,
                //because (0,0,-2) is not considered for a 4x4x4 lattice
                // but (0,0,2) is

            }

        }
        //after iteration over all vectors the grpoup should be complete
        //therefore adding them to our lists
        list_of_groups.emplace_back(correlation_group);
        multiplicity_of_groups.emplace_back(multiplicity);            


        //deleting starting vector for group from copy vectors
        connecting_vector_copy.erase(connecting_vector_copy.begin());
        distances_copy.erase(distances_copy.begin());
    }
    
    //cout << "Amount of groups " << list_of_groups.size() << endl;

    /*
    cout << "\n#########Printing the first different correlation types and their multiplicity#########" <<endl;
    cout << "Amount of groups " << list_of_groups.size() << endl;
    int multiplicity_sum = 0;
    for(int i = 0; i < multiplicity_of_groups.size(); i ++){
        cout << "\tVector \t\t\tMultiplicity" <<endl;

        cout << "\t" << list_of_groups[i][0].transpose() * 1./lattice_const << "\t" << multiplicity_of_groups[i] <<endl;

        multiplicity_sum += multiplicity_of_groups[i];
    }
    cout <<"\tSum over Multiplicities yields " << multiplicity_sum << endl;
    */
    
    //assigning values

    list_of_correlation_vector_groups = list_of_groups;
    multiplicity_of_correlation_groups = multiplicity_of_groups;

}









//used to check which vectors satisfy my conditions
//also function that uses indices of functions that satisfy that to calculate correlations
//returned vector includes all index pairs for that certain wanted correlation

vector< pair<uint, uint> > lattice_database::find_suiting_lattice_pairs(Eigen::Vector3cd exemplary_vector){
 
    
    double exemplary_length = exemplary_vector.norm();

    //cout << "FIND FOR " << exemplary_vector.transpose() * 1./lattice_const << " with length " << exemplary_length * 1./lattice_const << endl;

    //This case for autocorrelationis handled in the loop
    //if(exemplary_length == 0){
    //    cout << "\n\n\tERROR VECTOR FOR PAIR CORRELATION HAS LENGTH 0\n\n" << endl;
    //}

    //vector to store all index pairs
    vector< pair<uint, uint> > index_pairs{};

    
    if(amount_per_row == 1){
        //special case for three particle system
        if(exemplary_length == 0){
            pair<uint, uint> auto_pair(0,0);
            index_pairs.emplace_back(auto_pair);
            pair<uint, uint> auto_pair1(1,1);
            index_pairs.emplace_back(auto_pair1);
            pair<uint, uint> auto_pair2(2,2);
            index_pairs.emplace_back(auto_pair2);
            return index_pairs;
        }
        else if( abs(exemplary_length - lattice_const) < pow(10,-14) ){
            pair<uint, uint> pair01(0,1);
            index_pairs.emplace_back(pair01);
            pair<uint, uint> pair10(1,0);
            index_pairs.emplace_back(pair10);
            pair<uint, uint> pair12(1,2);
            index_pairs.emplace_back(pair12);
            pair<uint, uint> pair21(2,1);
            index_pairs.emplace_back(pair21);
            return index_pairs;
        }
        else if( abs(exemplary_length - sqrt(2.) * lattice_const) < pow(10,-14) ){
            pair<uint, uint> pair02(0,2);
            index_pairs.emplace_back(pair02);
            pair<uint, uint> pair20(2,0);
            index_pairs.emplace_back(pair20);
            return index_pairs;
        }
        else{
            cout << "\n\n\tERROR NO SUITING PAIR FOUND IN 3 PARTICLE SYSTEM\n\n" << endl;
            return index_pairs;
        }
    }
    
    //this loop checks if the distance between two lattice sites correspond to the lenth of the exemplary vector for the pari correlation we want to calculate
    // additionaly it checks if its parallel to the eexmeplary_vector, so that we only get the right correlations
    //calculated index pairs are being appended to the index pair vector and later returned
    for(int j = 0; j < system_size; j++){

        //special case for autocorrelation
        if (exemplary_length == 0){
            pair<uint, uint> auto_pair(j,j);
            index_pairs.emplace_back(auto_pair);
            continue;
        }

        //if( abs(exemplary_length - distances[0][j]) < pow(10,-14) && abs(connecting_vector[0][j].dot(exemplary_vector) - exemplary_vector.dot(exemplary_vector)) < pow(10,-14) ){
            for(int i = 0; i < system_size; i++){

                //just add brute force comparison here
                if (abs(connecting_vector[j][i].dot(exemplary_vector) - exemplary_vector.dot(exemplary_vector)) < pow(10,-8) && abs(distances[j][i] - exemplary_length) < pow(10,-8) ){

                    int new_first = j;
                    int new_second = i;

                    //cout << "new first index " << new_first << "\tnew second index " << new_second << endl;

                    pair<uint, uint> pair(new_first,new_second); 
                    index_pairs.emplace_back(pair);

                }


                
                //int second_index = i + j;
                //if(second_index >= system_size){
                //    second_index = i + j - (system_size); 
                //}

               //cout <<"first index " << i << "\tsecond index " <<  second_index  << endl;
               //cout << "new first index " << new_first << "\tnew second index " << new_second << endl;

                //pair<uint, uint> pair(i,second_index);  
                //pair<uint, uint> pair(new_first,new_second); 
                //index_pairs.emplace_back(pair);
            }
        //}

    }
    if(index_pairs.size() != system_size){
        cout << "\t\tAmount of suiting pairs found " << index_pairs.size() << " expected " << system_size << " for vector " << exemplary_vector.transpose() * 1./lattice_const <<  endl;
    }
    return index_pairs;
}

    //functions that calculates for every vector in evry correlation type group all suiting lattice pairs ans saves
    //This is later returned and saved to lattice_system

    void lattice_database::calc_all_lattice_pairs(){
        if (list_of_correlation_vector_groups.size() == 0 || multiplicity_of_correlation_groups.size() == 0){
            cout << "\n\n\tERROR TRYING TO SORT EMPTY VECTORS\n\n" << endl;
        }
        //for good measure sort again
        sort_correlation_groups();


        double amount_of_pairs = 0;

        //butt ugly variables
        vector<pair<uint, uint>> single_vector_lattice_pairs{};
        vector<vector<pair<uint, uint>>> group_lattice_pairs{};
        //now everything combined
        vector<vector<vector<pair<uint, uint>>>> all_group_pairs;

        for(int i = 0; i < list_of_correlation_vector_groups.size(); i++){
            group_lattice_pairs = {};

            for(int j = 0; j < list_of_correlation_vector_groups[i].size(); j++){
                single_vector_lattice_pairs = find_suiting_lattice_pairs(list_of_correlation_vector_groups[i][j]);
                group_lattice_pairs.emplace_back(single_vector_lattice_pairs);
                amount_of_pairs += single_vector_lattice_pairs.size();
            }
            all_group_pairs.emplace_back(group_lattice_pairs);
        }
        all_correlation_groups_lattice_pairs = all_group_pairs;

        //cout << "\nAmount of lattice pairs found " << amount_of_pairs << endl;
        //cout <<"amount expected " << system_size * system_size << endl;
        //cout << all_group_pairs.size() << " groups found " << endl;
        //cout << all_group_pairs[0].size() << " vectors in first group " << endl;
        //cout << all_group_pairs[1].size() << " vectors in second group " << endl;
        //cout << all_group_pairs[2].size() << " vectors in third group " << endl;
    }






    //function that sorts the correlations according to their connecting vector length
    //also sorts the multiplicitys the same way
    //mainly written by chatgpt
    //!!!! Needs a condition how to handle same distances but different interactions



void lattice_database::sort_correlation_groups() {
    if (list_of_correlation_vector_groups.empty() || multiplicity_of_correlation_groups.empty()) {
        std::cerr << "\n\n\tERROR: Trying to sort empty vectors\n\n" << std::endl;
        return;
    }

    coupling_const_sorted = coupling_const[0];
    // Build index list [0, 1, 2, ..., N-1]
    std::vector<size_t> indices(list_of_correlation_vector_groups.size());
    std::iota(indices.begin(), indices.end(), 0);

    const double tol = 1e-7; // tolerance for floating point comparisons

    auto comparator = [this, tol](size_t a, size_t b) {
        // --- Primary key: vector norm ---
        double normA = list_of_correlation_vector_groups[a].empty() ? 0.0 : list_of_correlation_vector_groups[a][0].norm();
        double normB = list_of_correlation_vector_groups[b].empty() ? 0.0 : list_of_correlation_vector_groups[b][0].norm();


        // --- Secondary key: sum of abs(coupling_const[row][col]) ---
        //fuck it lets go brute force
        double couplingA = 0.0, couplingB = 0.0;
            Eigen::Vector3cd vecA = list_of_correlation_vector_groups[a][0].normalized();
            Eigen::Vector3cd vecB = list_of_correlation_vector_groups[b][0].normalized();
            Eigen::Vector3cd B_field_direction_norm = B_field_direction.normalized();

            couplingA = (1. - 3. * (vecA.real().dot(B_field_direction_norm.real())) * (vecA.real().dot(B_field_direction_norm.real()))) * 1./(normA * normA * normA);
            couplingB = (1. - 3. * (vecB.real().dot(B_field_direction_norm.real())) * (vecB.real().dot(B_field_direction_norm.real()))) * 1./(normB * normB * normB);
            //couplingA = std::fabs(coupling_const_sorted[a]);
            //couplingB = std::fabs(coupling_const_sorted[b]);

            //cout << a <<" | " << list_of_correlation_vector_groups[a][0].real().transpose() * 1./lattice_const << " with norm " << normA * 1./lattice_const ;
            //cout << " and coupling " << couplingA / (1. / (lattice_const * lattice_const * lattice_const)) << endl;

        //cout << "Comparing groups " << a << " and " << b << " with norms " << normA << " and " << normB << " and couplings " << couplingA/coupling_const[0][1] << " and " << couplingB/coupling_const[0][1] << endl;
        if (std::fabs(normA - normB) < tol) {

            return abs(couplingA) > abs(couplingB); // descending
        }

        // --- Tertiary key: original index (deterministic tie-break) ---
        //true = tauschen
        return normA < normB;
    };

    std::sort(indices.begin(), indices.end(), comparator);

    // Reorder groups + multiplicities consistently
    std::vector<std::vector<Eigen::Vector3cd>> sorted_groups;
    std::vector<unsigned int> sorted_multiplicity;

    sorted_groups.reserve(indices.size());
    sorted_multiplicity.reserve(indices.size());

    for (size_t i : indices) {
        sorted_groups.push_back(std::move(list_of_correlation_vector_groups[i]));
        sorted_multiplicity.push_back(multiplicity_of_correlation_groups[i]);
        //cout << coupling_const_sorted[i]/coupling_const[0][1] << endl;
    }

    list_of_correlation_vector_groups = std::move(sorted_groups);
    multiplicity_of_correlation_groups = std::move(sorted_multiplicity);
}






//function that prints correlation groups and their multiplicity to the terminal

void lattice_database::print_correlation_groups(){

    cout << "\n#########Printing the different correlation types and their multiplicity for a " << dimension <<"D lattice#########" <<endl;
    cout << "\n\tAmount of vector correlation groups found is " << list_of_correlation_vector_groups.size() << endl <<endl;
    int multiplicity_sum = 0;
    
    vector<double> coupling = coupling_const[0];

    double coupling_norm = coupling[1];

    for(int i = 1; i < coupling.size(); i++){
        if( abs(coupling[i]) > abs(coupling_norm) ){
            coupling_norm = coupling[i];
        }
    }

    cout << "\nThe largest coupling constant has the value " << coupling_norm << endl <<endl;

    for(int i = 1; i < coupling.size(); i++){
	    coupling[i] /= coupling_norm;
        //to reproduce timos table
        coupling[i] = coupling[i];
    } 	

    if(dimension ==2){
        coupling = coupling_const_heisenberg[0];
    }



    cout << "\t\tVector  \t\tVector length \t\tCoupling const/C(0) \t\t\tMultiplicity" <<endl << endl;
    cout.precision(4);
    for(int i = 0; i < multiplicity_of_correlation_groups.size(); i ++){


        uint second_index = all_correlation_groups_lattice_pairs[i][0][0].second;
        cout << "\t" <<i <<"\t" << list_of_correlation_vector_groups[i][0].real().transpose() * 1./lattice_const << "\t\t\t"  ;
        cout << list_of_correlation_vector_groups[i][0].norm() * 1./lattice_const << "\t\t\t" << coupling[second_index] << "\t\t\t\t\t" << multiplicity_of_correlation_groups[i] << endl;

        multiplicity_sum += multiplicity_of_correlation_groups[i];

        if(i == multiplicity_of_correlation_groups.size() - 1 && multiplicity_of_correlation_groups[i] == 1){
            cout << "Special Case: Due to rotational and Inversion symmetrie this vector should have the multiplicity 8.";
            cout << "Due to the fact that this is the longest vector pointing to the middle of the cube and due to periodic boundarie conditions its multiplicity is only 1";
            cout << "This stems to the fact that all these transformations map to the same lattice site, therefore only one vector exists " <<endl;
        }
    }
    cout <<"\n\tSum over Multiplicities yields " << multiplicity_sum << endl;
}





