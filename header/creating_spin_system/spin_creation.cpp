
#include <iostream>
#include "../../header/creating_spin_system/semi_classical_spins.h"





// function that creates spins for each lattice site. Each component is normal distributed.

void spin_system::create_spins(){

    Eigen::Vector3cd spin;

    for(int x = 0; x < amount_per_row; x++){
        for(int y = 0; y < amount_per_row; y++){
            for(int z = 0; z < amount_per_row; z++){

                //test with second method to draw spin vector


                spin_vec.emplace_back(create_Spin_Vector2( ));
                //spin_vec.emplace_back(create_Spin_Vector_non_normalized()); //create_Spin_Vector2( ));
            }
        }
    }
}

//Randomly drawing a spcin vector by drawing each component independently following a normal distribution
//vectoors is not normalized

Eigen::Vector3cd spin_system::create_Spin_Vector_non_normalized(){
    // Vector to store the components of the spin vector
    Eigen::Vector3cd spinVector;
    double mean = 0.0;
    // variance is 1/4
    double stddev = 0.5; // Standard deviation
    // Create a random number generator
    
    //random_device rd() ;  // Seed generator nomally gen(rd())

    //seed &generator are located in matrices& consz.h

    // Create a normal distribution with the given mean and standard deviation
    normal_distribution<double> dist(mean, stddev);

    // Fill the vector with normally distributed values

    spinVector << dist(gen), dist(gen), dist(gen);
    //spinVector = sqrt(3.) / 2. * spinVector.normalized();

    return spinVector;
}


//alternative way to draw the spin vectors randomly. This time by utilizing a vector in spherical coordinates

Eigen::Vector3cd spin_system::create_Spin_Vector2() {
    // Vector to store the components of the spin vector
    Eigen::Vector3cd spinVector;

    // Create a random number generator
    
    //random_device rd() ;  // Seed generator nomally gen(rd())

    //seed &generator are located in matrices& consz.h

    // Create a normal distribution with the given mean and standard deviation
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> dist_theta(-1.0, 1.0);

    // Generate random spherical coordinates
    double phi = dist_phi(gen);
    double cos_theta = dist_theta(gen);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    // Convert to Cartesian coordinates
    double x = sin_theta * std::cos(phi);
    double y = sin_theta * std::sin(phi);
    double z = cos_theta;


    spinVector << x, y, z;

    return sqrt(3.) / 2. * spinVector;
}

