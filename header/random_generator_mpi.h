#include<iostream>
#include<random>
#include<mpi.h>

/* generates an overall seed from inserted seed string and broadcasts it to all cores */
inline uint generate_seed( std::string seed_str, const uint my_rank )
{
    if( seed_str.find("random") != std::string::npos ) // random seed
    {
        uint seed{};
        if( my_rank == 0 ) // draw seed on rank 0 
        {
            std::random_device my_seed{};
            seed = my_seed();
        }
        MPI_Bcast( &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD ); // broadcast seed to all the other ranks
        return seed;
    }  
    else // preset seed (transform from string to uint)
    {
        return std::stoi( seed_str );
    }    
}

/* generates a rank dependent seed for the Mersenne Twister Generator
the unsigned integer seed is build from the formula |seed - rank| ensuring different seeds on each rank */
inline uint throw_seed( const uint seed, const uint my_rank )
{
    return std::abs(static_cast<int>(seed) - static_cast<int>(my_rank));
}
