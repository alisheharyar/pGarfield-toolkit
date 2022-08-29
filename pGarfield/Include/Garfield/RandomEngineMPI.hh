// MPI random number generator

#ifndef G_RANDOM_ENGINE_MPI_H
#define G_RANDOM_ENGINE_MPI_H

#include <mpi.h>
#include "RandomEngine.hh"

#define RAND_BUFFER_SIZE 1000


namespace Garfield {

class RandomEngineMPI : public RandomEngine {

  public:
    enum MessageType { DRAW, SEED, KILL};

  public:
    // Constructor
    RandomEngineMPI();
    // Destructor    
    ~RandomEngineMPI();    
    // Call the random number generator
    double Draw();
    // Initialise the random number generator
    void Seed(unsigned int s);
    
  private:
    // Buffer where random numbers received from server will be stored
    double rand_numbers[RAND_BUFFER_SIZE];
    // index of the next random number
    int index;
    // MPI size
    int mpi_size;
    // MPI rank
    int mpi_rank;
};

}

#endif
