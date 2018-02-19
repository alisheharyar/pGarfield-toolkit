#ifndef G_RANDOM_ENGINE_MPI_SERVER_H
#define G_RANDOM_ENGINE_MPI_SERVER_H

#include <mpi.h>
#include <TRandom3.h>
#include "RandomEngineMPI.hh"

namespace Garfield {

class RandomEngineMPIServer {

  public:
    // Constructor
    RandomEngineMPIServer();
    // Destructor    
    ~RandomEngineMPIServer();    

    void Run();

private:
    TRandom3 rng;

    int mpi_rank;
    int mpi_size;
};

}

#endif