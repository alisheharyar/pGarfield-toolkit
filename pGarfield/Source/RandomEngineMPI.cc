/**
   @author Ali Sheharyar
   @organization Texas A&M University at Qatar
*/

#include <iostream>
#include "Garfield/RandomEngineMPI.hh"

namespace Garfield {

RandomEngineMPI randomEngine;

RandomEngineMPI::RandomEngineMPI() : index(-1) {

  std::cout << "RandomEngineMPI:\n";
  std::cout << "    Generator type: TRandom3\n";
  //std::cout << "    Seed: " << rng.GetSeed() << "\n";
}

RandomEngineMPI::~RandomEngineMPI() {

}

double
RandomEngineMPI::Draw() {

	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// Check if there is a need to retrieve the random numbers
	// from the server. This will happen when it is drawn for
	// the first time or when all are consumed
	if(index==-1 || index==RAND_BUFFER_SIZE) {
		int N = RAND_BUFFER_SIZE;
		MPI_Status recv_status;

		// The server will ignore the N and will send RAND_BUFFER_SIZE random numbers
		MPI_Ssend( &N, 1, MPI_INT, mpi_size-1, RandomEngineMPI::DRAW, MPI_COMM_WORLD );
		MPI_Recv( rand_numbers, N, MPI_DOUBLE, mpi_size-1, RandomEngineMPI::DRAW, 
			MPI_COMM_WORLD, &recv_status );

		// reset the index
		index=0;
	}

	return rand_numbers[index++];
}

void
RandomEngineMPI::Seed(unsigned int s) {

  //rng.SetSeed(s);
  std::cout << "RandomEngineMPI::Seed:\n";
  std::cout << "	NOT IMPLEMENTED YET\n";
  //std::cout << "    Seed: " << rng.GetSeed() << "\n";

}

}

