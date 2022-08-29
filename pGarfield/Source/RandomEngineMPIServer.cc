/**
   @author Ali Sheharyar
   @organization Texas A&M University at Qatar
*/

#include <iostream>
#include "Garfield/RandomEngineMPIServer.hh"

namespace Garfield {

RandomEngineMPIServer::RandomEngineMPIServer() : rng(0) {

	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  std::cout << "RandomEngineMPIServer:\n";
  std::cout << "    Generator type: MPI TRandom3\n";
  std::cout << "    Seed: " << rng.GetSeed() << "\n";
  std::cout << "	MPI size: " << mpi_size << "\n";
  std::cout << "	MPI rank: " << mpi_rank << "\n";

}

RandomEngineMPIServer::~RandomEngineMPIServer() {

}

void
RandomEngineMPIServer::Run() {

	int N;
	MPI_Status recv_status;
	double rand_numbers[RAND_BUFFER_SIZE];

	while(1) {
		MPI_Recv( &N, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
			MPI_COMM_WORLD, &recv_status );

		if(recv_status.MPI_TAG == RandomEngineMPI::DRAW) {
			rng.RndmArray(RAND_BUFFER_SIZE, rand_numbers);
			MPI_Ssend( rand_numbers, RAND_BUFFER_SIZE, MPI_DOUBLE, 
				recv_status.MPI_SOURCE, RandomEngineMPI::DRAW, MPI_COMM_WORLD );
		}
		else if(recv_status.MPI_TAG == RandomEngineMPI::KILL && recv_status.MPI_SOURCE==0) {
			std::cerr << "RandomEngineMPIServer: Received KILL message from the root process\n";
			break;
		}
		else {
			std::cerr << "RandomEngineMPIServer: Unknown MPI tag\n";
		}
	}
}

}
