# What is pGarfield-toolkit?
The pGarfield-toolkit (Parallel Garfield Toolkit) is an MPI-based version of the Garfield++ (https://garfieldpp.web.cern.ch/garfieldpp/). It supports simulation of multiple detector events in parallel over several machines in a distributed fashion. It consists of two components - core simulator and the simulation client. The core simulator is basically the original Garfield SDK from CERN with slight modifications to add the MPI support. The simulation client in the main MPI program where a user constructs the simulation scenarios and performs the simulation using the API provided by Garfield. 

The simulation client spawns multiple MPI processes. One of the processes acts as a master that assigns jobs (event simulations) to worker processes. Another MPI process generates random numbers and distributes to worker processes as they request. All other process, called workers, perform the actual simulation using the Garfield API.

# Dependencies
Boost v1.58.0 or later
ROOT v5.34.36 or later

# Compiling
Note: The makefile uses the Intel C++ compilers (icc/mpicc) and Intel MPI. Make sure the intel compiler is accessible before compiling. Make changes to use a different C++ compiler as desired.

1) Build pGarfield-sim
source <root_install_dir>/bin/thisroot.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<boost_install_dir>/lib
export GARFIELD_HOME=<pgarfield_sim_source_dir>
cd <pgarfield_sim_source_dir>
PARALLEL=1 make

Note that running just make command will produce the serial version of the Garfield. The flag PARALLEL=1 enables enables the MPI support and also links with the system MPI header files and libraries.


# Documentation
The repository includes a sample simulation client that supports the input of simulation parameters through a file (card file). The card file is specified through the command line argument. See pGarfield-client/example1/card.ini for a sample card file. 

To run the main program (let us call it 'program')
mpirun -np 16 ./program --card card.ini


# Information
Send questions and comments to: Ali Sheharyar ali.sheharyar@gmail.com
