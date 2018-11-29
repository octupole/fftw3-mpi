//============================================================================
// Name        : test-fftw3.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#define OMPI_SKIP_MPICXX
#ifdef HAVE_MPI
#include <mpi.h>
#endif
using namespace std;

class MPI{
private:
public:
	static void Init(int * argc, char **argv[]){
#ifdef HAVE_MPI
		MPI_Init(argc,argv);
#endif
		
	}
	static void Comm_size(MPI_Comm comm, int *world_size){
#ifdef HAVE_MPI
		MPI_Comm_size(comm,world_size);
#else
		*world_size=1;
#endif
	}
	static void Comm_rank(MPI_Comm comm, int * world_rank){
#ifdef HAVE_MPI
		MPI_Comm_rank(comm,world_rank);
#else
		*world_rank=1;
#endif
	}
	static void Get_processor_name(char * processor_name,int * length){
#ifdef HAVE_MPI
		MPI_Get_processor_name(processor_name,length);
#else
		
#endif

	}
};

int main(int argc, char *argv[]) {
	MPI::Init(&argc,&argv);
	int world_size{0},world_rank{0};
	MPI::Comm_size(MPI_COMM_WORLD,&world_size);

	MPI::Comm_rank(MPI_COMM_WORLD,&world_rank);
	char name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI::Get_processor_name(name,&name_len);
	cout << "Hello from processor "<< name << ",  rank "<< world_rank << " out of "<< world_size << " processors" <<endl;
	return 0;
}

