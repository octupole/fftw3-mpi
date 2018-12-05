//============================================================================
// Name        : test-fftw3.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <random>


//#include <mpi.h>
#include <fftw3-mpi.h>
using namespace std;
double Wtime(){
	timeval tim;
	double temp=tim.tv_sec+(tim.tv_usec/1000000.0);
	return temp;
	}


int main(int argc, char **argv)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);


    const ptrdiff_t L =256, M = 256, N = 256;
    fftw_plan plan;
    double *rin;
    fftw_complex *cout;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j, k;


    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_3d(L, M, N/2+1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    rin = fftw_alloc_real(2 * alloc_local);
    cout = fftw_alloc_complex(alloc_local);

    /* create plan for out-of-place r2c DFT */
    plan = fftw_mpi_plan_dft_r2c_3d(L, M, N, rin, cout, MPI_COMM_WORLD,
                                    FFTW_MEASURE);
    double time0=MPI_Wtime();
    /* initialize rin to some function my_func(x,y,z) */
    for (i = 0; i < local_n0; ++i)
       for (j = 0; j < M; ++j)
         for (k = 0; k < N; ++k)
       rin[(i*M + j) * (2*(N/2+1)) + k] = distribution(generator);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank == 2) {
    	std::cout << local_0_start <<endl;
    }
    /* compute transforms as many times as desired */
    fftw_execute(plan);

    double time1=MPI_Wtime();

    if(rank==0) std::cout << time1-time0 << std::endl;

    fftw_destroy_plan(plan);

    MPI_Finalize();
}
