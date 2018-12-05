//============================================================================
// Name        : test-fftw3.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <random>

#define OMPI_SKIP_MPICXX
#include <mpi.h>
#include <fftw3-mpi.h>
using namespace std;

double my_func(int i,int j, int k){
	static std::default_random_engine generator;
	static std::uniform_real_distribution<double> distribution(0.0,1.0);
	return distribution(generator);;
}
int main(int argc, char **argv)
{
    const ptrdiff_t L =128, M = 128, N = 128;
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

    /* initialize rin to some function my_func(x,y,z) */
    for (i = 0; i < local_n0; ++i)
       for (j = 0; j < M; ++j)
         for (k = 0; k < N; ++k)
       rin[(i*M + j) * (2*(N/2+1)) + k] = my_func(local_0_start+i, j, k);

    /* compute transforms as many times as desired */
    fftw_execute(plan);

    fftw_destroy_plan(plan);

    MPI_Finalize();
}
