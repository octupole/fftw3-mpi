//============================================================================
// Name        : test-fftw3.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <random>

#include "Array.h"
#include <mpi.h>
#include <fftw3-mpi.h>
#include "fftw++.h"

using namespace std;
using namespace Array;

double Wtime(){
	timeval tim;
	double temp=tim.tv_sec+(tim.tv_usec/1000000.0);
	return temp;
	}


int bamby(MPI_Comm & myComm){
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
    int rank,size;
    MPI_Comm_rank(myComm,&rank);
    MPI_Comm_size(myComm,&size);

    const ptrdiff_t L =256, M = 256, N = 256;
    unsigned int nx{L},ny{M},nz{N},nzp{(N/2+1)};
	array3<double> ri(nx,ny,nz);

	for (size_t i = 0; i < nx; ++i)
       for (size_t j = 0; j < ny; ++j)
         for (size_t k = 0; k < nz; ++k)
        	 ri[i][j][k] = distribution(generator);

	array3<double> rin;
	array3<Complex> c_out;

    ptrdiff_t alloc_local, local_n0, local_0_start;

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_3d(L, M, N/2+1, myComm,
                                         &local_n0, &local_0_start);

    rin.Allocate(local_n0,ny,nzp*2);
    c_out.Allocate(local_n0,ny,nzp);

    fftwpp::rcfft3d Forward3(nx,ny,nz,&rin[0][0][0],&c_out[0][0][0],myComm);
    fftwpp::crfft3d Backward3(nx,ny,nz,&c_out[0][0][0],&rin[0][0][0],myComm);

    for (size_t i = 0; i < local_n0; ++i)
       for (size_t j = 0; j < M; ++j)
         for (size_t k = 0; k < N; ++k)
        	 rin[i][j][k] = ri[local_0_start+i][j][k];
    if(rank == 0) {
    	cout << rin[4][129][20]<<endl;
    }

    double time0=MPI_Wtime();
    /* initialize rin to some function my_func(x,y,z) */

    /* compute transforms as many times as desired */
    Forward3.fft(&rin[0][0][0],&c_out[0][0][0]);
    Backward3.fftNormalized(&c_out[0][0][0],&rin[0][0][0]);

    if(size == 1) {
    	cout << rin[64][129][20]<<endl;
    }
    double time1=MPI_Wtime();

    if(rank == 0) {
    	cout << rin[4][129][20]<<endl;
    }
    return 0;
}
fftw_plan myPlan0(ptrdiff_t L, ptrdiff_t M, ptrdiff_t N, double * rin, Complex * c_out, MPI_Comm comm){
	  return fftw_mpi_plan_dft_r2c_3d(L, M, N, rin, (fftw_complex *) c_out, MPI_COMM_WORLD,FFTW_MEASURE);
}
fftw_plan myPlan1(ptrdiff_t L, ptrdiff_t M, ptrdiff_t N, double * rin, Complex * c_out, MPI_Comm comm){
	  return fftw_mpi_plan_dft_c2r_3d(L, M, N, (fftw_complex *) c_out, rin, MPI_COMM_WORLD,FFTW_MEASURE);
}

int bamba(){
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
    const ptrdiff_t L = 128, M = 128, N = 128;
    unsigned int nx{L},ny{M},nz{N},nzp{(N/2+1)};
    fftw_plan plan0,plan1;
    ptrdiff_t alloc_local, local_n0, local_0_start;

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_3d(L, M, N/2+1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
	array3<double> rin;
	array3<Complex> c_out;
    rin.Allocate(local_n0,ny,nzp*2);
    c_out.Allocate(local_n0,ny,nzp);
    double * p_rin=&rin[0][0][0];
    Complex * p_cout=&c_out[0][0][0];

    /* create plan for out-of-place r2c DFT */
    plan0 = myPlan0(L, M, N, p_rin, p_cout, MPI_COMM_WORLD);
    plan1 = myPlan1(L, M, N, p_rin, p_cout, MPI_COMM_WORLD);
	array3<double> ri(nx,ny,nz);

	for (size_t i = 0; i < nx; ++i)
       for (size_t j = 0; j < ny; ++j)
         for (size_t k = 0; k < nz; ++k)
        	 ri[i][j][k] = distribution(generator);

    /* initialize rin to some function my_func(x,y,z) */
    double time0=MPI_Wtime();
    for(int w{0};w<200;w++){
    	if(!rank) cout << "Step No. "<<w<<endl;
    	for (int i = 0; i < local_n0; ++i)
    		for (int j = 0; j < M; ++j)
    			for (int k = 0; k < N; ++k)
    				rin[i][j][k] = ri(local_0_start+i, j, k);

    	/* compute transforms as many times as desired */
    	fftw_execute(plan0);

    	fftw_execute(plan1);


    }
    double time1=MPI_Wtime();
    if(rank==0) cout << time1-time0<<endl;
    fftw_destroy_plan(plan0);
    fftw_destroy_plan(plan1);

    MPI_Finalize();
	return 0;
}
int threads_ok;
int main(int argc, char **argv)
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    threads_ok = provided >= MPI_THREAD_FUNNELED;

    if (threads_ok) threads_ok = fftw_init_threads();
    cout << threads_ok << endl;
    fftw_mpi_init();

    if (threads_ok) fftw_plan_with_nthreads(1);

    bamba();
    return 0;
}
