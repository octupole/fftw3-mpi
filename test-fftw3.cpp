//============================================================================
// Name        : test-fftw3.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <random>

#include <mpi.h>
#include "fftw3-mpi.h"

#include "fftw3pp/Array.h"
#include "fftw3pp/Pcrfft3d.h"
#include "fftw3pp/Prcfft3d.h"



using namespace std;
using namespace Array;

double Wtime(){
	timeval tim;
	double temp=tim.tv_sec+(tim.tv_usec/1000000.0);
	return temp;
	}


int bamby(){
	int rank,size;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);

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

     alloc_local = fftw_mpi_local_size_3d(L, M, N/2+1, MPI_COMM_WORLD,
                                            &local_n0, &local_0_start);


    c_out.Allocate(local_n0,ny,nzp);
    rin.Allocate(local_n0,ny,nzp*2);

	fftwpp::Prcfft3d Forward3(nx,ny,nz,rin,c_out);

    fftwpp::Pcrfft3d Backward3(nx,ny,nz,c_out,rin);

    for (size_t i = 0; i < local_n0; ++i)
       for (size_t j = 0; j < M; ++j)
         for (size_t k = 0; k < N; ++k)
        	 rin[i][j][k] = ri[local_0_start+i][j][k];

    /* compute transforms as many times as desired */

    double time0=MPI_Wtime();
    for(size_t w{0};w<1;w++){
    	Forward3.fft(rin,c_out);

    	Backward3.fftNormalized(c_out,rin);
    }
    double time1=MPI_Wtime();

    if(rank==0) cout << time1-time0<<endl;
    return 0;
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
    fftw_mpi_init();
    bamby();
    MPI_Finalize();
    return 0;
}
