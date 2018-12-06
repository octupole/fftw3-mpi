/*
 * Fftw3pp.h
 *
 *  Created on: Nov 28, 2018
 *      Author: marchi
 */

#ifndef FFTW3PP_FFTW3PP_H_
#define FFTW3PP_FFTW3PP_H_
#ifdef HAVE_FFTW_MPI
#include <fftw3-mpi.h>
#else
#include <fftw3.h>
#endif
namespace fftw3pp {

class Fftw3pp {
public:
	Fftw3pp();
	~Fftw3pp(){};
};

} /* namespace fftw3pp */

#endif /* FFTW3PP_FFTW3PP_H_ */
