/*
 * Support.h
 *
 *  Created on: Dec 28, 2017
 *      Author: diego
 */

#ifndef SUPPORT_H_
#define SUPPORT_H_

#include <cmath>
#include <complex>

namespace support{
	/*
	 * This functon returns whether the first argument s larger or not than the second one
	 */
	inline bool cmpcmp(std::complex<double> a, std::complex<double> b){
		return std::abs(a) < std::abs(b);
	}

	/*
	 * This funciton implements the BLAS level I operation of scaling its vector argument by a.
	 */
	inline void scal(std::complex<double>* v, unsigned int vsize, double a){
		for(unsigned int i = 0; i < vsize; i++) v[i] *= a;
	}

	inline void scal(std::complex<double>* v, unsigned int vsize, std::complex<double> a){
		for(unsigned int i = 0; i < vsize; i++) v[i] *= a;
	}

	inline void scal(double* v, unsigned int vsize, double a){
		for(unsigned int i = 0; i < vsize; i++) v[i] *= a;
	}

	inline double norm(std::complex<double>* v, unsigned int vsize){
		double norm = 0.0;
		for(unsigned int i = 0; i < vsize; i++) norm+=std::norm(v[i]);
		return std::sqrt(norm);
	}

	inline double norm(double* v, unsigned int vsize){
		double norm = 0.0;
		for(unsigned int i = 0; i < vsize; i++) norm+=(v[i]*v[i]);
		return std::sqrt(norm);
	}

	bool check_stop(std::complex<double>* v, std::complex<double>* vv, std::complex<double>* vvv, unsigned int vsize, double prec, double *rnormd);

	std::complex<double> max_mag_cmp(std::complex<double>* v, unsigned int vsize);
	double max_mag(double* v, unsigned int vsize);
}

#endif /* SUPPORT_H_ */
