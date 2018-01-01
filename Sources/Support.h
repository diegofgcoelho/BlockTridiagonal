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

/*
 * Structure used for storing the matrix data (DenseMatrix and SymMatrix) with a attached counter
 * that stores the number marices pointing to that location.
 */
typedef struct {
	//Array of actual data
	double* data;
	//Counter of all the matrices pointing to this structure data
	unsigned int counter;
} Data;

namespace support{
	/*
	 * This functon returns whether the first argument s larger or not than the second one
	 */
	inline bool cmpcmp(std::complex<double> a, std::complex<double> b){
		return std::abs(a) < std::abs(b);
	}

	/*
	 * This function implements the BLAS level I operation of scaling its vector argument by a.
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

	//Return the component with largest magnitude for array of complex numbers
	std::complex<double> max_mag_cmp(std::complex<double>* v, unsigned int vsize);
	//Return the component with largest magnitude for array of doubles
	double max_mag(double* v, unsigned int vsize);

	//Returns the variance
	inline double variance(double* v, unsigned int vsize, double mean){
		double var = 0.0;
		for(unsigned int i = 0; i < vsize; i++)	var = var + (v[i]-mean)*(v[i]-mean);
		var = std::sqrt(var); var /= vsize;
		return var;
	}


}

#endif /* SUPPORT_H_ */
