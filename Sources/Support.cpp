/*
 * Support.cpp
 *
 *  Created on: Dec 28, 2017
 *      Author: diego
 */

#include "Support.h"

namespace support {

bool check_stop(std::complex<double>* v, std::complex<double>* vv, std::complex<double>* vvv, unsigned int vsize, double prec, double *rnormd){
	/*Input:
	 * v is a vector
	 * vv is a vector
	 * vvv is a vector
	 * vsize is the vector size
	 * prec is a double
	 * rnormd is a pointer for double
	 */
	/*Output:
	 *flag is a bool
	 */
	/*Requirement:
	 * v and vv must be non-empty vectors and prec must be greater than 0
	 */
	/*Description:
	 * The function compare vv and vvv against v (the current eigenvector estimate) and
	 * determines if the power method have achieved convergence
	 */

	//Output variable
	bool flagstop;

	//Auxiliary variables representing the norms of v, v-vv and v-vvv, respectively
	double normv = 0.0, normvvd = 0.0, normvvvd = 0.0;

	for(unsigned i = 0; i < vsize; i++){
		std::complex<double> tempv = v[i];
		std::complex<double> tempvv = vv[i];
		std::complex<double> tempvvv = vvv[i];

		normv += std::norm(tempv);
		normvvd += std::norm(tempv-tempvv);
		normvvvd += std::norm(tempv-tempvvv);
	}

	normv = std::sqrt(normv);
	normvvd = std::sqrt(normvvd);
	normvvvd = std::sqrt(normvvvd);

	//Assigning the value of the output variable
	*rnormd = std::min(normvvd/normv, normvvvd/normv);

	(*rnormd <= prec)?flagstop=true:flagstop=false;

	return flagstop;
}

double max_mag(double* v, unsigned int vsize){
	/*Input:
	 * v is a vector
	 */
	/*Output:
	 * maxv is a double
	 */
	/*Requirement:
	 * v must be a non-empty vector
	 */
	/*Description: this function returns the absolute value component of
	 * the largest component
	 */
	//Value to be passed to the output
	double maxv = 0.0;

	for(unsigned i = 0; i < vsize; i++){
		//If v[i] is greater than maxv, update maxv
		if(std::abs(v[i]) > std::abs(maxv)){
			maxv = v[i];
		}
	}

	return maxv;
}


std::complex<double> max_mag_cmp(std::complex<double>* v, unsigned int vsize){
	/*Input:
	 * v is a complex vector
	 */
	/*Output:
	 * maxv is a complex
	 */
	/*Requirement:
	 * v must be a non-empty vector
	 */
	/*Description: this function returns the absolute value component of
	 * the largest component
	 */
	//Value to be passed to the output
	std::complex<double> maxv(0.0,0.0);

	for(unsigned i = 0; i < vsize; i++){
		//If v[i] is greater than maxv, update maxv
		if(std::norm(v[i]) > std::norm(maxv)){
			maxv = v[i];
		}
	}
	return maxv;
}



}  // namespace support


