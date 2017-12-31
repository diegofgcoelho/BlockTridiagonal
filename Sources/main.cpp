/*
 * main.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#include "DenseMatrix.h"
#include "SymMatrix.h"
#include "SBlockTrid.h"
#include "SBlockPent.h"
#include <omp.h>

int main(int agrc, char* argv[]){
//	std::cout << "*************************Testing the methods for DenseMatrix************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	matrix::testPrintOperatorDenseMatrix();
//	matrix::testAssignOperatorDenseMatrix();
//	matrix::testCloneToDenseMatrix();
//	matrix::testCloneFromDenseMatrix();
//	matrix::testMultiplyOperatorDenseMatrixXDenseMatrix();
//	matrix::testMultiplyOperatorDenseMatrixXSymMatrix();
//	matrix::testDenseMatrixXVector();
//	matrix::testDenseMatrixXComplexVector();
//	matrix::testDenseMatrixTransposedXComplexVector();
//	matrix::testDenseMatrixTransposedXVector();
//	matrix::testAAT();
//	matrix::testATA();
//	matrix::testAddOperatorDenseMatrixXDenseMatrix();
//	matrix::testCloneTransposeToDenseMatrix();
//	matrix::testCloneTransposeFromDenseMatrix();
//	std::cout << "*************************Testing the methods for SymMatrix************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	matrix::testPrintOperatorSymMatrix();
//	matrix::testAssignOperatorSymMatrix();
//	matrix::testCloneToSymMatrix();
//	matrix::testCloneFromSymMatrix();
	matrix::testMultiplyOperatorSymMatrixXDenseMatrix();
//	matrix::testMultiplyOperatorSymMatrixXSymMatrix();
//	matrix::testSymMatrixXVector();
//	matrix::testSymMatrixXComplexVector();
//	matrix::testAddOperatorSymMatrixXSymMatrix();
//	matrix::testSquareSymMatrix();
//	std::cout << "*************************Testing the methods for SBlockTrid************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	matrix::testPrintOperatorSBlockTrid();
//	matrix::testSetBlocksSBlockTrid();
//	matrix::testAssignOperatorSBlockTrid();
//	matrix::testSymMatrixXVectorPlusDenseMatrixXVector();
//	matrix::testSymvTrid();
//	matrix::testSymvTridComplex();
//	matrix::testSymvTridComplexAndScale();
//	matrix::testPowerMethodTrid();
//	matrix::testPowerMethodTrid2();
//	std::cout << "*************************Testing the methods for SBlockPent************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	std::cout << "************************************************************************************" << std::endl;
//	matrix::testPrintOperatorSBlockPent();
//	matrix::testSetBlocksSBlockPent();
//	matrix::testAssignOperatorSBlockPent();
//	matrix::testSquareSBlockTridMatlab();
//	matrix::testSquareSBlockTridEigen();
//	matrix::testSymvPent()
//	matrix::testSymvPentComplex();
//	matrix::testPowerMethodSBlockPent();
//	matrix::testPowerMethodSBlockPent2();


//	//Setting the number of blocks on the main diagonal and their size
//	const unsigned int nblocks = 2000;
//	const unsigned int bsize = 5;
//
//	//Lambda
//	std::complex<double> lambdaPent(1.0,0.0), lambdaTrid(1.0,0.0);
//	//Precision
//	double prec = 1e-3;
//	//Maximum number of iterations
//	unsigned int maxIter = 1000;
//	//Number of actual iterations
//	unsigned int iterPent = 0, iterTrid = 0;
//	//Norm of the difference between consecutive vectors
//	double dnormPent = 0.0, dnormTrid = 0.0;
//
//
//	//Forming the tridiagonal block matrix
//	matrix::SBlockTrid Trid;
//	Trid.random(nblocks, bsize);
//	matrix::SBlockPent Pent = Trid.square();
//
//	//Generating random vector
//	Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(bsize*nblocks);
//	std::complex<double>* vpent = new std::complex<double>[bsize*nblocks]();
//	std::complex<double>* Tvpent = new std::complex<double>[bsize*nblocks]();
//	std::complex<double>* avpent = new std::complex<double>[bsize*nblocks]();
//
//	std::complex<double>* vtrid = new std::complex<double>[bsize*nblocks]();
//	std::complex<double>* Tvtrid = new std::complex<double>[bsize*nblocks]();
//	std::complex<double>* avtrid = new std::complex<double>[bsize*nblocks]();
//
//
//	//Setting up the initial guesses
//	for(unsigned int i = 0; i < bsize*nblocks; i++) {
//		vpent[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
//		vtrid[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
//	}
//
//	double time_pent = omp_get_wtime();
//	//Executing power method
//	Pent.powerMethod(vpent, &lambdaPent, prec, maxIter, &iterPent, &dnormPent);
//	time_pent = omp_get_wtime()-time_pent;
//	lambdaPent = std::sqrt(lambdaPent);
//
//	std::cout << "Pent: The returned largest eigenvalue (after sqrt) is: " << lambdaPent << " with dnorm = " << dnormPent << " and " << iterPent << " iterations." << std::endl;
//
//	std::fill_n(Tvpent, bsize*nblocks, std::complex<double>(0.0,0.0));
//	Trid.symv(vpent, Tvpent);
//
//	std::memcpy(avpent, vpent, sizeof(std::complex<double>)*bsize*nblocks);
//
//	support::scal(avpent, bsize*nblocks, lambdaPent);
//
//	if (bsize*nblocks < 500) {
//		std::cout << "\n" << "The returned eigenvalue   vec_y[";
//		for(unsigned int j = 0; j < nblocks*bsize; j++){
//			std::cout  << vpent[j] << "  ";
//		}std::cout << "\n" << std::endl;
//		std::cout << "\n" << "The lambda*v[";
//		for(unsigned int j = 0; j < nblocks*bsize; j++){
//			std::cout  << avpent[j] << "  ";
//		}std::cout << "\n" << std::endl;
//		std::cout << "\n" << "The T*v     [";
//		for(unsigned int j = 0; j < nblocks*bsize; j++){
//			std::cout  << Tvpent[j] << "  ";
//		}std::cout << "\n" << std::endl;
//	}
//
//	for(unsigned int i = 0; i < nblocks*bsize; i++){
//		avpent[i] = avpent[i]-Tvpent[i];
//	}
//
//	std::cout << "\n\nThe norm of difference in SBlockPent is l*v-Tx:" << support::norm(avpent, bsize*nblocks)/support::norm(vpent, bsize*nblocks) << std::endl;
//
//	double time_trid = omp_get_wtime();
//	//Executing power method
//	Trid.powerMethod(vtrid, &lambdaTrid, prec, maxIter, &iterTrid, &dnormTrid);
//	time_trid = omp_get_wtime()-time_trid;
//
//	std::cout << "\n\nTrid: The returned largest eigenvalue is: " << lambdaTrid << " with dnorm = " << dnormTrid << " and " << iterTrid << " iterations." << std::endl;
//
//	std::fill_n(Tvtrid, bsize*nblocks, std::complex<double>(0.0,0.0));
//	Trid.symv(vtrid, Tvtrid);
//
//	std::memcpy(avtrid, vtrid, sizeof(std::complex<double>)*bsize*nblocks);
//
//	support::scal(avtrid, bsize*nblocks, lambdaTrid);
//
//	if (bsize*nblocks < 500) {
//		std::cout << "\n" << "The returned eigenvalue   vec_y[";
//		for(unsigned int j = 0; j < nblocks*bsize; j++){
//			std::cout  << vtrid[j] << "  ";
//		}std::cout << "\n" << std::endl;
//		std::cout << "\n" << "The lambda*v[";
//		for(unsigned int j = 0; j < nblocks*bsize; j++){
//			std::cout  << avtrid[j] << "  ";
//		}std::cout << "\n" << std::endl;
//		std::cout << "\n" << "The T*v     [";
//		for(unsigned int j = 0; j < nblocks*bsize; j++){
//			std::cout  << Tvtrid[j] << "  ";
//		}std::cout << "\n" << std::endl;
//	}
//
//	for(unsigned int i = 0; i < nblocks*bsize; i++){
//		avtrid[i] = avtrid[i]-Tvtrid[i];
//	}
//
//	std::cout << "\n\nThe norm of difference in SBlockTrid is l*v-Tx:" << support::norm(avtrid, bsize*nblocks)/support::norm(vtrid, bsize*nblocks) << std::endl;
//
//
//	std::cout << "Time for Pentadiagonal is " << time_pent << " and for Tridiagonal is " << time_trid << std::endl;
//
//
//	delete [] vpent;
//	delete [] avpent;
//	delete [] Tvpent;
//
//	delete [] vtrid;
//	delete [] avtrid;
//	delete [] Tvtrid;

}


