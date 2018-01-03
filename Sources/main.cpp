/*
 * main.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#include <numeric>
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
//	matrix::testMultiplyOperatorSymMatrixXDenseMatrix();
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
	matrix::testPrintOperatorSBlockPent();
	matrix::testSetBlocksSBlockPent();
	matrix::testAssignOperatorSBlockPent();
	matrix::testSquareSBlockTridMatlab();
	matrix::testSquareSBlockTridEigen();
	matrix::testSymvPent();
	matrix::testSymvPentComplex();
	matrix::testPowerMethodSBlockPent();
	matrix::testPowerMethodSBlockPent2();


/*	//Setting the number of blocks on the main diagonal and their size
	const unsigned int nblocks_v[] = {100, 500, 1000, 2000};
	const unsigned int nblocks_vn = 1;
	const unsigned int bsize_v[] = {5, 10, 20};
	const unsigned int bsize_vn = 2;

	//Number of replicates for each block size and number of blocks in the main diagonal
	const unsigned int total_reps = 100;

	//Times and error arrays
	double* ptimes_pent = new double[total_reps]();
	double* ptimes_trid = new double[total_reps]();
	double* perror_pent = new double[total_reps]();
	double* perror_trid = new double[total_reps]();


	//Table (Table) that will store all the metrics
	//Col1, Col2 ....
	//Number of blocks, size of each block, avg time(trid), time variance (trid), error (trid), avg time(pent), time variance (pent), error (pent)
	double* metrics_data = new double[nblocks_vn*bsize_vn*8]();
	matrix::DenseMatrix metrics(metrics_data, nblocks_vn*bsize_vn, 8);

	for (unsigned int nbs = 0; nbs < nblocks_vn; nbs++) {
		for (unsigned int bs = 0; bs < bsize_vn; bs++) {
			//Cleaning the array storing the time and error measures
			std::fill_n(ptimes_pent, total_reps, 0.0);
			std::fill_n(ptimes_trid, total_reps, 0.0);
			std::fill_n(perror_pent, total_reps, 0.0);
			std::fill_n(perror_trid, total_reps, 0.0);

			//Setting sizes
			unsigned int nblocks = nblocks_v[nbs];
			unsigned int bsize = bsize_v[bs];

			unsigned int metric_row = nbs*nblocks_vn+bs;
			metrics.setEle(metric_row, 0, nblocks);//Number of blocks
			metrics.setEle(metric_row, 1, bsize);//Size of blocks

			std::cout << "Starting the simulation for matrices with " << nblocks << " blocks on the main diagonal, each with size " << bsize << "." << std::endl;

			for (unsigned int reps = 0; reps < total_reps; reps++) {
				std::cout << "Reps = " << reps << std::endl;
				//Lambda
				std::complex<double> lambdaPent(1.0,0.0), lambdaTrid(1.0,0.0);
				//Precision
				double prec = 1e-3;
				//Maximum number of iterations
				unsigned int maxIter = 1000;
				//Number of actual iterations
				unsigned int iterPent = 0, iterTrid = 0;
				//Norm of the difference between consecutive vectors
				double dnormPent = 0.0, dnormTrid = 0.0;

				//Forming the tridiagonal block matrix
				matrix::SBlockTrid Trid;
				Trid.random(nblocks, bsize);
				matrix::SBlockPent Pent = Trid.square();

				//Generating random vector
				Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(bsize*nblocks);
				std::complex<double>* vpent = new std::complex<double>[bsize*nblocks]();
				std::complex<double>* Tvpent = new std::complex<double>[bsize*nblocks]();
				std::complex<double>* avpent = new std::complex<double>[bsize*nblocks]();

				std::complex<double>* vtrid = new std::complex<double>[bsize*nblocks]();
				std::complex<double>* Tvtrid = new std::complex<double>[bsize*nblocks]();
				std::complex<double>* avtrid = new std::complex<double>[bsize*nblocks]();


				//Setting up the initial guesses
				for(unsigned int i = 0; i < bsize*nblocks; i++) {
					vpent[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
					vtrid[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
				}

				ptimes_pent[reps] = omp_get_wtime();
				//Executing power method
				Pent.powerMethod(vpent, &lambdaPent, prec, maxIter, &iterPent, &dnormPent);
				ptimes_pent[reps] = omp_get_wtime()- ptimes_pent[reps];
				lambdaPent = std::sqrt(lambdaPent);

				std::fill_n(Tvpent, bsize*nblocks, std::complex<double>(0.0,0.0));
				Trid.symv(vpent, Tvpent);

				std::memcpy(avpent, vpent, sizeof(std::complex<double>)*bsize*nblocks);

				support::scal(avpent, bsize*nblocks, lambdaPent);

				if (bsize*nblocks < 500) {
					std::cout << "\n" << "The returned eigenvalue   vec_y[";
					for(unsigned int j = 0; j < nblocks*bsize; j++){
						std::cout  << vpent[j] << "  ";
					}std::cout << "\n" << std::endl;
					std::cout << "\n" << "The lambda*v[";
					for(unsigned int j = 0; j < nblocks*bsize; j++){
						std::cout  << avpent[j] << "  ";
					}std::cout << "\n" << std::endl;
					std::cout << "\n" << "The T*v     [";
					for(unsigned int j = 0; j < nblocks*bsize; j++){
						std::cout  << Tvpent[j] << "  ";
					}std::cout << "\n" << std::endl;
				}

				for(unsigned int i = 0; i < nblocks*bsize; i++){
					avpent[i] = avpent[i]-Tvpent[i];
				}
				//Error for Pent
				perror_pent[reps] = support::norm(avpent, bsize*nblocks)/support::norm(vpent, bsize*nblocks);

				ptimes_trid[reps] = omp_get_wtime();
				//Executing power method
				Trid.powerMethod(vtrid, &lambdaTrid, prec, maxIter, &iterTrid, &dnormTrid);
				ptimes_trid[reps] = omp_get_wtime()- ptimes_trid[reps];

				std::fill_n(Tvtrid, bsize*nblocks, std::complex<double>(0.0,0.0));
				Trid.symv(vtrid, Tvtrid);

				std::memcpy(avtrid, vtrid, sizeof(std::complex<double>)*bsize*nblocks);

				support::scal(avtrid, bsize*nblocks, lambdaTrid);

				if (bsize*nblocks < 500) {
					std::cout << "\n" << "The returned eigenvalue   vec_y[";
					for(unsigned int j = 0; j < nblocks*bsize; j++){
						std::cout  << vtrid[j] << "  ";
					}std::cout << "\n" << std::endl;
					std::cout << "\n" << "The lambda*v[";
					for(unsigned int j = 0; j < nblocks*bsize; j++){
						std::cout  << avtrid[j] << "  ";
					}std::cout << "\n" << std::endl;
					std::cout << "\n" << "The T*v     [";
					for(unsigned int j = 0; j < nblocks*bsize; j++){
						std::cout  << Tvtrid[j] << "  ";
					}std::cout << "\n" << std::endl;
				}

				for(unsigned int i = 0; i < nblocks*bsize; i++){
					avtrid[i] = avtrid[i]-Tvtrid[i];
				}
				//Error for Trid
				perror_trid[reps] = support::norm(avtrid, bsize*nblocks)/support::norm(vtrid, bsize*nblocks);

				delete [] vpent;
				delete [] avpent;
				delete [] Tvpent;

				delete [] vtrid;
				delete [] avtrid;
				delete [] Tvtrid;

				std::cout << "Reps = " << reps << std::endl;

			}

			std::cout << "Finished the simulation for matrices with " << nblocks << " blocks on the main diagonal, each with size " << bsize << "." << std::endl;
			//Setting metrics table values
			double avg_time_trid = std::accumulate(&ptimes_trid[0],&ptimes_trid[total_reps-1], 0.0)/total_reps;
			metrics.setEle(metric_row, 2, avg_time_trid);//Avg for trid
			metrics.setEle(metric_row, 3, support::variance(&ptimes_trid[0], total_reps, avg_time_trid));//Var for trid
			double avg_error_trid = std::accumulate(&perror_trid[0],&perror_trid[total_reps-1], 0.0)/total_reps;
			metrics.setEle(metric_row, 4, avg_error_trid);//Error for trid

			double avg_time_pent = std::accumulate(&ptimes_pent[0],&ptimes_pent[total_reps-1], 0.0)/total_reps;
			metrics.setEle(metric_row, 5, avg_time_pent);//Avg for pent
			metrics.setEle(metric_row, 6, support::variance(&ptimes_pent[0], total_reps, avg_time_pent));//Var for pent
			double avg_error_pent = std::accumulate(&perror_pent[0],&perror_pent[total_reps-1], 0.0)/total_reps;
			metrics.setEle(metric_row, 7, avg_error_pent);//Error for pent
		}
	}

	delete [] ptimes_pent;
	delete [] ptimes_trid;

	std::cout << "The metrics are:\n" << metrics << std::endl;*/
}


