/*
 * SBlockPent_test.cpp
 *
 *  Created on: Dec 20, 2017
 *      Author: diego
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include "../Sources/SymMatrix.h"
#include "../Sources/DenseMatrix.h"
#include "../Sources/SBlockTrid.h"
#include "../Sources/SBlockPent.h"
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

namespace matrix {

void testPrintOperatorSBlockPent(){

	std::cout << "Please, check the two tex files generated for testing the << operator and compare them with meld for checking the matrices blocks." << std::endl;

	std::ofstream eigen_file("Eigen_SBlockPent_print_file.txt");
	std::ofstream sblockpent_file("SBlockPent_print_file.txt");

	eigen_file << "Testing << operator for SBlockPent\n";
	sblockpent_file << "Testing << operator for SBlockPent\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[2*nblocks-3];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-2; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);
		Eigen::MatrixXd C = Eigen::MatrixXd::Random(nsize, nsize);

		eigen_file << "Blocks of row " << k << " :\n" << A << "\nand\n" << B << "\nand\n" << C << "\n";

		//Assigning matrix elements
		double* dataA = new double[nsize*(nsize+1)/2];
		double* dataB = new double[nsize*nsize];
		double* dataC = new double[nsize*nsize];
		unsigned int datapos = 0;
		for(unsigned int i = 0; i < A.rows(); i++){
			for(unsigned int j = 0; j < A.cols(); j++){
				dataB[i*nsize+j] = B(i,j);
				dataC[i*nsize+j] = C(i,j);
				if(j >= i) {
					dataA[datapos] = A(i,j);
					datapos++;
				}
			}
		}

		//Setting the data to the array of DenseMatrix and SymMatrix
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[2*k].setData(dataB, nsize, nsize);
		denseBlocks[2*k+1].setData(dataC, nsize, nsize);
	}

	//The before the last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	eigen_file << "Blocks of row " << nblocks-2 << " :\n" << A << "\nand\n" << B << "\n";

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataB[i*nsize+j] = B(i,j);
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-2].setData(dataA, nsize);
	denseBlocks[2*(nblocks-2)].setData(dataB, nsize, nsize);

	//The last row in separate
	A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	eigen_file << "Blocks of row " << nblocks-1 << " :\n" << A << "\n";

	//Assigning matrix elements
	dataA = new double[nsize*(nsize+1)/2];
	datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-1].setData(dataA, nsize);

	SBlockPent Pent(symBlocks, denseBlocks, nblocks);
	sblockpent_file << "Pent created ...\n";

	//End of data setting. Testing the operator <<

	sblockpent_file << "The tested operator output:" << Pent << "\n";

	sblockpent_file << "End of testing for << operator for SBlockPent\n";
	eigen_file << "End of testing for << operator for SBlockPent\n";

	eigen_file.close();
	sblockpent_file.close();

}

void testSetBlocksSBlockPent(){

	std::cout << "Please, check the two tex files generated for testing the setBlocks method and compare them with meld for checking the matrices blocks." << std::endl;

	std::ofstream eigen_file("Eigen_SBlockPent_setBlocks_file.txt");
	std::ofstream sblockpent_file("SBlockPent_setBlocks_file.txt");

	eigen_file << "Testing setBlocks method for SBlockPent\n";
	sblockpent_file << "Testing setBlocks method for SBlockPent\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[2*nblocks-3];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-2; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);
		Eigen::MatrixXd C = Eigen::MatrixXd::Random(nsize, nsize);

		eigen_file << "Blocks of row " << k << " :\n" << A << "\nand\n" << B << "\nand\n" << C << "\n";

		//Assigning matrix elements
		double* dataA = new double[nsize*(nsize+1)/2];
		double* dataB = new double[nsize*nsize];
		double* dataC = new double[nsize*nsize];
		unsigned int datapos = 0;
		for(unsigned int i = 0; i < A.rows(); i++){
			for(unsigned int j = 0; j < A.cols(); j++){
				dataB[i*nsize+j] = B(i,j);
				dataC[i*nsize+j] = C(i,j);
				if(j >= i) {
					dataA[datapos] = A(i,j);
					datapos++;
				}
			}
		}

		//Setting the data to the array of DenseMatrix and SymMatrix
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[2*k].setData(dataB, nsize, nsize);
		denseBlocks[2*k+1].setData(dataC, nsize, nsize);
	}

	//The before the last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	eigen_file << "Blocks of row " << nblocks-2 << " :\n" << A << "\nand\n" << B << "\n";

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataB[i*nsize+j] = B(i,j);
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-2].setData(dataA, nsize);
	denseBlocks[2*(nblocks-2)].setData(dataB, nsize, nsize);

	//The last row in separate
	A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	eigen_file << "Blocks of row " << nblocks-1 << " :\n" << A << "\n";

	//Assigning matrix elements
	dataA = new double[nsize*(nsize+1)/2];
	datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-1].setData(dataA, nsize);

	SBlockPent Pent;
	Pent.setBlocks(symBlocks, denseBlocks, nblocks);
	sblockpent_file << "Pent created ...\n";

	//End of data setting. Testing the operator <<

	sblockpent_file << "The tested operator output:" << Pent << "\n";

	sblockpent_file << "End of testing for << operator for SBlockPent\n";
	eigen_file << "End of testing for << operator for SBlockPent\n";

	eigen_file.close();
	sblockpent_file.close();

}

void testAssignOperatorSBlockPent(){

	std::cout << "Please, check the two tex files generated for testing the = operator and compare them with meld for checking the matrices blocks." << std::endl;

	std::ofstream sblockpent_file("SBlockPent_assign_file_original.txt");
	std::ofstream sblockpent_file_b("SBlockPent_assign_file_copied.txt");

	sblockpent_file_b << "Testing = operator for SBlockPent\n";
	sblockpent_file << "Testing = operator for SBlockPent\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[2*nblocks-3];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-2; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);
		Eigen::MatrixXd C = Eigen::MatrixXd::Random(nsize, nsize);

		//Assigning matrix elements
		double* dataA = new double[nsize*(nsize+1)/2];
		double* dataB = new double[nsize*nsize];
		double* dataC = new double[nsize*nsize];
		unsigned int datapos = 0;
		for(unsigned int i = 0; i < A.rows(); i++){
			for(unsigned int j = 0; j < A.cols(); j++){
				dataB[i*nsize+j] = B(i,j);
				dataC[i*nsize+j] = C(i,j);
				if(j >= i) {
					dataA[datapos] = A(i,j);
					datapos++;
				}
			}
		}

		//Setting the data to the array of DenseMatrix and SymMatrix
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[2*k].setData(dataB, nsize, nsize);
		denseBlocks[2*k+1].setData(dataC, nsize, nsize);
	}

	//The before the last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataB[i*nsize+j] = B(i,j);
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-2].setData(dataA, nsize);
	denseBlocks[2*(nblocks-2)].setData(dataB, nsize, nsize);

	//The last row in separate
	A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Assigning matrix elements
	dataA = new double[nsize*(nsize+1)/2];
	datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-1].setData(dataA, nsize);

	SBlockPent Pent;
	Pent.setBlocks(symBlocks, denseBlocks, nblocks);
	sblockpent_file << "Pent created ...\n";
	sblockpent_file_b << "Pent created ...\n";

	SBlockPent Pentb;
	//End of data setting. Testing the operator <<

	sblockpent_file << "The tested operator output:" << Pent << "\n";

	//Operator under test
	Pentb = Pent;

	sblockpent_file_b << "The tested operator output:" << Pentb << "\n";

	sblockpent_file << "End of testing for << operator for SBlockPent\n";
	sblockpent_file_b << "End of testing for << operator for SBlockPent\n";

	sblockpent_file_b.close();
	sblockpent_file.close();

}

void testSymvPent(){

	std::cout << "*****************************Testing the Pent::symv method*****************************" << std::endl;

	const unsigned int nsize = 5;
	const unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[2*nblocks-3];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	Eigen::MatrixXd EigenPent = Eigen::MatrixXd::Zero(nsize*nblocks, nsize*nblocks);

	for(unsigned int k = 0; k < nblocks-2; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);
		Eigen::MatrixXd C = Eigen::MatrixXd::Random(nsize, nsize);

		//Setting the Eigen tridiagonal block matrix
		EigenPent.block<nsize, nsize>(nsize*k, nsize*k) = A;
		EigenPent.block<nsize, nsize>(nsize*k, nsize*(k+1)) = B;
		EigenPent.block<nsize, nsize>(nsize*k, nsize*(k+2)) = C;
		EigenPent.block<nsize, nsize>(nsize*(k+1), nsize*k) = B.transpose();
		EigenPent.block<nsize, nsize>(nsize*(k+2), nsize*k) = C.transpose();


		//Assigning matrix elements
		double* dataA = new double[nsize*(nsize+1)/2];
		double* dataB = new double[nsize*nsize];
		double* dataC = new double[nsize*nsize];
		unsigned int datapos = 0;
		for(unsigned int i = 0; i < A.rows(); i++){
			for(unsigned int j = 0; j < A.cols(); j++){
				dataB[i*nsize+j] = B(i,j);
				dataC[i*nsize+j] = C(i,j);
				if(j >= i) {
					dataA[datapos] = A(i,j);
					datapos++;
				}
			}
		}

		//Setting the data to the array of DenseMatrix and SymMatrix
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[2*k].setData(dataB, nsize, nsize);
		denseBlocks[2*k+1].setData(dataC, nsize, nsize);
	}

	//The before the last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	//Setting the Eigen pentadiagonal block matrix
	EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-2)) = A;
	EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-1)) = B;
	EigenPent.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-2)) = B.transpose();

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataB[i*nsize+j] = B(i,j);
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-2].setData(dataA, nsize);
	denseBlocks[2*(nblocks-2)].setData(dataB, nsize, nsize);

	//The last row in separate
	A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Setting the last block in the main diagonal
	EigenPent.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1)) = A;

	//Assigning matrix elements
	dataA = new double[nsize*(nsize+1)/2];
	datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-1].setData(dataA, nsize);

	SBlockPent Pent;
	Pent.setBlocks(symBlocks, denseBlocks, nblocks);


	//Generateing random vector
	Eigen::VectorXd eigen_vec_x = Eigen::VectorXd::Random(nsize*nblocks);
	Eigen::VectorXd eigen_vec_y = Eigen::VectorXd::Random(nsize*nblocks);
	double* vec_x = new double[nsize*nblocks]();
	double* vec_y = new double[nsize*nblocks]();

//	std::cout << " vec_x = [";
	for(unsigned int i = 0; i < nsize*nblocks; i++) {
		vec_x[i] = eigen_vec_x(i);
//		std::cout << "   " << vec_x[i];
	}
	std::cout << std::endl;
//	std::cout << "eigen_vec_x = [" << eigen_vec_x.transpose() << std::endl;

	//Matrix vector product
	eigen_vec_y = EigenPent*eigen_vec_x;
	Pent.symv(vec_x, vec_y);

//	std::cout << "EigenTrid\n" << EigenTrid << std::endl;


	for(unsigned int i = 0; i < nsize*nblocks; i++) {
		double derr = std::abs(vec_y[i]-eigen_vec_y(i));

		if(derr > 1e-2){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << std::endl;
		}

	}

	delete [] vec_x;
	delete [] vec_y;


	std::cout << "End of Trid::symv method test. If nothing was print on the screen, everything is fine." << std::endl;

}

void testSymvPentComplex(){

	std::cout << "*****************************Testing the Pent::symv method*****************************" << std::endl;

	const unsigned int nsize = 5;
	const unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[2*nblocks-3];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	Eigen::MatrixXd EigenPent = Eigen::MatrixXd::Zero(nsize*nblocks, nsize*nblocks);

	for(unsigned int k = 0; k < nblocks-2; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);
		Eigen::MatrixXd C = Eigen::MatrixXd::Random(nsize, nsize);

		//Setting the Eigen tridiagonal block matrix
		EigenPent.block<nsize, nsize>(nsize*k, nsize*k) = A;
		EigenPent.block<nsize, nsize>(nsize*k, nsize*(k+1)) = B;
		EigenPent.block<nsize, nsize>(nsize*k, nsize*(k+2)) = C;
		EigenPent.block<nsize, nsize>(nsize*(k+1), nsize*k) = B.transpose();
		EigenPent.block<nsize, nsize>(nsize*(k+2), nsize*k) = C.transpose();


		//Assigning matrix elements
		double* dataA = new double[nsize*(nsize+1)/2];
		double* dataB = new double[nsize*nsize];
		double* dataC = new double[nsize*nsize];
		unsigned int datapos = 0;
		for(unsigned int i = 0; i < A.rows(); i++){
			for(unsigned int j = 0; j < A.cols(); j++){
				dataB[i*nsize+j] = B(i,j);
				dataC[i*nsize+j] = C(i,j);
				if(j >= i) {
					dataA[datapos] = A(i,j);
					datapos++;
				}
			}
		}

		//Setting the data to the array of DenseMatrix and SymMatrix
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[2*k].setData(dataB, nsize, nsize);
		denseBlocks[2*k+1].setData(dataC, nsize, nsize);
	}

	//The before the last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	//Setting the Eigen pentadiagonal block matrix
	EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-2)) = A;
	EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-1)) = B;
	EigenPent.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-2)) = B.transpose();

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataB[i*nsize+j] = B(i,j);
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-2].setData(dataA, nsize);
	denseBlocks[2*(nblocks-2)].setData(dataB, nsize, nsize);

	//The last row in separate
	A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Setting the last block in the main diagonal
	EigenPent.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1)) = A;

	//Assigning matrix elements
	dataA = new double[nsize*(nsize+1)/2];
	datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-1].setData(dataA, nsize);

	SBlockPent Pent;
	Pent.setBlocks(symBlocks, denseBlocks, nblocks);


	//Generateing random vector
	Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(nsize*nblocks);
	Eigen::VectorXcd eigen_vec_y = Eigen::VectorXcd::Random(nsize*nblocks);
	std::complex<double>* vec_x = new std::complex<double>[nsize*nblocks]();
	std::complex<double>* vec_y = new std::complex<double>[nsize*nblocks]();

//	std::cout << " vec_x = [";
	for(unsigned int i = 0; i < nsize*nblocks; i++) {
		vec_x[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
//		std::cout << "   " << vec_x[i];
	}
	std::cout << std::endl;
//	std::cout << "eigen_vec_x = [" << eigen_vec_x.transpose() << std::endl;

	//Matrix vector product
	eigen_vec_y = EigenPent*eigen_vec_x;
	Pent.symv(vec_x, vec_y);

//	std::cout << "EigenTrid\n" << EigenTrid << std::endl;


	for(unsigned int i = 0; i < nsize*nblocks; i++) {
		double derr = std::abs(vec_y[i]-eigen_vec_y(i));

		if(derr > 1e-2){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << std::endl;
		}
	}

	delete [] vec_x;
	delete [] vec_y;


	std::cout << "End of Trid::symv method test. If nothing was print on the screen, everything is fine." << std::endl;

}

void testPowerMethodSBlockPent(){
	srand((unsigned int) time(0));

	const unsigned int nblocks = 2000;
	const unsigned int bsize = 10;

	//Lambda
	std::complex<double> lambda(1.0,0.0), elambda(1.0,0.0);
	//Precision
	double prec = 1e-3;
	//Maximum number of iterations
	unsigned int maxIter = 1000;
	//Number of actual iterations
	unsigned int iter = 0;
	//Norm of the difference between consecutive vectors
	double dnorm = 0.0;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks-1];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	Eigen::MatrixXd EigenTrid = Eigen::MatrixXd::Zero(bsize*nblocks, bsize*nblocks);

	for(unsigned int k = 0; k < nblocks-1; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(bsize, bsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(bsize, bsize);

		//Setting the Eigen tridiagonal block matrix
		EigenTrid.block<bsize, bsize>(bsize*k, bsize*k) = A;
		EigenTrid.block<bsize, bsize>(bsize*k, bsize*(k+1)) = B;
		EigenTrid.block<bsize, bsize>(bsize*(k+1), bsize*k) = B.transpose();

		//Assigning matrix elements
		double* dataA = new double[bsize*(bsize+1)/2];
		double* dataB = new double[bsize*bsize];
		unsigned int datapos = 0;
		for(unsigned int i = 0; i < A.rows(); i++){
			for(unsigned int j = 0; j < A.cols(); j++){
				dataB[i*bsize+j] = B(i,j);
				if(j >= i) {
					dataA[datapos] = A(i,j);
					datapos++;
				}
			}
		}

		//Setting the data to the array of DenseMatrix and SymMatrix
		symBlocks[k].setData(dataA, bsize);
		denseBlocks[k].setData(dataB, bsize, bsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(bsize, bsize); A = A*A.transpose();

	//Setting the Eigen tridiagonal block matrix
	EigenTrid.block<bsize, bsize>(bsize*(nblocks-1), bsize*(nblocks-1)) = A;

	//Assigning matrix elements
	double* dataA = new double[bsize*(bsize+1)/2];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				datapos++;
			}
		}
	}

	//Setting the data to the array of DenseMatrix and SymMatrix
	symBlocks[nblocks-1].setData(dataA, bsize);

	//Forming the tridiagonal block matrix
	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);
	SBlockPent Pent = Trid.square();

	//End of matrix generation. From this point on, all the matrices are set.

	//Generating random vector
	Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(bsize*nblocks);
	Eigen::VectorXcd eigen_vec_y = Eigen::VectorXcd::Random(bsize*nblocks);
	Eigen::VectorXcd eigen_vec_z = Eigen::VectorXcd::Random(bsize*nblocks);
	std::complex<double>* vec_x = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_y = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_z = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_w = new std::complex<double>[bsize*nblocks]();

	//Setting up the initial guesses
	for(unsigned int i = 0; i < bsize*nblocks; i++) {
		vec_x[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
		vec_y[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
	}

	//Executing power method
	Pent.powerMethod(vec_y, &lambda, prec, maxIter, &iter, &dnorm);

	lambda = std::sqrt(lambda);

	std::cout << "The returned largest eigenvalue (after sqrt) is: " << lambda << " with dnorm = " << dnorm << " and " << iter << " iterations." << std::endl;

	std::fill_n(vec_z, bsize*nblocks, std::complex<double>(0.0,0.0));
	Trid.symv(vec_y, vec_z);

	std::memcpy(vec_w, vec_y, sizeof(std::complex<double>)*bsize*nblocks);

	support::scal(vec_w, bsize*nblocks, lambda);

	if (bsize*nblocks < 500) {
		std::cout << "\n" << "The returned eigenvalue   vec_y[";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << vec_y[j] << "  ";
		}std::cout << "\n" << std::endl;
		std::cout << "\n" << "The lambda*v[";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << vec_w[j] << "  ";
		}std::cout << "\n" << std::endl;
		std::cout << "\n" << "The T*v     [";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << vec_z[j] << "  ";
		}std::cout << "\n" << std::endl;
	}

	for(unsigned int i = 0; i < nblocks*bsize; i++){
		vec_w[i] = vec_w[i]-vec_z[i];
	}

	std::cout << "\n\nThe norm of difference in SBlockTrid is l*v-Tx:" << support::norm(vec_w, bsize*nblocks)/support::norm(vec_y, bsize*nblocks) << std::endl;

	delete [] vec_x;
	delete [] vec_y;
	delete [] vec_w;
	delete [] vec_z;

}

void testPowerMethodSBlockPent2(){

	const unsigned int nblocks = 2000;
	const unsigned int bsize = 10;

	//Lambda
	std::complex<double> lambda(1.0,0.0);
	//Precision
	double prec = 1e-3;
	//Maximum number of iterations
	unsigned int maxIter = 1000;
	//Number of actual iterations
	unsigned int iter = 0;
	//Norm of the difference between consecutive vectors
	double dnorm = 0.0;


	//Forming the tridiagonal block matrix
	SBlockTrid Trid;
	Trid.random(nblocks, bsize);
	SBlockPent Pent = Trid.square();

	//Generating random vector
	Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(bsize*nblocks);
	std::complex<double>* vec_y = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_z = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_w = new std::complex<double>[bsize*nblocks]();

	//Setting up the initial guesses
	for(unsigned int i = 0; i < bsize*nblocks; i++) {
		vec_y[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
	}

	//Executing power method
	Pent.powerMethod(vec_y, &lambda, prec, maxIter, &iter, &dnorm);

	lambda = std::sqrt(lambda);

	std::cout << "The returned largest eigenvalue (after sqrt) is: " << lambda << " with dnorm = " << dnorm << " and " << iter << " iterations." << std::endl;

	std::fill_n(vec_z, bsize*nblocks, std::complex<double>(0.0,0.0));
	Trid.symv(vec_y, vec_z);

	std::memcpy(vec_w, vec_y, sizeof(std::complex<double>)*bsize*nblocks);

	support::scal(vec_w, bsize*nblocks, lambda);

	if (bsize*nblocks < 500) {
		std::cout << "\n" << "The returned eigenvalue   vec_y[";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << vec_y[j] << "  ";
		}std::cout << "\n" << std::endl;
		std::cout << "\n" << "The lambda*v[";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << vec_w[j] << "  ";
		}std::cout << "\n" << std::endl;
		std::cout << "\n" << "The T*v     [";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << vec_z[j] << "  ";
		}std::cout << "\n" << std::endl;
	}

	for(unsigned int i = 0; i < nblocks*bsize; i++){
		vec_w[i] = vec_w[i]-vec_z[i];
	}

	std::cout << "\n\nThe norm of difference in SBlockTrid is l*v-Tx:" << support::norm(vec_w, bsize*nblocks)/support::norm(vec_y, bsize*nblocks) << std::endl;

	delete [] vec_y;
	delete [] vec_w;
	delete [] vec_z;

}


} /* namespace matrix */
