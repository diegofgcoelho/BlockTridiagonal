/*
 * SBlockTrid_test.cpp
 *
 *  Created on: Dec 20, 2017
 *      Author: diego
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include "../Sources/SymMatrix.h"
#include "../Sources/DenseMatrix.h"
#include "../Sources/SBlockTrid.h"
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

namespace matrix {

void testPrintOperatorSBlockTrid(){

	std::cout << "Please, check the two tex files generated for testing the << operator and compare them with meld for checking the matrices blocks." << std::endl;

	std::ofstream eigen_file("Eigen_print_file.txt");
	std::ofstream sblocktrid_file("SBlockTrid_print_file.txt");

	eigen_file << "Testing << operator for SBlockTrid\n";
	sblocktrid_file << "Testing << operator for SBlockTrid\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-1; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

		eigen_file << "Blocks of row " << k << " :\n" << A << "\nand\n" << B << "\n";

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
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[k].setData(dataB, nsize, nsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	eigen_file << "Blocks of row " << nblocks-1 << " :\n" << A << "\n";

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
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
	symBlocks[nblocks-1].setData(dataA, nsize);


	SBlockTrid Trid(symBlocks, denseBlocks, nblocks);
	sblocktrid_file << "Trid created ...\n";

	//End of data setting. Testing the operator <<

	sblocktrid_file << "The tested operator output:" << Trid << "\n";

	sblocktrid_file << "End of testing for << operator for SBlockTrid\n";
	eigen_file << "End of testing for << operator for SBlockTrid\n";

	eigen_file.close();
	sblocktrid_file.close();
}

void testSetBlocksSBlockTrid(){

	std::cout << "Please, check the two tex files generated for testing the setBlocks method and compare them with meld for checking the matrices blocks." << std::endl;

	std::ofstream eigen_file("Eigen_setBlocks_file.txt");
	std::ofstream sblocktrid_file("SBlockTrid_setBlocks_file.txt");


	eigen_file << "Testing setBlocks method for SBlockTrid\n";
	sblocktrid_file << "Testing setBlocks method for SBlockTrid\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-1; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

		eigen_file << "Blocks of row " << k << " :\n" << A << "\nand\n" << B << "\n";

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
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[k].setData(dataB, nsize, nsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	eigen_file << "Blocks of row " << nblocks-1 << " :\n" << A <<  "\n";

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
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
	symBlocks[nblocks-1].setData(dataA, nsize);


	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);
	sblocktrid_file << "Trid created ...\n";

	//End of data setting. Testing the operator <<

	sblocktrid_file << "Stored matrix:" << Trid <<  "\n";

	eigen_file << "End of test for setBlocks method for SBlockTrid";
	sblocktrid_file << "End of test for setBlocks method for SBlockTrid";

	eigen_file.close();
	sblocktrid_file.close();
}

void testAssignOperatorSBlockTrid(){

	std::cout << "Please, check the two tex files generated for testing the = operator method and compare them with meld for checking the matrices blocks." << std::endl;

	std::ofstream sblocktrid_file_b("SBlockTrid_assign_file_b.txt");
	std::ofstream sblocktrid_file("SBlockTrid_assign_file.txt");


	sblocktrid_file_b << "Testing setBlocks method for SBlockTrid\n";
	sblocktrid_file << "Testing setBlocks method for SBlockTrid\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-1; k++){
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
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[k].setData(dataB, nsize, nsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
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
	symBlocks[nblocks-1].setData(dataA, nsize);


	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);
	sblocktrid_file << "Trid created ...\n";
	sblocktrid_file_b << "Trid created ...\n";

	//End of data setting. Testing the operator <<

	sblocktrid_file << "Stored matrix:" << Trid <<  "\n";

	SBlockTrid Trid_b;
	Trid_b = Trid;

	sblocktrid_file_b << "Stored matrix:" << Trid_b <<  "\n";

	sblocktrid_file_b << "End of test for setBlocks method for SBlockTrid";
	sblocktrid_file << "End of test for setBlocks method for SBlockTrid";

	sblocktrid_file_b.close();
	sblocktrid_file.close();
}

void testSquareSBlockTridMatlab(){
	std::cout << "Please, check the text file generated for testing the square method and compare them matlab output (check the Matlab folder)." << std::endl;


	std::ofstream sblocktrid_file("SBlockTrid_square_Matlab_file.txt");


	sblocktrid_file << "The original SBlockTrid matrix is\n";

	unsigned int nsize = 5;
	unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-1; k++){
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
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[k].setData(dataB, nsize, nsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
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
	symBlocks[nblocks-1].setData(dataA, nsize);


	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);

	sblocktrid_file << Trid <<  "\n";

	SBlockPent Pent;
	Pent = Trid.square();

	sblocktrid_file << "The resulting SBlockPent is:\n";
	sblocktrid_file << Pent;

	sblocktrid_file.close();

}

void testSquareSBlockTridEigen(){
	std::cout << "Please, check the two text file generated for testing the square method and compare them using meld." << std::endl;


	std::ofstream sblocktrid_file("SBlockTrid_square_file.txt");
	std::ofstream eigen_file("Eigen_square_file.txt");


	sblocktrid_file << "The original SBlockTrid matrix is\n";
	eigen_file << "The original SBlockTrid matrix is\n";

	const unsigned int nsize = 5;
	const unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	Eigen::MatrixXd EigenTrid = Eigen::MatrixXd::Zero(nsize*nblocks, nsize*nblocks);

	for(unsigned int k = 0; k < nblocks-1; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

		//Setting the Eigen tridiagonal block matrix
		EigenTrid.block<nsize, nsize>(nsize*k, nsize*k) = A;
		EigenTrid.block<nsize, nsize>(nsize*k, nsize*(k+1)) = B;
		EigenTrid.block<nsize, nsize>(nsize*(k+1), nsize*k) = B.transpose();

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
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[k].setData(dataB, nsize, nsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Setting the Eigen tridiagonal block matrix
	EigenTrid.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1)) = A;

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
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
	symBlocks[nblocks-1].setData(dataA, nsize);


	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);

	sblocktrid_file << Trid <<  "\n";
	eigen_file << EigenTrid << "\n";

	SBlockPent Pent;
	Pent = Trid.square();

	Eigen::MatrixXd EigenPent = EigenTrid*EigenTrid;

	//Compute the norm of the difference of the resulting matrices
	{
		SymMatrix* symPent = NULL;
		DenseMatrix* denPent = NULL;
		Pent.getBlocks(&symPent, &denPent);
		double derr = 0.0;
		for(unsigned int k = 0; k < nblocks-2; k++){

			std::cout << "k = " << k << std::endl;

			for (unsigned int i = 0; i < nsize; i++) {
				for (unsigned int j = 0; j < nsize; j++) {
					//Getting the Eigen pentadiagonal block matrix
					derr = std::fabs(symPent[k].getEle(i,j) - EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j));
					if(derr > 1e-2){
						std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
						std::cout << "symPent[k].getEle(i,j) = " << symPent[k].getEle(i,j) << " EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) << std::endl;
					}

					derr = std::fabs(denPent[2*k].getEle(i,j) - EigenPent.block<nsize, nsize>(nsize*k, nsize*(k+1))(i,j));

					if(derr > 1e-2){
						std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
						std::cout << "denPent[2*k].getEle(i,j) = " << denPent[2*k].getEle(i,j) << " EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) << std::endl;
					}


					derr = std::fabs(denPent[2*k+1].getEle(i,j) - EigenPent.block<nsize, nsize>(nsize*k, nsize*(k+2))(i,j));

					if(derr > 1e-2){
						std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
						std::cout << "denPent[2*k].getEle(i,j) = " << denPent[2*k].getEle(i,j) << " EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) << std::endl;
					}

				}
			}

		}

		//k = nblocks-2
		std::cout << "k = " << nblocks-2 << std::endl;
		for (unsigned int i = 0; i < nsize; i++) {
			for (unsigned int j = 0; j < nsize; j++) {
				//Getting the Eigen pentadiagonal block matrix
				derr = std::fabs(symPent[nblocks-2].getEle(i,j) - EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-2))(i,j));

				if(derr > 1e-2){
					std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
					std::cout << "symPent[k].getEle(i,j) = " << symPent[nblocks-2].getEle(i,j) << " EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-2))(i,j) << std::endl;
				}


				derr = std::fabs(denPent[2*(nblocks-2)].getEle(i,j) - EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-1))(i,j));

				if(derr > 1e-2){
					std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
					std::cout << "denPent[2*k].getEle(i,j) = " << denPent[2*(nblocks-2)].getEle(i,j) << " EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-2))(i,j) << std::endl;
				}

			}
		}

		//k = nblocks-1
		std::cout << "k = " << nblocks-1 << std::endl;
		for (unsigned int i = 0; i < nsize; i++) {
			for (unsigned int j = 0; j < nsize; j++) {
				//Getting the Eigen pentadiagonal block matrix
				derr = std::fabs(symPent[nblocks-1].getEle(i,j) - EigenPent.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1))(i,j));

				if(derr > 1e-2){
					std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
					std::cout << "symPent[k].getEle(i,j) = " << symPent[nblocks-2].getEle(i,j) << " EigenPent.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenPent.block<nsize, nsize>(nsize*(nblocks-2), nsize*(nblocks-2))(i,j) << std::endl;
				}

			}
		}

		std::cout << "The magnitude of all the differences is " << derr << std::endl;

	}

	sblocktrid_file << "The resulting SBlockPent is:\n";
	sblocktrid_file << Pent;

	eigen_file << "The resulting SBlockPent is:\n";
	eigen_file << EigenPent;


	sblocktrid_file.close();
	eigen_file.close();

}

void testSymMatrixXVectorPlusDenseMatrixXVector(){

	const unsigned int nsize = 5;

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
	SymMatrix symBlock(dataA, nsize);
	DenseMatrix denseBlock(dataB, nsize, nsize);

	//Generateing random vector
	Eigen::VectorXd eigen_vec_x = Eigen::VectorXd::Random(2*nsize);
	double* vec_x = new double[2*nsize]();
	double* vec_y = new double[nsize]();

	for(unsigned int i = 0; i < 2*nsize; i++) {
		vec_x[i] = eigen_vec_x(i);
	}

	std::cout << "Eigen output for symv and gemv combined:\n" << (A*eigen_vec_x.segment(0, nsize)+B*eigen_vec_x.segment(nsize, nsize)).transpose() << std::endl;

	symBlock.symv(vec_x, vec_y);
	denseBlock.gemv(&vec_x[nsize], vec_y);

	std::cout << "SymMatrix and DenseMatrix output or symv and gemv combined:" << std::endl;
	for(unsigned int i = 0; i < nsize; i++) std::cout << "  " << vec_y[i];
	std::cout << std::endl;

	std::cout << "Eigen output for symv and gemtv combined:\n" << (A*eigen_vec_x.segment(0, nsize)+B.transpose()*eigen_vec_x.segment(nsize, nsize)).transpose() << std::endl;

	std::fill_n(vec_y, nsize, 0.0);
	symBlock.symv(vec_x, vec_y);
	denseBlock.gemtv(&vec_x[nsize], vec_y);

	std::cout << "SymMatrix and DenseMatrix output or symv and gemtv combined:" << std::endl;
	for(unsigned int i = 0; i < nsize; i++) std::cout << "  " << vec_y[i];
	std::cout << std::endl;



	delete [] vec_x;
	delete [] vec_y;

}


void testSymvTrid(){

	std::cout << "*****************************Testing the Trid::symv method*****************************" << std::endl;

	const unsigned int nsize = 5;
	const unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	Eigen::MatrixXd EigenTrid = Eigen::MatrixXd::Zero(nsize*nblocks, nsize*nblocks);

	for(unsigned int k = 0; k < nblocks-1; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

		//Setting the Eigen tridiagonal block matrix
		EigenTrid.block<nsize, nsize>(nsize*k, nsize*k) = A;
		EigenTrid.block<nsize, nsize>(nsize*k, nsize*(k+1)) = B;
		EigenTrid.block<nsize, nsize>(nsize*(k+1), nsize*k) = B.transpose();

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
		symBlocks[k].setData(dataA, nsize);
		denseBlocks[k].setData(dataB, nsize, nsize);
	}

	//The last row in separate
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	//Setting the Eigen tridiagonal block matrix
	EigenTrid.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1)) = A;

	//Assigning matrix elements
	double* dataA = new double[nsize*(nsize+1)/2];
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
	symBlocks[nblocks-1].setData(dataA, nsize);

	//Compute the norm of the difference of the resulting matrices
//	{
//		double derr = 0.0;
//		for(unsigned int k = 0; k < nblocks-1; k++){
//
//			std::cout << "k = " << k << std::endl;
//
//			for (unsigned int i = 0; i < nsize; i++) {
//				for (unsigned int j = 0; j < nsize; j++) {
//					//Getting the Eigen pentadiagonal block matrix
//					derr = std::fabs(symBlocks[k].getEle(i,j) - EigenTrid.block<nsize, nsize>(nsize*k, nsize*k)(i,j));
//					if(derr > 1e-2){
//						std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
//						std::cout << "symBlocks[k].getEle(i,j) = " << symBlocks[k].getEle(i,j) << " EigenTrid.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenTrid.block<nsize, nsize>(nsize*k, nsize*k)(i,j) << std::endl;
//					}
//
//					derr = std::fabs(denseBlocks[k].getEle(i,j) - EigenTrid.block<nsize, nsize>(nsize*k, nsize*(k+1))(i,j));
//
//					if(derr > 1e-2){
//						std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
//						std::cout << "denseBlocks[k].getEle(i,j) = " << denseBlocks[k].getEle(i,j) << " EigenTrid.block<nsize, nsize>(nsize*k, nsize*k)(i,j) = " << EigenTrid.block<nsize, nsize>(nsize*k, nsize*k)(i,j) << std::endl;
//					}
//
//				}
//			}
//
//		}
//
//		//k = nblocks-1
//		std::cout << "k = " << nblocks-1 << std::endl;
//		for (unsigned int i = 0; i < nsize; i++) {
//			for (unsigned int j = 0; j < nsize; j++) {
//				//Getting the Eigen pentadiagonal block matrix
//				derr = std::fabs(symBlocks[nblocks-1].getEle(i,j) - EigenTrid.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1))(i,j));
//
//				if(derr > 1e-2){
//					std::cout << "(i,j) = (" << i << "," << j << ") -> derr = " << derr << std::endl;
//					std::cout << "symBlocks[nblocks-1].getEle(i,j) = " << symBlocks[nblocks-1].getEle(i,j) << " EigenTrid.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1))(i,j) = " << EigenTrid.block<nsize, nsize>(nsize*(nblocks-1), nsize*(nblocks-1))(i,j) << std::endl;
//				}
//
//			}
//		}
//
//	}

	//Forming the tridiagonal block matrix
	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);


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
	eigen_vec_y = EigenTrid*eigen_vec_x;
	Trid.symv(vec_x, vec_y);

//	std::cout << "EigenTrid\n" << EigenTrid << std::endl;


	for(unsigned int i = 0; i < nsize*nblocks; i++) {
		double derr = std::fabs(vec_y[i]-eigen_vec_y(i));

		if(derr > 1e-2){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << std::endl;
		}
	}

	delete [] vec_x;
	delete [] vec_y;

	std::cout << "End of Trid::symv method test. If nothing was print on the screen, everything is fine." << std::endl;

}



} /* namespace matrix */
