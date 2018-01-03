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

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks-1];
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
		double derr = std::abs(static_cast<double>(vec_y[i]-eigen_vec_y(i)));

		if(derr > 1e-2){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << "\n";
		}
	}

	delete [] vec_x;
	delete [] vec_y;

	std::cout << "End of Trid::symv method test. If nothing was print on the screen, everything is fine." << std::endl;

}

void testSymvTridComplex(){

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

	//Forming the tridiagonal block matrix
	SBlockTrid Trid;
	Trid.setBlocks(symBlocks, denseBlocks, nblocks);


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
	eigen_vec_y = EigenTrid*eigen_vec_x;
	Trid.symv(vec_x, vec_y);

//	std::cout << "EigenTrid\n" << EigenTrid << std::endl;


	for(unsigned int i = 0; i < nsize*nblocks; i++) {
		double derr = std::abs(static_cast<std::complex<double> >(vec_y[i]-eigen_vec_y(i)));

		if(derr > 1e-3){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << "\n";
		}
	}

	delete [] vec_x;
	delete [] vec_y;

	std::cout << "End of Trid::symv method test. If nothing was print on the screen, everything is fine." << std::endl;

}

void testSymvTridComplexAndScale(){

	std::cout << "*****************************Testing the Trid::symv and support::scal method*****************************" << std::endl;

	const unsigned int bsize = 5;
	const unsigned int nblocks = 7;

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
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


	//Generateing random vector
	Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(bsize*nblocks);
	Eigen::VectorXcd eigen_vec_y = Eigen::VectorXcd::Random(bsize*nblocks);
	std::complex<double>* vec_x = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_y = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_z = new std::complex<double>[bsize*nblocks]();

	//Copying eigen_vec_x value to vec_x
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		vec_x[i] = eigen_vec_x(i);
	}

	//Testing the symv a few times

	std::fill_n(&vec_y[0], nblocks*bsize, std::complex<double>(0.0,0.0));
	Trid.symv(vec_x, vec_y);
	std::cout << "\n vec_y      [";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		std::cout  << vec_y[i] << "  ";
	}

	eigen_vec_y = EigenTrid*eigen_vec_x;
	double derrtest = 0.0;
	std::cout << "\n" << 0 << "eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_y(j) << "  ";
		derrtest+=std::norm(vec_y[j]-eigen_vec_y(j));
	}std::cout  << "\n\n1st difference acc = " << derrtest<< "\n\n\n" << std::endl;

	derrtest = 0.0;
	//std::memcpy(vec_z, vec_y, nblocks*bsize*sizeof(std::complex<double>));
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		vec_z[i] = vec_y[i];
	}
	std::fill_n(&vec_y[0], nblocks*bsize, std::complex<double>(0.0,0.0));
	Trid.symv(vec_z, vec_y);
	std::cout << "\n vec_z      [";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		std::cout  << vec_z[i] << "  ";
	}
	std::cout << "\n vec_y      [";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		std::cout  << vec_y[i] << "  ";
	}

	std::cout << "\n" << "eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_y(j) << "  ";
	}
	eigen_vec_y = EigenTrid*eigen_vec_y;
	std::cout << "\n" << "eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_y(j) << "  ";
		derrtest+=std::norm(vec_y[j]-eigen_vec_y(j));
	}std::cout  << "\n\n2nd difference acc = " << derrtest<< "\n\n\n" << std::endl;

	derrtest = 0.0;
	//std::memcpy(vec_z, vec_y, nblocks*bsize*sizeof(std::complex<double>));
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		vec_z[i] = vec_y[i];
	}
	std::fill_n(&vec_y[0], nblocks*bsize, std::complex<double>(0.0,0.0));
	Trid.symv(vec_z, vec_y);
	std::cout << "\n vec_z      [";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		std::cout  << vec_z[i] << "  ";
	}
	std::cout << "\n vec_y      [";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		std::cout  << vec_y[i] << "  ";
	}

	std::cout << "\n" << "eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_y(j) << "  ";
	}
	eigen_vec_y = EigenTrid*eigen_vec_y;
	std::cout << "\n" << "eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_y(j) << "  ";
		derrtest+=std::norm(vec_y[j]-eigen_vec_y(j));
	}std::cout  << "\n\n3rd difference acc = " << derrtest<< "\n\n\n" << std::endl;


	//End of testing the symv a few times

	for(unsigned int i = 0; i < bsize*nblocks; i++) {
		double derr = std::abs(static_cast<std::complex<double> >(vec_y[i]-eigen_vec_y(i)));

		if(derr > 1e-3){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << "\n";
		}else {
			std::cout << "Same at i = " << i << " with eigen_vec_y = " << eigen_vec_y[i] << " and vec_y = " << vec_y[i] << "\n";
		}
	}

	delete [] vec_x;
	delete [] vec_y;
	delete [] vec_z;

	std::cout << "End of Trid::symv method test. If nothing was print on the screen, everything is fine." << std::endl;

}


void testPowerMethodTrid(){

	srand((unsigned int) time(0));

	const unsigned int nblocks = 7;
	const unsigned int bsize = 5;

	//Lambda
	std::complex<double> lambda(1.0,0.0), elambda(1.0,0.0);
	//Precision
	double prec = 1e-3;
	//Maximum number of iterations
	unsigned int maxIter = 100;
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

	//End of matrix generation. From this point on, all the matrices are set.

	//Generating random vector
	Eigen::VectorXcd eigen_vec_x = Eigen::VectorXcd::Random(bsize*nblocks);
	Eigen::VectorXcd eigen_vec_y = Eigen::VectorXcd::Random(bsize*nblocks);
	Eigen::VectorXcd eigen_vec_z = Eigen::VectorXcd::Random(bsize*nblocks);
	std::complex<double>* vec_x = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_y = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_z = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* vec_w = new std::complex<double>[bsize*nblocks]();
	std::complex<double>* data_eigen_vec_y = &eigen_vec_y[0];

	//Setting up the initial guesses
	for(unsigned int i = 0; i < bsize*nblocks; i++) {
		vec_x[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
		vec_y[i] = static_cast<std::complex<double> >(eigen_vec_x(i));
	}


	eigen_vec_y = EigenTrid*eigen_vec_x;
	eigen_vec_y = eigen_vec_y*(1.0/elambda);
	std::cout << "\n" << 0 << "th eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_y(j) << "  ";
	}
	std::cout << "\n" << 0 << "th elambda = " << elambda << std::endl;
	elambda = support::max_mag_cmp(data_eigen_vec_y, bsize*nblocks);


	for(unsigned int i = 1; i < maxIter; i++){

		eigen_vec_y = EigenTrid*eigen_vec_y;
		eigen_vec_y = eigen_vec_y/elambda;
		std::cout << "\n" << i << "th eigen_vec_y[";
		for(unsigned int j = 0; j < nblocks*bsize; j++){
			std::cout  << eigen_vec_y(j) << "  ";
		}
		std::cout << "\n" << i << "th elambda = " << elambda << std::endl;
		elambda = support::max_mag_cmp(data_eigen_vec_y, bsize*nblocks);
	}

	std::cout << "\n" << " elambda*eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << elambda*eigen_vec_y(j) << "  ";
	}std::cout << "\n" << std::endl;

	eigen_vec_z = EigenTrid*eigen_vec_y;
	std::cout << "\n" << " T*eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_z(j) << "  ";
	}std::cout << "\n" << std::endl;

	std::cout << "The norm of difference in eigen is l*v-Tx:" << (elambda*eigen_vec_y-eigen_vec_z).norm()/eigen_vec_y.norm() << std::endl;

	//Executing power method
	Trid.powerMethod(vec_y, &lambda, prec, maxIter, &iter, &dnorm);

	std::cout << "The lambda returned is: " << lambda << std::endl;

	std::cout << "\n" << "The returned eigenvalue   vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << vec_y[j] << "  ";
	}std::cout << "\n" << std::endl;

	std::fill_n(vec_z, bsize*nblocks, std::complex<double>(0.0,0.0));
	Trid.symv(vec_y, vec_z);

	std::memcpy(vec_w, vec_y, sizeof(std::complex<double>)*bsize*nblocks);

	support::scal(vec_w, bsize*nblocks, lambda);

	std::cout << "\n" << " elambda*eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << elambda*eigen_vec_y(j) << "  ";
	}std::cout << "\n" << std::endl;

	std::cout << "\n" << " T*eigen_vec_y[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << eigen_vec_z(j) << "  ";
	}std::cout << "\n" << std::endl;


	std::cout << "\n" << "The lambda*v[";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << vec_w[j] << "  ";
	}std::cout << "\n" << std::endl;
	std::cout << "\n" << "The T*v     [";
	for(unsigned int j = 0; j < nblocks*bsize; j++){
		std::cout  << vec_z[j] << "  ";
	}std::cout << "\n" << std::endl;

	std::cout << "\n" << "diff ours lambda*v-T*v [";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		vec_w[i] = vec_w[i]-vec_z[i];
		std::cout  << vec_w[i] << "  ";
	}std::cout << "\n" << std::endl;

	std::cout << "\n" << "diff Eigen lambda*v-T*v[";
	for(unsigned int i = 0; i < nblocks*bsize; i++){
		std::cout  << elambda*eigen_vec_y(i)-eigen_vec_z(i) << "  ";
	}std::cout << "\n" << std::endl;

	std::cout << "\n\nThe norm of difference in SBlockTrid is l*v-Tx:" << support::norm(vec_w, bsize*nblocks)/support::norm(vec_y, bsize*nblocks) << std::endl;

	double derr = 0.0;
	for(unsigned int i = 0; i < bsize*nblocks; i++) {
		derr += std::abs(static_cast<std::complex<double> >(vec_w[i]-(elambda*eigen_vec_y(i)-eigen_vec_z(i))));

		if(derr > 1e-2){
			std::cout << "Different at i = " << i << " with eigen_vec_y = " << eigen_vec_y(i) << " and vec_y = " << vec_y[i] << "\n";
		}
	}
	std::cout << "derr = " << derr << std::endl;
//
//	std::cout << "The returned largest eigenvalue is: " << lambda << " with dnorm = " << dnorm << " and " << iter << " iterations." << std::endl;
//
//	std::cout << "\n vec_x      [";
//	for(unsigned int i = 0; i < nblocks*bsize; i++){
//		std::cout  << vec_x[i] << "  ";
//	}
//	std::cout << "\n vec_y      [";
//	for(unsigned int i = 0; i < nblocks*bsize; i++){
//		std::cout  << vec_y[i] << "  ";
//	}
//	std::cout << "\n eigen_vec_y[";
//	for(unsigned int i = 0; i < nblocks*bsize; i++){
//		std::cout  << eigen_vec_y(i) << "  ";
//	}



	delete [] vec_x;
	delete [] vec_y;
	delete [] vec_w;
	delete [] vec_z;

}

void testPowerMethodTrid2(){

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
	Trid.powerMethod(vec_y, &lambda, prec, maxIter, &iter, &dnorm);

	std::cout << "The returned largest eigenvalue is: " << lambda << " with dnorm = " << dnorm << " and " << iter << " iterations." << std::endl;

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


} /* namespace matrix */
