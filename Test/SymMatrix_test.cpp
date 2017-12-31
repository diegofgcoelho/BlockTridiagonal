/*
 * SymMatrix_test.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#include <iostream>
#include <cstdlib>
#include "../Sources/SymMatrix.h"
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

namespace matrix {

void testPrintOperatorSymMatrix(){
	std::cout << "*****************************Testing the << operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	SymMatrix matrix;
	double* data = new double[nsize*(nsize+1)/2];

	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*(nsize+1)/2; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize);

	//Displaying the array elements
	for(unsigned int i = 0; i < nsize*(nsize+1)/2; i++) std::cout << " data[" << i << "] = " << data[i] << std::endl;

	//Displaying the matrix
	std::cout << matrix << std::endl;
	matrix.free();
	std::cout << "End of << operator test." << std::endl;
}

void testAssignOperatorSymMatrix(){
	std::cout << "*****************************Testing the = operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	SymMatrix matrix;
	double* data = new double[nsize*(nsize+1)/2];

	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*(nsize+1)/2; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	SymMatrix matrixb;
	matrixb = matrix;
	std::cout << "Copyed matrix:\n" << matrixb << std::endl;
	matrix.free();
	matrixb.free();
	std::cout << "End of = operator test." << std::endl;
}

void testCloneToSymMatrix(){
	std::cout << "*****************************Testing the cloneTo() method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	SymMatrix matrix;
	double* data = new double[nsize*(nsize+1)/2];


	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*(nsize+1)/2; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	SymMatrix matrixb;
	matrix.cloneTo(matrixb);
	std::cout << "Copyed matrix:\n" << matrixb << std::endl;
	matrix.free();
	matrixb.free();
	std::cout << "End of cloneTo test." << std::endl;
}

void testCloneFromSymMatrix(){
	std::cout << "*****************************Testing the cloneFrom() method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	SymMatrix matrix;
	double* data = new double[nsize*nsize];


	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*(nsize+1)/2; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	SymMatrix matrixb;
	matrixb.cloneFrom(matrix);
	std::cout << "Copyed matrix:\n" << matrixb << std::endl;
	matrix.free();
	matrixb.free();
	std::cout << "End of cloneFrom test." << std::endl;
}

void testMultiplyOperatorSymMatrixXDenseMatrix(){

	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	std::cout << "A :\n" << A << std::endl;
	std::cout << "B :\n" << B << std::endl;
	std::cout << "Eigen A*B :\n" << A*B << std::endl;

	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int dataApos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[dataApos] = A(i,j);
				dataApos++;
			}
			dataB[i*A.rows()+j] = B(i,j);
		}
	}

	//DenseMatrix objects
	SymMatrix symA(dataA, nsize);
	DenseMatrix denseB(dataB, nsize, nsize);
	DenseMatrix denseAB;
	denseAB = symA*denseB;
	std::cout << "Tested operator A*B :" << denseAB << std::endl;
	symA.free();
	denseB.free();
	denseAB.free();
	std::cout << "End of * operator test." << std::endl;
}

void testMultiplyOperatorSymMatrixXSymMatrix(){

	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize); B = B*B.transpose();

	std::cout << "A :\n" << A << std::endl;
	std::cout << "B :\n" << B << std::endl;
	std::cout << "Eigen A*B :\n" << A*B << std::endl;

	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				dataB[datapos] = B(i,j);
				datapos++;
			}
		}
	}

	//DenseMatrix objects
	SymMatrix symA(dataA, nsize);
	SymMatrix symB(dataB, nsize);
	DenseMatrix denseAB;
	denseAB = symA*symB;
	std::cout << "Tested operator A*B :" << denseAB << std::endl;
	symA.free();
	symB.free();
	denseAB.free();
	std::cout << "End of * operator test." << std::endl;
}

void testSymMatrixXVector(){

	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A= A*A.transpose();
	Eigen::VectorXd v = Eigen::VectorXd::Random(nsize);
	std::cout << "A :\n" << A << std::endl;

	std::cout << "Eigen A*v :\n" << A*v << std::endl;

	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataV = new double[nsize]();
	double* dataU = new double[nsize]();
	unsigned int dataApos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		dataV[i] = v(i);
		for(unsigned int j = i; j < A.cols(); j++){
			dataA[dataApos] = A(i,j);
			dataApos++;
		}
	}

	//DenseMatrix objects
	SymMatrix symA(dataA, nsize);
	std::cout << "symA: " << symA << std::endl;
	symA.symv(dataV, dataU);
	std::cout << "Tested method output :" << std::endl;
	for(unsigned int i = 0; i < symA.getNRows(); i++){
		std::cout << " v[" << i << "] = " << dataU[i] << std::endl;
	}
	symA.free();
	delete [] dataU;
	delete [] dataV;
	std::cout << "End of gemv test." << std::endl;
}

void testSymMatrixXComplexVector(){
	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A= A*A.transpose();
	Eigen::VectorXcd v = Eigen::VectorXcd::Random(nsize);
	std::cout << "A :\n" << A << std::endl;

	std::cout << "Eigen A*v :\n" << A*v << std::endl;

	double* dataA = new double[nsize*(nsize+1)/2];
	std::complex<double>* dataV = new std::complex<double>[nsize]();
	std::complex<double>* dataU = new std::complex<double>[nsize]();
	unsigned int dataApos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		dataV[i] = v(i);
		for(unsigned int j = i; j < A.cols(); j++){
			dataA[dataApos] = A(i,j);
			dataApos++;
		}
	}

	//DenseMatrix objects
	SymMatrix symA(dataA, nsize);
	std::cout << "symA: " << symA << std::endl;
	symA.symv(dataV, dataU);
	std::cout << "Tested method output :" << std::endl;
	for(unsigned int i = 0; i < symA.getNRows(); i++){
		std::cout << " v[" << i << "] = " << dataU[i] << std::endl;
	}
	symA.free();
	delete [] dataU;
	delete [] dataV;
	std::cout << "End of gemv test." << std::endl;

}

void testAddOperatorSymMatrixXSymMatrix(){

	std::cout << "*****************************Testing the + operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize); B = B*B.transpose();

	std::cout << "A :\n" << A << std::endl;
	std::cout << "B :\n" << B << std::endl;
	std::cout << "Eigen A+B :\n" << A+B << std::endl;

	double* dataA = new double[nsize*(nsize+1)/2];
	double* dataB = new double[nsize*nsize];
	unsigned int datapos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataA[datapos] = A(i,j);
				dataB[datapos] = B(i,j);
				datapos++;
			}
		}
	}

	//DenseMatrix objects
	SymMatrix symA(dataA, nsize);
	SymMatrix symB(dataB, nsize);
	SymMatrix symAB;
	symAB = symA+symB;
	std::cout << "Tested operator A+B :" << symAB << std::endl;
	symA.free();
	symB.free();
	symAB.free();
	std::cout << "End of + operator test." << std::endl;
}


void testSquareSymMatrix(){
	std::cout << "*****************************Testing the square method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize); A = A*A.transpose();

	std::cout << "A :\n" << A << std::endl;
	std::cout << "Eigen A*A :\n" << A*A << std::endl;

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

	//DenseMatrix objects
	SymMatrix symA(dataA, nsize);
	SymMatrix symAA = symA.square();
	std::cout << "Tested square method :" << symAA << std::endl;
	symA.free();
	symAA.free();
	std::cout << "End of square method test." << std::endl;

}

} /* namespace matrix */
