/*
 * DenseMatrix_test.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#include <iostream>
#include <cstdlib>
#include "../Sources/DenseMatrix.h"
#include "/home/diego/softwares/eigen3.3/Eigen/Dense"

#include <omp.h>

namespace matrix {

void testPrintOperatorDenseMatrix(){
	std::cout << "*****************************Testing the << operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	DenseMatrix matrix;
	double* data = new double[nsize*nsize];

	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*nsize; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize, nsize);

	//Displaying the array elements
	for(unsigned int i = 0; i < nsize*nsize; i++) std::cout << " data[" << i << "] = " << data[i] << std::endl;

	//Displaying the matrix
	std::cout << matrix << std::endl;
	std::cout << "End of << operator test." << std::endl;
}

void testAssignOperatorDenseMatrix(){
	std::cout << "*****************************Testing the = operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	DenseMatrix matrix;
	double* data = new double[nsize*nsize];

	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*nsize; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	DenseMatrix matrixb;
	std::cout << "The counter = " << matrix.checkCounter() << "." << std::endl;
	matrixb = matrix;
	std::cout << "Copyed matrix:\n" << matrixb << std::endl;
	std::cout << "End of = operator test." << std::endl;
}

void testCloneToDenseMatrix(){
	std::cout << "*****************************Testing the cloneTo() method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	DenseMatrix matrix;
	double* data = new double[nsize*nsize];


	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*nsize; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	DenseMatrix matrixb;
	matrix.cloneTo(matrixb);
	std::cout << "Copyed matrix:\n" << matrixb << std::endl;
	std::cout << "End of cloneTo test." << std::endl;
}

void testCloneTransposeToDenseMatrix(){
	std::cout << "*****************************Testing the cloneTransposeTo() method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	DenseMatrix matrix;
	double* data = new double[nsize*nsize];


	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*nsize; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	DenseMatrix matrixb;
	matrix.cloneTransposeTo(matrixb);
	std::cout << "Copyed transpose matrix:\n" << matrixb << std::endl;
	std::cout << "End of cloneTransposeTo test." << std::endl;
}


void testCloneFromDenseMatrix(){
	std::cout << "*****************************Testing the cloneFrom() method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	DenseMatrix matrix;
	double* data = new double[nsize*nsize];


	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*nsize; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	DenseMatrix matrixb;
	matrixb.cloneFrom(matrix);
	std::cout << "Copyed matrix:\n" << matrixb << std::endl;
	std::cout << "End of cloneFrom test." << std::endl;
}

void testCloneTransposeFromDenseMatrix(){
	std::cout << "*****************************Testing the cloneTransposeFrom() method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	DenseMatrix matrix;
	double* data = new double[nsize*nsize];


	//Generating matrix elements
	for(unsigned int i = 0; i < nsize*nsize; i++) data[i] = static_cast<double>(rand())/RAND_MAX;

	matrix.setData(data, nsize, nsize);

	//Displaying the original matrix
	std::cout << "Original matrix:\n" << matrix << std::endl;

	//Displaying the copyed matrix
	DenseMatrix matrixb;
	matrixb.cloneTransposeFrom(matrix);
	std::cout << "Copyed transpose matrix:\n" << matrixb << std::endl;
	std::cout << "End of cloneTransposeFrom test." << std::endl;
}

void testMultiplyOperatorDenseMatrixXDenseMatrix(){

	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	std::cout << "A :\n" << A << std::endl;
	std::cout << "B :\n" << A << std::endl;
	std::cout << "Eigen A*B :\n" << A*B << std::endl;

	double* dataA = new double[nsize*nsize];
	double* dataB = new double[nsize*nsize];
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
			dataB[i*A.rows()+j] = B(i,j);
		}
	}

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	DenseMatrix denseB(dataB, nsize, nsize);
	DenseMatrix denseAB;
	denseAB = denseA*denseB;
	std::cout << "Tested operator A*B :" << denseAB << std::endl;
	std::cout << "End of * operator test." << std::endl;
}

void testMultiplyOperatorDenseMatrixXSymMatrix(){

	std::cout << "*****************************Testing the * operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize); B = B*B.transpose();

	std::cout << "A :\n" << A << std::endl;
	std::cout << "B :\n" << B << std::endl;
	std::cout << "Eigen A*B :\n" << A*B << std::endl;

	double* dataA = new double[nsize*nsize];
	double* dataB = new double[nsize*(nsize+1)/2];
	unsigned int dataApos = 0;
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			if(j >= i) {
				dataB[dataApos] = B(i,j);
				dataApos++;
			}
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	//DenseMatrix objects
	SymMatrix symB(dataB, nsize);
	DenseMatrix denseA(dataA, nsize, nsize);
	DenseMatrix denseAB;
	denseAB = denseA*symB;
	std::cout << "Tested operator A*B :" << denseAB << std::endl;
	std::cout << "End of * operator test." << std::endl;
}

void testDenseMatrixXVector(){

	std::cout << "*****************************Testing the gemv method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::VectorXd v = Eigen::VectorXd::Random(nsize);
	std::cout << "A :\n" << A << std::endl;

	std::cout << "Eigen A*v :\n" << A*v << std::endl;

	double* dataA = new double[nsize*nsize];
	double* dataV = new double[nsize]();
	double* dataU = new double[nsize]();

	for(unsigned int i = 0; i < A.rows(); i++){
		dataV[i] = v(i);
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	denseA.gemv(dataV, dataU);
	std::cout << "Tested method output :" << std::endl;
	for(unsigned int i = 0; i < denseA.getNRows(); i++){
		std::cout << " u[" << i << "] = " << dataU[i] << std::endl;
	}
	delete [] dataU;
	delete [] dataV;
	std::cout << "End of gemv test." << std::endl;
}

void testDenseMatrixXComplexVector(){

	std::cout << "*****************************Testing the gemv method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::VectorXcd v = Eigen::VectorXcd::Random(nsize);
	std::cout << "A :\n" << A << std::endl;

	std::cout << "Eigen A*v :\n" << A*v << std::endl;

	double* dataA = new double[nsize*nsize];
	std::complex<double>* dataV = new std::complex<double>[nsize]();
	std::complex<double>* dataU = new std::complex<double>[nsize]();

	for(unsigned int i = 0; i < A.rows(); i++){
		dataV[i] = static_cast<std::complex<double> >(v(i));
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	denseA.gemv(dataV, dataU);
	std::cout << "Tested method output :" << std::endl;
	for(unsigned int i = 0; i < denseA.getNRows(); i++){
		std::cout << " u[" << i << "] = " << dataU[i] << std::endl;
	}
	delete [] dataU;
	delete [] dataV;
	std::cout << "End of gemv test." << std::endl;
}


void testDenseMatrixTransposedXVector(){
	std::cout << "*****************************Testing the gemtv method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::VectorXd v = Eigen::VectorXd::Random(nsize);
	std::cout << "A :\n" << A << std::endl;

	std::cout << "Eigen A^T*v :\n" << A.transpose()*v << std::endl;

	double* dataA = new double[nsize*nsize];
	double* dataV = new double[nsize]();
	double* dataU = new double[nsize]();

	for(unsigned int i = 0; i < A.rows(); i++){
		dataV[i] = v(i);
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	denseA.gemtv(dataV, dataU);
	std::cout << "Tested method output :" << std::endl;
	for(unsigned int i = 0; i < denseA.getNRows(); i++){
		std::cout << " u[" << i << "] = " << dataU[i] << std::endl;
	}
	delete [] dataU;
	delete [] dataV;
	std::cout << "End of gemtv test." << std::endl;
}

void testDenseMatrixTransposedXComplexVector(){

	std::cout << "*****************************Testing the gemv method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::VectorXcd v = Eigen::VectorXcd::Random(nsize);
	std::cout << "A :\n" << A << std::endl;

	std::cout << "Eigen A^T*v :\n" << A.transpose()*v << std::endl;

	double* dataA = new double[nsize*nsize];
	std::complex<double>* dataV = new std::complex<double>[nsize]();
	std::complex<double>* dataU = new std::complex<double>[nsize]();

	for(unsigned int i = 0; i < A.rows(); i++){
		dataV[i] = static_cast<std::complex<double> >(v(i));
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	denseA.gemtv(dataV, dataU);
	std::cout << "Tested method output :" << std::endl;
	for(unsigned int i = 0; i < denseA.getNRows(); i++){
		std::cout << " u[" << i << "] = " << dataU[i] << std::endl;
	}
	delete [] dataU;
	delete [] dataV;
	std::cout << "End of gemv test." << std::endl;
}

void testAAT(){

	std::cout << "*****************************Testing the AAT method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);

	std::cout << "A :\n" << A << std::endl;

	double* dataA = new double[nsize*nsize];
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	A = A*A.transpose();
	std::cout << "A*A^T :\n" << A << std::endl;

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	std::cout << "DenseMatrix A :" << denseA << std::endl;
	std::cout << "Resulting SymMatrix A*A^T :" << denseA.AAT() << std::endl;
	std::cout << "End of AAT method test." << std::endl;
}

void testATA(){

	std::cout << "*****************************Testing the ATA method*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);

	std::cout << "A :\n" << A << std::endl;

	double* dataA = new double[nsize*nsize];
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
		}
	}

	A = A.transpose()*A;
	std::cout << "A^T*A :\n" << A << std::endl;

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	std::cout << "DenseMatrix A :" << denseA << std::endl;
	std::cout << "Resulting SymMatrix A^T*A :" << denseA.ATA() << std::endl;
	std::cout << "End of ATA method test." << std::endl;
}


void testAddOperatorDenseMatrixXDenseMatrix(){

	std::cout << "*****************************Testing the + operator*****************************" << std::endl;
	srand((unsigned int) time(0));

	unsigned int nsize = 5;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(nsize, nsize);
	Eigen::MatrixXd B = Eigen::MatrixXd::Random(nsize, nsize);

	std::cout << "A :\n" << A << std::endl;
	std::cout << "B :\n" << B << std::endl;
	std::cout << "Eigen A+B :\n" << A+B << std::endl;

	double* dataA = new double[nsize*nsize];
	double* dataB = new double[nsize*nsize];
	for(unsigned int i = 0; i < A.rows(); i++){
		for(unsigned int j = 0; j < A.cols(); j++){
			dataA[i*A.rows()+j] = A(i,j);
			dataB[i*A.rows()+j] = B(i,j);
		}
	}

	//DenseMatrix objects
	DenseMatrix denseA(dataA, nsize, nsize);
	DenseMatrix denseB(dataB, nsize, nsize);
	DenseMatrix denseAB;
	denseAB = denseA+denseB;
	std::cout << "Tested operator A+B :" << denseAB << std::endl;
	std::cout << "End of + operator test." << std::endl;
}



} /* namespace matrix */
