/*
 * SBlockTrid.cpp
 *
 *  Created on: Dec 20, 2017
 *      Author: diego
 */

#include "SBlockTrid.h"

namespace matrix {

SBlockTrid::SBlockTrid() {
	this->symBlocks = NULL;
	this->denseBlocks = NULL;
	this->nblocks = 0;
}

SBlockTrid::SBlockTrid(SymMatrix* symBlocks, DenseMatrix* denseBlocks, unsigned int nblocks) {
	this->nblocks = nblocks;
	this->symBlocks = symBlocks;
	this->denseBlocks = denseBlocks;
}

SBlockTrid::~SBlockTrid() {
	this->free();
}

void SBlockTrid::setBlocks(SymMatrix* symBlocks,
		DenseMatrix* denseBlocks,
		unsigned int nblocks) {
	this->symBlocks = symBlocks;
	this->denseBlocks = denseBlocks;
	this->nblocks = nblocks;

}

void SBlockTrid::setBlockByClone(SymMatrix* symBlock, unsigned int i) {
	this->symBlocks[i].cloneFrom(*symBlock);
}

void SBlockTrid::setBlockByClone(DenseMatrix* denseBlock,
		unsigned int i) {
	this->denseBlocks[i].cloneFrom(*denseBlock);
}

void SBlockTrid::getBlocks(SymMatrix** symBlocks,
		DenseMatrix** denseMatrix) {
	*symBlocks = this->symBlocks;
	*denseMatrix = this->denseBlocks;
}

SBlockPent SBlockTrid::square() {
	SymMatrix* psymBlocks = new SymMatrix[this->nblocks];
	DenseMatrix* pdenseMatrix = new DenseMatrix[2*this->nblocks-3];

	//Block row 0
	//00
	psymBlocks[0] = this->symBlocks[0].square()+this->denseBlocks[0].AAT();

	//01
	pdenseMatrix[0] = this->symBlocks[0]*this->denseBlocks[0]+this->denseBlocks[0]*this->symBlocks[1];

	//02
	pdenseMatrix[1] = this->denseBlocks[0]*this->denseBlocks[1];

	for(unsigned int i = 1; i < this->nblocks-2; i++){

		//ii
		psymBlocks[i] = this->denseBlocks[i-1].ATA()
					+this->symBlocks[i].square()
					+this->denseBlocks[i].AAT();
		//i(i+1)
		pdenseMatrix[2*i] = this->symBlocks[i]*this->denseBlocks[i]+this->denseBlocks[i]*this->symBlocks[i+1];

		//i(i+2)
		pdenseMatrix[2*i+1] = this->denseBlocks[i]*this->denseBlocks[i+1];
	}

	//(this->nblocks-2)(this->nblocks-2)
	psymBlocks[this->nblocks-2] = this->denseBlocks[this->nblocks-3].ATA()
				+this->symBlocks[this->nblocks-2].square()
				+this->denseBlocks[this->nblocks-2].AAT();
	//(this->nblocks-2)(this->nblocks-1)
	pdenseMatrix[2*(this->nblocks-2)] = this->symBlocks[this->nblocks-2]*this->denseBlocks[this->nblocks-2]
					+this->denseBlocks[this->nblocks-2]*this->symBlocks[this->nblocks-1];

	//(this->nblocks-1)(this->nblocks-1)
	psymBlocks[this->nblocks-1] = this->denseBlocks[this->nblocks-2].ATA()
				+this->symBlocks[this->nblocks-1].square();

	SBlockPent C(psymBlocks, pdenseMatrix, this->nblocks);
	return C;
}

void SBlockTrid::symv(double* src, double* dest) {

	const unsigned int bsize = this->symBlocks[0].getNRows();

	//Running over all the blocks in the main diagonal
	for(unsigned int i = 0; i < this->nblocks-1; i++){
		//dest[i*this->nblocks]
		this->symBlocks[i].symv(&src[i*bsize], &dest[i*bsize]);
		this->denseBlocks[i].gemv(&src[(i+1)*bsize],&dest[i*bsize]);
		this->denseBlocks[i].gemtv(&src[i*bsize], &dest[(i+1)*bsize]);

	}

	//For the last block
	this->symBlocks[this->nblocks-1].symv(&src[(this->nblocks-1)*bsize], &dest[(this->nblocks-1)*bsize]);

}

void SBlockTrid::symv(std::complex<double>* src, std::complex<double>* dest) {

	const unsigned int bsize = this->symBlocks[0].getNRows();

	//Running over all the blocks in the main diagonal
	for(unsigned int i = 0; i < this->nblocks-1; i++){
		//dest[i*this->nblocks]
		this->symBlocks[i].symv(&src[i*bsize], &dest[i*bsize]);
		this->denseBlocks[i].gemv(&src[(i+1)*bsize],&dest[i*bsize]);
		this->denseBlocks[i].gemtv(&src[i*bsize], &dest[(i+1)*bsize]);

	}

	//For the last block
	this->symBlocks[this->nblocks-1].symv(&src[(this->nblocks-1)*bsize], &dest[(this->nblocks-1)*bsize]);

}


void SBlockTrid::free(){
	if(this->symBlocks != NULL && this->denseBlocks != NULL){
		for(unsigned int i = 0; i < this->nblocks-1; i++){
			this->symBlocks[i].~SymMatrix();
			this->denseBlocks[i].~DenseMatrix();
		}

		this->symBlocks[this->nblocks-1].~SymMatrix();

		delete [] this->symBlocks;
		delete [] this->denseBlocks;

		this->symBlocks = NULL;
		this->denseBlocks = NULL;
	}
}

std::ostream& operator<<(std::ostream& os, const SBlockTrid& obj){
	os << "\nThe Matrix Blocks are";
	for(unsigned int i = 0; i < obj.nblocks-1; i++){
		os << "blocks at row ";
		os << i;
		os << " are : ";
		os << obj.symBlocks[i];
		os << "\nand\n";
		os << obj.denseBlocks[i];
		os << "\n";
	}
	os << "blocks at row ";
	os << obj.nblocks-1;
	os << " are : ";
	os << obj.symBlocks[obj.nblocks-1];
	os << "\n";

	return os;
}

SBlockTrid& SBlockTrid::operator=(const SBlockTrid& obj){

	this->free();

	this->symBlocks = new SymMatrix[obj.nblocks];
	this->denseBlocks = new DenseMatrix[obj.nblocks-1];
	for(unsigned int i = 0; i < obj.nblocks-1; i++){
		this->symBlocks[i] = obj.symBlocks[i];
		this->denseBlocks[i] = obj.denseBlocks[i];
	}

	this->symBlocks[obj.nblocks-1] = obj.symBlocks[obj.nblocks-1];

	this->nblocks = obj.nblocks;

	return *this;
}

void SBlockTrid::powerMethod(std::complex<double>* v, std::complex<double>* lambda,
		double prec, unsigned int maxIter, unsigned int* iter, double* dnorm) {
	/*Input:
	 * v is the initial guess for the eigenvector associated with the
	 * largest eigenvalue
	 * lambda is the initial guess for the eigenvalue
	 * maxIter is the maximum number of iterations
	 * prec is the precision using for the stopping criteria for the vector difference norm
	 */
	/*Output:
	 * v represents the eigenvector associated with the largest eigenvalue
	 * lambda represents the largest eigenvalue
	 * iter is the number of iterations used for computing the largest eigenvalue
	 * dnorm is the norm of the difference between the last consecutive eigenvector estimates
	 */
	/*Requirement:
	 * this object must be initialized (with non null symBlocks and denseBlocks)
	 * the largest eigenvalue and its eigenvector must be complex
	 * v must be non-null
	 * maxIter must be greater than 1
	 */
	/*Description: this function returns the largest eigenvalue of this matrix object using
	 * power method
	 */

	if(maxIter <= 1){
		std::cout << "Error: the maximum number of iterations must be at least 1." << std::endl;
		return;
	}

	//Vectors size
	const unsigned int vsize = this->nblocks*this->symBlocks[0].getNRows();

	//Auxiliary vectors
	std::complex<double>* vv = new std::complex<double>[vsize]();
	std::complex<double>* vvv = new std::complex<double>[vsize]();

	for(*iter = 0; *iter < maxIter; (*iter)++){

		//Update the past vectors
		std::memcpy(vvv, vv, sizeof(std::complex<double>)*vsize);
		std::memcpy(vv, v, sizeof(std::complex<double>)*vsize);

		//Cleaning up v
		std::fill_n(v, vsize, std::complex<double>(0.0,0.0));
		//Perform the matrix vector multiplication
		this->symv(vv, v);

		//Perform the vector scaling
		support::scal(v, vsize, 1.0/(*lambda));

		//Find the largest coefficient of resulting vector
//		*lambda = *std::max_element(v, v+vsize, support::cmpcmp);
		*lambda = support::max_mag_cmp(v, vsize);

		//Stopping criteria
		if(support::check_stop(v, vv, vvv, vsize, prec, dnorm)){
			break;
		}

	}

	delete [] vv;
	delete [] vvv;

}

void SBlockTrid::random(unsigned int nblocks, unsigned int bsize){

	srand((unsigned int) time(0));

	DenseMatrix* denseBlocks = new DenseMatrix[nblocks];
	SymMatrix* symBlocks = new SymMatrix[nblocks];

	for(unsigned int k = 0; k < nblocks-1; k++){
		Eigen::MatrixXd A = Eigen::MatrixXd::Random(bsize, bsize); A = A*A.transpose();
		Eigen::MatrixXd B = Eigen::MatrixXd::Random(bsize, bsize);

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

	this->symBlocks = symBlocks;
	this->denseBlocks = denseBlocks;
	this->nblocks = nblocks;

}

} /* namespace matrix */
