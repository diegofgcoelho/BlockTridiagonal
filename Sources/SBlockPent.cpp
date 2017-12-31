/*
 * SBlockPent.cpp
 *
 *  Created on: Dec 20, 2017
 *      Author: diego
 */

#include "SBlockPent.h"

namespace matrix {

SBlockPent::SBlockPent() {
	this->symBlocks = NULL;
	this->denseBlocks = NULL;
	this->nblocks = 0;
}

SBlockPent::SBlockPent(SymMatrix* symBlocks, DenseMatrix* denseBlocks, unsigned int nblocks){
	this->symBlocks = symBlocks;
	this->denseBlocks = denseBlocks;
	this->nblocks = nblocks;
}

SBlockPent::~SBlockPent() {
	this->free();
}

void SBlockPent::free(){
		if(this->symBlocks != NULL && this->denseBlocks != NULL){
			for(unsigned int i = 0; i < this->nblocks-2; i++){
				this->symBlocks[i].~SymMatrix();
				this->denseBlocks[2*i].~DenseMatrix();
				this->denseBlocks[2*i+1].~DenseMatrix();
			}

			this->symBlocks[this->nblocks-2].~SymMatrix();
			this->denseBlocks[2*(this->nblocks-2)].~DenseMatrix();

			this->symBlocks[this->nblocks-1].~SymMatrix();

			delete [] this->symBlocks;
			delete [] this->denseBlocks;

			this->symBlocks = NULL;
			this->denseBlocks = NULL;
		}
	}

void SBlockPent::setBlocks(SymMatrix* symBlocks, DenseMatrix* denseBlocks, unsigned int nblocks) {
	this->symBlocks = symBlocks;
	this->denseBlocks = denseBlocks;
	this->nblocks = nblocks;
}

void SBlockPent::setBlockByClone(SymMatrix* symBlock, unsigned int i) {
	symBlock->cloneTo(this->symBlocks[i]);
}

void SBlockPent::setBlockByClone(DenseMatrix* denseBlock, unsigned int i) {
	denseBlock->cloneTo(this->denseBlocks[i]);
}

void SBlockPent::getBlocks(SymMatrix** symBlocks, DenseMatrix** denseBlocks) {
	*symBlocks = this->symBlocks;
	*denseBlocks = this->denseBlocks;
}

SBlockPent& SBlockPent::operator=(const SBlockPent& obj){

	this->free();

	this->symBlocks = new SymMatrix[obj.nblocks];
	this->denseBlocks = new DenseMatrix[2*obj.nblocks-3];
	for(unsigned int i = 0; i < obj.nblocks-2; i++){
		this->symBlocks[i] = obj.symBlocks[i];
		this->denseBlocks[2*i] = obj.denseBlocks[2*i];
		this->denseBlocks[2*i+1] = obj.denseBlocks[2*i+1];
	}

	this->symBlocks[obj.nblocks-2] = obj.symBlocks[obj.nblocks-2];
	this->denseBlocks[2*(obj.nblocks-2)] = obj.denseBlocks[2*(obj.nblocks-2)];

	this->symBlocks[obj.nblocks-1] = obj.symBlocks[obj.nblocks-1];

	this->nblocks = obj.nblocks;

	return *this;
}

void SBlockPent::symv(double* src, double* dest) {
	//Size of the blocks (for organization)

	const unsigned int bsize = this->symBlocks[0].getNRows();

	//Running over all the blocks in the main diagonal
	for(unsigned int i = 0; i < this->nblocks-2; i++){
		//dest[i*this->nblocks]
		this->symBlocks[i].symv(&src[i*bsize], &dest[i*bsize]);
		this->denseBlocks[2*i].gemv(&src[(i+1)*bsize], &dest[i*bsize]);
		this->denseBlocks[2*i+1].gemv(&src[(i+2)*bsize], &dest[i*bsize]);
		this->denseBlocks[2*i].gemtv(&src[i*bsize], &dest[(i+1)*bsize]);
		this->denseBlocks[2*i+1].gemtv(&src[i*bsize], &dest[(i+2)*bsize]);
	}

	//For the before last block
	this->symBlocks[this->nblocks-2].symv(&src[(this->nblocks-2)*bsize], &dest[(this->nblocks-2)*bsize]);
	this->denseBlocks[2*(this->nblocks-2)].gemv(&src[(this->nblocks-1)*bsize], &dest[(this->nblocks-2)*bsize]);
	this->denseBlocks[2*(this->nblocks-2)].gemtv(&src[(this->nblocks-2)*bsize], &dest[(this->nblocks-1)*bsize]);

	//For the last block
	this->symBlocks[this->nblocks-1].symv(&src[(this->nblocks-1)*bsize], &dest[(this->nblocks-1)*bsize]);

}

void SBlockPent::symv(std::complex<double>* src, std::complex<double>* dest) {
	//Size of the blocks (for organization)

	const unsigned int bsize = this->symBlocks[0].getNRows();

	//Running over all the blocks in the main diagonal
	for(unsigned int i = 0; i < this->nblocks-2; i++){
		//dest[i*this->nblocks]
		this->symBlocks[i].symv(&src[i*bsize], &dest[i*bsize]);
		this->denseBlocks[2*i].gemv(&src[(i+1)*bsize], &dest[i*bsize]);
		this->denseBlocks[2*i+1].gemv(&src[(i+2)*bsize], &dest[i*bsize]);
		this->denseBlocks[2*i].gemtv(&src[i*bsize], &dest[(i+1)*bsize]);
		this->denseBlocks[2*i+1].gemtv(&src[i*bsize], &dest[(i+2)*bsize]);
	}

	//For the before last block
	this->symBlocks[this->nblocks-2].symv(&src[(this->nblocks-2)*bsize], &dest[(this->nblocks-2)*bsize]);
	this->denseBlocks[2*(this->nblocks-2)].gemv(&src[(this->nblocks-1)*bsize], &dest[(this->nblocks-2)*bsize]);
	this->denseBlocks[2*(this->nblocks-2)].gemtv(&src[(this->nblocks-2)*bsize], &dest[(this->nblocks-1)*bsize]);

	//For the last block
	this->symBlocks[this->nblocks-1].symv(&src[(this->nblocks-1)*bsize], &dest[(this->nblocks-1)*bsize]);

}


std::ostream& operator<<(std::ostream& os, const SBlockPent& obj){
	os << "\nThe Matrix Blocks are:\n";
	for(unsigned int i = 0; i < obj.nblocks-2; i++){
		os << "blocks at row ";
		os << i;
		os << " are : ";
		os << obj.symBlocks[i];
		os << "\nand\n";
		os << obj.denseBlocks[2*i];
		os << "\nand\n";
		os << obj.denseBlocks[2*i+1];
		os << "\n";
	}

	os << "blocks at row " << obj.nblocks-2;
	os << " are : ";
	os << obj.symBlocks[obj.nblocks-2];
	os << "\nand\n";
	os << obj.denseBlocks[2*(obj.nblocks-2)];

	os << "blocks at row " << obj.nblocks-1;
	os << " are : ";
	os << obj.symBlocks[obj.nblocks-1];
	os << "\n";

	return os;
}

void SBlockPent::powerMethod(std::complex<double>* v, std::complex<double>* lambda,
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



} /* namespace matrix */
