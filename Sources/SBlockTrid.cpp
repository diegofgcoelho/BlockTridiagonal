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

void SBlockTrid::getBlocks(SymMatrix* symBlocks,
		DenseMatrix* denseMatrix) {
	symBlocks = this->symBlocks;
	denseMatrix = this->denseBlocks;
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

void SBlockTrid::free(){
	if(this->symBlocks != NULL && this->denseBlocks != NULL){
		for(unsigned int i = 0; i < this->nblocks-1; i++){
			this->symBlocks[i].~SymMatrix();
			this->denseBlocks[i].~DenseMatrix();
		}
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

} /* namespace matrix */
