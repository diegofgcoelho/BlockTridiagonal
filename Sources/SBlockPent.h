/*
 * SBlockPent.h
 *
 *  Created on: Dec 20, 2017
 *      Author: diego
 */

#ifndef SBLOCKPENT_H_
#define SBLOCKPENT_H_

#include <iostream>
#include <cstring>
#include <complex>
#include "SymMatrix.h"
#include "DenseMatrix.h"
#include "Support.h"

namespace matrix {

/*
 * This class represents a symmetric block pentadiagonal matrix that is the result of the squaring
 * of a block tridigonal matrix. The size of the matrix is represented by
 * the size of the first block (number of rows) multiplied by the attribute nblocks (the number of blocks in
 * the main diagonal).
 */

class SBlockPent {
public:
	SBlockPent();
	SBlockPent(SymMatrix* symBlocks, DenseMatrix* denseBlocks, unsigned int nblocks);
	virtual ~SBlockPent();
	/*
	 * Set the internal blocks by pointing the denseBlocks to the DenseMatrix pointer passed as
	 * argument and the symBlocks to SymMatrix passed as argument.
	 */
	void setBlocks(SymMatrix* symBlocks, DenseMatrix* denseBlocks, unsigned int nblocks);
	/*
	 * These methods set the symBlocks and denseBlocks by deep copy of the blocks in the array of
	 * the input arguments.
	 */
	void setBlockByClone(SymMatrix* symBlock, unsigned int i);
	void setBlockByClone(DenseMatrix* denseBlock, unsigned int i);
	/*
	 * This method copy the memory position for the attributes of this object to the input arguments.
	 * Considering to remove this method for possibly violating OO principles.
	 */
	void getBlocks(SymMatrix** symBlocks, DenseMatrix** denseBlocks);
	inline unsigned int getNRows(){
		return this->nblocks*this->symBlocks->getNRows();
	}
	SBlockPent& operator=(const SBlockPent& obj);
	friend std::ostream& operator<<(std::ostream& os, const SBlockPent& obj);
	void symv(double* src, double* dest);
	void symv(std::complex<double>* src, std::complex<double>* dest);
	void free();
	/*
	 * This method performs the power iteration and returns the largest eignevalue associated with the matrix
	 * represented by this object. It also returns its associated eigenvector and norm of the difference of the
	 * two last eigenvector estimates.
	 */
	void powerMethod(std::complex<double>* v, std::complex<double>* lambda, double prec, unsigned int maxIter, unsigned int* iter, double* dnorm);
private:
	unsigned int nblocks;
	SymMatrix* symBlocks;
	DenseMatrix* denseBlocks;
};

void testPrintOperatorSBlockPent();
void testSetBlocksSBlockPent();
void testAssignOperatorSBlockPent();
void testSymvPent();
void testSymvPentComplex();
void testPowerMethodSBlockPent();
void testPowerMethodSBlockPent2();

} /* namespace matrix */

#endif /* SBLOCKPENT_H_ */
