/*
 * SBlockTrid.h
 *
 *  Created on: Dec 20, 2017
 *      Author: diego
 */

#ifndef SBLOCKTRID_H_
#define SBLOCKTRID_H_
#include <cstring>
#include "DenseMatrix.h"
#include "SymMatrix.h"
#include "SBlockPent.h"

namespace matrix {

//class SBlockPent;

/*
 * This class represents a symmetric block tridiagonal matrix. The size of the matrix is represented by
 * the size of the first block (number of rows) multiplied by the attribute nblocks (the number of blocks in
 * the main diagonal).
 */

class SBlockTrid {
public:
	SBlockTrid();
	SBlockTrid(SymMatrix* symBlocks, DenseMatrix* DenseMatrix, unsigned int nblocks);
	virtual ~SBlockTrid();
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
	void getBlocks(SymMatrix* symBlocks, DenseMatrix* DenseMatrix);
	inline unsigned int getNRows(){
		return this->nblocks*this->symBlocks->getNRows();
	}
	SBlockPent square();
	void symv(double* src, double* dest);
	void free();
	friend std::ostream& operator<<(std::ostream& os, const SBlockTrid& obj);
	SBlockTrid& operator=(const SBlockTrid& obj);
private:
	unsigned int nblocks;
	SymMatrix* symBlocks;
	DenseMatrix* denseBlocks;
};

void testPrintOperatorSBlockTrid();
void testSetBlocksSBlockTrid();
void testAssignOperatorSBlockTrid();
void testSquareSBlockTridMatlab();
void testSquareSBlockTridEigen();
void testSymMatrixXVectorPlusDenseMatrixXVector();
void testSymvTrid();

} /* namespace matrix */

#endif /* SBLOCKTRID_H_ */
