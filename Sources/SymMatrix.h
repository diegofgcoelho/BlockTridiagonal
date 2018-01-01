/*
 * SymMatrix.h
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#ifndef SYMMATRIX_H_
#define SYMMATRIX_H_

#include <cstring>
#include <complex>
#include "Support.h"
#include "DenseMatrix.h"

namespace matrix {

class DenseMatrix;

/*
 * SymMatrix class represents a dense symetric matrix whose entries are in double format and are
 * organized in row major format. Just the upper triangular part is stored.
 *
 * The matrix coeffcieints are stores as an linear array of doubles that is dynamically allocated. The object has
 * hold of its matrix coeffcients by storing the position of the first matrix coefficient.
 */


class SymMatrix {
public:
	friend class DenseMatrix;
	SymMatrix();
	SymMatrix(double* data, unsigned int nrows);
	virtual ~SymMatrix();
	/*
	 * The methods cloneTo and cloneFrom provide deep copy to and from the input argument.
	 */
	void cloneTo(SymMatrix& obj);
	void cloneFrom(const SymMatrix& obj);
	/*
	 * This method smartClone provides soft copy: the nrows and ncols attributes are copied,
	 * but the data attribute will point to the data of the input argument.
	 */
	void smartClone(const SymMatrix& obj);

	inline unsigned int getNRows(){
		return this->nrows;
	}
	/*
	 * Set the data stored by the object. The array data must be allocated by the calling method.
	 */
	void setData(double* data, unsigned int nrows);
	/*
	 * Set and get methods for a particular matrix element.
	 */
	void setEle(unsigned int i, unsigned int j, double ele);
	double getEle(unsigned int i, unsigned int j);
	/*
	 * The = operator is similar to the CloneFrom method. It provides deep copy.
	 */
	SymMatrix& operator=(const SymMatrix& obj);
	/*
	 * Overloaded operator performs the multiplication of DenseMatrix objects
	 */
	DenseMatrix operator*(DenseMatrix& rightMatrix);
	/*
	 * Overloaded operator performs the multiplication of SymMatrix objects
	 */
	DenseMatrix operator*(SymMatrix& rightMatrix);
	/*
	 * Overloaded operator that performs the addition of two SymMatrix objects
	 */
	SymMatrix operator+(const SymMatrix& rightMatrix);
	/*
	 * This method computes the matrix square of this SymMatrix and return another SymMatrix object.
	 */
	SymMatrix square();
	/*
	 * This method performs the matrix-vector multiplication. The space for the output argument,
	 * dest, must be already allocated.
	 */
	void symv(double* src, double* dest);
	void symv(std::complex<double>* src, std::complex<double>* dest);
	/*
	 * Overloaded operator prints the matrix elements
	 */
	friend std::ostream& operator<<(std::ostream& os, const SymMatrix& obj);
	inline bool isEmpty(){
		if(this->data == NULL) return true; else return false;
	}
	inline unsigned int checkCounter(){
		return this->data->counter;
	}
protected:
	inline void free(){
		//If the structure is not empty
		if(this->data != NULL){
			//Decrement the counter by one
			this->data->counter--;
			//If the data inside the structure is null, do nothing, if not, deallocate
			if(this->data->data != NULL){
				if (this->data->counter == 0) {
					delete [] this->data->data;
				}
				delete this->data;
				this->data = NULL;
			}
		}
	}
private:
	/*
	 * Pointer to Data structure
	 */
	Data* data;
	/*
	 * The number of rows and columns (because the matrix is symmetric, just one is needed).
	 */
	unsigned int nrows;
};

void testPrintOperatorSymMatrix();
void testAssignOperatorSymMatrix();
void testSetSymMatrix();
void testGetSymMatrix();
void testCloneToSymMatrix();
void testCloneFromSymMatrix();
void testMultiplyOperatorSymMatrixXDenseMatrix();
void testMultiplyOperatorSymMatrixXSymMatrix();
void testSymMatrixXVector();
void testSymMatrixXComplexVector();
void testAddOperatorSymMatrixXSymMatrix();
void testSquareSymMatrix();

} /* namespace matrix */

#endif /* SYMMATRIX_H_ */
