/*
 * DenseMatrix.h
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <iostream>
#include <cstring>
#include <complex>
#include "SymMatrix.h"

namespace matrix {

class SymMatrix;

/*
 * DenseMatrix class represents a dense general matrix with no symmetries whose entries are in double format and are
 * organized in row major format.
 *
 * The matrix coeffcieints are stores as an linear array of doubles that is dynamically allocated. The object has
 * hold of its matrix coeffcients by storing the position of the first matrix coefficient.
 */

class DenseMatrix {
public:
	friend class SymMatrix;
	DenseMatrix();
	DenseMatrix(double* data, unsigned int nrows, unsigned int ncols);
	virtual ~DenseMatrix();
	/*
	 * The methods cloneTo and cloneFrom provide deep copy to and from the input argument.
	 */
	void cloneTo(DenseMatrix& obj);
	void cloneTransposeTo(DenseMatrix& obj);
	void cloneFrom(const DenseMatrix& obj);
	void cloneTransposeFrom(const DenseMatrix& obj);
	/*
	 * This method smartClone provides soft copy: the nrows and ncols attributes are copied,
	 * but the data attribute will point to the data of the input argument.
	 */
	void smartClone(const DenseMatrix& obj);
	inline unsigned int getNRows(){
		return this->nrows;
	}
	inline unsigned int getNCols(){
		return this->ncols;
	}
	/*
	 * Set the data stored by the object. The array data must be allocated by the calling method.
	 */
	void setData(double* data, unsigned int nrows, unsigned int ncols);
	/*
	 * Set and get methods for a particular matrix element.
	 */
	inline void setEle(unsigned int i, unsigned int j, double ele){
		this->data->data[i*ncols+j] = ele;
	}
	inline double getEle(unsigned int i, unsigned int j){
		return this->data->data[i*ncols+j];
	}
	/*
	 * The = operator is similar to the CloneFrom method. It provides deep copy.
	 */
	DenseMatrix& operator=(const DenseMatrix& obj);
	/*
	 * Overloaded operator that performs the multiplication of DesenMatrix objects
	 */
	DenseMatrix operator*(DenseMatrix& rightMatrix);
	/*
	 * Overloaded operator that performs the multiplication of DesenMatrix by SymMatrix object
	 */
	DenseMatrix operator*(SymMatrix& rightMatrix);
	/*
	 * Overloaded operator that performs the addition of two DenseMatrix objects
	 */
	DenseMatrix operator+(const DenseMatrix& rightMatrix);
	/*
	 * This method performs the matrix multicplication A*A^T and returns a symMatrix object
	 */
	SymMatrix AAT();
	/*
	 * This method performs the matrix multicplication A^T*A and returns a symMatrix object
	 */
	SymMatrix ATA();

	/*
	 * This method performs the matrix-vector multiplication (and the transposed matrix). The space for the output argument,
	 * dest, must be already allocated.
	 */
	void gemv(double* src, double* dest);
	void gemv(std::complex<double>* src, std::complex<double>* dest);
	void gemtv(double* src, double* dest);
	void gemtv(std::complex<double>* src, std::complex<double>* dest);
	/*
	 * Overloaded operator prints the matrix elements
	 */
	friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& obj);
	inline unsigned int checkCounter(){
		return this->data->counter;
	}
protected:
	inline void free(){
		//If the structu is not empty
		if(this->data != NULL){
			//Decrement the counter by one
			this->data->counter--;
			//If the data inside the structure is null, do nothing, if not, deallocate
			if(this->data->data != NULL){
				if (this->data->counter == 0) {
					delete [] this->data->data;
					delete this->data;
				}
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
	 * The number of rows and columns, respectively.
	 */
	unsigned int nrows;
	unsigned int ncols;
};

void testPrintOperatorDenseMatrix();
void testAssignOperatorDenseMatrix();
void testCloneToDenseMatrix();
void testCloneFromDenseMatrix();
void testMultiplyOperatorDenseMatrixXDenseMatrix();
void testMultiplyOperatorDenseMatrixXSymMatrix();
void testDenseMatrixXVector();
void testDenseMatrixXComplexVector();
void testDenseMatrixTransposedXVector();
void testDenseMatrixTransposedXComplexVector();
void testAAT();
void testATA();
void testAddOperatorDenseMatrixXDenseMatrix();
void testCloneTransposeToDenseMatrix();
void testCloneTransposeFromDenseMatrix();

} /* namespace matrix */



#endif /* DENSEMATRIX_H_ */
