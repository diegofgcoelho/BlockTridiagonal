/*
 * DenseMatrix.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#include "DenseMatrix.h"

namespace matrix {

DenseMatrix::DenseMatrix() {
	this->data = NULL;
	this->nrows = 0;
	this->ncols = 0;
}

DenseMatrix::DenseMatrix(double* data, unsigned int nrows, unsigned int ncols) {
	this->data = data;
	this->nrows = nrows;
	this->ncols = ncols;
}

DenseMatrix::~DenseMatrix() {
	this->free();
}

void DenseMatrix::cloneTo(DenseMatrix& obj){

	if(obj.data != NULL)	{
		delete [] obj.data;
		obj.data = NULL;
	}

	obj.data = new double[this->nrows*this->ncols];
	std::memcpy(&obj.data[0], &this->data[0], this->nrows*this->ncols*sizeof(double));
	obj.nrows = this->nrows;
	obj.ncols = this->ncols;
}

void DenseMatrix::cloneTransposeTo(DenseMatrix& obj){

	if(obj.data != NULL)	{
		delete [] obj.data;
		obj.data = NULL;
	}

	obj.data = new double[this->nrows*this->ncols];
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			obj.data[i*this->nrows+j] = this->data[j*this->nrows+i];
		}
	}
	obj.nrows = this->nrows;
	obj.ncols = this->ncols;
}

void DenseMatrix::cloneFrom(const DenseMatrix& obj){

	if(this->data != NULL)	{
		delete [] this->data;
		this->data = NULL;
	}

	this->data = new double[obj.nrows*obj.ncols];
	std::memcpy(&this->data[0], &obj.data[0], obj.nrows*obj.ncols*sizeof(double));
	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
}

void DenseMatrix::cloneTransposeFrom(const DenseMatrix& obj){

	if(this->data != NULL)	{
		delete [] this->data;
		this->data = NULL;
	}

	this->data = new double[obj.nrows*obj.ncols];
	for(unsigned int i = 0; i < obj.nrows; i++){
		for(unsigned int j = 0; j < obj.ncols; j++){
			this->data[i*obj.nrows+j] = obj.data[j*obj.nrows+i];
		}
	}
	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
}


void DenseMatrix::smartClone(const DenseMatrix& obj){
	this->data = obj.data;
	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
}

void DenseMatrix::setData(double* data, unsigned int nrows,
		unsigned int ncols) {
	this->data = data;
	this->nrows = nrows;
	this->ncols = ncols;
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& obj){
	this->free();

	this->data = new double[obj.nrows*obj.ncols];
	std::memcpy(&this->data[0], &obj.data[0], obj.nrows*obj.ncols*sizeof(double));

	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
	return *this;
}

DenseMatrix DenseMatrix::operator *(DenseMatrix& rightMatrix) {
	//Data for the resulting matrix
	double *result = new double[this->nrows*rightMatrix.ncols];

	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < rightMatrix.ncols; j++){
			result[j+this->nrows*i] = 0.0;
			for(unsigned int k = 0; k < this->nrows; k++){
				result[j+this->nrows*i] += this->data[this->nrows*i+k]*rightMatrix.data[j+this->nrows*k];
			}
		}
	}

	//Output matrix
	DenseMatrix C(result, this->nrows, rightMatrix.ncols);
	return C;
}

DenseMatrix DenseMatrix::operator *(SymMatrix& rightMatrix) {
	//Data for the resulting matrix
	double *result = new double[this->nrows*rightMatrix.nrows];

	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < rightMatrix.nrows; j++){
			result[j+this->nrows*i] = 0.0;
			unsigned int BColpos = j;
			for(unsigned int k = 0; k < j; k++){
				result[j+this->nrows*i] += this->data[this->nrows*i+k]*rightMatrix.data[BColpos];
				BColpos = BColpos+(rightMatrix.nrows-1-k);
			}
			for(unsigned int k = j; k < rightMatrix.nrows; k++){
				result[j+this->nrows*i] += this->data[this->nrows*i+k]*rightMatrix.data[BColpos];
				BColpos++;
			}
		}
	}

	//Output matrix
	DenseMatrix C(result, this->nrows, rightMatrix.nrows);
	return C;
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix& rightMatrix){
	double* result = new double[this->nrows*this->ncols]();
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			result[i*this->nrows+j] = this->data[i*this->nrows+j]+rightMatrix.data[i*this->nrows+j];
		}
	}
	DenseMatrix C(result, this->nrows, this->ncols);
	return C;
}

SymMatrix DenseMatrix::AAT(){
	double* result = new double[this->nrows*(this->nrows+1)/2]();

	unsigned int resultPos = 0;
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = i; j < this->ncols; j++){

			for(unsigned int k = 0; k < this->ncols; k++){
				result[resultPos] += this->data[i*this->nrows+k]*this->data[j*this->nrows+k];
			}
			resultPos++;
		}
	}

	SymMatrix C(result, this->nrows);
	return C;
}

SymMatrix DenseMatrix::ATA(){
	double* result = new double[this->nrows*(this->nrows+1)/2]();

	unsigned int resultPos = 0;
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = i; j < this->ncols; j++){

			for(unsigned int k = 0; k < this->ncols; k++){
				result[resultPos] += this->data[i+this->nrows*k]*this->data[j+this->nrows*k];
			}
			resultPos++;
		}
	}

	SymMatrix C(result, this->nrows);
	return C;
}

void DenseMatrix::gemv(double* src, double* dest){
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			dest[i]+=this->data[i*this->nrows+j]*src[j];
		}
	}
}

void DenseMatrix::gemtv(double* src, double* dest){
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			dest[i]+=this->data[j*this->nrows+i]*src[j];
		}
	}
}


std::ostream& operator<<(std::ostream& os, const DenseMatrix& obj){
	os << "\n[";
	for(unsigned int i = 0; i <  obj.nrows; i++){
		for(unsigned int j = 0; j <  obj.ncols; j++){
			os << obj.data[i*obj.ncols+j];
			if(j != obj.ncols-1) os << " ";
		}
		if(i != obj.nrows-1) os << "\n";
	}
	os << "]\n";
	return os;
}

} /* namespace matrix */

