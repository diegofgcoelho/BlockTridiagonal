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
	this->data = new Data;
	this->data->data = data;
	this->data->counter = 1;
	this->nrows = nrows;
	this->ncols = ncols;
}

DenseMatrix::~DenseMatrix() {
	this->free();
}

void DenseMatrix::cloneTo(DenseMatrix& obj){
	//If obj matrix is empty, just return nothing
	if(this->data == NULL) return;
	//Free the memory of the rights-side object
	obj.free();
	//Allocate new memory
	obj.data = new Data;
	obj.data->counter = 1;
	obj.data->data = new double[this->nrows*this->ncols];
	std::memcpy(&obj.data->data[0], &this->data->data[0], this->nrows*this->ncols*sizeof(double));
	obj.nrows = this->nrows;
	obj.ncols = this->ncols;
}

void DenseMatrix::cloneTransposeTo(DenseMatrix& obj){
	//If this matrix is empty, just return nothing
	if(this->data == NULL) return;
	//Free the memory of the rights-side object
	obj.free();
	//Allocate new memory
	obj.data = new Data;
	obj.data->counter = 1;
	obj.data->data = new double[this->nrows*this->ncols];
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			obj.data->data[i*this->nrows+j] = this->data->data[j*this->nrows+i];
		}
	}
	obj.nrows = this->nrows;
	obj.ncols = this->ncols;
}

void DenseMatrix::cloneFrom(const DenseMatrix& obj){
	//If obj matrix is empty, just return nothing
	if(obj.data == NULL) return;
	//Free the memory of the rights-side object
	this->free();
	//Allocate new memory
	this->data = new Data;
	this->data->counter = 1;
	this->data->data = new double[obj.nrows*obj.ncols];
	std::memcpy(&this->data->data[0], &obj.data->data[0], obj.nrows*obj.ncols*sizeof(double));
	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
}

void DenseMatrix::cloneTransposeFrom(const DenseMatrix& obj){
	//If obj matrix is empty, just return nothing
	if(obj.data == NULL) return;
	//Free the memory of the rights-side object
	this->free();
	//Allocate new memory
	this->data = new Data;
	this->data->counter = 1;
	this->data->data = new double[obj.nrows*obj.ncols];
	for(unsigned int i = 0; i < obj.nrows; i++){
		for(unsigned int j = 0; j < obj.ncols; j++){
			this->data->data[i*obj.nrows+j] = obj.data->data[j*obj.nrows+i];
		}
	}
	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
}


void DenseMatrix::smartClone(const DenseMatrix& obj){
	this->data = obj.data;
	this->data->counter++;
	this->nrows = obj.nrows;
	this->ncols = obj.ncols;
}

void DenseMatrix::setData(double* data, unsigned int nrows,
		unsigned int ncols) {
	this->free();
	//Allocate new memory
	this->data = new Data;
	this->data->counter = 1;
	this->data->data = data;
	this->nrows = nrows;
	this->ncols = ncols;
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& obj){
	this->free();

	this->data = obj.data;
	obj.data->counter++;
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
				result[j+this->nrows*i] += this->data->data[this->nrows*i+k]*rightMatrix.data->data[j+this->nrows*k];
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
				result[j+this->nrows*i] += this->data->data[this->nrows*i+k]*rightMatrix.data->data[BColpos];
				BColpos = BColpos+(rightMatrix.nrows-1-k);
			}
			for(unsigned int k = j; k < rightMatrix.nrows; k++){
				result[j+this->nrows*i] += this->data->data[this->nrows*i+k]*rightMatrix.data->data[BColpos];
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
			result[i*this->nrows+j] = this->data->data[i*this->nrows+j]+rightMatrix.data->data[i*this->nrows+j];
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
				result[resultPos] += this->data->data[i*this->nrows+k]*this->data->data[j*this->nrows+k];
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
				result[resultPos] += this->data->data[i+this->nrows*k]*this->data->data[j+this->nrows*k];
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
			dest[i]+=this->data->data[i*this->nrows+j]*src[j];
		}
	}
}

void DenseMatrix::gemv(std::complex<double>* src, std::complex<double>* dest){
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			dest[i]+=this->data->data[i*this->nrows+j]*src[j];
		}
	}
}

void DenseMatrix::gemtv(double* src, double* dest){
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			dest[i]+=this->data->data[j*this->nrows+i]*src[j];
		}
	}
}

void DenseMatrix::gemtv(std::complex<double>* src, std::complex<double>* dest){
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < this->ncols; j++){
			dest[i]+=this->data->data[j*this->nrows+i]*src[j];
		}
	}
}


std::ostream& operator<<(std::ostream& os, const DenseMatrix& obj){
	if (obj.data != NULL) {
		os << "\n[";
		for(unsigned int i = 0; i <  obj.nrows; i++){
			for(unsigned int j = 0; j <  obj.ncols; j++){
				os << obj.data->data[i*obj.ncols+j];
				if(j != obj.ncols-1) os << " ";
			}
			if(i != obj.nrows-1) os << "\n";
		}
		os << "]\n";
	} else os << "This matrix is empty\n";
	return os;
}

} /* namespace matrix */

