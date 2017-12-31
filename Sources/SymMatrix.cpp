/*
 * SymMatrix.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: diego
 */

#include "SymMatrix.h"

namespace matrix {

SymMatrix::SymMatrix() {
	this->data = NULL;
	this->nrows = 0;
}

SymMatrix::SymMatrix(double* data, unsigned int nrows) {
	this->data = data;
	this->nrows = nrows;
}

SymMatrix::~SymMatrix() {
	this->free();
}

void SymMatrix::cloneTo(SymMatrix& obj) {
	if(obj.data != NULL)	{
		delete [] obj.data;
		obj.data = NULL;
	}

	obj.data = new double[(this->nrows*(this->nrows+1))/2];
	std::memcpy(&obj.data[0], &this->data[0], ((this->nrows*(this->nrows+1))/2)*sizeof(double));
	obj.nrows = this->nrows;
}

void SymMatrix::cloneFrom(const SymMatrix& obj) {
	if(this->data != NULL)	{
		delete [] this->data;
		this->data = NULL;
	}

	this->data = new double[obj.nrows*(obj.nrows+1)/2];
	std::memcpy(&this->data[0], &obj.data[0], obj.nrows*(obj.nrows+1)/2*sizeof(double));
	this->nrows = obj.nrows;
}

void SymMatrix::smartClone(const SymMatrix& obj) {
	this->data = obj.data;
	this->nrows = obj.nrows;
}

void SymMatrix::setData(double* data, unsigned int nrows) {
	this->data = data;
	this->nrows = nrows;
}

void SymMatrix::setEle(unsigned int i, unsigned int j, double ele){
	if(i < j){
		unsigned int temp = i;
		i = j;
		j = temp;
	}
	//Assuming that i > j
	unsigned int datapos = (2*nrows-i)*(i+1)/2+i-j;
	this->data[datapos] = ele;
}

double SymMatrix::getEle(unsigned int i, unsigned int j){
	if(i > j){
		unsigned int temp = i;
		i = j;
		j = temp;
	}
	//Assuming that j > i
	unsigned int datapos = ((2*this->nrows-i+2)*(i+1))/2+j-i-(this->nrows+1);
	return this->data[datapos];
}

SymMatrix& SymMatrix::operator =(const SymMatrix& obj) {
	this->free();

	this->data = new double[obj.nrows*(obj.nrows+1)/2];
	std::memcpy(&this->data[0], &obj.data[0], obj.nrows*(obj.nrows+1)/2*sizeof(double));

	this->nrows = obj.nrows;
	return *this;
}

DenseMatrix SymMatrix::operator *(DenseMatrix& rightMatrix) {

	double *result = new double[this->nrows*rightMatrix.ncols];

	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = 0; j < rightMatrix.ncols; j++){
			result[j+this->nrows*i] = 0.0;
			unsigned int BRowpos = i;
			for(unsigned int k = 0; k < i; k++){
				result[j+this->nrows*i] += this->data[BRowpos]*rightMatrix.data[this->nrows*k+j];
				BRowpos = BRowpos+(this->nrows-1-k);
			}
			for(unsigned int k = i; k < rightMatrix.nrows; k++){
				result[j+this->nrows*i] += this->data[BRowpos]*rightMatrix.data[this->nrows*k+j];
				BRowpos++;
			}
		}
	}

	DenseMatrix C(result, this->nrows, rightMatrix.ncols);
	return C;
}

DenseMatrix SymMatrix::operator*(SymMatrix& rightMatrix){

	double* result = new double[this->nrows*this->nrows];

	unsigned int CEleLower = 0;
	unsigned int CEleUpper = 0;
	for(unsigned int i = 0; i < this->nrows-1; i++){

		CEleLower+=i;

		//Equivalent to j = i
		result[CEleLower] = 0.0;
		unsigned int ARowpos = i;
		unsigned int BColpos = i;
		for(unsigned int k = 0; k < i; k++){
			result[CEleLower] += rightMatrix.data[ARowpos]*this->data[ARowpos];
			ARowpos = ARowpos+(this->nrows-1-k);
			BColpos = BColpos+(this->nrows-1-k);
		}
		for(unsigned int k = i; k < this->nrows; k++){
			result[CEleLower] += rightMatrix.data[ARowpos]*this->data[ARowpos];
			ARowpos++;
			BColpos++;
		}
		CEleLower++;
		CEleUpper = CEleLower+this->nrows-1;
		result[CEleUpper] = 0.0;
		for(unsigned int j = i+1; j < this->nrows; j++){
			result[CEleLower] = 0.0;
			result[CEleUpper] = 0.0;
			unsigned int ARowpos = i;
			unsigned int BColpos = j;
			for(unsigned int k = 0; k < i; k++){
				result[CEleUpper] += rightMatrix.data[ARowpos]*this->data[BColpos];
				result[CEleLower] += rightMatrix.data[BColpos]*this->data[ARowpos];
				ARowpos = ARowpos+(this->nrows-1-k);
				BColpos = BColpos+(this->nrows-1-k);
			}
			for(unsigned int k = i; k < j; k++){
				result[CEleUpper] += rightMatrix.data[ARowpos]*this->data[BColpos];
				result[CEleLower] += rightMatrix.data[BColpos]*this->data[ARowpos];
				ARowpos++;
				BColpos = BColpos+(this->nrows-1-k);
			}
			for(unsigned int k = j; k < this->nrows; k++){
				result[CEleUpper] += rightMatrix.data[ARowpos]*this->data[BColpos];
				result[CEleLower] += rightMatrix.data[BColpos]*this->data[ARowpos];
				ARowpos++;
				BColpos++;
			}
			CEleLower++;
			CEleUpper+=this->nrows;
		}
	}

	//Case i = this->nrows-1
	CEleLower = CEleLower+this->nrows-1;

	//Equivalent to j = i
	result[CEleLower] = 0.0;
	unsigned int ARowpos = this->nrows-1;
	unsigned int BColpos = this->nrows-1;
	for(unsigned int k = 0; k < this->nrows-1; k++){
		result[CEleLower] += rightMatrix.data[ARowpos]*this->data[ARowpos];
		ARowpos = ARowpos+(this->nrows-1-k);
		BColpos = BColpos+(this->nrows-1-k);
	}
	for(unsigned int k = this->nrows-1; k < this->nrows; k++){
		result[CEleLower] += rightMatrix.data[ARowpos]*this->data[ARowpos];
		ARowpos++;
		BColpos++;
	}

	DenseMatrix C(result, this->nrows, this->nrows);
	return C;
}

SymMatrix SymMatrix::operator+(const SymMatrix& rightMatrix){
	double* result = new double[this->nrows*(this->nrows+1)/2]();
	unsigned int resultPos = 0;
	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = i; j < this->nrows; j++){
			result[resultPos] = this->data[resultPos]+rightMatrix.data[resultPos];
			resultPos++;
		}
	}
	SymMatrix C(result, this->nrows);
	return C;
}

SymMatrix SymMatrix::square(){
	double* result = new double[this->nrows*(this->nrows+1)/2]();
//
//	unsigned int CEleLower = 0;
//	unsigned int CEleUpper = 0;
//	for(unsigned int i = 0; i < this->nrows-1; i++){
//
//		CEleLower+=i;
//
//		//Equivalent to j = i
//		result[CEleLower] = 0.0;
//		unsigned int ARowpos = i;
//		unsigned int BColpos = i;
//		for(unsigned int k = 0; k < i; k++){
//			result[CEleLower] += this->data[ARowpos]*this->data[ARowpos];
//			ARowpos = ARowpos+(this->nrows-1-k);
//			BColpos = BColpos+(this->nrows-1-k);
//		}
//		for(unsigned int k = i; k < this->nrows; k++){
//			result[CEleLower] += this->data[ARowpos]*this->data[ARowpos];
//			ARowpos++;
//			BColpos++;
//		}
//		CEleLower++;
//		CEleUpper = CEleLower+this->nrows-1;
//		for(unsigned int j = i+1; j < this->nrows; j++){
//			result[CEleLower] = 0.0;
//			unsigned int ARowpos = i;
//			unsigned int BColpos = j;
//			for(unsigned int k = 0; k < i; k++){
//				result[CEleLower] += this->data[BColpos]*this->data[ARowpos];
//				ARowpos = ARowpos+(this->nrows-1-k);
//				BColpos = BColpos+(this->nrows-1-k);
//			}
//			for(unsigned int k = i; k < j; k++){
//				result[CEleLower] += this->data[BColpos]*this->data[ARowpos];
//				ARowpos++;
//				BColpos = BColpos+(this->nrows-1-k);
//			}
//			for(unsigned int k = j; k < this->nrows; k++){
//				result[CEleLower] += this->data[BColpos]*this->data[ARowpos];
//				ARowpos++;
//				BColpos++;
//			}
//			CEleLower++;
//			CEleUpper+=this->nrows;
//		}
//	}
//
//	//Case i = this->nrows-1
//	CEleLower = CEleLower+this->nrows-1;
//
//	//Equivalent to j = i
//	result[CEleLower] = 0.0;
//	unsigned int ARowpos = this->nrows-1;
//	unsigned int BColpos = this->nrows-1;
//	for(unsigned int k = 0; k < this->nrows-1; k++){
//		result[CEleLower] += this->data[ARowpos]*this->data[ARowpos];
//		ARowpos = ARowpos+(this->nrows-1-k);
//		BColpos = BColpos+(this->nrows-1-k);
//	}
//	for(unsigned int k = this->nrows-1; k < this->nrows; k++){
//		result[CEleLower] += this->data[ARowpos]*this->data[ARowpos];
//		ARowpos++;
//		BColpos++;
//	}


	unsigned int rPos = 0;
	unsigned int RowPos = 0;
	unsigned int ColPos = 0;

	for(unsigned int i = 0; i < this->nrows; i++){
		for(unsigned int j = i; j < this->nrows; j++){

//			std::cout << "i = " << i << " j = " << j << std::endl;
			RowPos = i;
			ColPos = j;
			for(unsigned int k = 0; k < i; k++){

				result[rPos] += this->data[RowPos]*this->data[ColPos];

				RowPos = RowPos+this->nrows-k-1;
				ColPos = ColPos+this->nrows-k-1;
			}

			for(unsigned int k = i; k < j; k++){

				result[rPos] += this->data[RowPos]*this->data[ColPos];

				RowPos++;
				ColPos = ColPos+this->nrows-k-1;
			}

			for(unsigned int k = j; k < this->nrows; k++){

				result[rPos] += this->data[RowPos]*this->data[ColPos];

				RowPos++;
				ColPos++;
			}


			rPos++;
		}
	}


	SymMatrix C(result, this->nrows);
	return C;
}

void SymMatrix::symv(double* src, double* dest){
	unsigned int AEleUpper = 0;
	for(unsigned int i = 0; i < this->nrows; i++){

		dest[i] += this->data[AEleUpper]*src[i];
		AEleUpper++;

		for(unsigned int j = i+1; j < this->nrows; j++){
			dest[i] += this->data[AEleUpper]*src[j];
			dest[j] += this->data[AEleUpper]*src[i];
			AEleUpper++;
		}
	}
}

void SymMatrix::symv(std::complex<double>* src, std::complex<double>* dest){
	unsigned int AEleUpper = 0;
	for(unsigned int i = 0; i < this->nrows; i++){

		dest[i] += this->data[AEleUpper]*src[i];
		AEleUpper++;

		for(unsigned int j = i+1; j < this->nrows; j++){
			dest[i] += this->data[AEleUpper]*src[j];
			dest[j] += this->data[AEleUpper]*src[i];
			AEleUpper++;
		}
	}
}


std::ostream& operator<<(std::ostream& os, const SymMatrix& obj){
	unsigned int pos = 0;
	os << "\n[";
	for(unsigned int i = 0; i <  obj.nrows; i++){
		for(unsigned int j = 0; j <  i; j++){
			os << "x.xxxxx ";
		}

		for(unsigned int j = i; j <  obj.nrows; j++){
			os << obj.data[pos];
			pos++;
			if(j != obj.nrows-1) os << " ";
		}
		if(i != obj.nrows-1) os << "\n";
	}
	os << "]\n";
	return os;
}

} /* namespace matrix */
