#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
	double **data;

    // Parameterized constructor
	Matrix();
	Matrix(const int v_size);
    Matrix(const int n_row, const int n_column);
	
	// Member operators
	double& operator()(const int n);
    double& operator()(const int row, const int column);
	Matrix& operator+(Matrix &m);
    Matrix& operator+(const double val);
    Matrix& operator-(Matrix &m);
    Matrix& operator-(const double val);
    Matrix& operator*(Matrix &m);
    Matrix& operator*(const double val);
    Matrix& operator/(Matrix &m);
    Matrix& operator/(const double val);
    Matrix& operator=(Matrix m);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
	
	double norm();
    double dot(Matrix &m);
    Matrix extract_vector(const int index, bool row);
    Matrix extract_row(const int r);
    Matrix extract_column(const int c);
    void assign_row(const int r, Matrix &rowMat);
    void assign_column(const int c, Matrix &colMat);
    Matrix union_vector(Matrix &v, bool horizontal);
    Matrix cross(Matrix &m);
	Matrix eye(const int n);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix transpose(Matrix &m);
Matrix inv(Matrix &m);
Matrix& zeros(const int n);
Matrix& zeros(const int n_row, const int n_column);

#endif