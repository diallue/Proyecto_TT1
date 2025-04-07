#include "..\include\matrix.h"

Matrix::Matrix(const int v_size) {
	if (v_size < 0) {
		cout << "Matrix: Error in v_size\n";
		exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = v_size;
	this->data = (double**)malloc(n_row*sizeof(double));
	
	if (this->data == NULL) {
		cout << "Matrix: Error in data\n";
		exit(EXIT_FAILURE);
	}
	
	this->data[v_size] = (double*)calloc(v_size, sizeof(double));
}

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

double& Matrix::operator() (const int n) {
	if (n <= 0 || n >this->n_column*this->n_row) {
		cout << "Matrix get: Error in n\n";
		exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n - 1)%this->n_column];
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& Matrix::operator * (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix mult: error, columnas de A != filas de B\n";
		exit(EXIT_FAILURE);
	}

	Matrix* result = new Matrix(this->n_row, m.n_column);

	for (int i = 1; i <= this->n_row; i++) {
		for (int j = 1; j <= m.n_column; j++) {
			double sum = 0.0;
			for (int k = 1; k <= this->n_column; k++) {
				sum += (*this)(i, k) * m(k, j);
			}
			(*result)(i, j) = sum;
		}
	}

	return *result;
}

Matrix& Matrix::operator / (Matrix &m) {
	if (m.n_row != m.n_column) {
		cout << "Matrix div: error, matriz no cuadrada\n";
		exit(EXIT_FAILURE);
	}

	Matrix inverse = m.inversa();

	if (this->n_column != inverse.n_row) {
		cout << "Matrix div: error, dimensiones incompatibles\n";
		exit(EXIT_FAILURE);
	}

	Matrix& result = (*this) * inverse;

	return result;
}


Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n, n_column);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}