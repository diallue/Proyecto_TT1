#include "..\include\matrix.hpp"

Matrix::Matrix() {
	this->n_row = 0;
	this->n_column = 0;
	this->data = nullptr;
}

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

double& Matrix::operator () (const int n) {
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

Matrix Matrix::operator + (Matrix &m) {
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

Matrix& Matrix::operator + (const double val) {
    Matrix* res = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            (*res)(i,j) = (*this)(i,j) + val;
        }
    }
    return *res;
}

Matrix Matrix::operator - (Matrix &m) {
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

Matrix& Matrix::operator - (const double val) {
    Matrix* res = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            (*res)(i,j) = (*this)(i,j) - val;
        }
    }
    return *res;
}

ostream& operator << (ostream &o, Matrix m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix Matrix::operator * (const Matrix &m) {
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

Matrix& Matrix::operator * (const double val) {
    Matrix* res = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            (*res)(i,j) = (*this)(i,j) * val;
        }
    }
    return *res;
}

Matrix Matrix::operator / (Matrix &m) {
	if (m.n_row != m.n_column) {
		cout << "Matrix div: error, matriz no cuadrada\n";
		exit(EXIT_FAILURE);
	}

	Matrix inverse = inv(m);

	if (this->n_column != inverse.n_row) {
		cout << "Matrix div: error, dimensiones incompatibles\n";
		exit(EXIT_FAILURE);
	}

	Matrix& result = (*this) * inverse;

	return result;
}

Matrix& Matrix::operator / (const double val) {
    if (val == 0) {
        cout << "Matrix div: división por cero\n";
        exit(EXIT_FAILURE);
    }
    Matrix* res = new Matrix(n_row, n_column);
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            (*res)(i,j) = (*this)(i,j) / val;
        }
    }
    return *res;
}

Matrix& Matrix::operator = (Matrix m) {
    if (this == &m) return *this;

    n_row = m.n_row;
    n_column = m.n_column;
    data = (double**) malloc(n_row * sizeof(double*));
    for (int i = 0; i < n_row; i++) {
        data[i] = (double*) malloc(n_column * sizeof(double));
        for (int j = 0; j < n_column; j++) {
            data[i][j] = m.data[i][j];
        }
    }
    return *this;
}

Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n, n);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix zeros(const int n_row, const int n_column) {
    Matrix m(n_row, n_column);
    for(int i = 1; i <= n_row; i++) {
        for(int j = 1; j <= n_column; j++) {
            m(i,j) = 0.0;
        }
    }
    return m;
}

double Matrix::norm() {
    double sum = 0.0;
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            sum += (*this)(i,j) * (*this)(i,j);
        }
    }
    return sqrt(sum);
}

double Matrix::dot(Matrix &m) {
    if (n_row != m.n_row || n_column != m.n_column) {
        cout << "Matrix dot: dimensiones incompatibles\n";
        exit(EXIT_FAILURE);
    }
    double res = 0.0;
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            res += (*this)(i,j) * m(i,j);
        }
    }
    return res;
}

Matrix Matrix::extract_vector(const int index, bool row) {
    if (row) return extract_row(index);
    else return extract_column(index);
}

Matrix Matrix::extract_row(const int r) {
    if (r <= 0 || r > n_row) {
        cout << "Matrix extract_row: fila inválida\n";
        exit(EXIT_FAILURE);
    }
    Matrix rowMat(1, n_column);
    for (int j = 1; j <= n_column; j++) rowMat(1,j) = (*this)(r,j);
    return rowMat;
}

Matrix Matrix::extract_column(const int c) {
    if (c <= 0 || c > n_column) {
        cout << "Matrix extract_column: columna inválida\n";
        exit(EXIT_FAILURE);
    }
    Matrix colMat(n_row, 1);
    for (int i = 1; i <= n_row; i++) colMat(i,1) = (*this)(i,c);
    return colMat;
}

void Matrix::assign_row(const int r, Matrix &rowMat) {
    if (r <= 0 || r > n_row || rowMat.n_row != 1 || rowMat.n_column != n_column) {
        cout << "Matrix assign_row: dimensiones incompatibles\n";
        exit(EXIT_FAILURE);
    }
    for (int j = 1; j <= n_column; j++) (*this)(r,j) = rowMat(1,j);
}

void Matrix::assign_column(const int c, Matrix &colMat) {
    if (c <= 0 || c > n_column || colMat.n_column != 1 || colMat.n_row != n_row) {
        cout << "Matrix assign_column: dimensiones incompatibles\n";
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= n_row; i++) (*this)(i,c) = colMat(i,1);
}

Matrix Matrix::union_vector(Matrix &v, bool horizontal) {
    if (horizontal) {
        if (n_row != v.n_row) {
            cout << "Matrix union_vector: filas incompatibles\n";
            exit(EXIT_FAILURE);
        }
        Matrix result(n_row, n_column + v.n_column);
        for (int i = 1; i <= n_row; i++) {
            for (int j = 1; j <= n_column; j++) result(i,j) = (*this)(i,j);
            for (int j = 1; j <= v.n_column; j++) result(i,n_column+j) = v(i,j);
        }
        return result;
    } else {
        if (n_column != v.n_column) {
            cout << "Matrix union_vector: columnas incompatibles\n";
            exit(EXIT_FAILURE);
        }
        Matrix result(n_row + v.n_row, n_column);
        for (int j = 1; j <= n_column; j++) {
            for (int i = 1; i <= n_row; i++) result(i,j) = (*this)(i,j);
            for (int i = 1; i <= v.n_row; i++) result(n_row+i,j) = v(i,j);
        }
        return result;
    }
}

Matrix Matrix::eye(const int n) {
    if (n <= 0) {
        cout << "Matrix eye: tamaño inválido\n";
        exit(EXIT_FAILURE);
    }

    Matrix result(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            result(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }
    return result;
}

/**
 * Función para calcular la matriz traspuesta
 * @param m Matriz de entrada
 * @return Matriz traspuesta
 */
Matrix transpose(Matrix& m) {
    Matrix result(m.n_column, m.n_row);
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

/**
 * Función para calcular la matriz inversa usando eliminación gaussiana
 * @param m Matriz de entrada (debe ser cuadrada)
 * @return Matriz inversa
 */
Matrix inv(Matrix& m) {
    if (m.n_row != m.n_column) {
        cout << "inv: Error - La matriz no es cuadrada\n";
        exit(EXIT_FAILURE);
    }

    int n = m.n_row;
    Matrix augmented(n, 2*n);
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            augmented(i, j) = m(i, j);
            augmented(i, j + n) = (i == j) ? 1.0 : 0.0;
        }
    }
    
    for (int col = 1; col <= n; col++) {
        int max_row = col;
        for (int row = col + 1; row <= n; row++) {
            if (fabs(augmented(row, col)) > fabs(augmented(max_row, col))) {
                max_row = row;
            }
        }
        
        if (max_row != col) {
            for (int j = 1; j <= 2*n; j++) {
                double temp = augmented(col, j);
                augmented(col, j) = augmented(max_row, j);
                augmented(max_row, j) = temp;
            }
        }
        
        if (fabs(augmented(col, col)) < 1e-12) {
            cout << "inv: Error - Matriz singular (no invertible)\n";
            exit(EXIT_FAILURE);
        }
        
        double pivot = augmented(col, col);
        for (int j = 1; j <= 2*n; j++) {
            augmented(col, j) /= pivot;
        }
        
        for (int row = 1; row <= n; row++) {
            if (row != col && augmented(row, col) != 0.0) {
                double factor = augmented(row, col);
                for (int j = 1; j <= 2*n; j++) {
                    augmented(row, j) -= factor * augmented(col, j);
                }
            }
        }
    }
    
    Matrix inverse(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            inverse(i, j) = augmented(i, j + n);
        }
    }
    
    return inverse;
}

Matrix Matrix::cross(Matrix& m) {
    if (n_row != 3 || n_column != 1 || m.n_row != 3 || m.n_column != 1) {
        cout << "Matrix cross: solo válido para vectores 3D\n";
        exit(EXIT_FAILURE);
    }

    Matrix result(3, 1);
    result(1, 1) = (*this)(2, 1) * m(3, 1) - (*this)(3, 1) * m(2, 1);
    result(2, 1) = (*this)(3, 1) * m(1, 1) - (*this)(1, 1) * m(3, 1);
    result(3, 1) = (*this)(1, 1) * m(2, 1) - (*this)(2, 1) * m(1, 1);
    
    return result;
}
