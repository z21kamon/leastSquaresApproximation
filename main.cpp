#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

class Matrix {
public:
    int rows{}, columns{};
    vector<vector<double>> matrix;

    Matrix(int rows, int columns) : rows(rows), columns(columns) {
        for (int i = 0; i < rows; ++i) {
            vector<double> row;
            for (int j = 0; j < columns; ++j) row.push_back(0);
            matrix.push_back(row);
        }
    }

    Matrix() : rows(0), columns(0) {}

    friend istream& operator>> (istream& in, Matrix& m);
    friend ostream& operator<< (ostream& out, Matrix& m);

    Matrix operator+ (const Matrix& other) {
        if (this->rows == other.rows && this->columns == other.columns) {
            Matrix result(rows, columns);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    result.matrix[i][j] = this->matrix[i][j] + other.matrix[i][j];
                }
            }
            return result;
        }
        else cout << "Error: the dimensional problem occurred\n";
        return {0, 0};
    }

    Matrix operator- (const Matrix& other) {
        if (this->rows == other.rows && this->columns == other.columns) {
            Matrix result(rows, columns);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    result.matrix[i][j] = this->matrix[i][j] - other.matrix[i][j];
                }
            }
            return result;
        }
        else cout << "Error: the dimensional problem occurred\n";
        return {0, 0};
    }

    Matrix operator* (const Matrix& other) {
        if (this->columns == other.rows) {
            Matrix result(this->rows, other.columns);
            for (int i = 0; i < this->rows; ++i) {
                for (int j = 0; j < other.columns; ++j) {
                    result.matrix[i][j] = 0;
                    for (int r = 0; r < other.rows; ++r) result.matrix[i][j] += this->matrix[i][r] * other.matrix[r][j];
                }
            }
            return result;
        }
        else cout << "Error: the dimensional problem occurred\n";
        return {0, 0};
    }

    Matrix operator= (const Matrix& other) {
        for (int i = 0; i < other.rows; ++i) {
            for (int j = 0; j < other.columns; j++) {
                matrix[i][j] = other.matrix[i][j];
            }
        }
        return *this;
    }

    vector<double>* operator[](int index) {
        return &(this->matrix[index]);
    }

    Matrix transpose() {
        Matrix result(this->columns, this->rows);
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) result.matrix[i][j] = this->matrix[j][i];
        }
        return result;
    }
};

istream &operator>>(istream &in, Matrix &m) {
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.columns; ++j) {
            double value;
            in >> value;
            m.matrix[i][j] = value;
        }
    }
    return in;
}

ostream &operator<<(ostream &out, Matrix &m) {
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.columns; ++j) {
            out << m.matrix[i][j];
            if (j < m.columns - 1) out << " ";
        }
        cout << "\n";
    }
    return out;
}


bool operator==(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.rows != rhs.rows && lhs.columns != rhs.columns) return false;
    for (int i = 0; i < lhs.rows; ++i) {
        for (int j = 0; j < lhs.columns; ++j) {
            if (lhs.matrix[i][j] != rhs.matrix[i][j]) return false;
        }
    }
    return true;
}

Matrix inverse(Matrix m) {
    int dimension = m.rows;
    Matrix augm = Matrix(dimension, 2 * dimension);
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            augm.matrix[i][j] = m.matrix[i][j];
            if (i == j) augm.matrix[i][j + dimension] = 1;
        }
    }

    for (int i = dimension - 1; i > 0; --i) {
        if (augm.matrix[i - 1][0] < augm.matrix[i][0]) {
            vector<double> temp = augm.matrix[i];
            augm.matrix[i] = augm.matrix[i - 1];
            augm.matrix[i - 1] = temp;
        }
    }

    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            if (j != i) {
                double temp = augm.matrix[j][i] / augm.matrix[i][i];
                for (int k = 0; k < 2 * dimension; ++k) {
                    augm.matrix[j][k] -= augm.matrix[i][k] * temp;
                }
            }
        }
    }

    for (int i = 0; i < dimension; ++i) {
        double temp = augm.matrix[i][i];
        for (int j = 0; j < 2 * dimension; j++) {
            augm.matrix[i][j] = augm.matrix[i][j] / temp;
        }
    }

    Matrix res(dimension,dimension);
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            res.matrix[i][j] = augm.matrix[i][dimension + j];
        }
    }
    return res;
}


int main() {
    cout << setprecision(4) << fixed;
    int m, n;
    vector<double> t_set;
    cin >> m;
    Matrix b_vector = Matrix(m, 1);
    for (int i = 0; i < m; ++i) {
        double t, b;
        cin >> t >> b;
        t_set.push_back(t);
        b_vector.matrix[i][0] = b;
    }
    cin >> n;

    Matrix a = Matrix(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            a.matrix[i][j] = pow(t_set[i], j);
        }
    }

    cout << "A:\n" << a;
    Matrix a_transposed = a.transpose();
    Matrix a_t_a = a_transposed * a;
    cout << "A_T*A:\n" << a_t_a;

    Matrix a_t_a_inverse = inverse(a_t_a);
    cout << "(A_T*A)^-1:\n" << a_t_a_inverse;

    Matrix a_t_b = a_transposed * b_vector;
    cout << "A_T*b:\n" << a_t_b;

    Matrix result = a_t_a_inverse * a_t_b;
    cout << "x~:\n" << result;
    return 0;
}
