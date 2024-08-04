#include "s21_matrix_oop.h"

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  S21Matrix res(rows, cols_);
  if (rows > rows_) {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        if (i >= rows_) {
          res.matrix_[i][j] = 0;
        } else {
          res.matrix_[i][j] = matrix_[i][j];
        }
      }
    }
  } else {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        res.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  *this = res;
}

void S21Matrix::SetCols(int cols) {
  S21Matrix res(rows_, cols);
  if (cols > cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols; j++) {
        if (j >= cols_) {
          res.matrix_[i][j] = 0;
        } else {
          res.matrix_[i][j] = matrix_[i][j];
        }
      }
    }
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols; j++) {
        res.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  *this = res;
}

S21Matrix::S21Matrix() : rows_(1), cols_(1) { CreateMatrix(rows_, cols_); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  CreateMatrix(rows_, cols_);
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  matrix_ = other.matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;
  other.matrix_ = nullptr;
  other.cols_ = 0;
  other.rows_ = 0;
}

void S21Matrix::CreateMatrix(int rows, int cols) {
  matrix_ = new double*[rows];
  for (int i = 0; i < rows; i++) {
    matrix_[i] = new double[cols]();
  }
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    for (int i = 0; i < rows_; ++i) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(other.matrix_[i][j] - matrix_[i][j]) > 1e-6) return false;
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::invalid_argument("Mistake - different size of matrixs");
  } else {
    for (int i = 0; i != rows_; i++) {
      for (int j = 0; j != cols_; j++) {
        matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::invalid_argument("Mistake - different size of matrixs");
  } else {
    for (int i = 0; i != rows_; i++) {
      for (int j = 0; j != cols_; j++) {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (other.rows_ != cols_) {
    throw std::invalid_argument("Mistake - different size of matrixs");
  }
  S21Matrix res(rows_, other.cols_);
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = res;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(cols_, rows_);
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      res.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Mistake - not square matrix");
  }
  double result = 0;
  if (rows_ == 1) {
    result = matrix_[0][0];
  } else if (rows_ == 2) {
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int j = 0; j < cols_; j++) {
      result += matrix_[0][j] * Minor(0, j).Determinant() * pow(-1, j);
    }
  }
  return result;
}

S21Matrix S21Matrix::Minor(int iRow, int jCol) {
  int rows = rows_ - 1, cols = cols_ - 1;
  if (rows_ == 1) rows = 1;
  if (cols_ == 1) cols = 1;
  S21Matrix minor(rows, cols);

  int row_no_skip = 0, column_no_skip;
  for (int r = 0; r < rows_; r++) {
    if (iRow == r && rows_ != 1) continue;
    column_no_skip = 0;
    for (int c = 0; c < cols_; c++) {
      if (jCol == c && cols_ != 1) continue;
      minor.matrix_[row_no_skip][column_no_skip] = matrix_[r][c];
      if (cols_ != 1) column_no_skip++;
    }
    if (rows_ != 1) row_no_skip++;
  }
  return minor;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Mistake - not square matrix");
  }
  S21Matrix calc(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      calc.matrix_[i][j] = Minor(i, j).Determinant() * pow(-1, i + j);
    }
  }
  return calc;
}

S21Matrix S21Matrix::InverseMatrix() {
  const double det = Determinant();
  if (det == 0) {
    throw std::invalid_argument("Error: determinant is 0");
  }
  S21Matrix inverse(cols_, rows_);
  inverse = CalcComplements().Transpose();
  inverse.MulNumber(1 / det);
  return inverse;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  cols_ = other.cols_;
  rows_ = other.rows_;
  matrix_ = new double*[other.rows_];
  for (int i = 0; i < other.rows_; i++) {
    matrix_[i] = new double[other.cols_];
  }
  for (int i = 0; i < other.rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res = *this;
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res = *this;
  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res = *this;
  res.MulNumber(num);
  return res;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

void S21Matrix::operator+=(const S21Matrix& other) { SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix& other) { SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix& other) { MulMatrix(other); }

void S21Matrix::operator*=(const double& num) { MulNumber(num); }