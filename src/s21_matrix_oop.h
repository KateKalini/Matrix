#ifndef _S21_MATRIX_OOP_H_
#define _S21_MATRIX_OOP_H_

#include <math.h>

#include <iostream>

class S21Matrix {
 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated

 public:
  S21Matrix();                        // Default constructor
  ~S21Matrix();                       // Destructor
  S21Matrix(int rows, int cols);      // parameterized constructor
  S21Matrix(const S21Matrix& other);  // copy cnstructor
  S21Matrix(S21Matrix&& other);       // move

  // ACCESSOR AND MUTATOR:
  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);
  void FillMatrix(double* tmp_matrix, S21Matrix& other);
  void Matrix_pr(S21Matrix& mtrx);
  S21Matrix Minor(int iRow, int jCol);

  // METHODS
  void CreateMatrix(int rows, int cols);
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  // Operators
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);

  bool operator==(const S21Matrix& other) const;
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(const double& num);
};

#endif