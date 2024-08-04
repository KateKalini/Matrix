#include <gtest/gtest.h>
#include "../s21_matrix_oop.h"

void S21Matrix::FillMatrix(double* tmp_matrix, S21Matrix& other) {
  int n = 0;
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other.matrix_[i][j] = tmp_matrix[n++];
    }
  }
}

void S21Matrix::Matrix_pr(S21Matrix& mtrx) {
    for (int i = 0; i < mtrx.GetRows(); i++) {
    for (int j = 0; j < mtrx.GetCols(); j++) {
      printf("%lf ", mtrx.matrix_[i][j]);
    }
    printf("\n");
  }
}

TEST(CreateMatrix, basic) {
    S21Matrix mtrx;
    ASSERT_EQ(mtrx.GetRows(), 1);
    ASSERT_EQ(mtrx.GetCols(), 1);
}

TEST(CreateMatrix, params) {
    S21Matrix mtrx(3,3);
    ASSERT_EQ(mtrx.GetRows(), 3);
    ASSERT_EQ(mtrx.GetCols(), 3);
}

TEST(EqMatrix, equal) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    mtrx.FillMatrix(arr, mtrx);
    mtrx2.FillMatrix(arr, mtrx2);
    ASSERT_EQ(mtrx.EqMatrix(mtrx2), true);
}

TEST(Copy, cop) {
    S21Matrix mtrx(2,2);
    S21Matrix mtrx2(mtrx);
    ASSERT_EQ(mtrx == mtrx2, true);
}

TEST(SumMatrix, sum) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    mtrx2.FillMatrix(arr, mtrx2);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx.SumMatrix(mtrx2);
    ASSERT_EQ(mtrx.EqMatrix(mtrx_res), true);
}

TEST(SumMatrix, OutOfRangeException) {
    S21Matrix mtrx(3,3);
    S21Matrix other(3,2);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double other_arr[6] = {2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    other.FillMatrix(other_arr, other);

    EXPECT_THROW(mtrx.SumMatrix(other), std::invalid_argument);
}

TEST(SubMatrix, sub) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {0,0,0,0,0,0,0,0,0};
    mtrx.FillMatrix(arr, mtrx);
    mtrx2.FillMatrix(arr, mtrx2);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx.SubMatrix(mtrx2);
    ASSERT_EQ(mtrx.EqMatrix(mtrx_res), true);
}

TEST(SubMatrix, OutOfRangeException) {
    S21Matrix mtrx(3,3);
    S21Matrix other(3,2);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double other_arr[6] = {2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    other.FillMatrix(other_arr, other);

    EXPECT_THROW(mtrx.SubMatrix(other), std::invalid_argument);
}

TEST(MulNumber, number) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double number = 2;
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx.MulNumber(number);
    ASSERT_EQ(mtrx.EqMatrix(mtrx_res), true);
}

TEST(MulMatrix, mul) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    double arr_3[9] = {6,6,6,6,6,6,6,6,6};
    mtrx.FillMatrix(arr, mtrx);
    mtrx2.FillMatrix(arr_2, mtrx2);
    mtrx_res.FillMatrix(arr_3, mtrx_res);
    mtrx.MulMatrix(mtrx2);
    ASSERT_EQ(mtrx.EqMatrix(mtrx_res), true);
}

TEST(MulMatrix, OutOfRangeException) {
    S21Matrix mtrx(3,3);
    S21Matrix other(2,2);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double other_arr[6] = {2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    other.FillMatrix(other_arr, other);
    EXPECT_THROW(mtrx.MulMatrix(other), std::invalid_argument);
}

TEST(Transpose, tr) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,2,3,1,2,3,1,2,3};
    double arr_2[9] = {1,1,1,2,2,2,3,3,3};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx2 = mtrx.Transpose();
    ASSERT_EQ(mtrx2.EqMatrix(mtrx_res), true);
}

TEST(Minor, min) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(2,2);
    S21Matrix mtrx_res(2,2);
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    double arr_res[4] = {5,6,8,9};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_res, mtrx_res);
    mtrx2 = mtrx.Minor(0,0);
    ASSERT_EQ(mtrx2.EqMatrix(mtrx_res), true);
}

TEST(Minor, min2) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(2,2);
    S21Matrix mtrx_res(2,2);
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    double arr_res[4] = {1,3,7,9};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_res, mtrx_res);
    mtrx2 = mtrx.Minor(1,1);
    ASSERT_EQ(mtrx2.EqMatrix(mtrx_res), true);
}

TEST(Determinant, det1) {
    S21Matrix mtrx(3,3);
    double res = 121, res_mtrx;
    double arr[9] = {1,4,2,5,3,7,6,2,1};
    mtrx.FillMatrix(arr, mtrx);
    res_mtrx = mtrx.Determinant();
    ASSERT_EQ(res_mtrx, res);
}

TEST(Determinant, det2) {
    S21Matrix mtrx(3,3);
    double res = 0, res_mtrx;
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    mtrx.FillMatrix(arr, mtrx);
    res_mtrx = mtrx.Determinant();
    ASSERT_EQ(res_mtrx, res);
}

TEST(Determinant, det3) {
    S21Matrix mtrx(1,1);
    double res = 1, res_mtrx;
    double arr[4] = {1,2,3,4};
    mtrx.FillMatrix(arr, mtrx);
    res_mtrx = mtrx.Determinant();
    ASSERT_EQ(res_mtrx, res);
}

TEST(Determinant, OutOfRangeException) {
    S21Matrix mtrx(3,2);
    double arr[6] = {1,1,1,1,1,1};
    mtrx.FillMatrix(arr, mtrx);
    EXPECT_THROW(mtrx.Determinant(), std::invalid_argument);
}

TEST(CalcComplements, calc) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    double arr_2[9] = {-3,6,-3,6,-12,6,-3,6,-3};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx2 = mtrx.CalcComplements();
    ASSERT_EQ(mtrx2.EqMatrix(mtrx_res), true);
}

TEST(CalcComplements, OutOfRangeException) {
    S21Matrix mtrx(3,2);
    double arr[6] = {1,1,1,1,1,1};
    mtrx.FillMatrix(arr, mtrx);
    EXPECT_THROW(mtrx.CalcComplements(), std::invalid_argument);
}

TEST(InverseMatrix, inv) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {2,5,7,6,3,4,5,-2,-3};
    double arr_2[9] = {1,-1,1,-38,41,-34,27,-29,24};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx2 = mtrx.InverseMatrix();
    ASSERT_EQ(mtrx2.EqMatrix(mtrx_res), true);
}

TEST(InverseMatrix, OutOfRangeException) {
    S21Matrix mtrx(1,1);
    double arr[4] = {0,0,0,0};
    mtrx.FillMatrix(arr, mtrx);
    EXPECT_THROW(mtrx.InverseMatrix(), std::invalid_argument);
}

TEST(SetRows, srow) {
    S21Matrix mtrx(3,3);
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    mtrx.FillMatrix(arr, mtrx);
    S21Matrix mtrx2(5,3);
    double arr2[15] = {1,2,3,4,5,6,7,8,9,0,0,0,0,0,0};
    mtrx2.FillMatrix(arr2, mtrx2);
    mtrx.SetRows(5);
    ASSERT_EQ(mtrx2.EqMatrix(mtrx), true);
}

TEST(SetRows, smaller) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx_res(2,3);
    double arr[9] = {1,2,3,4,5,6,7,8,9};
    double arr2[6] = {1,2,3,4,5,6};
    mtrx.FillMatrix(arr,mtrx);
    mtrx_res.FillMatrix(arr2,mtrx_res);
    mtrx.SetRows(2);
    EXPECT_EQ(mtrx == mtrx_res , true);
}

TEST(SetCols, smaller) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx_res(3,2);
    double arr[9] = {9,8,7,6,5,4,3,2,1};
    double arr2[6] = {9,8,6,5,3,2};
    mtrx.FillMatrix(arr,mtrx);
    mtrx_res.FillMatrix(arr2,mtrx_res);
    mtrx.SetCols(2);
    EXPECT_EQ(mtrx == mtrx_res , true);
}

TEST(SetCols, bigger) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx_res(3,4);
    double arr[9] = {9,8,7,6,5,4,3,2,1};
    double arr2[12] = {9,8,7,0,6,5,4,0,3,2,1,0};
    mtrx.FillMatrix(arr,mtrx);
    mtrx_res.FillMatrix(arr2,mtrx_res);
    mtrx.SetCols(4);
    EXPECT_EQ(mtrx == mtrx_res , true);
}

TEST(Operator_plus, plus) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx2 = mtrx + mtrx;
    ASSERT_EQ(mtrx2 == mtrx_res, true);
}

TEST(Operator_plus, OperatorPr) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx += mtrx;
    ASSERT_EQ(mtrx == mtrx_res, true);
}

TEST(Operator_minus, min) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {0,0,0,0,0,0,0,0,0};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx2 = mtrx - mtrx;
    ASSERT_EQ(mtrx2 == mtrx_res, true);
}

// TEST(Operator_minus, OperatorPr2) {
//     S21Matrix mtrx(3,3);
//     S21Matrix mtrx_res(3,3);
//     double arr[9] = {1,1,1,1,1,1,1,1,1};
//     double arr_2[9] = {0,0,0,0,0,0,0,0,0};
//     mtrx.FillMatrix(arr, mtrx);
//     mtrx_res.FillMatrix(arr_2, mtrx_res);
//     mtrx -= mtrx;
//     ASSERT_EQ(mtrx == mtrx_res, true);
// }

TEST(Operator_mul, mul2) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx3(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    double arr_3[9] = {6,6,6,6,6,6,6,6,6};
    mtrx.FillMatrix(arr, mtrx);
    mtrx2.FillMatrix(arr_2, mtrx2);
    mtrx_res.FillMatrix(arr_3, mtrx_res);
    mtrx3 = mtrx * mtrx2;
    ASSERT_EQ(mtrx3 == mtrx_res, true);
}

TEST(Operator_mul, OperatorPr2) {
    S21Matrix mtrx(2,2);
    S21Matrix mtrx_res(2,2);
    double arr[4] = {1,1,1,1};
    double arr_2[4] = {2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx *= mtrx;
    ASSERT_EQ(mtrx == mtrx_res, true);
}

TEST(Operator_mul, number) {
    S21Matrix mtrx(3,3);
    S21Matrix mtrx2(3,3);
    S21Matrix mtrx_res(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double arr_2[9] = {2,2,2,2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    mtrx_res.FillMatrix(arr_2, mtrx_res);
    mtrx2 = mtrx * 2;
    ASSERT_EQ(mtrx2 == mtrx_res, true);
}

TEST(Operator_mul, OperatorPr) {
    S21Matrix mtrx(3,3);
    double number = 2.0;
    S21Matrix check_mtrx(3,3);
    double arr[9] = {1,1,1,1,1,1,1,1,1};
    double check_arr[9] = {2,2,2,2,2,2,2,2,2};
    mtrx.FillMatrix(arr, mtrx);
    check_mtrx.FillMatrix(check_arr,check_mtrx);
    mtrx *= number;
    EXPECT_EQ(mtrx.EqMatrix(check_mtrx), true);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS(); 
}