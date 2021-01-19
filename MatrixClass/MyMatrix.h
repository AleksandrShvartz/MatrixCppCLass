#pragma once
#include <iostream>
class MyMatrix
{
public:
  MyMatrix(int s);
  MyMatrix(const MyMatrix&);
  ~MyMatrix();
  MyMatrix(int, int);
  MyMatrix& operator= (const MyMatrix&);
  double& operator()(const int row, const int col);
  double operator()(const int row, const int col)const;
  MyMatrix operator*(const double f)const;
  friend MyMatrix operator* (const double f, const MyMatrix& M);
  MyMatrix operator* (const MyMatrix& M)const;
  friend std::ostream& operator<<(std::ostream& out, const MyMatrix& M);
  friend std::istream& operator>>(std::istream& in, MyMatrix& M);
  MyMatrix preMinor(const int row, const int col) const;
  double determinant() const;
  MyMatrix allied()const;
  MyMatrix inverse()const;
  MyMatrix operator!()const;
  double norm(int value = 2)const;
  double cond()const;
  MyMatrix makeCond(double cond);
  MyMatrix toIdentity();
  MyMatrix LUanalysis(const MyMatrix& B)const;
  MyMatrix operator+(const MyMatrix& M)const;
  MyMatrix operator-(const MyMatrix& M)const;
  MyMatrix inverseGauss()const;
  MyMatrix makeCustomDet(double det);
  MyMatrix row(int index)const;
  MyMatrix col(int index)const;
  MyMatrix swapRows(int firstRowIndex, int secondRowIndex);
  MyMatrix swapCols(int firstColIndex, int secondColIndex);
 

private:

  double** matrix;
  int key;
  int size, rows, cols;
};

MyMatrix eye(int n);
