#include "MyMatrix.h"
#define KEY_NUMBER_MAT       79834
#define INCORRECT_SIZE       1
#define UNINITIALIZED_MATIRX 2
#define INCORRECT_DIMENSIONS 3
#define NOT_SQUARE_MATRIX    4
#define INCORRECT_B_MATRIX   5
#define IRREVERSIBLE_MATRIX  6

MyMatrix::MyMatrix(int s)
{
  if (s < 1)
    throw INCORRECT_SIZE;

  key = KEY_NUMBER_MAT;

  cols = rows = size = s;

  matrix = (double**)new double[s];
  for (int i = 0; i < s; i++)
  {
    matrix[i] = new double[s];
    for (int j = 0; j < s; j++)
      matrix[i][j] = 0;
  }
}

MyMatrix eye(int n)
{
  MyMatrix M(n);

  for (int i = 0; i < n; i++)
    M(i, i) = 1;

  return M;
}

MyMatrix::MyMatrix(int r, int c)
{
  if (r < 0 || c < 0)
    throw INCORRECT_SIZE;

  rows = r;
  cols = c;
  size = (rows == cols) ? rows : 0;
  key = KEY_NUMBER_MAT;

  matrix = (double**)new double[rows];

  for (int i = 0; i < rows; i++)
  {
    matrix[i] = new double[cols];
    for (int j = 0; j < cols; j++)
      (*this)(i, j) = 0;
  }

}

double& MyMatrix::operator()(int row, int col)
{
  return (double&)matrix[row][col];
}

double MyMatrix::operator()(int row, int col)const
{
  return (double&)matrix[row][col];
}

MyMatrix::MyMatrix(const MyMatrix& M)
{
  if (M.key != KEY_NUMBER_MAT)
    throw UNINITIALIZED_MATIRX;

  key = KEY_NUMBER_MAT;

  rows = M.rows;
  cols = M.cols;
  size = M.size;

  matrix = (double**)new double[rows];
  for (int i = 0; i < rows; i++)
  {
    matrix[i] = new double[cols];
    for (int j = 0; j < cols; j++)
      (*this)(i, j) = M(i, j);
  }
}

MyMatrix MyMatrix::row(int index)const
{
  MyMatrix M(1, cols);

  memcpy(M.matrix[0], matrix[index], cols * sizeof(double));

  return M;
}

MyMatrix MyMatrix::col(int index)const
{
  MyMatrix M(rows, 1);

  for (int i = 0; i < rows; i++)
    M(i, 0) = (*this)(i, index);

  return M;
}

MyMatrix MyMatrix::swapRows(int firstRowIndex, int secondRowIndex)
{
  double* buf = matrix[firstRowIndex];
  matrix[firstRowIndex] = matrix[secondRowIndex];
  matrix[secondRowIndex] = buf;
  return *this;
}

MyMatrix MyMatrix::swapCols(int firstColIndex, int secondColIndex)
{
  double buf;

  for (int i = 0; i < rows; i++)
  {
    buf = (*this)(i, firstColIndex);
    (*this)(i, firstColIndex) = (*this)(i, secondColIndex);
    (*this)(i, secondColIndex) = buf;
  }

  return *this;
}

MyMatrix& MyMatrix::operator=(const MyMatrix& M)
{
  if (M.key != KEY_NUMBER_MAT)
    throw UNINITIALIZED_MATIRX;


  for (int i = 0; i < rows; i++)
    delete matrix[i];

  delete matrix;

  rows = M.rows;
  cols = M.cols;
  size = M.size;

  matrix = (double**)new double[rows];
  for (int i = 0; i < rows; i++)
  {
    matrix[i] = new double[cols];
    memcpy(matrix[i], M.matrix[i], sizeof(double) * cols);
  }

  return *this;
}

MyMatrix::~MyMatrix()
{
  key = 0;
  for (int i = 0; i < rows; i++)
    delete matrix[i];
  delete matrix;
}

MyMatrix MyMatrix::operator*(const double f)const
{
  MyMatrix M2(rows, cols);

  for (int i = 0; i < rows; i++)
    for (int k = 0; k < cols; k++)
      M2(i, k) = (*this)(i, k) * f;

  return M2;
}

MyMatrix operator*(const double f, const MyMatrix& M)
{
  return M * f;
}

MyMatrix MyMatrix::operator*(const MyMatrix& M)const
{
  MyMatrix M3(rows, M.cols);

  if (cols != M.rows)
  {
    throw INCORRECT_DIMENSIONS;
  }

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < M.cols; j++)
      for (int k = 0; k < cols; k++)
        M3(i, j) += (*this)(i, k) * M(k, j);

  return M3;
}

std::ostream& operator<<(std::ostream& out, const MyMatrix& M)
{
  for (int i = 0; i < M.rows; i++)
  {
    for (int j = 0; j < M.cols; j++)
      out << M(i, j) << " ";
    out << "\n";
  }

  return out;
}

std::istream& operator>>(std::istream& in, MyMatrix& M)
{
  for (int i = 0; i < M.rows; i++)
    for (int j = 0; j < M.cols; j++)
      in >> M(i, j);

  return in;
}

void swap(double* p1, double* p2, int size)
{
  double a;
  for (int i = 0; i < size; i++)
  {
    a = p1[i];
    p1[i] = p2[i];
    p2[i] = a;
  }
}

MyMatrix MyMatrix::operator!()const
{
  MyMatrix newMatrix(cols, rows);

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      newMatrix(j, i) = (*this)(i, j);
  return newMatrix;
}

double MyMatrix::norm(int value)const
{
  double s = 0, max = 0;
  if (value == 1)
  {
    for (int j = 0; j < cols; j++)
    {
      for (int i = 0; i < rows; i++)
        s += abs((*this)(i, j));
      if (s > max)
        max = s;
      s = 0;
    }
  }
  else if (value == 2)
  {
    if (cols != 1)
    {
      MyMatrix V(rows, 1);
      double k = 1, k0 = 2;
      MyMatrix M = (*this) * !(*this);
      V(0, 0) = 1;
      while (abs(k - k0) > 1e-12)
      {
        V = M * V;
        k0 = k;
        k = V.norm(3);
        V = V * (1 / k);
        k = sqrt(k);
      }
      return k;
    }
    else
    {
      for (int i = 0; i < rows; i++)
        s += (*this)(i, 0) * (*this)(i, 0);
      return sqrt(s);
    }
  }
  else if (value == 3)
  {
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols; j++)
        s += abs((*this)(i, j));
      if (s > max)
        max = s;
      s = 0;
    }
  }
  return max;
}

double MyMatrix::cond()const
{
  return (*this).inverseGauss().norm(2) * (*this).norm(2);
}

MyMatrix MyMatrix::makeCond(double cond)
{
  if (rows != cols)
    throw NOT_SQUARE_MATRIX;

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
    {
      if (i > j)
        (*this)(i, j) = rand() / (double)RAND_MAX;
      else if (i < j)
        (*this)(i, j) = 0;
      else
      {
        if (i == 0)
          (*this)(i, j) = 1;
        else if (i == 1)
          (*this)(i, j) = cond;
        else
          (*this)(i, j) = 1 + rand() % (int)cond;
      }
    }

  MyMatrix Q(rows, cols), E(rows, cols), W(rows, 1);
  E.toIdentity();

  for (int i = 0; i < rows; i++)
    W(i, 0) = 1 + rand() % 10;

  Q = E - (2 / (W.norm(1) * W.norm(1))) * (W * (!W));

  (*this) = (Q * (*this));
  (*this) = (*this) * !Q;

  return *this;
}

MyMatrix MyMatrix::toIdentity()
{
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      if (i == j)
        (*this)(i, j) = 1;
      else
        (*this)(i, j) = 0;
  return *this;
}

MyMatrix MyMatrix::LUanalysis(const MyMatrix& B)const
{
  if (B.rows != rows || B.cols != 1)
    throw INCORRECT_B_MATRIX;

  if (size < 1)
    throw INCORRECT_SIZE;

  if (rows != cols)
    throw INCORRECT_DIMENSIONS;

  MyMatrix L(size), U(size), Y(size, 1), X(size, 1);

  for (int m = 0; m < size; m++)
  {
    L(m, m) = 1;
    for (int j = m; j < size; j++)
    {
      U(m, j) = (*this)(m, j);

      for (int k = 0; k < m; k++)
        U(m, j) -= L(m, k) * U(k, j);
    }

    for (int i = m + 1; i < size; i++)
    {
      L(i, m) = (*this)(i, m);

      for (int k = 0; k < m; k++)
        L(i, m) -= L(i, k) * U(k, m);

      L(i, m) /= U(m, m);
      L(m, i) = 0;
    }
  }

  for (int i = 0; i < size; i++)
  {
    Y(i, 0) = B(i, 0);
    for (int j = 0; j < i; j++)
      Y(i, 0) -= Y(j, 0) * L(i, j);
  }
  for (int i = size - 1; i >= 0; i--)
  {
    X(i, 0) = Y(i, 0);
    for (int j = i + 1; j < size; j++)
      X(i, 0) -= U(i, j) * X(j, 0);
    X(i, 0) /= U(i, i);
  }
  return X;
}

MyMatrix MyMatrix::operator+(const MyMatrix& M)const
{
  if (cols != M.cols || rows != M.rows)
    throw INCORRECT_DIMENSIONS;

  MyMatrix Result(rows, cols);

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Result(i, j) = (*this)(i, j) + M(i, j);

  return Result;
}

MyMatrix MyMatrix::operator-(const MyMatrix& M)const
{
  if (cols != M.cols || rows != M.rows)
    throw INCORRECT_DIMENSIONS;

  MyMatrix Result(rows, cols);

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Result(i, j) = (*this)(i, j) - M(i, j);

  return Result;
}

MyMatrix MyMatrix::inverseGauss()const
{
  if (rows != cols)
    throw INCORRECT_DIMENSIONS;

  MyMatrix right = eye(size);
  MyMatrix M = *this;

  double maxElem, mult;
  int maxI;

  for (int j = 0; j < cols; j++)
  {
    //выбор ведущего элемента
    maxI = j;
    for (int i = j + 1; i < rows; i++)
      if (abs(M(i, j)) > abs(M(maxI, j)))
        maxI = i;

    if (M(maxI, j) == 0)
      throw IRREVERSIBLE_MATRIX;

    if (maxI != j)
    {
      M.swapRows(maxI, j);
      right.swapRows(maxI, j);
    }

    //прямой и обратный ход

    for (int i = 0; i < rows; i++)
    {
      if (i == j)
        continue;

      mult = M(i, j) / M(j, j);

      for (int k = 0; k < cols; k++)
      {
        if (k >= j)
          M(i, k) -= M(j, k) * mult;
        right(i, k) -= right(j, k) * mult;
      }
    }
    mult = M(j, j);
    for (int i = 0; i < cols; i++)
    {
      if (i >= j)
        M(j, i) /= mult;
      right(j, i) /= mult;
    }

  }
  return right;
}

MyMatrix MyMatrix::makeCustomDet(double det)
{
  if (rows != cols)
    throw 1;
  MyMatrix B(rows);
  double a = 1;
  double g = pow(det, (1. / rows)), gmin = 0.95 * g;
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
    {
      if (i == j)
      {

        (*this)(i, j) = g * (1 + (1. / (95 + (int)rand() % 6)) * pow(-1, rand() % 2));
        a *= (*this)(i, j);
      }
      else
        (*this)(i, j) = 0;
    }

  (*this)(rows - 1, cols - 1) = (*this)(rows - 1, cols - 1) * det / a;

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      if (i == j)
        B(i, j) = 0.001 + rand() / (double)RAND_MAX;
      else
        B(i, j) = 10 + rand() % 11;

  *this = B * (*this) * B.inverseGauss();

  std::cout << (*this).cond();

  return *this;

}