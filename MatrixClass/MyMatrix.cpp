#include "MyMatrix.h"
#define KEY_NUMBER_MAT       79834
#define INCORRECT_SIZE       1
#define UNINITIALIZED_MATIRX 2
#define INCORRECT_DIMENSIONS 3
#define NOT_SQUARE_MATRIX    4
#define INCORRECT_B_MATRIX   5

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

MyMatrix MyMatrix::row(int index)
{
  MyMatrix M(1, cols);

  memcpy(M.matrix[0], matrix[index], cols*sizeof(double));

  return M;
}

MyMatrix MyMatrix::col(int index)
{
  MyMatrix M(rows, 1);

  for (int i = 0; i < rows; i++)
    M(i, 0) = (*this)(i, index);

  return M;
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

MyMatrix MyMatrix::preMinor(int row, int col) const
{
  MyMatrix newMatrix(size - 1);

  for (int i = 0, in = 0; i < size; i++)
    if (i != row)
    {
      for (int j = 0, jn = 0; j < size; j++)
        if (j != col)
          newMatrix(in, jn++) = (*this)(i, j);
      in++;
    }

  return newMatrix;
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

double MyMatrix::determinant() const
{
  const double EPS = pow(0.1, 9);
  int n = rows;

  MyMatrix a(*this);
  double det = 1;
  for (int i = 0; i < n; ++i)
  {
    int k = i;
    for (int j = i + 1; j < n; ++j)
      if (abs(a(j, i)) > abs(a(k, i)))
        k = j;
    if (abs(a(k, i)) < EPS)
    {
      det = 0;
      break;
    }
    swap(a.matrix[i], a.matrix[k], rows);
    if (i != k)
      det = -det;
    det *= a(i, i);
    for (int j = i + 1; j < n; ++j)
      a(i, j) = a(i, j) / a(i, i);
    for (int j = 0; j < n; ++j)
      if (j != i && abs(a(j, i)) > EPS)
        for (int k = i + 1; k < n; ++k)
          a(j, k) -= a(i, k) * a(j, k);
  }
  return det;
}//меняет матрицу

MyMatrix MyMatrix::allied()const
{
  MyMatrix newMatrix(size);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      newMatrix(i, j) = this->preMinor(i, j).determinant() * ((i + j) % 2 ? -1 : 1);

  return newMatrix;
}

//по определению
MyMatrix MyMatrix::inverse()const
{
  double det = (*this).determinant();

  if (abs(det) < 0.000000001)
  {
    throw 10;
  }
  static MyMatrix newMatrix(size);

  newMatrix = *this;
  return !((1 / det) * newMatrix.allied());
}

MyMatrix MyMatrix::operator!()const
{
  MyMatrix newMatrix(cols, rows);

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      newMatrix(j, i) = (*this)(i, j);
  return newMatrix;
}

double MyMatrix::norm()const
{
  double s = 0, max = 0;

  for (size_t j = 0; j < cols; ++j) {
    for (size_t i = 0; i < rows; ++i) {
      s += ((*this)(i, j)) * ((*this)(i, j));
    }
    if (max < s)
      max = s;
    s = 0;
  }
  return sqrt(max);
}

double MyMatrix::cond()const
{
  return (*this).inverse().norm() * (*this).norm();
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

  Q = E - (2 / (W.norm() * W.norm())) * (W * (!W));

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

MyMatrix MyMatrix::inverseGauss()
{
  if (rows != cols)
    throw INCORRECT_DIMENSIONS;

  MyMatrix Right(rows, cols), Left(*this);
  double d = 0, r = 0, tmp = 0;
  Right.toIdentity();//единичная

  for (int i = 0; i < rows; i++)
  {
    if (Left(i, i) == 0)
    {
      for (int j = i + 1; j < rows; j++)//меняем столбцы местами если на диагонали 0
        if (Left(j, i) != 0)
        {
          for (int k = i; k < cols; k++)
          {
            tmp = Left(j, k);
            Left(j, k) = Left(i, k);
            Left(i, k) = tmp;

            tmp = Right(j, k);
            Right(j, k) = Right(i, k);
            Right(i, k) = tmp;
          }
          break;
        }
    }

    r = Left(i, i);
    for (int l = 0; l < cols; l++)//делим всю строку на диагональный элемент
    {
      Left(i, l) = Left(i, l) / r;
      Right(i, l) = Right(i, l) / r;
    }

    for (int j = 0; j < cols; j++)
      if (j != i)
      {
        d = Left(j, i);
        for (int k = 0; k < rows; k++)
        {
          Left(j, k) = Left(j, k) - d * Left(i, k);
          Right(j, k) = Right(j, k) - d * Right(i, k);
        }
      }
  }

  return Right;
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