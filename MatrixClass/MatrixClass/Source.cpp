#include "../MyMatrix.h"
#include <iostream>

using namespace std;

int main()
{
  MyMatrix A(2);

  A(0, 0) = 1; A(0, 1) = 2;
  A(1, 0) = 3; A(1, 1) = 4;

  MyMatrix B(2);

  B = A;

  cout <<B.col(0);
}