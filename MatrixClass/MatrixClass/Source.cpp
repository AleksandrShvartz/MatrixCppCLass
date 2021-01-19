#include "../MyMatrix.h"
#include <iostream>

using namespace std;

int main()
{
  MyMatrix A(3,1);

  /*A(0, 0) = 1; A(0, 1) = 2;
  A(1, 0) = 3; A(1, 1) = 4;*/
  cin >> A;

  MyMatrix B(3,1);

  B = A;


  cout <<"cond(A) = "<<A.norm();
}