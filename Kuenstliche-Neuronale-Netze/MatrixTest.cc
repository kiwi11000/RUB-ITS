#include "Matrix.h"
#include <iostream>

// main()-Funktion, no parameter
int main(int argc, char* argv[]){
  // initialize random number generator
  initializeRNG();

  // generate a Matrix with 3 rows and 4 columns
  Matrix M(3, 4);
  // fill it randomly
  M.fillRandom(-5., 5.);

  // copy the matrix and its contents
  Matrix M_inv(M);
  // invert Matrix
  M_inv.invert();

  // calculate the product M*M^{-1} (=E)
  Matrix I1 = M * M_inv;
  // calculate M^{-1}*M*M^{-1} ( should be equal to M^{-1})
  Matrix M_inv_test = M_inv * M * M_inv;

  // generate a Matrix with 5x5 entries
  Matrix M2(5, 5);
  // fill the Matrix randomly with values between -5 and 5
  M2.fillRandom(-5., 5.);
  // generate M2^{-1} as above
  Matrix M2_inv(M2);
  M2_inv.invert();
  // multiply M2 * M2^{-1}
  Matrix I2 = M2 * M2_inv;
  Matrix I3 = M2_inv * M2;

  // write all results to std::cout
  std::cout << "\nrandom matrix:\n" << M << "\ninverted:\n" << M_inv << "\nunity matrices: \n" << I1 << "\ntest of M^{-1}:\n" << M_inv_test << "\nrandom matrix:\n" << M2 << "\ninverted:\n" << M2_inv << "\nunity matrix:\n"  << I2 << "\nand:\n" << I3 << "\n\n";

  return 0;
}
