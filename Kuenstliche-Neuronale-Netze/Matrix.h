#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>

// matrix class for rectangular matrices
class Matrix{
public:
  // creates Matrix with given number of rows (xSizeA) and columns (ySizeA)
  Matrix(unsigned xSizeA = 0, unsigned ySizeA = 0);
  // copy constructor; copies size and matrix entries
  Matrix(const Matrix& srcA);
  // default destructor; frees the memory
  ~Matrix();

  // assignement operator, copies a Matrix
  Matrix& operator = (const Matrix& srcA);

  // initializes all entries of this Matrix with valueA
  void initialize(double valueA);

  // returnes the number of rows of the matrix
  unsigned xSize() const;

  // returnes the number of columns of the matrix
  unsigned ySize() const;

  // returnes the pointer to the row xA;
  double* operator [](unsigned xA);

  // returnes the read-only pointer to the row xA; 
  const double* const operator[](unsigned xA) const;

  // Matrix addition: A = B + C
  Matrix operator +(const Matrix& mA) const;

  // Matrix multiplication: A = B * C
  Matrix operator *(const Matrix& mA) const;

  // inverts this matrix using singular value decomposition
  // all eigenvalues < epsilonA are set to 0
  void invert(double epsilonA = 1e-12);

  // fills the matrix with quasi-random values between lowerBoundA and upperBoundA
  void fillRandom(double lowerBoundA, double upperBoundA);

protected:
  // number of rows and columns
  unsigned xSizeE, ySizeE;

  // data array
  double* dataE;

private:
  // private function, that computes the SVD
  void computeSVD(Matrix& u, Matrix& v, double* w);

  // updates the element of the matrix at position [xA][yA]
  void element(unsigned xA, unsigned yA, double valueA);

  // returnes the element of the matrix at position [xA][yA]
  const double element(unsigned xA, unsigned yA) const;
};

// generates quasi-random real value between lowerBoundA and upperBoundA
double randomFromInterval(double lowerBoundA, double upperBoundA);

// initializes the random number generator with the given value or the current time (if no value is specified)
void initializeRNG(unsigned valueA = 0);

// overload ostream operator << for class Matrix
std::ostream& operator <<(std::ostream& streamA, const Matrix& matrixA);




// --------------------------------------
//          Implementations
// --------------------------------------
#include <math.h>
#include <algorithm>
#include <functional>

#include <time.h>

#ifdef WIN32
#define _CRT_SECURE_NO_WARNING 1
#endif

// global random function; generates quasi-random values
inline double randomFromInterval(double lowerBoundA, double upperBoundA){
  assert(lowerBoundA <= upperBoundA);
  // get a random value
  double valueL = rand();
  // generate random value between 0 and (upperBoundA - lowerBoundA)
  valueL *= (upperBoundA - lowerBoundA) / RAND_MAX;

  //return shifted random value
  return valueL + lowerBoundA;
}

// initializes the random function with the current time; only possible on Linux Systems
inline void initializeRNG(unsigned valueA){
  // check if the value is greater than 0
  if (valueA){
    // take the given value to initialize the RNG
    srand(valueA);
  }else{
    // take the current time to initialize the RNG
    srand(time(0));
  }
}


// construct a new matrix with the given size
inline Matrix::Matrix(unsigned xSizeA, unsigned ySizeA)
:xSizeE(xSizeA), ySizeE(ySizeA){
  // generate the new data array
  if (xSizeE*ySizeE > 0)
    dataE = new double[xSizeE * ySizeE];
  else
    dataE = 0;
}

// copy constructor
inline Matrix::Matrix(const Matrix& srcA)
:xSizeE(srcA.xSizeE), ySizeE(srcA.ySizeE){
  //generate new data array
  if (xSizeE*ySizeE > 0)
    dataE = new double[xSizeE * ySizeE];
  else
    dataE = 0;
  // copy data from the source Matrix
  for (unsigned i = 0; i < xSizeE * ySizeE; ++i){
    dataE[i] = srcA.dataE[i];
  }
}

// destructor
inline Matrix::~Matrix(){
  // delete data
  if (dataE)
    delete[] dataE;
}

// initializes all entries to the same value
inline void Matrix::initialize(double valueA){
  // init the data array with the given value
  std::fill(dataE, dataE + xSizeE * ySizeE, valueA);
}

// return the number of rows
inline unsigned Matrix::xSize() const{
  return xSizeE;
}

// return the number of columns
inline unsigned Matrix::ySize() const{
  return ySizeE;
}

// assignment operator
inline Matrix& Matrix::operator =(const Matrix& srcA){
  // check if we assign me to me.
  if (&srcA == this)
    return *this;

  // assign new size
  xSizeE = srcA.xSizeE;
  ySizeE = srcA.ySizeE;

  // delete old data array
  if (dataE)
    delete[] dataE;
  // generate new data array
  if (xSizeE*ySizeE > 0)
    dataE = new double[xSizeE * ySizeE];
  else
    dataE = 0;
  // copy data from the source Matrix
  std::copy(srcA.dataE, srcA.dataE + xSizeE * ySizeE, dataE);

  // return *this as usual in assignment operators
  return *this;
}


// return matrix data at [xA][yA]
inline const double Matrix::element(unsigned xA, unsigned yA) const{
  // check, whether row and column is ok
  assert(xA < xSizeE && yA < ySizeE);
  // return the desired element
  return dataE[xA * ySizeE + yA];
}

// set  matrix data at [xA][yA] to valueA
inline void Matrix::element(unsigned xA, unsigned yA, double valueA) {
  // check, whether row and column is ok
  assert(xA < xSizeE && yA < ySizeE);
  // set the element
  dataE[xA * ySizeE + yA] = valueA;
}


// overloaded version of ostream::operator << for Matrix class
inline std::ostream& operator <<(std::ostream& streamA, const Matrix& matrixA){
  // write all rows
  for (unsigned xL = 0; xL < matrixA.xSize(); xL++){
    // write all entries of the row
    for (unsigned yL = 0; yL < matrixA.ySize(); ++yL){
      streamA << matrixA[xL][yL]<< '\t';
    }
    streamA << '\n';
  }
  // return the new position in the stream, as usual for overloads of ostream::operator <<
  return streamA;
}

// overload of operator + for two matrices
inline Matrix Matrix::operator+(const Matrix& mA) const{
  assert(mA.xSizeE == xSizeE && mA.ySizeE == ySizeE);
  // generate matrix with needed size
  Matrix resL(xSizeE, ySizeE);

  // do the elementwise addition
  std::transform(dataE, dataE + xSizeE * ySizeE, mA.dataE, resL.dataE, std::plus<double>());

  // return the summed matrices
  return resL;
}

// overload of operator * for matrices; implements matrix-multiplication; doesn't change *this Matrix
inline Matrix Matrix::operator*(const Matrix& mA) const{
  assert(mA.xSizeE == ySizeE);
  // generate Matrix with the needed size
  Matrix resL(xSizeE, mA.ySizeE);
  // fill the result matrix element by element
  for (unsigned rowL = 0; rowL < resL.xSizeE; ++rowL){
    for (unsigned colL = 0; colL < resL.ySizeE; ++colL){
      // init counter
      double resValueL = 0.;
      // iterate over the current row/column of the source matrices
      for (unsigned iL = 0; iL < ySizeE; ++iL){
        // add up the multiplied source matrices
        resValueL += element(rowL, iL) * mA.element(iL, colL);
      }
      // set the value
      resL.element(rowL, colL, resValueL);
    }
  }
  // return the resulting Matrix
  return resL;
}

// fills the matrix randomly with uniformly distributed values between lowerBoundA and upperBoundA
inline void Matrix::fillRandom(double lowerBoundA, double upperBoundA){
  // iterate over all elements
  for (unsigned iL = 0; iL < xSizeE * ySizeE; ++iL){
    // fill in a random element
    dataE[iL] = randomFromInterval(lowerBoundA, upperBoundA);
  }
}

// invert this matrix using singular value decomposition
inline void Matrix::invert(double epsilonA){
  // copy this matrix to U
  Matrix U(*this);
  Matrix V(ySizeE, ySizeE);
  double* W = new double[ySizeE];
  // fill U, V and W with the singular value decomposition matrices
  computeSVD(U, V, W);

  // delete all singular values (i.e. the eigenvalues less than epsilonA)
  for (unsigned i = 0; i < ySizeE; ++i){
    if (W[i] < epsilonA){
      W[i] = 0.;
    }else{
      W[i] = 1./ W[i];
    }
  }

  // fill this matrix with the new values
  unsigned t = xSizeE;
  xSizeE = ySizeE;
  ySizeE = t;
  // no need to get new memory, because the size (xSizeE * ySizeE) haven't changed

  // add up A = V^T W U
  for (unsigned i = 0; i < xSizeE; ++i){
    for (unsigned j = 0; j < ySizeE; ++j){
      double sumL = 0.;
      for (unsigned k = 0; k < xSizeE; ++k){
        sumL += V[i][k] * W[k] * U[j][k];
      }
      element(i, j, sumL);
    }
  }

  delete W;
  // done! :)
}



//////////////////////////////////////////////////////
// SVD stuff copied from Shark
// no more good documentation from here on
static double SIGN (double a, double b){
    return b > 0 ? fabs (a) : -fabs (a);
}

// operator [] for usual [row][col]-access to matrix data, inaccassable from outside
inline double* Matrix::operator [](unsigned xA){
  assert(xA < xSizeE);
  return &dataE[xA * ySizeE];
}

// operator [] for usual [row][col]-read-only-access to matrix data
inline const double* const Matrix::operator[](unsigned xA) const{
  assert(xA < xSizeE);
  return &dataE[xA * ySizeE];
}


// computes the singular value decomposition
inline void Matrix::computeSVD(Matrix& u, Matrix& v, double* w){

  double* rv1 = new double[ySizeE];
  int m = xSizeE, n = ySizeE;

  int flag;
  int i, its, j, jj, k, l, nm(0);
  double anorm, c, f, g, h, p, s, scale, x, y, z;

  // householder reduction to bidiagonal form
  g = scale = anorm = 0.0;

  for (i = 0; i < n; i++){
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m){
      for (k = i; k < m; k++){
        scale += fabs (u[k][i]);
      }

      if (scale != 0.0){
        for (k = i; k < m; k++){
          u[k][i] /= scale;
          s += u[k][i] * u[k][i];
        }

        f = u[i][i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        u[i][i] = f - g;

        for (j = l; j < n; j++){
          s = 0.0;
          for (k = i; k < m; k++){
            s += u[k][i] * u[k][j];
          }

          f = s / h;
          for (k = i; k < m; k++){
            u[k][j] += f * u[k][i];
          }
        }

        for (k = i; k < m; k++){
          u[k][i] *= scale;
        }
      }
    }

    w[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m && i != n-1){
      for (k = l; k < n; k++){
        scale += fabs(u[i][k]);
      }

      if (scale != 0.0){
        for (k = l; k < n; k++){
          u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }

        f = u[i][l];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        u[i][l] = f - g;

        for (k = l; k < n; k++){
          rv1[k] = u[i][k] / h;
        }

        for (j = l; j < m; j++){
          s = 0.0;
          for (k = l; k < n; k++){
            s += u[j][k] * u[i][k];
          }

          for (k = l; k < n; k++){
            u[j][k] += s * rv1[k];
          }
        }

        for (k = l; k < n; k++){
          u[i][k] *= scale;
        }
      }
    }

    anorm = std::max(anorm, fabs(w[i]) + fabs(rv1[i]));
  }

  // accumulation of right-hand transformations
  for (l = i = n; i--; l--){
    if (l < n){
      if (g != 0.0){
        for (j = l; j < n; j++){
          // double division avoids possible underflow
          v[j][i] = (u[i][j] / u[i][l]) / g;
        }

        for (j = l; j < n; j++){
          s = 0.0;
          for (k = l; k < n; k++){
            s += u[i][k] * v[k][j];
          }

          for (k = l; k < n; k++){
            v[k][j] += s * v[k][i];
          }
        }
      }

      for (j = l; j < n; j++){
        v[i][j] = v[j][i] = 0.0;
      }
    }

    v[i][i] = 1.0;
    g = rv1[i];
  }

  // accumulation of left-hand transformations
  for (l = i = std::min(m, n); i--; l--){
    g = w[i];

    for (j = l; j < n; j++){
      u[i][j] = 0.0;
    }

    if (g != 0.0){
      g = 1.0 / g;

      for (j = l; j < n; j++){
        s = 0.0;
        for (k = l; k < m; k++){
          s += u[k][i] * u[k][j];
        }

        // double division avoids possible underflow
        f = (s / u[i][i]) * g;

        for (k = i; k < m; k++){
          u[k][j] += f * u[k][i];
        }
      }

      for (j = i; j < m; j++){
        u[j][i] *= g;
      }
    }else{
      for (j = i; j < m; j++){
        u[j][i] = 0.0;
      }
    }

    u[i][i]++;
  }

  // diagonalization of the bidiagonal form
  for (k = n; k--; ){
    for (its = 1; its <= 30; its++){
      flag = 1;

      // test for splitting
      for (l = k + 1; l--; ){
        // rv1 [0] is always zero, so there is no exit
        nm = l - 1;

        if (fabs(rv1[l]) + anorm == anorm){
          flag = 0;
          break;
        }

        if (fabs(w[nm]) + anorm == anorm){
          break;
        }
      }

      if (flag){
        // cancellation of rv1 [l] if l greater than 0
        c = 0.0;
        s = 1.0;

        for (i = l; i <= k; i++){
          f = s * rv1[i];
          rv1[i] *= c;

          if (fabs (f) + anorm == anorm){
            break;
          }

          g = w[i];
          h = hypot(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g * h;
          s = -f * h;

          for (j = 0; j < m; j++){
            y = u[j][nm];
            z = u[j][i];
            u[j][nm] = y * c + z * s;
            u[j][i] = z * c - y * s;
          }
        }
      }

      // test for convergence
      z = w[k];

      if (l == k){
        if (z < 0.0){
          w[k] = -z;
          for (j = 0; j < n; j++){
            v[j][k] = -v[j][k];
          }
        }
        break;
      }

      if (its == 30){
        throw k;
      }

      // shift from bottom 2 by 2 minor
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = hypot(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

      // next qr transformation
      c = s = 1.0;

      for (j = l; j < k; j++){
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g *= c;
        z = hypot(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;

        for (jj = 0; jj < n; jj++)
        {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }

        z = hypot(f, h);
        w[j] = z;

        // rotation can be arbitrary if z is zero
        if (z != 0.0){
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }

        f = c * g + s * y;
        x = c * y - s * g;

        for (jj = 0; jj < m; jj++){
          y = u[jj][j];
          z = u[jj][i];
          u[jj][j] = y * c + z * s;
          u[jj][i] = z * c - y * s;
        }
      }

      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  //////////////////////////////////////////////
  //sort the eigenvalues in descending order
  for (i = 0; i < n - 1; i++){
    p = w[k = i];

    for (j = i + 1; j < n; j++){
      if (w[j] >= p){
        p = w[k = j];
      }
    }

    if (k != i){
      w[k] = w[i];
      w[i] = p;

      for (j = 0; j < n; j++){
        p = v[j][i];
        v[j][i] = v[j][k];
        v[j][k] = p;
      }

      for (j = 0; j < m; j++){
        p = u[j][i];
        u[j][i] = u[j][k];
        u[j][k] = p;
      }
    }
  }

  //free the acquired memory
  delete[] rv1;
}

#endif //MATRIX_H_INCLUDED
