/*
* Copyright (c) 2008, Carnegie Mellon University
* All rights reserved.
*
* Code author: Jernej Barbic
* Research: Jernej Barbic, Doug L. James
* Funding: NSF, Link Foundation
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Carnegie Mellon University, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE* DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  This class implements a 3x3 matrix, with elementary algebra,
  including the computation of eigenvalues and eigenvectors of a *symmetric* 3x3 matrix.
  (using a public domain external routine; see mat3d.cpp and eig3.h)

  Note: this code was inspired by Andrew Willmott's VL and SVL Libraries:
    http://www.cs.cmu.edu/afs/cs/user/ajw/www/software/index.html#VL
    (these two libraries contain a lot of useful functionality and
     are highly recommended)
    My library offers just the basic functionality, hence "minivector".
    It was written from scratch for a course project at CMU.

  Version: 1.2
*/

#ifndef _MINIMAT3D_H_
#define _MINIMAT3D_H_

#include "vec3d.h"

class Mat3d {
public:
  
  inline Mat3d() {}

  /*
          [ x0  x1  x2 ]
    M  =  [ x3  x4  x5 ]
          [ x6  x7  x8 ]
  */
  inline Mat3d(double x0, double x1, double x2,
               double x3, double x4, double x5,
               double x6, double x7, double x8);
  inline Mat3d(double mat[9]); // "mat" must be given in row-major order
  inline Mat3d(Vec3d rows[3]); 
  inline Mat3d(Vec3d row0, Vec3d row1, Vec3d row2);
  inline Mat3d(double diag); // create a diagonal matrix with all entries "diag" (can create zero matrix by passing 0.0)

  inline void set(double x0, double x1, double x2,
		  double x3, double x4, double x5,
		  double x6, double x7, double x8);
  inline void set(double value); // set all elements to value

  inline Mat3d & operator=(const Mat3d & source); 

  inline Mat3d operator+ (const Mat3d & );
  inline Mat3d & operator+= (const Mat3d & );

  inline Mat3d operator- (const Mat3d & );
  inline Mat3d & operator-= (const Mat3d & );

  inline Mat3d & operator*= (double scalar);
  inline Mat3d & operator/= (double scalar);

  friend inline Mat3d operator* (double scalar, const Mat3d & mat2);
  friend inline Mat3d operator/ (double scalar, const Mat3d & mat2);

  inline Mat3d & multiplyDiagRight(Vec3d & v); // M = M * diag(v)
  inline Mat3d & multiplyDiagLeft(Vec3d & v); // M = diag(v) * M

  friend inline Mat3d tensorProduct(const Vec3d & vec1, const Vec3d & vec2);

  friend inline Mat3d inv(const Mat3d &); // inverse matrix
  friend inline double det(const Mat3d & mat); // determinant
  friend inline Mat3d trans(const Mat3d & mat); // transpose

  friend inline std::ostream &operator << (std::ostream &s, const Mat3d &v);
  inline void print();

  inline Vec3d & operator[] (int index); // M[i] returns i-th row
  inline const Vec3d & operator[] (int index) const; // M[i] returns i-th row

  // converts the matrix into a C array of length 9
  inline void convertToArray(double * array) const; // in row-major order

  // matrix-vector multiply
  inline const Vec3d operator* (const Vec3d & vec) const;

  // matrix-matrix multiply
  inline const Mat3d operator* (const Mat3d & mat2) const;

  // Computes eigenvalues and eigenvectors of a 3x3 matrix
  // Assumes symmetric matrix; contents of matrix "a" are not modified by the routine
  // Eigenvalues are sorted in decreasing order (not decreasing absolute order)
  // Returned eigenvectors are unit length
  friend void eigen_sym(Mat3d & a, Vec3d & eig_val, Vec3d eig_vec[3]);

  // NOTE: This particular routine is not publicly released, as its implementation
  // was taken from Numerical Recipes, which is not free software.
  // (you can use the eigen_sym routine above which is public domain)
  // Computes eigenvalues and eigenvectors of a 3x3 matrix, using Jacobi Iteration
  // Assumes symmetric matrix; contents of matrix "a" are overwritten (destroyed)
  // Returns true if iteration succeeded in making the sum of abs values of non-diagonal values below epsilon, and false otherwise
  // Default epsilon is machine precision (which always converged with our matrices)
  // Eigenvalues are sorted in decreasing absolute order
  // Returned eigenvectors are unit length
  friend bool eigen_sym_NR(Mat3d & a, Vec3d & eig_val, Vec3d eig_vec[3],
                           int maxIterations=50, double epsilon=0.0);

protected:
  Vec3d elt[3]; // the three rows of the matrix
};

// ===== below is the implementation =====

inline Mat3d::Mat3d(double x0_g, double x1_g, double x2_g,
                    double x3_g, double x4_g, double x5_g,
                    double x6_g, double x7_g, double x8_g)
{
  elt[0] = Vec3d(x0_g,x1_g,x2_g);
  elt[1] = Vec3d(x3_g,x4_g,x5_g);
  elt[2] = Vec3d(x6_g,x7_g,x8_g);
}

inline Mat3d::Mat3d(double mat[9])
{
  elt[0] = Vec3d(mat[0],mat[1],mat[2]);
  elt[1] = Vec3d(mat[3],mat[4],mat[5]);
  elt[2] = Vec3d(mat[6],mat[7],mat[8]);
}

inline Mat3d::Mat3d(Vec3d rows[3]) 
{
  elt[0] = rows[0];
  elt[1] = rows[1];
  elt[2] = rows[2];
}

inline Mat3d::Mat3d(Vec3d row0, Vec3d row1, Vec3d row2)
{
  elt[0] = row0;
  elt[1] = row1;
  elt[2] = row2;
}

inline Mat3d::Mat3d(double diag)
{
  elt[0] = Vec3d(diag,0,0);
  elt[1] = Vec3d(0,diag,0);
  elt[2] = Vec3d(0,0,diag);
}

inline void Mat3d::set(double x0, double x1, double x2,
		       double x3, double x4, double x5,
		       double x6, double x7, double x8)
{
  elt[0][0] = x0;
  elt[0][1] = x1;
  elt[0][2] = x2;

  elt[1][0] = x3;
  elt[1][1] = x4;
  elt[1][2] = x5;

  elt[2][0] = x6;
  elt[2][1] = x7;
  elt[2][2] = x8;
}

inline void Mat3d::set(double value)
{
  elt[0][0] = elt[0][1] = elt[0][2] = 
  elt[1][0] = elt[1][1] = elt[1][2] =
  elt[2][0] = elt[2][1] = elt[2][2] = value;
}

inline Mat3d & Mat3d::operator=(const Mat3d & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];
  elt[2] = source.elt[2];

  return *this;
}

inline Mat3d Mat3d::operator+ (const Mat3d & mat2)
{
  Mat3d sum = *this;
  sum.elt[0] += mat2.elt[0];
  sum.elt[1] += mat2.elt[1];
  sum.elt[2] += mat2.elt[2];

  return sum;
}

inline Mat3d & Mat3d::operator+= (const Mat3d & mat2)
{
  elt[0] += mat2.elt[0];
  elt[1] += mat2.elt[1];
  elt[2] += mat2.elt[2];
  return *this;
}

inline Mat3d Mat3d::operator- (const Mat3d & mat2)
{
  Mat3d sum = *this;
  sum.elt[0] -= mat2.elt[0];
  sum.elt[1] -= mat2.elt[1];
  sum.elt[2] -= mat2.elt[2];

  return sum;
}

inline Mat3d & Mat3d::operator-= (const Mat3d & mat2)
{
  elt[0] -= mat2.elt[0];
  elt[1] -= mat2.elt[1];
  elt[2] -= mat2.elt[2];

  return *this;
}

inline Vec3d & Mat3d::operator[] (int index)
{
  return elt[index];
}

inline const Vec3d & Mat3d::operator[] (int index) const
{
  return elt[index];
}

inline Mat3d & Mat3d::operator*= (double scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  elt[2] *= scalar;

  return *this;
}

inline Mat3d & Mat3d::multiplyDiagRight(Vec3d & v) // M = M * diag(v)
{
  elt[0][0] *= v[0];
  elt[1][0] *= v[0];
  elt[2][0] *= v[0];

  elt[0][1] *= v[1];
  elt[1][1] *= v[1];
  elt[2][1] *= v[1];

  elt[0][2] *= v[2];
  elt[1][2] *= v[2];
  elt[2][2] *= v[2];

  return *this;
}

inline Mat3d & Mat3d::multiplyDiagLeft(Vec3d & v) // M = diag(v) * M
{
  elt[0] *= v[0];
  elt[1] *= v[1];
  elt[2] *= v[2];

  return *this;
}

inline Mat3d operator* (double scalar, const Mat3d & mat2)
{
  Mat3d result = mat2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;
  result.elt[2] *= scalar;

  return result;
}

inline Mat3d operator/ (double scalar, const Mat3d & mat2)
{
  Mat3d result = mat2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;
  result.elt[2] /= scalar;

  return result;
}
 
inline Mat3d tensorProduct(Vec3d & vecA, Vec3d & vecB)
{
  Mat3d result(vecA[0]*vecB[0],vecA[0]*vecB[1],vecA[0]*vecB[2],
	       vecA[1]*vecB[0],vecA[1]*vecB[1],vecA[1]*vecB[2],
	       vecA[2]*vecB[0],vecA[2]*vecB[1],vecA[2]*vecB[2]);

  return result;
}

inline Mat3d & Mat3d::operator/= (double scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  elt[2] /= scalar;

  return *this;
}

inline const Vec3d Mat3d::operator* (const Vec3d & vec) const
{
  return(Vec3d(
    dot(elt[0],vec),
    dot(elt[1],vec),
    dot(elt[2],vec)));
}

inline const Mat3d Mat3d::operator* (const Mat3d & mat2) const
{
  return(Mat3d(
    dot(elt[0],Vec3d(mat2.elt[0][0],mat2.elt[1][0],mat2.elt[2][0])),
    dot(elt[0],Vec3d(mat2.elt[0][1],mat2.elt[1][1],mat2.elt[2][1])),
    dot(elt[0],Vec3d(mat2.elt[0][2],mat2.elt[1][2],mat2.elt[2][2])),

    dot(elt[1],Vec3d(mat2.elt[0][0],mat2.elt[1][0],mat2.elt[2][0])),
    dot(elt[1],Vec3d(mat2.elt[0][1],mat2.elt[1][1],mat2.elt[2][1])),
    dot(elt[1],Vec3d(mat2.elt[0][2],mat2.elt[1][2],mat2.elt[2][2])),

    dot(elt[2],Vec3d(mat2.elt[0][0],mat2.elt[1][0],mat2.elt[2][0])),
    dot(elt[2],Vec3d(mat2.elt[0][1],mat2.elt[1][1],mat2.elt[2][1])),
    dot(elt[2],Vec3d(mat2.elt[0][2],mat2.elt[1][2],mat2.elt[2][2])) ));
}

inline Mat3d inv(const Mat3d & mat)
{
  double invDeterminant = 1.0 / 
    (-mat[0][2] * mat[1][1] * mat[2][0] + 
      mat[0][1] * mat[1][2] * mat[2][0] + 
      mat[0][2] * mat[1][0] * mat[2][1] - 
      mat[0][0] * mat[1][2] * mat[2][1] - 
      mat[0][1] * mat[1][0] * mat[2][2] + 
      mat[0][0] * mat[1][1] * mat[2][2] );

  return Mat3d(
    invDeterminant * (-mat[1][2] * mat[2][1] + mat[1][1] * mat[2][2]),
    invDeterminant * (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]),
    invDeterminant * (-mat[0][2] * mat[1][1] + mat[0][1] * mat[1][2]),
    invDeterminant * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]),
    invDeterminant * (-mat[0][2] * mat[2][0] + mat[0][0] * mat[2][2]),
    invDeterminant * (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]),
    invDeterminant * (-mat[1][1] * mat[2][0] + mat[1][0] * mat[2][1]),
    invDeterminant * (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]),
    invDeterminant * (-mat[0][1] * mat[1][0] + mat[0][0] * mat[1][1])
  );
}

inline double det(const Mat3d & mat) 
{
  return
   (-mat[0][2] * mat[1][1] * mat[2][0] + 
     mat[0][1] * mat[1][2] * mat[2][0] + 
     mat[0][2] * mat[1][0] * mat[2][1] - 
     mat[0][0] * mat[1][2] * mat[2][1] - 
     mat[0][1] * mat[1][0] * mat[2][2] + 
     mat[0][0] * mat[1][1] * mat[2][2] );
}

inline Mat3d trans(const Mat3d & mat)
{
  return
       Mat3d( mat[0][0], mat[1][0], mat[2][0],
              mat[0][1], mat[1][1], mat[2][1],
              mat[0][2], mat[1][2], mat[2][2] );
}

inline void Mat3d::convertToArray(double * array) const // in row-major order
{
  array[0] = elt[0][0]; array[1] = elt[0][1]; array[2] = elt[0][2];
  array[3] = elt[1][0]; array[4] = elt[1][1]; array[5] = elt[1][2];
  array[6] = elt[2][0]; array[7] = elt[2][1]; array[8] = elt[2][2];
}

inline std::ostream &operator << (std::ostream &s, const Mat3d &v)
{
  double a00 = v.elt[0][0]; double a01 = v.elt[0][1]; double a02 = v.elt[0][2];
  double a10 = v.elt[1][0]; double a11 = v.elt[1][1]; double a12 = v.elt[1][2];
  double a20 = v.elt[2][0]; double a21 = v.elt[2][1]; double a22 = v.elt[2][2];

  return(
    s << '[' << a00 << ' ' << a01 << ' ' << a02 << ']' << std::endl <<
    s << '[' << a10 << ' ' << a11 << ' ' << a12 << ']' << std::endl <<
    s << '[' << a20 << ' ' << a21 << ' ' << a22 << ']'
  );
}

inline void Mat3d::print()
{
  double a00 = elt[0][0]; double a01 = elt[0][1]; double a02 = elt[0][2];
  double a10 = elt[1][0]; double a11 = elt[1][1]; double a12 = elt[1][2];
  double a20 = elt[2][0]; double a21 = elt[2][1]; double a22 = elt[2][2];
  
  printf("[%G %G %G;\n", a00, a01, a02);
  printf(" %G %G %G;\n", a10, a11, a12);
  printf(" %G %G %G]\n", a20, a21, a22);
}

#endif

