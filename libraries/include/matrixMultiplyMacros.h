/*

* Copyright (c) 2011, Jernej Barbic, Yili Zhao, University of Southern California
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of University of Southern California, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY JERNEJ BARBIC, YILI ZHAO AND UNIVERSITY OF SOUTHERN CALIFORNIA 
* ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
* IN NO EVENT SHALL JERNEJ BARBIC, YILI ZHAO OR UNIVERSITY OF SOUTHERN CALIFORNIA BE 
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  3x3 and 4x4 matrix multiplication macros

  This code was used in the following publication:

  Jernej Barbic, Yili Zhao:
  Real-time Large-deformation Substructuring, ACM Transactions on Graphics (SIGGRAPH 2011), 30(4), Aug 2011

  @article{Barbic:2011:RLS,
    author =  {Jernej Barbi\v{c} and Yili Zhao},
    journal = {ACM Trans. on Graphics (SIGGRAPH 2011)},
    number =  "4",
    title =   "Real-time Large-deformation Substructuring",
    volume =  "30",
    year =    "2011",
    pages =   "91:1--91:7",
  }

  Authors of this code: Jernej Barbic, Yili Zhao
*/

#ifndef _MATRIX_MULTIPLY_MACROS_H_
#define _MATRIX_MULTIPLY_MACROS_H_

#ifndef MIN
  #define MIN(x,y) ((x)<(y) ? (x) : (y))
#endif

#ifndef MAX
  #define MAX(x,y) ((x)>(y) ? (x) : (y))
#endif

#ifndef ELT
  #define ELT(rows,i,j) (((long)j)*((long)rows)+((long)i))
#endif

// all macros assume that the 3x3 and 4x4 matrices are stored row-major

// a = b
#define VECTOR_SET3(a,b)\
(a)[0] = (b)[0];\
(a)[1] = (b)[1];\
(a)[2] = (b)[2];

// a = b
#define VECTOR_SET4(a,b)\
(a)[0] = (b)[0];\
(a)[1] = (b)[1];\
(a)[2] = (b)[2];\
(a)[3] = (b)[3];

// c = a + b
#define VECTOR_ADD3(a,b,c)\
(c)[0] = (a)[0] + (b)[0];\
(c)[1] = (a)[1] + (b)[1];\
(c)[2] = (a)[2] + (b)[2];

// c = a - b
#define VECTOR_SUBTRACT3(a,b,c)\
(c)[0] = (a)[0] - (b)[0];\
(c)[1] = (a)[1] - (b)[1];\
(c)[2] = (a)[2] - (b)[2];

// c = a + scalar * b
#define VECTOR_SCALE_ADD3(a,scalar,b,c)\
(c)[0] = (a)[0] + (scalar) * (b)[0];\
(c)[1] = (a)[1] + (scalar) * (b)[1];\
(c)[2] = (a)[2] + (scalar) * (b)[2];

// a += b
#define VECTOR_ADDEQUAL3(a,b)\
(a)[0] += (b)[0];\
(a)[1] += (b)[1];\
(a)[2] += (b)[2];

// a -= b
#define VECTOR_SUBTRACTEQUAL3(a,b)\
(a)[0] -= (b)[0];\
(a)[1] -= (b)[1];\
(a)[2] -= (b)[2];

// a *= scalar
#define VECTOR_SCALE3(a,scalar)\
(a)[0] *= (scalar);\
(a)[1] *= (scalar);\
(a)[2] *= (scalar);

// c = a + b
#define VECTOR_ADD4(a,b,c)\
(c)[0] = (a)[0] + (b)[0];\
(c)[1] = (a)[1] + (b)[1];\
(c)[2] = (a)[2] + (b)[2];\
(c)[3] = (a)[3] + (b)[3];

// c = a - b
#define VECTOR_SUBTRACT4(a,b,c)\
(c)[0] = (a)[0] - (b)[0];\
(c)[1] = (a)[1] - (b)[1];\
(c)[2] = (a)[2] - (b)[2];\
(c)[3] = (a)[3] - (b)[3];

// c = a + scalar * b
#define VECTOR_SCALE_ADD4(a,scalar,b,c)\
(c)[0] = (a)[0] + (scalar) * (b)[0];\
(c)[1] = (a)[1] + (scalar) * (b)[1];\
(c)[2] = (a)[2] + (scalar) * (b)[2];\
(c)[3] = (a)[3] + (scalar) * (b)[3];

// a += b
#define VECTOR_ADDEQUAL4(a,b)\
(a)[0] += (b)[0];\
(a)[1] += (b)[1];\
(a)[2] += (b)[2];\
(a)[3] += (b)[3];

// a -= b
#define VECTOR_SUBTRACTEQUAL4(a,b)\
(a)[0] -= (b)[0];\
(a)[1] -= (b)[1];\
(a)[2] -= (b)[2];\
(a)[3] -= (b)[3];

// a *= scalar
#define VECTOR_SCALE4(a,scalar)\
(a)[0] *= (scalar);\
(a)[1] *= (scalar);\
(a)[2] *= (scalar);\
(a)[3] *= (scalar);

// a = b
#define MATRIX_SET3X3(a,b)\
  (a)[0] = (b)[0];\
  (a)[1] = (b)[1];\
  (a)[2] = (b)[2];\
  (a)[3] = (b)[3];\
  (a)[4] = (b)[4];\
  (a)[5] = (b)[5];\
  (a)[6] = (b)[6];\
  (a)[7] = (b)[7];\
  (a)[8] = (b)[8];

// a = b
#define MATRIX_SET4X4(a,b)\
  (a)[0] = (b)[0];\
  (a)[1] = (b)[1];\
  (a)[2] = (b)[2];\
  (a)[3] = (b)[3];\
  (a)[4] = (b)[4];\
  (a)[5] = (b)[5];\
  (a)[6] = (b)[6];\
  (a)[7] = (b)[7];\
  (a)[8] = (b)[8];\
  (a)[9] = (b)[9];\
  (a)[10] = (b)[10];\
  (a)[11] = (b)[11];\
  (a)[12] = (b)[12];\
  (a)[13] = (b)[13];\
  (a)[14] = (b)[14];\
  (a)[15] = (b)[15];

// a *= scalar
#define MATRIX_SCALE3X3(a,scalar)\
  (a)[0] *= (scalar);\
  (a)[1] *= (scalar);\
  (a)[2] *= (scalar);\
  (a)[3] *= (scalar);\
  (a)[4] *= (scalar);\
  (a)[5] *= (scalar);\
  (a)[6] *= (scalar);\
  (a)[7] *= (scalar);\
  (a)[8] *= (scalar);

// a += b
#define MATRIX_ADDEQUAL3X3(a,b)\
  (a)[0] += (b)[0];\
  (a)[1] += (b)[1];\
  (a)[2] += (b)[2];\
  (a)[3] += (b)[3];\
  (a)[4] += (b)[4];\
  (a)[5] += (b)[5];\
  (a)[6] += (b)[6];\
  (a)[7] += (b)[7];\
  (a)[8] += (b)[8];

// c = a + b
#define MATRIX_ADD3X3(a,b,c)\
  (c)[0] = (a)[0] + (b)[0];\
  (c)[1] = (a)[1] + (b)[1];\
  (c)[2] = (a)[2] + (b)[2];\
  (c)[3] = (a)[3] + (b)[3];\
  (c)[4] = (a)[4] + (b)[4];\
  (c)[5] = (a)[5] + (b)[5];\
  (c)[6] = (a)[6] + (b)[6];\
  (c)[7] = (a)[7] + (b)[7];\
  (c)[8] = (a)[8] + (b)[8];

// c = a * b
#define MATRIX_MULTIPLY3X3(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[3]+(a)[2]*(b)[6];\
(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[4]+(a)[2]*(b)[7];\
(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[5]+(a)[2]*(b)[8];\
(c)[3]=(a)[3]*(b)[0]+(a)[4]*(b)[3]+(a)[5]*(b)[6];\
(c)[4]=(a)[3]*(b)[1]+(a)[4]*(b)[4]+(a)[5]*(b)[7];\
(c)[5]=(a)[3]*(b)[2]+(a)[4]*(b)[5]+(a)[5]*(b)[8];\
(c)[6]=(a)[6]*(b)[0]+(a)[7]*(b)[3]+(a)[8]*(b)[6];\
(c)[7]=(a)[6]*(b)[1]+(a)[7]*(b)[4]+(a)[8]*(b)[7];\
(c)[8]=(a)[6]*(b)[2]+(a)[7]*(b)[5]+(a)[8]*(b)[8];

//c = a * bT
#define MATRIX_MULTIPLY3X3ABT(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2];\
(c)[1]=(a)[0]*(b)[3]+(a)[1]*(b)[4]+(a)[2]*(b)[5];\
(c)[2]=(a)[0]*(b)[6]+(a)[1]*(b)[7]+(a)[2]*(b)[8];\
(c)[3]=(a)[3]*(b)[0]+(a)[4]*(b)[1]+(a)[5]*(b)[2];\
(c)[4]=(a)[3]*(b)[3]+(a)[4]*(b)[4]+(a)[5]*(b)[5];\
(c)[5]=(a)[3]*(b)[6]+(a)[4]*(b)[7]+(a)[5]*(b)[8];\
(c)[6]=(a)[6]*(b)[0]+(a)[7]*(b)[1]+(a)[8]*(b)[2];\
(c)[7]=(a)[6]*(b)[3]+(a)[7]*(b)[4]+(a)[8]*(b)[5];\
(c)[8]=(a)[6]*(b)[6]+(a)[7]*(b)[7]+(a)[8]*(b)[8];

//c = aT * b
#define MATRIX_MULTIPLY3X3ATB(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[3]*(b)[3]+(a)[6]*(b)[6];\
(c)[1]=(a)[0]*(b)[1]+(a)[3]*(b)[4]+(a)[6]*(b)[7];\
(c)[2]=(a)[0]*(b)[2]+(a)[3]*(b)[5]+(a)[6]*(b)[8];\
(c)[3]=(a)[1]*(b)[0]+(a)[4]*(b)[3]+(a)[7]*(b)[6];\
(c)[4]=(a)[1]*(b)[1]+(a)[4]*(b)[4]+(a)[7]*(b)[7];\
(c)[5]=(a)[1]*(b)[2]+(a)[4]*(b)[5]+(a)[7]*(b)[8];\
(c)[6]=(a)[2]*(b)[0]+(a)[5]*(b)[3]+(a)[8]*(b)[6];\
(c)[7]=(a)[2]*(b)[1]+(a)[5]*(b)[4]+(a)[8]*(b)[7];\
(c)[8]=(a)[2]*(b)[2]+(a)[5]*(b)[5]+(a)[8]*(b)[8];

// w = A * v
#define MATRIX_VECTOR_MULTIPLY3X3(A,v,w)\
(w)[0] = (A)[0]*(v)[0]+(A)[1]*(v)[1]+(A)[2]*(v)[2];\
(w)[1] = (A)[3]*(v)[0]+(A)[4]*(v)[1]+(A)[5]*(v)[2];\
(w)[2] = (A)[6]*(v)[0]+(A)[7]*(v)[1]+(A)[8]*(v)[2];

// w = A^T * v
#define MATRIX_VECTOR_MULTIPLY3X3T(A,v,w)\
(w)[0] = (A)[0]*(v)[0]+(A)[3]*(v)[1]+(A)[6]*(v)[2];\
(w)[1] = (A)[1]*(v)[0]+(A)[4]*(v)[1]+(A)[7]*(v)[2];\
(w)[2] = (A)[2]*(v)[0]+(A)[5]*(v)[1]+(A)[8]*(v)[2];

// w = A * v
#define MATRIX_VECTOR_MULTIPLY4X4(A,v,w)\
(w)[0] = (A)[0]*(v)[0]+(A)[1]*(v)[1]+(A)[2]*(v)[2]+(A)[3]*(v)[3];\
(w)[1] = (A)[4]*(v)[0]+(A)[5]*(v)[1]+(A)[6]*(v)[2]+(A)[7]*(v)[3];\
(w)[2] = (A)[8]*(v)[0]+(A)[9]*(v)[1]+(A)[10]*(v)[2]+(A)[11]*(v)[3];\
(w)[3] = (A)[12]*(v)[0]+(A)[13]*(v)[1]+(A)[14]*(v)[2]+(A)[15]*(v)[3];

// w = A^T * v
#define MATRIX_VECTOR_MULTIPLY4X4T(A,v,w)\
(w)[0] = (A)[0]*(v)[0]+(A)[4]*(v)[1]+(A)[8]*(v)[2]+(A)[12]*(v)[3];\
(w)[1] = (A)[1]*(v)[0]+(A)[5]*(v)[1]+(A)[9]*(v)[2]+(A)[13]*(v)[3];\
(w)[2] = (A)[2]*(v)[0]+(A)[6]*(v)[1]+(A)[10]*(v)[2]+(A)[14]*(v)[3];\
(w)[3] = (A)[3]*(v)[0]+(A)[7]*(v)[1]+(A)[11]*(v)[2]+(A)[15]*(v)[3];

// B = A^T
#define MATRIX_TRANSPOSE3X3(A,B)\
(B)[0] = (A)[0]; (B)[1] = (A)[3]; (B)[2] = (A)[6];\
(B)[3] = (A)[1]; (B)[4] = (A)[4]; (B)[5] = (A)[7];\
(B)[6] = (A)[2]; (B)[7] = (A)[5]; (B)[8] = (A)[8];

// A = u * v^T
#define VECTOR_TENSOR_PRODUCT3X3(u,v,A)\
(A)[0] = (u)[0] * (v)[0]; (A)[1] = (u)[0] * (v)[1]; (A)[2] = (u)[0] * (v)[2];\
(A)[3] = (u)[1] * (v)[0]; (A)[4] = (u)[1] * (v)[1]; (A)[5] = (u)[1] * (v)[2];\
(A)[6] = (u)[2] * (v)[0]; (A)[7] = (u)[2] * (v)[1]; (A)[8] = (u)[2] * (v)[2];

// A += u * v^T
#define VECTOR_TENSOR_PRODUCT_ADD3X3(u,v,A)\
(A)[0] += (u)[0] * (v)[0]; (A)[1] += (u)[0] * (v)[1]; (A)[2] += (u)[0] * (v)[2];\
(A)[3] += (u)[1] * (v)[0]; (A)[4] += (u)[1] * (v)[1]; (A)[5] += (u)[1] * (v)[2];\
(A)[6] += (u)[2] * (v)[0]; (A)[7] += (u)[2] * (v)[1]; (A)[8] += (u)[2] * (v)[2];

// scalar = A : B
#define MATRIX_DOT_PRODUCT3X3(A,B)\
  ( (A)[0] * (B)[0] + (A)[1] * (B)[1] + (A)[2] * (B)[2] +\
    (A)[3] * (B)[3] + (A)[4] * (B)[4] + (A)[5] * (B)[5] +\
    (A)[6] * (B)[6] + (A)[7] * (B)[7] + (A)[8] * (B)[8] )

#define SKEW_MATRIX(a, A)\
(A)[0] =       0;  (A)[1] = -(a)[2]; (A)[2] =  (a)[1];\
(A)[3] =  (a)[2];  (A)[4] =       0; (A)[5] = -(a)[0];\
(A)[6] = -(a)[1];  (A)[7] =  (a)[0]; (A)[8] =       0;

#define SYM_MATRIX(a, A)\
(A)[0] =  (a)[0];  (A)[1] =  (a)[1]; (A)[2] =  (a)[2];\
(A)[3] =  (a)[1];  (A)[4] =  (a)[3]; (A)[5] =  (a)[4];\
(A)[6] =  (a)[2];  (A)[7] =  (a)[4]; (A)[8] =  (a)[5];

#define MATRIX_DETERMINANT3X3(A)\
  ( (A)[0] * ( (A)[4] * (A)[8] - (A)[5] * (A)[7] ) +\
    (A)[1] * ( (A)[5] * (A)[6] - (A)[3] * (A)[8] ) +\
    (A)[2] * ( (A)[3] * (A)[7] - (A)[4] * (A)[6] ) )

// a = skew( 0.5 * (A - A^T) )
#define SKEW_PART(A, a)\
(a)[0] = 0.5 * ((A)[7] - (A)[5]);\
(a)[1] = 0.5 * ((A)[2] - (A)[6]);\
(a)[2] = 0.5 * ((A)[3] - (A)[1]);

// a = upper-triangle( 0.5 * (A + A^T) )
#define SYM_PART(A, a)\
(a)[0] = ((A)[0]);\
(a)[1] = 0.5 * ((A)[1] + (A)[3]);\
(a)[2] = 0.5 * ((A)[2] + (A)[6]);\
(a)[3] = ((A)[4]);\
(a)[4] = 0.5 * ((A)[5] + (A)[7]);\
(a)[5] = ((A)[8]);\

//  a . b
#define VECTOR_DOT_PRODUCT3(a,b)\
  ( (a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2] )

//  a . b
#define VECTOR_DOT_PRODUCT4(a,b)\
  ( (a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2] + (a)[3] * (b)[3] )

// c = a x b
#define VECTOR_CROSS_PRODUCT(a,b,c)\
  (c)[0] = ((a)[1]) * ((b)[2]) - ((b)[1]) * ((a)[2]);\
  (c)[1] = ((b)[0]) * ((a)[2]) - ((a)[0]) * ((b)[2]);\
  (c)[2] = ((a)[0]) * ((b)[1]) - ((b)[0]) * ((a)[1])

// c += a x b
#define VECTOR_CROSS_PRODUCT_ADD(a,b,c)\
  (c)[0] += ((a)[1]) * ((b)[2]) - ((b)[1]) * ((a)[2]);\
  (c)[1] += ((b)[0]) * ((a)[2]) - ((a)[0]) * ((b)[2]);\
  (c)[2] += ((a)[0]) * ((b)[1]) - ((b)[0]) * ((a)[1])

#endif

