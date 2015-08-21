/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "sparseMatrix" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC   *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

/*
  The "SparseMatrix" class implements double-precision sparse matrices 
  with common algebraic operations such as incremental construction, 
  addition, mtx-vec multiplication, row-column deletion. 
  The matrices can be loaded from and saved to a file.

  The matrix can be rectangular (it need not be square). 
  In memory, the matrix is stored in a row-based format. 
  For each matrix row, the class stores the integer
  indices of the columns containing non-zero entries, 
  together with the corresponding double precision values. 
  All quantities (rows, columns, etc.) in this class are 0-indexed.

  Also included is a Conjugate Gradient iterative linear system solver 
  (for positive-definite large sparse symmetric matrices).
  The solver can be used without preconditioning, or with diagonal (Jacobi) preconditoning.
  The solver either uses an explicitly provided sparse matrix, or you can only
  provide the matrix-vector multiplication routine (without explicitly giving the matrix).
  The CG Solver was implemented by following Jonathan Shewchuk's 
  "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain":
  http://www.cs.cmu.edu/~jrs/jrspapers.html#cg

  The class also includes a Gauss-Seidel iterative linear system solver.

  There are two classes available to the user: SparseMatrixOutline and SparseMatrix.
  The SparseMatrixOutline class should be used to construct the non-zero 
  locations in the matrix, and (optionally) assign values to them.
  You should then transform it into a SparseMatrix class, via SparseMatrix's
  class constructor. The SparseMatrix class has all the computational routines,
  but you can only add new non-zero entries to SparseMatrixOutline, not
  to SparseMatrix.  The reason for this separation is that SparseMatrixOutline 
  does not know the number of matrix entries ahead of time, so it uses STL's 
  "map" datastructure to store the matrix data. In the SparseMatrix class, however, the number 
  of sparse entries and their locations are fixed, so all operations can use 
  known-length C arrays, which is faster.

  So: you should first create an instance of SparseMatrixOutline, then create 
  an instance of SparseMatrix by passing the SparseMatrixOutline object to 
  SparseMatrix's constructor.
  If your matrix is a text file on disk, you can load it to SparseMatrixOutline, 
  or directly load it into SparseMatrix (which will, however, internally still 
  proceed via SparseMatrixOutline).

  The text disk file format is as follows:
  <number of matrix rows>
  <number of matrix columns>
  one or more data lines, each giving one matrix entry, in the format:
  <row index> <column index> <data value> 
  (indices are 0-indexed)

  Example:
  
    [0 17 -1 0]
  A=[0  5  0 0]
    [3  8  6 0]  

  would be given as:

  3
  4
  0 1 17 
  0 2 -1
  1 1 5
  2 0 3
  2 1 8
  2 2 6
*/

#include <string.h>
#include <vector>
using std::vector;
using std::shared_ptr;
#include <map>

#include "sparseMatrixOutline.h"

class SparseSubMatrixLinkage;
class SparseSuperMatrixLinkage;
class SparseMatrixIndexRemapper;

class SparseMatrix: public std::enable_shared_from_this<SparseMatrix>
{
public:
    
    SparseMatrix(const char * filename); // load from text file (same text file format as SparseMatrixOutline)
    SparseMatrix(const SparseMatrixOutline * sparseMatrixOutline); // create it from the outline
    SparseMatrix(const SparseMatrixOutline& sparseMatrixOutline); // create it from the outline
    SparseMatrix(const SparseMatrix & source); // copy constructor
    virtual ~SparseMatrix();
    
    SparseMatrixOutline GetTopology() const; // create an outline for the topology of this matrix
    
    int Save(const char * filename, int oneIndexed=0) const; // save matrix to a disk text file 
    
    int SaveToMatlabFormat(const char * filename) const; // save matrix to a text file that can be imported into Matlab
    
    // set/add value to the j-th sparse entry in the given row (NOT to matrix element at (row,j))
    inline void SetEntry(int row, int j, double value) { columnEntries[row][j] = value; }
    inline void AddEntry(int row, int j, double value) { columnEntries[row][j] += value; }
    void ResetToZero(); // reset all entries to zero
    void ResetRowToZero(int row); // reset all entries in the row to zero
    
    inline int Getn() const { return GetNumRows(); } // get the number of rows
    inline int GetNumRows() const { return (int)columnEntries.size(); }
    inline int GetRowLength(int row) const { return (int)columnEntries[row].size(); }
    int GetNumColumns() const; // get the number of columns (i.e., search for max column index)
    // returns the j-th sparse entry in row i (NOT matrix element at (row, j))
    inline double GetEntry(int row, int j) const { return columnEntries[row][j]; }
    // returns the column index of the j-th sparse entry in the given row
    inline int GetColumnIndex(int row, int j) const { return columnIndices[row][j]; } 
    inline const vector<vector<double> >& GetEntries() const { return columnEntries; }
    inline const vector<vector<int> >& GetColumnIndices() const { return columnIndices; }
    vector<int> GetRowLengths() const;
    
    // finds the compressed column index of element at location (row, jDense)
    // returns -1 if column not found
    int GetInverseIndex(int row, int jDense) const;
    
    int GetNumEntries() const; // returns the total number of non-zero entries
    double SumEntries() const; // returns the sum of all matrix entries
    void SumRowEntries(double * rowSums) const; // returns the sum of all entries in each row
    double GetMaxAbsEntry() const; // max abs value of a matrix entry
    double GetInfinityNorm() const; // matrix infinity norm
    void Print() const; // prints the matrix out to standard output
    void PrintSparse() const;
    void PrintPartial(int startRow, int startDenseColumn, int endRow, int endDenseColumn);
    
    double GetRowNorm2(int row) const;
    
    // matrix algebra (all involved matrices must have the same pattern of non-zero entries)
    SparseMatrix operator+(const SparseMatrix & mat2) const;
    SparseMatrix operator-(const SparseMatrix & mat2) const;
    friend SparseMatrix operator* (const double alpha, const SparseMatrix & mat2); // warning: this function makes a local copy; "ScalarMultiply" is more efficient
    SparseMatrix & operator=(const SparseMatrix & source); // matrices must have same size and locations of non-zero entries
    SparseMatrix & operator*=(const double alpha);
    SparseMatrix & operator+=(const SparseMatrix & mat2);
    SparseMatrix & operator-=(const SparseMatrix & mat2);
    void ScalarMultiply(const double alpha, SparseMatrix * dest=NULL); // dest = alpha * dest (if dest=NULL, operation is applied to this object)
    void ScalarMultiplyAdd(const double alpha, SparseMatrix * dest=NULL); // dest += alpha * dest (if dest=NULL, operation is applied to this object)
    void MultiplyRow(int row, double scalar); // multiplies all elements in row 'row' with scalar 'scalar'
    
    // multiplies the sparse matrix with the given vector/matrix
    void MultiplyVector(const double * vector, double * result) const; // result = A * vector
    void MultiplyVectorAdd(const double * vector, double * result) const; // result += A * vector
    void MultiplyVector(int startRow, int endRow, const double * vector, double * result) const; // result = A(startRow:endRow-1,:) * vector
    void TransposeMultiplyVector(const double * vector, int resultLength, double * result) const; // result = trans(A) * vector
    void TransposeMultiplyVectorAdd(const double * vector, double * result) const; // result += trans(A) * vector
    void MultiplyMatrix(int numDenseRows, int numDenseColumns, const double * denseMatrix, double * result) const; // result = A * denseMatrix (denseMatrix is a numDenseRows x numDenseColumns dense matrix, result is a numRows x numDenseColumns dense matrix)
    void MultiplyMatrixAdd(int numDenseRows, int numDenseColumns, const double * denseMatrix, double * result) const; // result += A * denseMatrix (denseMatrix is a numDenseRows x numDenseColumns dense matrix, result is a numDenseRows x numDenseColumns dense matrix)
    void MultiplyMatrixTranspose(int numDenseColumns, const double * denseMatrix, double * result) const; // result = A * trans(denseMatrix) (trans(denseMatrix) is a dense matrix with 'numDenseColumns' columns, result is a numRows x numDenseColumns dense matrix)
    
    // computes <M * vector, vector> (assumes symmetric M)
    double QuadraticForm(const double * vector) const;
    // normalizes vector in the M-norm: vector := vector / sqrt(<M * vector, vector>)  (assumes symmetric M)
    void NormalizeVector(double * vector) const;
    void ConjugateMatrix(double * U, int r, double * MTilde); // computes MTilde = U^T M U (M can be a general square matrix, U need not be a square matrix; number of columns of U is r; sizes of M and U must be such that product is defined; output matrix will have size r x r, stored column-major)
    SparseMatrix ConjugateMatrix(SparseMatrix & U, int verbose=0); // computes U^T M U (M is this matrix, and can be a general square matrix, U need not be a square matrix; sizes of M and U must be such that product is defined)
    
    // builds indices for subsequent faster product computation (below)
    // input: U, MTilde; MTilde must equal U^T M U, computed using the "ConjugateMatrix" routine above
    // output: precomputedIndices
    typedef int *** precomputedIndicesType;
    void BuildConjugationIndices(SparseMatrix & U, SparseMatrix & MTilde, precomputedIndicesType * precomputedIndices); // note: must be debugged
    // input: precomputedIndices, U
    // output: MTilde
    void ConjugateMatrix(precomputedIndicesType precomputedIndices, SparseMatrix & U, SparseMatrix & MTilde); // note: must be debugged
    
    // writes all entries into the space provided by 'data'
    // space must be pre-allocated
    // data is written row after row, and by non-zero columns within each row
    void MakeLinearDataArray(double * data) const;
    // writes row indices of non-zero entries into array "indices"
    // same order as for data
    void MakeLinearRowIndexArray(int * indices) const;
    // indices in this function version are double to ensure compatibility with Matlab
    void MakeLinearRowIndexArray(double * indices) const;
    // writes column indices
    void MakeLinearColumnIndexArray(int * indices) const;
    void MakeLinearColumnIndexArray(double * indices) const;
    
    // make a dense matrix (column-major LAPACK style storage)
    // (this can be a huge matrix for large sparse matrices)
    // storage in denseMatrix must be pre-allocated
    void MakeDenseMatrix(double * denseMatrix) const;
    // also transposes the matrix:
    void MakeDenseMatrixTranspose(int numColumns, double * denseMatrix) const;
    
    // removes row(s) and column(s) from the matrix
    void RemoveRowColumn(int rowColumn); // 0-indexed
    void RemoveRowsColumns(int numRemovedRowColumns, int * removedRowColumns, int oneIndexed=0); // the rowColumns must be sorted (ascending)
    void RemoveRowsColumnsSlow(int numRemovedRowColumns, int * removedRowColumns, int oneIndexed=0); // the rowColumns need not be sorted
    
    // removes row(s) from the matrix
    void RemoveRow(int row); // 0-indexed
    void RemoveRows(int numRemovedRows, int * removedRows, int oneIndexed=0); // rows must be sorted (ascending)
    void RemoveRowsSlow(int numRemovedRows, int * removedRows, int oneIndexed=0); // the rows need not be sorted
    
    // removes column(s) from the matrix
    void RemoveColumn(int column); // 0-indexed
    void RemoveColumns(int numRemovedColumns, int * removedColumns, int oneIndexed=0); // columns must be sorted (ascending)
    void RemoveColumnsSlow(int numRemovedColumns, int * removedColumns, int oneIndexed=0); // columns need not be sorted 
    
    void IncreaseNumRows(int numAddedRows); // increases the number of matrix rows (new rows are added at the bottom of the matrix, and are all empty)
    void SetRows(SparseMatrix * source, int startRow, int startColumn=0); // starting with startRow, overwrites the rows with those of matrix "source"; data is written into columns starting at startColumn
    void AppendRowsColumns(SparseMatrix * source); // appends the matrix "source" at the bottom of matrix, and trans(source) to the right of the matrix
    
    void Append(SparseMatrix *source); // appends source at the bottom right of this, where the start of 'bottom right' is defined by (GetNumRows(),GetNumRows())
    int InsertNewEntry(int row, int denseColumn); // insert a new entry at (row,denseColumn), return sparseColumn
    
    void CreateEntriesIfNecessary(const SparseMatrixOutline& outline, unsigned int rowsColumnsOffset);
    
    // transposition (note: the matrix need not be symmetric)
    void BuildTranspositionIndices();
    void FreeTranspositionIndices();
    // returns the list position of the transposed element (row, list position j)
    inline int TransposedIndex(int row, int j) const { return transposedIndices[row][j]; }
    
    // returns the transposed matrix
    // numColumns is the number of columns in the original matrix;
    // this is important in case there are zero columns at the end of the matrix
    // if numColumns=-1 (default), GetNumColumns() will be called; however, this will lead to a transposed matrix with a fewer number of rows in case of empty columns at the end of the original matrix
    SparseMatrix * Transpose(int numColumns=-1);
    
    // checks if the matrix is skew-symmetric
    // the non-zero entry locations must form a symmetric pattern
    // returns max ( abs ( A^T + A ) ) = || A^T + A ||_{\infty}
    double SkewSymmetricCheck(); 
    // makes matrix symmetric by copying upper triangle + diagonal into the lower triangle
    // the non-zero entry locations must form a symmetric pattern
    void SymmetrizeMatrix();
    
    // pre-computes the sparse columns of diagonal matrix entries
    // this routine will accelerate subsequent GetDiagonal or AddDiagonalMatrix calls, but is not necessary for GetDiagonal or AddDiagonalMatrix
    void BuildDiagonalIndices();
    void FreeDiagonalIndices();
    void GetDiagonal(double * diagonal);
    void AddDiagonalMatrix(double * diagonalMatrix);
    void AddDiagonalMatrix(double constDiagonalElement);
    
    
    
    // Build submatrix indices is used for pair of matrices where the sparsity of one matrix is a subset of another matrix (for example, mass matrix and stiffness matrix).
    // The two matrices need not have the same number of entries but the subMatrix should be contained within the bounds of the super matrix after adding denseRowColumnOffset to
    // row and dense column indices in the sub matrix.
    
    // Call this once to establish the correspondence:
    shared_ptr<SparseSubMatrixLinkage> AttachSubMatrix(shared_ptr<SparseMatrix> subMatrix, int denseRowColumnOffset=0);
    void DetachSubMatrix(shared_ptr<SparseSubMatrixLinkage> linkage);
    shared_ptr<SparseSubMatrixLinkage> GetExistingSubMatrixLinkage(shared_ptr<SparseMatrix> subMatrix);
    

    //void BuildSubMatrixIndices(SparseMatrix & submatrix, int subMatrixID=0);
    //void FreeSubMatrixIndices(int subMatrixID=0);
    // add a matrix to the current matrix, whose elements are a subset of the elements of the current matrix
    // += factor * mat2 
    void AddFromSubMatrix(double factor, shared_ptr<SparseSubMatrixLinkage> subMatrixLinkage);
    // subMatrix must have been attached via AttachSubMatrix
    void AddFromSubMatrix(double factor, shared_ptr<SparseMatrix> subMatrix);
    
    /*
    // Build supermatrix indices is used for pair of matrices with rows/columns removed.
    // oneIndexed: tells whether the fixed rows and columns are specified 1-indexed or 0-indexed
    // First, call BuildSuperMatrixIndices once to inialize (all fixed rows and columns are indexed with respect the superMatrix):
    void BuildSuperMatrixIndices(int numFixedRowColumns, int * fixedRowColumns, SparseMatrix * superMatrix, int oneIndexed=0); // use this version if the indices of removed rows and columns are the same
    void BuildSuperMatrixIndices(int numFixedRows, int * fixedRows, int numFixedColumns, int * fixedColumns, SparseMatrix * superMatrix, int oneIndexed=0); // allows arbitrary row and column indices*/
    // Then, call this (potentially many times) to quickly assign the values at the appropriate places in the submatrix.
    // For example, you can use this to copy data from a matrix into a submatrix obtained by a previous call to RemoveRowColumns.
    
    shared_ptr<SparseSuperMatrixLinkage> AttachSuperMatrix(shared_ptr<SparseMatrix> superMatrix);
    
    void AssignFromSuperMatrix(std::shared_ptr<SparseMatrix> superMatrix);

    
    // returns the total number of non-zero entries in the lower triangle (including diagonal)
    int GetNumLowerTriangleEntries() const;
    int GetNumUpperTriangleEntries() const;
    // exports the matrix to format for NAG library
    int GenerateNAGFormat(double * a,int * irow,int * icol, int * istr) const;
    
    void GenerateCompressedRowMajorFormat(double * a, int * ia, int * ja, int upperTriangleOnly=0, int oneIndexed=0) const; 
    void GenerateCompressedRowMajorFormat_four_array(double * values, int * columns, int * pointerB, int * pointerE, int upperTriangleOnly=0, int oneIndexed=0) const; 
    
    // diagonal solve M * x = b
    // ASSUMES the sparse matrix is diagonal !
    // result is overwritten into rhs
    // (to solve non-diagonal linear systems, you need to use an external library; or you can use the CGSolver class, or you can use the Gauss-Seidel iteration below)
    void DiagonalSolve(double * rhs) const;
    
    // performs one Gauss-Seidel iteration of solving the system A * x = b
    // updates vector x in place, b is not modified
    // (A can be a general matrix)
    // assumes that diagonal entries of the matrix are set and are non-zero
    void DoOneGaussSeidelIteration(double * x, const double * b) const;
    void ComputeResidual(const double * x, const double * b, double * residual) const;
    // checks if A * x - b is close to zero and prints out the findings
    // returns ||A * x - b|| / ||b|| (all norms are infinity)
    // passing a buffer (length of n) will avoid a malloc/free pair to generate scratch space for the residual
    double CheckLinearSystemSolution(const double * x, const double * b, int verbose=1, double * buffer=NULL) const;
    
    // below are low-level routines which are rarely used
    inline std::vector<std::vector<double> >& GetDataHandle() { return columnEntries; }
    inline const std::vector<double>& GetRowHandle(int row) const { return columnEntries[row]; }
    
    // create a nxn identity matrix
    static SparseMatrix * CreateIdentityMatrix(int n);
    
protected:
    
    // compressed row storage
    vector<vector<int> > columnIndices; // indices of columns of non-zero entries in each row
    vector<vector<double> > columnEntries; // values of non-zero entries in each row
    
    bool HasCachedDiagonalIndices() { return diagonalIndices.size() > 0; }
    vector<int> diagonalIndices;
    bool HasCachedTransposedIndices() { return transposedIndices.size() > 0; }
    vector<vector<int> > transposedIndices;
    

    void AttachSubMatrix(shared_ptr<SparseSubMatrixLinkage> linkage);
    vector<shared_ptr<SparseSubMatrixLinkage> > subMatrixLinkages;
    
    shared_ptr<SparseSuperMatrixLinkage> superMatrixLinkage;
    
    void InitFromOutline(const SparseMatrixOutline& sparseMatrixOutline);
    void Allocate(size_t numRows);
};

#endif

