//
//  sparseMatrixOutline.h
//  SparseMatrixImproved
//
//  Created by Damian Stewart on 18/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#ifndef __SparseMatrixImproved__sparseMatrixOutline__
#define __SparseMatrixImproved__sparseMatrixOutline__

#include <iostream>
#include <string.h>
#include <vector>
using std::vector;
#include <map>

class SparseMatrix;

class SparseMatrixOutline
{
public:
  // makes an empty sparse matrix with numRows rows
  SparseMatrixOutline(int numRows);
  ~SparseMatrixOutline();

  // makes a diagonal numRows x numRows sparse matrix; with a constant diagonal
  SparseMatrixOutline(int numRows, double diagonal);
  // makes a diagonal numRows x numRows sparse matrix; diagonal is a vector of n numbers
  SparseMatrixOutline(int numRows, double * diagonal);

  // loads the sparse matrix from a text file
  // if expand is greater than 1, the routine also expands each element into a diagonal block of size expand x expand... 
  //   (expand option is useful for loading the mass matrix in structural mechanics (with expand=3 in 3D))
  SparseMatrixOutline(const char * filename, int expand=1); 

  // save matrix to a text file
  int Save(const char * filename, int oneIndexed=0) const;

  // add entry at location (i,j) in the matrix
  void AddEntry(int i, int j, double value=0.0);
  void AddBlock3x3Entry(int i, int j, double * matrix3x3); // matrix3x3 should be given in row-major order
  // add a block (sparse) matrix (optionally multiplied with "scalarFactor"), starting at row i, and column j
  void AddBlockMatrix(int i, int j, const SparseMatrix * block, double scalarFactor=1.0);
  void IncreaseNumRows(int numAddedRows); // increases the number of matrix rows (new rows are added at the bottom of the matrix, and are all empty)

  void MultiplyRow(int row, double scalar); // multiplies all elements in row 'row' with scalar 'scalar'

  inline int Getn() const { return numRows; } // get number of rows
  inline int GetNumRows() const { return numRows; } // get number of rows
  int GetNumColumns() const; // get the number of columns (i.e., search for max column index)
  int GetNumEntries() const; // get total number of non-zero matrix elements
	bool HasEntry(int i, int j) const;
  double GetEntry(int i, int j) const; // returns the matrix entry at location (i,j) in the matrix (or zero if entry has not been assigned)
  void Print() const;
	
	
  void PrintTopology(int clusterSize=3) const; // clusterSize is the size of the square blocks that are assumed to compose the matrix; eg for 3x3 blocks (3d FEM matrix), use clusteSize=3

  // low-level routine which is rarely used
  inline const std::map<int,double> & GetRow(int i) const { return columnEntries[i]; }
  friend class SparseMatrix;

protected:
  int numRows;
  std::vector< std::map<int,double> > columnEntries;
  void Allocate();
};

#endif /* defined(__SparseMatrixImproved__sparseMatrixOutline__) */
