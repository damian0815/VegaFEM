//
//  sparseMatrixIndexRemapper.h
//  VegaFEM
//
//  Created by Damian Stewart on 21/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#ifndef __VegaFEM__sparseMatrixIndexRemapper__
#define __VegaFEM__sparseMatrixIndexRemapper__

#include <iostream>
#include <map>
#include <vector>

class SparseMatrix;

using std::vector;
using std::map;
using std::shared_ptr;

class SparseMatrixIndexRemapper
{
public:
    /// for each entry in sourceMatrix, map to the matching entry in targetMatrix.
    /// Matrices must have the same number of rows.
    SparseMatrixIndexRemapper(shared_ptr<SparseMatrix> superMatrix, shared_ptr<SparseMatrix> subMatrix);
    
    inline int GetSuperMatrixSparseColumnForSubMatrixSparseColumn(int subMatrixRow, int subMatrixSparseColumn) const;
    
    inline int GetSubMatrixRowForSuperMatrixRow(int superMatrixRow) const { return superMatrixToSubMatrixRowMap.at(superMatrixRow); }
    int GetSuperMatrixRowForSubMatrixRow(int subMatrixRow) const;
    
    void RemoveSuperRowFromSubMatrix(int whichSuperMatrixRow);
    void RemoveSuperColumnFromSubMatrix(int whichSuperMatrixDenseColumn);
    
    void AssignSubMatrixFromSuperMatrix();
    
    void Print();
    
private:
    void SetupMapping();
    void PopulateRowMap();
    void PopulateColumnIndexMaps();
    
     /*
     length(subMatrixSparseToSuperMatrixSparseColumnMaps) == number of rows
     */
    // maps from subMatrix sparse column indices to superMatrix sparse column indices
    
    
    map<int, int> superMatrixToSubMatrixRowMap;
    
    // A list of the super matrix sparse indices, one for each of the sub matrix sparse indices
    typedef vector<int> SubMatrixSparseColumnToSuperMatrixSparseColumnMap;
    // One subMatirxSparseColumnToSuperMatrixSparseColumnMap for each row in the sub matrix
    vector<SubMatrixSparseColumnToSuperMatrixSparseColumnMap> subMatrixSparseToSuperMatrixSparseColumnMaps;
    
    shared_ptr<SparseMatrix> superMatrix;
    shared_ptr<SparseMatrix> subMatrix;
    
};

int SparseMatrixIndexRemapper::GetSuperMatrixSparseColumnForSubMatrixSparseColumn(int superMatrixRow, int subMatrixSparseColumn) const
{
    int subMatrixRow = GetSubMatrixRowForSuperMatrixRow(superMatrixRow);
    const auto& subMatrixToSuperMatrixSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow];
    return subMatrixToSuperMatrixSparseColumnMap.at(subMatrixSparseColumn);
}

#endif /* defined(__VegaFEM__sparseMatrixIndexRemapper__) */
