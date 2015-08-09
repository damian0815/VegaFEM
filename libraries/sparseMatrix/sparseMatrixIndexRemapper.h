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
    /// for each entry in sourceMatrix, map to the matching entry in targetMatrix, taking rowColumnOffset into account.
    SparseMatrixIndexRemapper(shared_ptr<SparseMatrix> superMatrix, shared_ptr<SparseMatrix> subMatrix, int denseRowColumnOffset);
    
    inline int GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(int subMatrixRow, int subMatrixSparseColumn) const;
    inline int GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SuperMatrixRow(int superMatrixRow, int subMatrixSparseColumn) const;
    
    inline bool HasSubMatrixRowForSuperMatrixRow(int superMatrixRow) const { return superMatrixToSubMatrixRowMap.count(superMatrixRow); }
    inline int GetSubMatrixRowForSuperMatrixRow(int superMatrixRow) const { return superMatrixToSubMatrixRowMap.at(superMatrixRow); }
    int GetSuperMatrixRowForSubMatrixRow(int subMatrixRow) const;
    
    void RemoveSuperRowFromSubMatrix(int whichSuperMatrixRow);
    void RemoveSuperColumnFromSubMatrix(int whichSuperMatrixDenseColumn);
    
    void AssignSubMatrixFromSuperMatrix();
    
    void Print();
    
    shared_ptr<SparseMatrix> GetSubMatrix() { return subMatrix; }
    shared_ptr<SparseMatrix> GetSuperMatrix() { return superMatrix; }
    
private:
    void SetupMapping();
    void PopulateRowMap();
    void PopulateColumnIndexMaps();
    
    map<int, int> superMatrixToSubMatrixRowMap;
    
    // A list of the super matrix sparse indices, one for each of the sub matrix sparse indices
    typedef vector<int> SubMatrixSparseColumnToSuperMatrixSparseColumnMap;
    // One subMatirxSparseColumnToSuperMatrixSparseColumnMap for each row in the sub matrix
    vector<SubMatrixSparseColumnToSuperMatrixSparseColumnMap> subMatrixSparseToSuperMatrixSparseColumnMaps;
    
    shared_ptr<SparseMatrix> superMatrix;
    shared_ptr<SparseMatrix> subMatrix;
    
    int denseRowColumnOffset;
    
};

int SparseMatrixIndexRemapper::GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(int subMatrixRow, int subMatrixSparseColumn) const
{
    const auto& subMatrixToSuperMatrixSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow];
    return subMatrixToSuperMatrixSparseColumnMap.at(subMatrixSparseColumn);
}

int SparseMatrixIndexRemapper::GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SuperMatrixRow(int superMatrixRow, int subMatrixSparseColumn) const
{
    int subMatrixRow = GetSubMatrixRowForSuperMatrixRow(superMatrixRow);
    return GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(subMatrixRow, subMatrixSparseColumn);
}

#endif /* defined(__VegaFEM__sparseMatrixIndexRemapper__) */
