//
//  sparseMatrixIndexRemapper.h
//  VegaFEM
//
//  Created by Damian Stewart on 21/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#ifndef __VegaFEM__sparseMatrixIndexRemapper__
#define __VegaFEM__sparseMatrixIndexRemapper__

#include <set>
#include <iostream>
#include <map>
#include <vector>

class SparseMatrix;

using std::vector;
using std::map;
using std::set;
using std::shared_ptr;

class SparseMatrixIndexRemapper
{
public:
    /// for each entry in sourceMatrix, map to the matching entry in targetMatrix, taking rowColumnOffset into account.
    SparseMatrixIndexRemapper(shared_ptr<SparseMatrix> superMatrix, shared_ptr<SparseMatrix> subMatrix, int denseRowColumnOffset);
       
    void AssignSubMatrixFromSuperMatrix();
    void AddSubMatrixToSuperMatrix(double factor=1.0);
    
    void Print();
    
    inline int GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(int subMatrixRow, int subMatrixSparseColumn) const;
    inline int GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SuperMatrixRow(int superMatrixRow, int subMatrixSparseColumn) const;
    
    inline bool HasSubMatrixRowForSuperMatrixRow(int superMatrixRow) const { return superMatrixToSubMatrixRowMap.count(superMatrixRow); }
    inline int GetSubMatrixRowForSuperMatrixRow(int superMatrixRow) const { return superMatrixToSubMatrixRowMap.at(superMatrixRow); }
    bool HasSuperMatrixRowForSubMatrixRow(int subMatrixRow) const;
    int GetSuperMatrixRowForSubMatrixRow(int subMatrixRow) const;
    
    bool HasSubMatrixSparseColumnForSuperMatrixSparseColumn(int superMatrixRow, int superMatrixSparseColumn);
    int GetSubMatrixSparseColumnForSuperMatrixSparseColumn(int superMatrixRow, int superMatrixSparseColumn);
    
    void RemoveSuperRowFromSubMatrix(int whichSuperMatrixRow);
    void RemoveSuperColumnFromSubMatrix(int whichSuperMatrixDenseColumn);
    
    bool IsSuperRowRemovedFromSubMatrix(int whichSuperRow) const { return removedSuperMatrixRows.count(whichSuperRow) != 0;  }
    bool IsSuperColumnRemovedFromSubMatrix(int whichSuperDenseColumn) const { return removedSuperMatrixDenseColumns.count(whichSuperDenseColumn) != 0; }
    
    void OnEntryWasInsertedIntoSuperMatrix(int superMatrixRow, int insertedSuperMatrixDenseColumn);
    void OnEntryWasInsertedIntoSubMatrix(int subMatrixRow, int insertedSubMatrixDenseColumn);
    
    void AddEntryToMap(int superRow, int superDenseColumn);
    
    shared_ptr<SparseMatrix> GetSubMatrix() { return subMatrix; }
    shared_ptr<SparseMatrix> GetSuperMatrix() { return superMatrix; }
    
private:
    void SetupMapping();
    void PopulateRowMap();
    void PopulateColumnIndexMaps();
    
    map<int, int> superMatrixToSubMatrixRowMap;
    
    // A list of the super matrix sparse indices, one for each of the sub matrix sparse indices
    typedef vector<int> SubMatrixSparseColumnToSuperMatrixSparseColumnMap;
    // One subMatrixSparseColumnToSuperMatrixSparseColumnMap for each row in the sub matrix
    vector<SubMatrixSparseColumnToSuperMatrixSparseColumnMap> subMatrixSparseToSuperMatrixSparseColumnMaps;
    
    shared_ptr<SparseMatrix> superMatrix;
    shared_ptr<SparseMatrix> subMatrix;
    
    set<int> removedSuperMatrixRows;
    set<int> removedSuperMatrixDenseColumns;
    
    int denseRowColumnOffset;

    void PrepareAddSubMatrixRowCaches(
            const vector<int> &thisSubMatrixToSuperMatrixSparseColumnMap, const vector<double> &thisSubColumnEntries,
            int subMatrixRowLength, int *sparseSuperJCache, double *subEntriesCache) const;

    void DoAddSubMatrixRow(double factor, const int *sparseSuperJCache, const double *subEntriesCache,
                           int subMatrixRowLength, vector<double> &thisSuperColumnEntries);
};

int SparseMatrixIndexRemapper::GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(int subMatrixRow, int subMatrixSparseColumn) const
{
    const auto& subMatrixToSuperMatrixSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow];
    return subMatrixToSuperMatrixSparseColumnMap[subMatrixSparseColumn];
}

int SparseMatrixIndexRemapper::GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SuperMatrixRow(int superMatrixRow, int subMatrixSparseColumn) const
{
    int subMatrixRow = GetSubMatrixRowForSuperMatrixRow(superMatrixRow);
    return GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(subMatrixRow, subMatrixSparseColumn);
}

#endif /* defined(__VegaFEM__sparseMatrixIndexRemapper__) */
