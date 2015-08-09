//
//  sparseMatrixIndexRemapper.cpp
//  VegaFEM
//
//  Created by Damian Stewart on 21/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include <assert.h>
#include "sparseMatrixIndexRemapper.h"
#include "sparseMatrix.h"

SparseMatrixIndexRemapper::SparseMatrixIndexRemapper(shared_ptr<SparseMatrix> superMatrix_, shared_ptr<SparseMatrix> subMatrix_)
:   superMatrix(superMatrix_),
    subMatrix(subMatrix_)
{
    SetupMapping();
}


void SparseMatrixIndexRemapper::SetupMapping()
{
    assert(superMatrix->GetNumRows() == subMatrix->GetNumRows() && "row count must match");
    
    PopulateRowMap();
    PopulateColumnIndexMaps();
}

void SparseMatrixIndexRemapper::PopulateRowMap()
{
    superMatrixToSubMatrixRowMap.clear();
    for (int superMatrixRow=0; superMatrixRow<superMatrix->GetNumRows(); superMatrixRow++)
    {
        int subMatrixRow = superMatrixRow;
        superMatrixToSubMatrixRowMap[superMatrixRow] = subMatrixRow;
    }
}

void SparseMatrixIndexRemapper::PopulateColumnIndexMaps()
{
    subMatrixSparseToSuperMatrixSparseColumnMaps.clear();
    subMatrixSparseToSuperMatrixSparseColumnMaps.resize(superMatrix->GetNumRows());
    
    for(int row=0; row<superMatrix->GetNumRows(); row++)
    {
        int submatrixRowLength = subMatrix->GetRowLength(row);
        subMatrixSparseToSuperMatrixSparseColumnMaps[row].resize(submatrixRowLength);
        const vector<int>& subMatrixDenseColumnIndices = subMatrix->GetColumnIndices()[row];
        
        for(int sparseColumn=0; sparseColumn < submatrixRowLength; sparseColumn++)
        {
            // finds the position in row i of element with column index jDense
            // int GetInverseIndex(int i, int jDense);
            int denseColumn = subMatrixDenseColumnIndices[sparseColumn];
            int superMatrixSparseIndex = superMatrix->GetInverseIndex(row, denseColumn);
            if (superMatrixSparseIndex == -1)
            {
                printf("Error (BuildSubMatrixIndices): given matrix is not a submatrix of this matrix. The following index does not exist in this matrix: (%d,%d)\n", row, denseColumn);
                assert(false);
            }
            subMatrixSparseToSuperMatrixSparseColumnMaps[row][sparseColumn] = superMatrixSparseIndex;
        }
    }
}


void SparseMatrixIndexRemapper::RemoveSuperRowFromSubMatrix(int whichSuperMatrixRow)
{
    auto it = superMatrixToSubMatrixRowMap.find(whichSuperMatrixRow);
    assert(it != superMatrixToSubMatrixRowMap.end() && "super matrix row does not exist in the sub matrix");
    int whichSubMatrixRow = (*it).second;
    
    it = superMatrixToSubMatrixRowMap.erase(it);
    for (; it != superMatrixToSubMatrixRowMap.end(); ++it)
    {
        --(*it).second;
    }
    
    assert(whichSubMatrixRow < subMatrixSparseToSuperMatrixSparseColumnMaps.size() && "sub matrix row has no column entries - probable data corruption");
    subMatrixSparseToSuperMatrixSparseColumnMaps.erase(subMatrixSparseToSuperMatrixSparseColumnMaps.begin() + whichSubMatrixRow);
}


void SparseMatrixIndexRemapper::RemoveSuperColumnFromSubMatrix(int whichSuperMatrixDenseColumn)
{
    for (auto& kvp: superMatrixToSubMatrixRowMap)
    {
        int superMatrixRow = kvp.first;
        int subMatrixRow = kvp.second;
        int superMatrixSparseColumn = superMatrix->GetInverseIndex(superMatrixRow, whichSuperMatrixDenseColumn);
        if (superMatrixSparseColumn != -1)
        {
            auto& subMatrixToSuperMatrixSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow];
            auto it = std::find(subMatrixToSuperMatrixSparseColumnMap.begin(), subMatrixToSuperMatrixSparseColumnMap.end(), superMatrixSparseColumn);
            if (it != subMatrixToSuperMatrixSparseColumnMap.end())
            {
                //int subMatrixSparseIndex = *it;
                //printf("removing sub matrix column sparse:%i dense:%i in row %i\n", subMatrixSparseIndex, whichSuperMatrixDenseColumn, subMatrixRow);
                subMatrixToSuperMatrixSparseColumnMap.erase(it);
            }
            else
            {
                //printf("super matrix dense column %i is not in row %i\n", whichSuperMatrixDenseColumn, subMatrixRow);
            }
        }
        
    }
    
}


void SparseMatrixIndexRemapper::Print()
{
    printf("SparseMatrixIndexRemapper\n");
    printf("subMatrix Row index -> superMatrix row index: superMatrixColumnIndices\n");
    for (auto& kvp: superMatrixToSubMatrixRowMap)
    {
        int superMatrixRow = kvp.first;
        int subMatrixRow = kvp.second;
        printf("%2i -> %2i: ", subMatrixRow, superMatrixRow);
        for (int j=0; j<subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow].size(); j++)
        {
            auto superMatrixSparseColumn = subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow][j];
            printf("%2i ", superMatrixSparseColumn);
        }
        printf("\n");
    }
}

int SparseMatrixIndexRemapper::GetSuperMatrixRowForSubMatrixRow(int subMatrixRow) const
{
    auto it = std::find_if(superMatrixToSubMatrixRowMap.begin(), superMatrixToSubMatrixRowMap.end(), [subMatrixRow](const std::pair<int,int>& kvp){ return kvp.second == subMatrixRow; });
    return (*it).first;
}

void SparseMatrixIndexRemapper::AssignSubMatrixFromSuperMatrix()
{
    auto& subColumnEntries = subMatrix->GetDataHandle();
    const auto& superMatrixEntries = superMatrix->GetDataHandle();
    for(int subRow=0; subRow<subMatrix->GetNumRows(); subRow++)
    {
        int superRow = GetSuperMatrixRowForSubMatrixRow(subRow);
        const vector<double>& superColumnEntries = superMatrixEntries[superRow];
        const vector<int>& subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subRow];
        int subRowLength = subMatrix->GetRowLength(subRow);
        for(int subSparseColumn=0; subSparseColumn < subRowLength; subSparseColumn++)
        {
            int superSparseColumn = subSparseToSuperSparseColumnMap[subSparseColumn];
            subColumnEntries[subRow][subSparseColumn] = superColumnEntries[superSparseColumn];
        }
    }
}