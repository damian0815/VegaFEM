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

SparseMatrixIndexRemapper::SparseMatrixIndexRemapper(shared_ptr<SparseMatrix> superMatrix_, shared_ptr<SparseMatrix> subMatrix_, int denseRowColumnOffset_)
:   superMatrix(superMatrix_),
    subMatrix(subMatrix_),
    denseRowColumnOffset(denseRowColumnOffset_)
{
    SetupMapping();
}


void SparseMatrixIndexRemapper::SetupMapping()
{
    assert(superMatrix->GetNumRows() >= denseRowColumnOffset + subMatrix->GetNumRows() && "not enough rows");
    
    PopulateRowMap();
    PopulateColumnIndexMaps();
}

void SparseMatrixIndexRemapper::PopulateRowMap()
{
    superMatrixToSubMatrixRowMap.clear();
    for (int subMatrixRow=0; subMatrixRow<subMatrix->GetNumRows(); subMatrixRow++) {
        auto superMatrixRow = subMatrixRow + denseRowColumnOffset;
        superMatrixToSubMatrixRowMap[superMatrixRow] = subMatrixRow;
    }
}

void SparseMatrixIndexRemapper::PopulateColumnIndexMaps()
{
    subMatrixSparseToSuperMatrixSparseColumnMaps.clear();
    subMatrixSparseToSuperMatrixSparseColumnMaps.resize(subMatrix->GetNumRows());
    
    for(int subMatrixRow=0; subMatrixRow<subMatrix->GetNumRows(); subMatrixRow++)
    {
        int submatrixRowLength = subMatrix->GetRowLength(subMatrixRow);
        subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow].resize(submatrixRowLength);
        const vector<int>& subMatrixDenseColumnIndices = subMatrix->GetColumnIndices()[subMatrixRow];
        
        int superMatrixRow = subMatrixRow + denseRowColumnOffset;
        
        for(int subMatrixSparseColumn=0; subMatrixSparseColumn < submatrixRowLength; subMatrixSparseColumn++)
        {
            // finds the position in row i of element with column index jDense
            // int GetInverseIndex(int i, int jDense);
            int subMatrixDenseColumn = subMatrixDenseColumnIndices[subMatrixSparseColumn];
            int superMatrixDenseColumn = subMatrixDenseColumn + denseRowColumnOffset;
            int superMatrixSparseIndex = superMatrix->GetInverseIndex(superMatrixRow, superMatrixDenseColumn);
            if (superMatrixSparseIndex == -1)
            {
                printf("Error (BuildSubMatrixIndices): given matrix is not a submatrix of this matrix. The following index does not exist in this matrix: (%d,%d)\n", superMatrixRow, superMatrixDenseColumn);
                assert(false);
            }
            subMatrixSparseToSuperMatrixSparseColumnMaps[subMatrixRow][subMatrixSparseColumn] = superMatrixSparseIndex;
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

static map<int,int>::const_iterator GetRowMapIterator(const map<int,int>& rowMap, int subMatrixRow)
{
    auto it = std::find_if(rowMap.begin(), rowMap.end(), [subMatrixRow](const std::pair<int,int>& kvp){ return kvp.second == subMatrixRow; });
    return it;
}

int SparseMatrixIndexRemapper::GetSuperMatrixRowForSubMatrixRow(int subMatrixRow) const
{
    auto it = GetRowMapIterator(superMatrixToSubMatrixRowMap, subMatrixRow);
    return (*it).first;
}

bool SparseMatrixIndexRemapper::HasSuperMatrixRowForSubMatrixRow(int subMatrixRow) const
{
    auto it = GetRowMapIterator(superMatrixToSubMatrixRowMap, subMatrixRow);
    return (it != superMatrixToSubMatrixRowMap.end());
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
        for(int subSparseColumn=0; subSparseColumn < subRowLength; subSparseColumn++) {
            int superSparseColumn = subSparseToSuperSparseColumnMap[subSparseColumn];
            subColumnEntries[subRow][subSparseColumn] = superColumnEntries[superSparseColumn];
        }
    }
}

void SparseMatrixIndexRemapper::OnEntryWasInsertedIntoSuperMatrix(int superRow, int insertedSuperDenseColumn)
{
    if (HasSubMatrixRowForSuperMatrixRow(superRow)) {

        int subRow = GetSubMatrixRowForSuperMatrixRow(superRow);
        vector<int>& subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps.at(subRow);
        
        int insertedSuperSparseColumn = superMatrix->GetInverseIndex(superRow, insertedSuperDenseColumn);
        assert(insertedSuperSparseColumn != -1 && "sparse column not found, something is broken");
        for (int i=0; i<subSparseToSuperSparseColumnMap.size(); i++) {
            //int subSparseColumn = i;
            int& superSparseColumn = subSparseToSuperSparseColumnMap[i];
            if (superSparseColumn >= insertedSuperSparseColumn) {
                ++superSparseColumn;
            }
        }
    }
}

void SparseMatrixIndexRemapper::OnEntryWasInsertedIntoSubMatrix(int subRow, int insertedSubDenseColumn)
{
    /*
    if (HasSuperMatrixRowForSubMatrixRow(subRow)) {
        int superRow = GetSuperMatrixRowForSubMatrixRow(subRow);
        vector<int> subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps.at(subRow);
        
        
        for (int i=0; i<subSparseToSuperSparseColumnMap.size(); i++) {
            int subSparseColumn = i;
            if (subSparseColumn )
    }*/
}

bool SparseMatrixIndexRemapper::HasSubMatrixSparseColumnForSuperMatrixSparseColumn(int superRow, int superSparseColumn)
{
	return GetSubMatrixSparseColumnForSuperMatrixSparseColumn(superRow, superSparseColumn) != -1;
}

int SparseMatrixIndexRemapper::GetSubMatrixSparseColumnForSuperMatrixSparseColumn(int superRow, int superSparseColumn)
{
	int subRow = GetSubMatrixRowForSuperMatrixRow(superRow);
    const vector<int>& subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps.at(subRow);
	for (int mapIndex=0; mapIndex<subSparseToSuperSparseColumnMap.size(); mapIndex++) {
		if (subSparseToSuperSparseColumnMap[mapIndex]==superSparseColumn) {
			return mapIndex;
		}
	}
	
	return -1;
	
}

