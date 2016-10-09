//
//  sparseMatrixIndexRemapper.cpp
//  VegaFEM
//
//  Created by Damian Stewart on 21/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include <set>
#include <assert.h>
#include <thread>
#include "sparseMatrixIndexRemapper.h"
#include "sparseMatrix.h"
using std::pair;

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
    superMatrixToSubMatrixRowMap.resize(subMatrix->GetNumRows() + denseRowColumnOffset);
    std::fill(superMatrixToSubMatrixRowMap.begin(), superMatrixToSubMatrixRowMap.end(), -1);

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
    auto it = std::find(superMatrixToSubMatrixRowMap.begin(), superMatrixToSubMatrixRowMap.end(), whichSuperMatrixRow);
    assert(it != superMatrixToSubMatrixRowMap.end() && "super matrix row does not exist in the sub matrix");
    assert((*it) >= 0 && "row has already been removed");
    int whichSubMatrixRow = *it;

    *it = -1;
    ++it;
    for (; it != superMatrixToSubMatrixRowMap.end(); ++it)
    {
        int& subMatrixRow = *it;
        if (subMatrixRow != -1) {
            assert(subMatrixRow >= 0);
            --subMatrixRow;
        }
    }
    
    assert(whichSubMatrixRow < subMatrixSparseToSuperMatrixSparseColumnMaps.size() && "sub matrix row has no column entries - probable data corruption");
    subMatrixSparseToSuperMatrixSparseColumnMaps.erase(subMatrixSparseToSuperMatrixSparseColumnMaps.begin() + whichSubMatrixRow);
    
    removedSuperMatrixRows.insert(whichSuperMatrixRow);
}


void SparseMatrixIndexRemapper::RemoveSuperColumnFromSubMatrix(int whichSuperMatrixDenseColumn)
{
    for (int superMatrixRow=0; superMatrixRow<superMatrixToSubMatrixRowMap.size(); ++superMatrixRow)
    {
        int subMatrixRow = superMatrixToSubMatrixRowMap[superMatrixRow];
        if (subMatrixRow == -1) {
            continue;
        }

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
    
    removedSuperMatrixDenseColumns.insert(whichSuperMatrixDenseColumn);
}


void SparseMatrixIndexRemapper::Print()
{
    for (int superRow=0; superRow<superMatrixToSubMatrixRowMap.size(); ++superRow) {
        int subRow = superMatrixToSubMatrixRowMap[superRow];
        if (subRow == -1) {
            continue;
        }

        const auto& subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps.at(subRow);
        if (subSparseToSuperSparseColumnMap.size()>0) {
            printf("sub row %2i -> super row %2i: ", subRow, superRow);
            for (int subSparseColumn=0; subSparseColumn<subSparseToSuperSparseColumnMap.size(); subSparseColumn++)
            {
                auto superSparseColumn = subSparseToSuperSparseColumnMap[subSparseColumn];
                auto subDenseColumn = subMatrix->GetColumnIndex(subRow, subSparseColumn);
                auto superDenseColumn = superMatrix->GetColumnIndex(superRow, superSparseColumn);
                printf("%2i:%2i, ", subDenseColumn, superDenseColumn);
            }
            printf("\n");
        }
    }
}

int SparseMatrixIndexRemapper::GetSuperMatrixRowForSubMatrixRow(int subMatrixRow) const
{
    auto it = std::find(superMatrixToSubMatrixRowMap.begin(), superMatrixToSubMatrixRowMap.end(), subMatrixRow);
    assert(it != superMatrixToSubMatrixRowMap.end());
    return (it - superMatrixToSubMatrixRowMap.begin());
}

bool SparseMatrixIndexRemapper::HasSuperMatrixRowForSubMatrixRow(int subMatrixRow) const
{
    auto it = std::find(superMatrixToSubMatrixRowMap.begin(), superMatrixToSubMatrixRowMap.end(), subMatrixRow);
    return (it != superMatrixToSubMatrixRowMap.end());
}

void SparseMatrixIndexRemapper::AssignSubMatrixFromSuperMatrix()
{
    const auto& superMatrixEntries = superMatrix->GetDataHandle();
    auto& subMatrixEntries = subMatrix->GetDataHandle();

    for (int superRow = 0; superRow<superMatrixToSubMatrixRowMap.size(); ++superRow)
    {
        int subRow = superMatrixToSubMatrixRowMap[superRow];
        if (subRow == -1) {
            continue;
        }
        //printf("assigning sub row %i from super row %i\n", subRow, superRow);
        const vector<double>& superColumnEntries = superMatrixEntries[superRow];
        vector<double>& subColumnEntries = subMatrixEntries[subRow];
        
        const vector<int>& subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subRow];
        int subRowLength = subMatrix->GetRowLength(subRow);
        for(int subSparseColumn=0; subSparseColumn < subRowLength; subSparseColumn++) {
            int superSparseColumn = subSparseToSuperSparseColumnMap[subSparseColumn];
            subColumnEntries[subSparseColumn] = superColumnEntries[superSparseColumn];
        }
    }
    /*
    
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
    }*/
}

static map<SparseMatrixIndexRemapper*, pair<std::chrono::nanoseconds, int>> doItDurations;

void SparseMatrixIndexRemapper::AddSubMatrixToSuperMatrix(double factor) {
    auto &superColumnEntries = superMatrix->GetDataHandle();
    const auto &subColumnEntries = subMatrix->GetDataHandle();

    const auto subRowLengths = subMatrix->GetRowLengths();

    const auto superRowLengths = superMatrix->GetRowLengths();
    int maxSuperRowLength = *(std::max_element(superRowLengths.begin(), superRowLengths.end()));
    double superEntriesCache[maxSuperRowLength];




    std::chrono::nanoseconds doItDuration(0);

    auto loop = [&](int startRow, int rowCount, bool storeDoIt) {
        for (int subRow = startRow; subRow < startRow + rowCount; subRow++) {
        //for (int subRow = 0; subRow < subMatrix->GetNumRows(); subRow++) {
            int superRow = GetSuperMatrixRowForSubMatrixRow(subRow);
            int subMatrixRowLength = subRowLengths[subRow];
            auto &thisSuperColumnEntries = superColumnEntries[superRow];
            const auto &thisSubColumnEntries = subColumnEntries[subRow];
            const auto &thisSubMatrixToSuperMatrixSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps[subRow];

            //auto doItStart = std::chrono::high_resolution_clock::now();


            for (int sparseSubJ = 0; sparseSubJ < subMatrixRowLength; sparseSubJ++) {
                int sparseSuperJ = thisSubMatrixToSuperMatrixSparseColumnMap[sparseSubJ];
                thisSuperColumnEntries[sparseSuperJ] += factor * thisSubColumnEntries[sparseSubJ];
            }

/*

        std::fill(superEntriesCache, superEntriesCache+maxSuperRowLength, 0);
        int superMatrixRowLength = superRowLengths[superRow];
        for(int sparseSubJ=0; sparseSubJ < subMatrixRowLength; sparseSubJ++) {
            int sparseSuperJ = thisSubMatrixToSuperMatrixSparseColumnMap[sparseSubJ];
            superEntriesCache[sparseSuperJ] = thisSubColumnEntries[sparseSubJ] * factor;
        }

        for (int sparseSuperJ=0; sparseSuperJ < superMatrixRowLength; sparseSuperJ++) {
            thisSuperColumnEntries[sparseSuperJ] += superEntriesCache[sparseSuperJ];
        }
*/

            //auto doItEnd = std::chrono::high_resolution_clock::now();
            //if (storeDoIt) {
            //    doItDuration += doItEnd - doItStart;
            //}
        }
    };

    const int numThreads = 1;
    if (numThreads == 1) {
        loop(0, subMatrix->GetNumRows(), true);
    } else {
        vector<std::thread> threads;
        int numRows = subMatrix->GetNumRows();
        for (int i = 0; i < numThreads; i++) {
            int rowCount = numRows / numThreads;
            int rowStart = i * rowCount;
            int actualRowCount = (i == (numThreads - 1) ? numRows - rowStart : rowCount);
            threads.push_back(std::thread(loop, rowStart, actualRowCount, i == 0));
        }

        for (auto& thread: threads) {
            thread.join();
        }
    }

    doItDurations[this].first += doItDuration;
    ++doItDurations[this].second;

    if ((doItDurations[this].second % 300) == 299) {
        double durationMillis = 1e-6 * doItDurations[this].first.count();
        //std::cout << "avg doItDuration (ms) " << this << ": " << (durationMillis/doItDurations[this].second) << std::endl;
        doItDurations[this] = std::make_pair(std::chrono::nanoseconds(0), 0);
    }
}

void SparseMatrixIndexRemapper::PrepareAddSubMatrixRowCaches(
        const vector<int> &thisSubMatrixToSuperMatrixSparseColumnMap, const vector<double> &thisSubColumnEntries,
        int subMatrixRowLength, int *sparseSuperJCache, double *subEntriesCache, float multiplicationFactor) const
{
    for(int sparseSubJ=0; sparseSubJ < subMatrixRowLength; sparseSubJ++) {
            int sparseSuperJ = thisSubMatrixToSuperMatrixSparseColumnMap[sparseSubJ];
            sparseSuperJCache[sparseSubJ] = sparseSuperJ;
            subEntriesCache[sparseSubJ] = multiplicationFactor * thisSubColumnEntries[sparseSubJ];
        }
}

void SparseMatrixIndexRemapper::DoAddSubMatrixRow(const int *sparseSuperJCache, const double *subEntriesCache,
                                                  int subMatrixRowLength, vector<double> &thisSuperColumnEntries) {
    for (int sparseSubJ=0; sparseSubJ < subMatrixRowLength; sparseSubJ++) {
            int sparseSuperJ = sparseSuperJCache[sparseSubJ];
            thisSuperColumnEntries[sparseSuperJ] += subEntriesCache[sparseSubJ];
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


void SparseMatrixIndexRemapper::AddEntryToMap(int superRow, int superDenseColumnToAdd)
{
    assert(HasSubMatrixRowForSuperMatrixRow(superRow) && "Must already have a row in sub matrix for this row in super matrix");
    
    auto superSparseColumnToAdd = superMatrix->GetInverseIndex(superRow, superDenseColumnToAdd);
    assert(superSparseColumnToAdd != -1 && "super entry does not exist");
    assert(!HasSubMatrixSparseColumnForSuperMatrixSparseColumn(superRow, superSparseColumnToAdd) && "Must not already have this entry");
    
    auto subRow = GetSubMatrixRowForSuperMatrixRow(superRow);
    auto subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps.at(subRow);
    
    auto addedSubSparseColumn = subMatrix->InsertNewEntry(subRow, superDenseColumnToAdd - denseRowColumnOffset);
    subSparseToSuperSparseColumnMap.insert(subSparseToSuperSparseColumnMap.begin()+addedSubSparseColumn, superSparseColumnToAdd);
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

bool SparseMatrixIndexRemapper::HasSubMatrixSparseColumnForSuperMatrixSparseColumn(int superMatrixRow, int superMatrixSparseColumn)
{
    return GetSubMatrixSparseColumnForSuperMatrixSparseColumn(superMatrixRow, superMatrixSparseColumn) != -1;
}

int SparseMatrixIndexRemapper::GetSubMatrixSparseColumnForSuperMatrixSparseColumn(int superMatrixRow, int superMatrixSparseColumn)
{
    assert(superMatrixRow < superMatrixToSubMatrixRowMap.size());
    int subMatrixRow = superMatrixToSubMatrixRowMap[superMatrixRow];
    assert(subMatrixRow != -1);
    const auto& subSparseToSuperSparseColumnMap = subMatrixSparseToSuperMatrixSparseColumnMaps.at(subMatrixRow);
    for (int i=0; i<subSparseToSuperSparseColumnMap.size(); i++) {
        int thisSuperSparseColumn = subSparseToSuperSparseColumnMap[i];
        if (thisSuperSparseColumn == superMatrixSparseColumn) {
            return i;
        }
    }
    
    return -1;
}

