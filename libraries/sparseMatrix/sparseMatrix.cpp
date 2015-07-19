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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "sparseMatrix.h"
#include "sparseSubMatrixLinkage.h"
using namespace std;

SparseMatrix::SparseMatrix(const char * filename)
{
    SparseMatrixOutline sparseMatrixOutline(filename);
    InitFromOutline(&sparseMatrixOutline);
}

SparseMatrix::SparseMatrix(SparseMatrixOutline * sparseMatrixOutline)
{
    InitFromOutline(sparseMatrixOutline);
}

// construct matrix from the outline
void SparseMatrix::InitFromOutline(SparseMatrixOutline * sparseMatrixOutline)
{
    Allocate(sparseMatrixOutline->GetNumRows());
    
    for(size_t i=0; i<GetNumRows(); i++)
    {
        int rowLength = (int)sparseMatrixOutline->columnEntries[i].size();
        columnEntries[i].resize(rowLength);
        columnIndices[i].resize(rowLength);
        
        map<int,double>::iterator pos;
        int j = 0;
        int prev = -1;
        for(pos = sparseMatrixOutline->columnEntries[i].begin(); pos != sparseMatrixOutline->columnEntries[i].end(); pos++)
        {
            columnIndices[i][j] = pos->first;
            if (columnIndices[i][j] <= prev)
                printf("Warning: entries not sorted in a row in a sparse matrix.\n");
            prev = columnIndices[i][j];
            columnEntries[i][j] = pos->second;
            j++;
        }
    }
}

// allocator
void SparseMatrix::Allocate(size_t numRows)
{
    // compressed row storage
    columnIndices.resize(numRows);
    columnEntries.resize(numRows);
    /*
    numSubMatrixIDs = 0;
    subMatrixIndices = NULL;
    subMatrixIndexLengths = NULL;*/
    superMatrixIndices = NULL;
    superRows = NULL;
    diagonalIndices.resize(0);
    transposedIndices.resize(0);
}

// destructor
SparseMatrix::~SparseMatrix()
{
    subMatrixLinkages.clear();
    
    /*
    if (subMatrixIndices != NULL)
    {
        for(int i=numSubMatrixIDs-1; i>=0; i--)
            FreeSubMatrixIndices(i);
    }*/
    
     
    if (superRows != NULL)
    {
        for(int i=0; i<GetNumRows(); i++)
        {
            free(superMatrixIndices[i]);
        }
        free(superMatrixIndices);
        free(superRows);
    }
    
    FreeTranspositionIndices();
}

// copy constructor
SparseMatrix::SparseMatrix(const SparseMatrix & source) 
{
    //printf("Copy constructor\n");fflush(NULL);
    int numRows = source.GetNumRows();
    Allocate(numRows);
    
    for(int i=0; i<numRows; i++)
    {
        int rowLength = source.GetRowLength(i);
        columnIndices[i].resize(rowLength);
        columnEntries[i].resize(rowLength);
        
        for(int j=0; j < rowLength; j++)
        {
            columnIndices[i][j] = source.columnIndices[i][j];
            columnEntries[i][j] = source.columnEntries[i][j];
        }
    }
    
    for (int i=0; i<source.subMatrixLinkages.size(); i++)
    {
        AttachSubMatrix(source.subMatrixLinkages[i]);
    }
    
    /*
    subMatrixIndices = NULL;
    subMatrixIndexLengths = NULL;
    numSubMatrixIDs = source.numSubMatrixIDs;
    if (source.subMatrixIndices != NULL)
    {
        subMatrixIndices = (int***) malloc (sizeof(int**) * numSubMatrixIDs);
        memcpy(subMatrixIndices, source.subMatrixIndices, sizeof(int**) * numSubMatrixIDs);
        subMatrixIndexLengths = (int**) malloc (sizeof(int*) * numSubMatrixIDs);
        memcpy(subMatrixIndexLengths, source.subMatrixIndexLengths, sizeof(int*) * numSubMatrixIDs);
        
        for(int matrixID=0; matrixID < numSubMatrixIDs; matrixID++)
        {
            if (source.subMatrixIndices[matrixID] == NULL)
            {
                subMatrixIndices[matrixID] = NULL;
                subMatrixIndexLengths[matrixID] = NULL;
                continue;
            }
            
            subMatrixIndices[matrixID] = (int**) malloc(sizeof(int*) * numRows);
            subMatrixIndexLengths[matrixID] = (int*) malloc(sizeof(int) * numRows);
            
            for(int i=0; i<numRows; i++)
            {
                subMatrixIndexLengths[matrixID][i] = source.subMatrixIndexLengths[matrixID][i];
                subMatrixIndices[matrixID][i] = (int*) malloc(sizeof(int) * subMatrixIndexLengths[matrixID][i]);
                for(int j=0; j < subMatrixIndexLengths[matrixID][i]; j++)
                {
                    subMatrixIndices[matrixID][i][j] = source.subMatrixIndices[matrixID][i][j];
                }
            }
        }
    }*/
    
    superRows = NULL;
    superMatrixIndices = NULL;
    if (source.superRows != NULL)
    {
        superRows = (int*) malloc(sizeof(int) * numRows);
        superMatrixIndices = (int**) malloc(sizeof(int*) * numRows);
        for(int i=0; i<numRows; i++)
        {
            superRows[i] = source.superRows[i];
            int rowLength = GetRowLength(i);
            superMatrixIndices[i] = (int*) malloc(sizeof(int) * rowLength);
            for(int j=0; j < rowLength; j++)
            {
                superMatrixIndices[i][j] = source.superMatrixIndices[i][j];
            }
        }
    }
    
    if (source.diagonalIndices.size())
    {
        diagonalIndices.resize(source.diagonalIndices.size());
        memcpy(&diagonalIndices[0], &source.diagonalIndices[0], sizeof(diagonalIndices[0]) * numRows);
    }
    
    if (source.transposedIndices.size())
    {
        transposedIndices.resize(numRows);
        for(int i=0; i<numRows; i++)
        {
            int rowLength = GetRowLength(i);
            transposedIndices[i].resize(rowLength);
            memcpy(&transposedIndices[i][0], &source.transposedIndices[i][0], sizeof(transposedIndices[i][0]) * rowLength);
        }
    }
}

void SparseMatrix::MultiplyVector(int startRow, int endRow, const double * vector, double * result) const // result = A(startRow:endRow-1,:) * vector
{
    for(int i=startRow; i<endRow; i++)
    {
        int rowLength = GetRowLength(i);
        result[i-startRow] = 0;
        for(int j=0; j < rowLength; j++)
            result[i-startRow] += vector[columnIndices[i][j]] * columnEntries[i][j];
    }
}

void SparseMatrix::MultiplyVector(const double * vector, double * result) const
{
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        result[i] = 0;
        for(int j=0; j < rowLength; j++)
            result[i] += vector[columnIndices[i][j]] * columnEntries[i][j];
    }
}

void SparseMatrix::MultiplyVectorAdd(const double * vector, double * result) const
{
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            result[i] += vector[columnIndices[i][j]] * columnEntries[i][j];
    }
}

void SparseMatrix::TransposeMultiplyVector(const double * vector, int resultLength, double * result) const
{
    for(int i=0; i<resultLength; i++)
        result[i] = 0;
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            result[columnIndices[i][j]] += vector[i] * columnEntries[i][j];
    }
}

void SparseMatrix::TransposeMultiplyVectorAdd(const double * vector, double * result) const
{
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
        {
            result[columnIndices[i][j]] += vector[i] * columnEntries[i][j];
        }
    }
}

void SparseMatrix::MultiplyMatrix(int numDenseRows, int numDenseColumns, const double * denseMatrix, double * result) const
{
    for(int column=0; column<numDenseColumns; column++)
        MultiplyVector(&denseMatrix[numDenseRows * column], &result[GetNumRows() * column]);
}

void SparseMatrix::MultiplyMatrixAdd(int numDenseRows, int numDenseColumns, const double * denseMatrix, double * result) const
{
    for(int column=0; column<numDenseColumns; column++)
        MultiplyVectorAdd(&denseMatrix[numDenseRows * column], &result[GetNumRows() * column]);
}

// result = A * trans(denseMatrix) 
// trans(denseMatrix) is a dense matrix with 'numDenseColumns' columns, result is a numRows x numDenseColumns dense matrix
void SparseMatrix::MultiplyMatrixTranspose(int numDenseColumns, const double * denseMatrix, double * result) const
{
    memset(result, 0, sizeof(double) * GetNumRows() * numDenseColumns);
    for(int column=0; column<numDenseColumns; column++)
    {
        for(int i=0; i<GetNumRows(); i++)
        {
            int rowLength = GetRowLength(i);
            for(int j=0; j < rowLength; j++)
                result[GetNumRows() * column + i] += denseMatrix[numDenseColumns * columnIndices[i][j] + column] * columnEntries[i][j];
        }
    }
}

double SparseMatrix::QuadraticForm(const double * vector) const
{
    double result = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
        {
            int index = columnIndices[i][j];
            if (index < i)
                continue;
            if (index == i)
                result += columnEntries[i][j] * vector[i] * vector[index];
            else
                result += 2.0 * columnEntries[i][j] * vector[i] * vector[index];
        }
    }
    
    return result;
} 

void SparseMatrix::NormalizeVector(double * vector) const
{
    double norm = sqrt(QuadraticForm(vector));
    for(int i=0; i<GetNumRows(); i++)
        vector[i] /= norm;
}

SparseMatrix SparseMatrix::operator+(const SparseMatrix & mat2) const
{
    SparseMatrix result(*this);
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            result.columnEntries[i][j] += mat2.columnEntries[i][j];
    }
    return result;
}

SparseMatrix SparseMatrix::operator-(const SparseMatrix & mat2) const
{
    SparseMatrix result(*this);
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            result.columnEntries[i][j] -= mat2.columnEntries[i][j];
    }
    return result;
}

SparseMatrix operator* (const double alpha, const SparseMatrix & mat2)
{
    SparseMatrix result(mat2);
    for(int i=0; i<result.GetNumRows(); i++)
    {
        int rowLength = result.GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            result.columnEntries[i][j] *= alpha;
    }
    return result;
}

SparseMatrix & SparseMatrix::operator*=(const double alpha)
{
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            columnEntries[i][j] *= alpha;
    }
    return *this;
}

SparseMatrix & SparseMatrix::operator+=(const SparseMatrix & mat2)
{   
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            columnEntries[i][j] += mat2.columnEntries[i][j];
    }
    return *this;
}

SparseMatrix & SparseMatrix::operator-=(const SparseMatrix & mat2)
{  
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            columnEntries[i][j] -= mat2.columnEntries[i][j];
    }
    return *this;
}

SparseMatrix & SparseMatrix::operator=(const SparseMatrix & source)
{
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            columnEntries[i][j] = source.columnEntries[i][j];
    }
    
    return *this;
}

void SparseMatrix::ScalarMultiply(const double alpha, SparseMatrix * dest)
{
    if (dest == NULL)
        dest = this;
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            dest->columnEntries[i][j] = columnEntries[i][j] * alpha;
    }
}

void SparseMatrix::ScalarMultiplyAdd(const double alpha, SparseMatrix * dest)
{
    if (dest == NULL)
        dest = this;
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            dest->columnEntries[i][j] += columnEntries[i][j] * alpha;
    }
}

void SparseMatrix::ResetToZero()
{
    for(int i=0; i<GetNumRows(); i++)
        ResetRowToZero(i);
}

void SparseMatrix::ResetRowToZero(int row)
{
    int rowLength = GetRowLength(row);
    columnEntries[row].assign(rowLength, 0);
}

void SparseMatrix::Print(bool sparsePrint) const
{
    if (sparsePrint)
    {
        for (int i=0; i<GetNumRows(); i++)
        {
            int rowLength = GetRowLength(i);
            for(int j=0; j< rowLength; j++)
                printf("%d %d %G\n", i, columnIndices[i][j], columnEntries[i][j]);
        }
    }
    else
    {
        int numColumns = GetNumColumns();
        for (int i=0; i<GetNumRows(); i++)
        {
            int index = 0;
            int rowLength = GetRowLength(i);
            for(int j=0; j< rowLength; j++)
            {
                while (index < columnIndices[i][j])
                {
                    index++;
                    printf("%f,",0.0);
                }
                printf("%f,",columnEntries[i][j]);
                index++;
            }
            
            while (index < numColumns)
            {
                index++;
                printf("%f,",0.0);
            }
            
            printf("\n");
        } 
    }
}

// seeks the element in column jDense (in row "row")
// if not found, returns -1
int SparseMatrix::GetInverseIndex(int row, int jDense) const
{
    int rowLength = GetRowLength(row);
    for(int j=0; j < rowLength; j++)
        if (columnIndices[row][j] == jDense)
            return j;
    
    return -1;
}

void SparseMatrix::BuildDiagonalIndices()
{
    if (HasCachedDiagonalIndices())
        return;
    
    diagonalIndices.resize(GetNumRows());
    for(int i=0; i<GetNumRows(); i++)
        diagonalIndices[i] = GetInverseIndex(i,i);
}

void SparseMatrix::FreeDiagonalIndices()
{
    diagonalIndices.clear();
}

void SparseMatrix::GetDiagonal(double * diagonal)
{
    if (HasCachedDiagonalIndices())
    {
        for(int i=0; i<GetNumRows(); i++)
            diagonal[i] = columnEntries[i][diagonalIndices[i]];
    }
    else
    {
        for(int i=0; i<GetNumRows(); i++)
            for(int j=0; j<GetRowLength(i); j++)
            {
                if (GetColumnIndex(i, j) == i)
                    diagonal[i] = columnEntries[i][j];
            }
    }
}

void SparseMatrix::AddDiagonalMatrix(double * diagonalMatrix)
{
    if (HasCachedDiagonalIndices())
    {
        for(int i=0; i<GetNumRows(); i++)
            columnEntries[i][diagonalIndices[i]] += diagonalMatrix[i];
    }
    else
    {
        for(int i=0; i<GetNumRows(); i++)
            for(int j=0; j<GetRowLength(i); j++)
            {
                if (GetColumnIndex(i, j) == i)
                    columnEntries[i][j] += diagonalMatrix[i];
            }
    }
}

void SparseMatrix::AddDiagonalMatrix(double constDiagonalElement)
{
    if (HasCachedDiagonalIndices())
    {
        for(int i=0; i<GetNumRows(); i++)
            columnEntries[i][diagonalIndices[i]] += constDiagonalElement;
    }
    else
    {
        for(int i=0; i<GetNumRows(); i++)
            for(int j=0; j<GetRowLength(i); j++)
            {
                if (GetColumnIndex(i, j) == i)
                    columnEntries[i][j] += constDiagonalElement;
            }
    }
}

int SparseMatrix::GetNumEntries() const
{
    int num = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        num += GetRowLength(i);
    }
    
    return num;
}

double SparseMatrix::SumEntries() const
{
    double sum=0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
            sum += columnEntries[i][j];
    }
    
    return sum;
}

void SparseMatrix::SumRowEntries(double * rowSums) const
{
    for(int i=0; i<GetNumRows(); i++)
    {
        double sum=0;
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
            sum += columnEntries[i][j];
        rowSums[i] = sum;
    }
}

void SparseMatrix::MakeLinearDataArray(double * data) const
{
    int count=0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            data[count] = columnEntries[i][j];
            count++;
        }
    }   
}

void SparseMatrix::MakeLinearRowIndexArray(double * indices) const
{
    int count=0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            indices[count] = i;
            count++;
        }
    }   
}

void SparseMatrix::MakeLinearColumnIndexArray(double * indices) const
{
    int count=0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            indices[count] = columnIndices[i][j];
            count++;
        }
    }   
}

void SparseMatrix::MakeLinearRowIndexArray(int * indices) const
{
    int count=0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            indices[count] = i;
            count++;
        }
    }   
}

void SparseMatrix::MakeLinearColumnIndexArray(int * indices) const
{
    int count=0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            indices[count] = columnIndices[i][j];
            count++;
        }
    }   
}

void SparseMatrix::FreeTranspositionIndices()
{
    transposedIndices.clear();
}

void SparseMatrix::BuildTranspositionIndices()
{
    if (HasCachedTransposedIndices())
        return;
    
    transposedIndices.resize(GetNumRows());
    
    vector<int> buffer(GetNumColumns());
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        transposedIndices[i].resize(rowLength);
        for(int j=0; j<rowLength; j++)
        {
            transposedIndices[i][j] = buffer[columnIndices[i][j]];
            buffer[columnIndices[i][j]]++;
        }   
    }  
}

double SparseMatrix::SkewSymmetricCheck()
{
    double maxEntry = 0;  
    
    BuildTranspositionIndices();  
    
    for(int i=0; i<GetNumRows(); i++)
    {    
        for(int j=0; j<GetRowLength(i); j++)    
        {      
            double entry1 = GetEntry(i, j);      
            int tindex = TransposedIndex(i, j);      
            double entry2 = GetEntry(GetColumnIndex(i,j), tindex);      
            
            // entry1 + entry2 should be zero          
            if (fabs(entry1 + entry2) > maxEntry)
                maxEntry = fabs(entry1 + entry2);
        }  
    }  
    
    FreeTranspositionIndices();
    
    return maxEntry;
}

void SparseMatrix::SymmetrizeMatrix()
{
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            int jAbs = columnIndices[i][j];
            
            if (jAbs >= i)
                break; 
            
            // copy elt (jAbs,i) into position (i,jAbs)
            columnEntries[i][j] = columnEntries[jAbs][TransposedIndex(i,j)];
        }
    }
}

double SparseMatrix::GetMaxAbsEntry() const
{
    double max = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            double el = fabs(GetEntry(i,j));
            if (el > max)
                max = el;
        }
    }
    
    return max;
}

double SparseMatrix::GetRowNorm2(int row) const
{
    double norm2 = 0;
    int rowLength = GetRowLength(row);
    for(int j=0; j<rowLength; j++)
    {
        double el = columnEntries[row][j];
        norm2 += el * el;
    }
    return norm2;
}

// solve M * x = b
// ASSUMES the sparse matrix is diagonal !!!
void SparseMatrix::DiagonalSolve(double * rhs) const
{
    for(int i=0; i<GetNumRows(); i++)
        rhs[i] /= columnEntries[i][0]; // the diagonal element
}

void SparseMatrix::BuildRenumberingVector(int nConstrained, int nSuper, int numFixedDOFs, int * fixedDOFs, int ** superDOFs, int oneIndexed)
{
    // superRows[i] is the row index in the super matrix corresponsing to row i of constrained matrix
    (*superDOFs) = (int*) malloc (sizeof(int) * nConstrained);
    int constrainedDOF = 0;
    int superDOF = 0;
    for(int i=0; i<numFixedDOFs; i++)
    {
        int nextSuperDOF = fixedDOFs[i];
        nextSuperDOF -= oneIndexed;
        if ( (nextSuperDOF >= nSuper) || (nextSuperDOF < 0) )
        {
            printf("Error: invalid fixed super DOF %d specified.\n", nextSuperDOF);
            exit(1);
        }
        
        while (superDOF < nextSuperDOF)
        {
            if (constrainedDOF >= nConstrained)
            {
                printf("Error: too many DOFs specified.\n");
                exit(1);
            }
            (*superDOFs)[constrainedDOF] = superDOF; 
            constrainedDOF++;
            superDOF++;
        }
        
        superDOF++; // skip the deselected DOF
    }
    while (superDOF < nSuper)
    {
        if (constrainedDOF >= nConstrained)
        {
            printf("Error: too many DOFs specified.\n");
            exit(1);
        }
        (*superDOFs)[constrainedDOF] = superDOF; 
        
        constrainedDOF++;
        superDOF++;
    }
}

void SparseMatrix::BuildSuperMatrixIndices(int numFixedRowColumns, int * fixedRowColumns, SparseMatrix * superMatrix, int oneIndexed)
{
    BuildSuperMatrixIndices(numFixedRowColumns, fixedRowColumns, numFixedRowColumns, fixedRowColumns, superMatrix, oneIndexed); 
}

void SparseMatrix::BuildSuperMatrixIndices(int numFixedRows, int * fixedRows, int numFixedColumns, int * fixedColumns, SparseMatrix * superMatrix, int oneIndexed)
{
    int numSuperColumns = superMatrix->GetNumColumns();
    int numColumns = numSuperColumns - numFixedColumns;
    
    if ((GetNumRows() + numFixedRows != superMatrix->GetNumRows()) || (GetNumColumns() + numFixedColumns > numSuperColumns) )
    {
        printf("Error in BuildSuperMatrixIndices: number of constrained DOFs does not match the size of the two matrices.\n");
        printf("num rows: %d num fixed rows in super matrix: %d num rows in super matrix: %d\n", GetNumRows(), numFixedRows, superMatrix->GetNumRows());
        printf("num columns: %d num fixed columns in super matrix: %d num columns in super matrix: %d\n", numColumns, numFixedColumns, numSuperColumns);
        exit(1);
    }
    
    // build row renumbering function:
    BuildRenumberingVector(GetNumRows(), superMatrix->GetNumRows(), numFixedRows, fixedRows, &superRows, oneIndexed);
    // build column renumbering function:
    int * superColumns_;
    BuildRenumberingVector(numColumns, numSuperColumns, numFixedColumns, fixedColumns, &superColumns_, oneIndexed);
    
    // superRows[i] is the row index in the super matrix corresponsing to row i of constrained matrix
    // superColumns_[i] is the dense column index in the super matrix corresponsing to the dense column i of constrained matrix
    
    // build column indices
    superMatrixIndices = (int**) malloc (sizeof(int*) * GetNumRows());
    for(int i=0; i < GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        superMatrixIndices[i] = (int*) malloc (sizeof(int) *  rowLength);
        for(int j=0; j < rowLength; j++)
        {
            int iConstrained = i;
            int jConstrainedDense = columnIndices[iConstrained][j];
            int iSuper = superRows[iConstrained];
            int jSuperDense = superColumns_[jConstrainedDense];
            int jSuper = superMatrix->GetInverseIndex(iSuper, jSuperDense);
            if (jSuper < 0)
            {
                printf("Error in BuildSuperMatrixIndices: failed to compute inverse index.\n");
                printf("i=%d j=%d iConstrained=%d jConstrainedDense=%d iSuper=%d jSuperDense=%d jSuper=%d\n", i, j, iConstrained, jConstrainedDense, iSuper, jSuperDense, jSuper);
                fflush(NULL);
                exit(1);
            }
            superMatrixIndices[i][j] = jSuper;
        }
    } 
    
    free(superColumns_);
}

void SparseMatrix::AssignSuperMatrix(SparseMatrix * superMatrix)
{
    for(int i=0; i<GetNumRows(); i++)
    {
        const vector<double>& row = superMatrix->columnEntries[superRows[i]];
        int * indices = superMatrixIndices[i];
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            columnEntries[i][j] = row[indices[j]];
    }
}

shared_ptr<SparseSubMatrixLinkage> SparseMatrix::AttachSubMatrix(shared_ptr<SparseMatrix> subMatrix)
{
    assert(GetExistingSubMatrixLinkage(subMatrix) == nullptr && "already have a linkage for this submatrix");
    
    shared_ptr<SparseMatrix> super = shared_from_this();
    auto linkage = std::make_shared<SparseSubMatrixLinkage>(super, subMatrix);
    AttachSubMatrix(linkage);
    return linkage;
}

void SparseMatrix::AttachSubMatrix(shared_ptr<SparseSubMatrixLinkage> link)
{
    subMatrixLinkages.push_back(link);
}

void SparseMatrix::DetachSubMatrix(shared_ptr<SparseSubMatrixLinkage> link)
{
    auto it = find(subMatrixLinkages.begin(), subMatrixLinkages.end(), link);
    assert(it != subMatrixLinkages.end() && "submatrix linkage not found");
    if (it != subMatrixLinkages.end())
    {
        subMatrixLinkages.erase(it);
    }
}

shared_ptr<SparseSubMatrixLinkage> SparseMatrix::GetExistingSubMatrixLinkage(shared_ptr<SparseMatrix> subMatrix)
{
    for (auto link: subMatrixLinkages)
    {
        if (link->GetSubMatrix() == subMatrix)
        {
            return link;
        }
    }
    
    return nullptr;
}

/*
void SparseMatrix::BuildSubMatrixIndices(SparseMatrix & submatrix, int subMatrixID)
{
    if (subMatrixID >= numSubMatrixIDs)
    {
        subMatrixIndices = (int***) realloc (subMatrixIndices, sizeof(int**) * (subMatrixID + 1));
        subMatrixIndexLengths = (int**) realloc (subMatrixIndexLengths, sizeof(int*) * (subMatrixID + 1));
        for(int i=numSubMatrixIDs; i <= subMatrixID; i++)
        {
            subMatrixIndices[i] = NULL;
            subMatrixIndexLengths[i] = NULL;
        }
        numSubMatrixIDs = subMatrixID + 1;
    }
    
    if ((subMatrixIndices[subMatrixID] != NULL) || (subMatrixIndexLengths[subMatrixID] != NULL))
    {
        free(subMatrixIndices[subMatrixID]);
        free(subMatrixIndexLengths[subMatrixID]);
        subMatrixIndices[subMatrixID] = NULL;
        subMatrixIndexLengths[subMatrixID] = NULL;
        //printf("Warning: old submatrix indices (matrixID %d) have not been de-allocated.\n", subMatrixID);
    }
    
    subMatrixIndices[subMatrixID] = (int**) malloc (sizeof(int*) * GetNumRows());
    subMatrixIndexLengths[subMatrixID] = (int*) malloc (sizeof(int) * GetNumRows());
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int submatrixRowLength = submatrix.GetRowLength(i);
        subMatrixIndices[subMatrixID][i] = (int*) malloc (sizeof(int) * submatrixRowLength);
        subMatrixIndexLengths[subMatrixID][i] = submatrixRowLength;
        const vector<int>& indices = submatrix.columnIndices[i];
        for(int j=0; j < submatrixRowLength; j++)
        {
            // finds the position in row i of element with column index jDense
            // int inverseIndex(int i, int jDense);
            subMatrixIndices[subMatrixID][i][j] = GetInverseIndex(i, indices[j]);
            if (subMatrixIndices[subMatrixID][i][j] == -1)
            {
                printf("Error (BuildSubMatrixIndices): given matrix is not a submatrix of this matrix. The following index does not exist in this matrix: (%d,%d)\n", i, indices[j]);
                exit(1);
            }
        }
    }
}

void SparseMatrix::FreeSubMatrixIndices(int subMatrixID)
{
    if (subMatrixID >= numSubMatrixIDs)
    {
        printf("Warning: attempted to free submatrix index that does not exist.\n");
        return;
    }
    
    if (subMatrixIndices[subMatrixID] != NULL)
    {
        for(int i=0; i<GetNumRows(); i++)
            free(subMatrixIndices[subMatrixID][i]); 
        free(subMatrixIndices[subMatrixID]);
        free(subMatrixIndexLengths[subMatrixID]);
        subMatrixIndices[subMatrixID] = NULL;
        subMatrixIndexLengths[subMatrixID] = NULL;
    }
    
    // check if this was the largest index
    for(int i=numSubMatrixIDs-1; i>=0; i--)
    {
        if (subMatrixIndices[i] != NULL)
        {
            numSubMatrixIDs = i + 1;
            subMatrixIndices = (int***) realloc (subMatrixIndices, sizeof(int**) * numSubMatrixIDs);
            subMatrixIndexLengths = (int**) realloc (subMatrixIndexLengths, sizeof(int*) * numSubMatrixIDs);
            break;
        }
        
        numSubMatrixIDs = i;
    }
    
    if (numSubMatrixIDs == 0)
    {
        free(subMatrixIndices);
        free(subMatrixIndexLengths);
        subMatrixIndices = NULL;
        subMatrixIndexLengths = NULL;
    }
}
 */


void SparseMatrix::AddSubMatrix(double factor, shared_ptr<SparseMatrix> subMatrix)
{
    auto linkage = GetExistingSubMatrixLinkage(subMatrix);
    assert(nullptr != linkage && "No linkage to the submatrix exists");
    AddSubMatrix(factor, linkage);
}

void SparseMatrix::AddSubMatrix(double factor, shared_ptr<SparseSubMatrixLinkage> link)
{
    assert(link->GetSuperMatrix().get() == this && "link is for a different super matrix");
    link->AddSubMatrixToSuperMatrix(factor);
    
    /*
    for(int i=0; i<GetNumRows(); i++)
    {
        int * indices = subMatrixIndices[subMatrixID][i];
        int submatrixRowLength = submatrix.GetRowLength(i);
        for(int j=0; j < submatrixRowLength; j++)
            columnEntries[i][indices[j]] += factor * submatrix.columnEntries[i][j];
    }
     */
}

int SparseMatrix::GetNumLowerTriangleEntries() const
{
    int num = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
        {
            if (columnIndices[i][j] <= i)
                num++;
        }
    }
    return num;
}

int SparseMatrix::GetNumUpperTriangleEntries() const
{
    int num = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
        {
            if (columnIndices[i][j] >= i)
                num++;
        }
    }
    return num;
}

int SparseMatrix::GenerateNAGFormat(double * a, int * irow, int * icol, int * istr) const
{
    int num = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        istr[i] = num; // starting address of row i
        for(int j=0; j < rowLength; j++)
        {
            if (columnIndices[i][j] <= i) // over lower triangle
            {
                a[num] = columnEntries[i][j];
                irow[num] = i+1; // NAG is 1-indexed
                icol[num] = columnIndices[i][j]+1; // NAG is 1-indexed
                num++;
            }
        }
    }
    
    istr[GetNumRows()] = num;
    
    return num;
} 

void SparseMatrix::GenerateCompressedRowMajorFormat(double * a, int * ia, int * ja, int upperTriangleOnly, int oneIndexed) const
{
    int count = 0;
    for(int row=0; row<GetNumRows(); row++)
    {
        if (ia != NULL)
            ia[row] = count + oneIndexed;
        
        int rowLength = GetRowLength(row);
        for(int j=0; j< rowLength; j++)
        {
            if ((!upperTriangleOnly) || (columnIndices[row][j] >= row))
            {
                if (a != NULL)
                    a[count] = columnEntries[row][j];         
                if (ja != NULL)
                    ja[count] = columnIndices[row][j] + oneIndexed; 
                count++;
            }
        }
    }
    
    if (ia != NULL)
        ia[GetNumRows()] = count + oneIndexed;
}

void SparseMatrix::GenerateCompressedRowMajorFormat_four_array(double * values, int * columns, int * pointerB, int * pointerE, int upperTriangleOnly, int oneIndexed) const
{
    int count = 0;
    for(int row=0; row<GetNumRows(); row++)
    {
        if (pointerB != NULL)
            pointerB[row] = count + oneIndexed;
        
        int rowLength = GetRowLength(row);
        for(int j=0; j< rowLength; j++)
        {
            if ((!upperTriangleOnly) || (columnIndices[row][j] >= row))
            {
                if (values != NULL)
                    values[count] = columnEntries[row][j];         
                if (columns != NULL)
                    columns[count] = columnIndices[row][j] + oneIndexed; 
                count++;
            }
        }
        
        if (pointerE != NULL)
            pointerE[row] = count + oneIndexed;
    }
}

int SparseMatrix::Save(const char * filename, int oneIndexed) const
{
    FILE * fout = fopen(filename,"w");
    if (!fout)
        return 1;
    
    fprintf(fout,"%d\n%d\n", GetNumRows(), GetNumColumns());
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
        {
            int index = columnIndices[i][j]; 
            double entry = columnEntries[i][j];
            fprintf(fout,"%d %d %.15G\n",i + oneIndexed, index + oneIndexed, entry);
        }
    }
    fclose(fout);
    
    return 0;
}

int SparseMatrix::SaveToMatlabFormat(const char * filename) const
{
    FILE * fout = fopen(filename,"w");
    if (!fout)
        return 1;
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j < rowLength; j++)
        {
            int index = columnIndices[i][j]; 
            double entry = columnEntries[i][j];
            fprintf(fout,"%d %d %.15G\n",i + 1, index + 1, entry);
        }
    }
    fclose(fout);
    
    return 0;
}

void SparseMatrix::RemoveRowColumn(int index)
{
    // remove row 'index'
    columnEntries.erase(columnEntries.begin()+index);
    columnIndices.erase(columnIndices.begin()+index);
    
    // remove column 'index'
    for(int i=0; i<GetNumRows(); i++)
    {
        // traverse all rows
        int cachedRowLength = GetRowLength(i);
        for(int j=0; j<cachedRowLength; j++)
        {
            // seek for entry 'index'
            if (columnIndices[i][j] == index) // found
            {
                // shift all elements ahead one step back
                columnIndices[i].erase(columnIndices[i].begin()+j);
                columnEntries[i].erase(columnEntries[i].begin()+j);
                --cachedRowLength;
            }
        }
        
        // decrease indices for DOFs above index
        for(int j=0; j<cachedRowLength; j++)
        {
            if(columnIndices[i][j] > index)     
            {
                // decrease index
                --columnIndices[i][j];
            }
        }   
    }
}

void SparseMatrix::RemoveRowsColumnsSlow(int numRemovedRowsColumns, int * removedRowsColumns, int oneIndexed)
{
    for(int i=0; i<numRemovedRowsColumns; i++)
        RemoveRowColumn(removedRowsColumns[i]-i-oneIndexed);
}

void SparseMatrix::RemoveRowsColumns(int numRemovedRowsColumns, int * removedRowsColumns, int oneIndexed)
{
    // the removed dofs must be pre-sorted
    // build a map from old dofs to new ones
    vector<int> oldToNew(GetNumRows());
    int dof = 0;
    int dofCount = 0;
    for(int i=0; i<numRemovedRowsColumns; i++)
    {
        while (dof < removedRowsColumns[i] - oneIndexed)
        {
            oldToNew[dof] = dofCount;
            dofCount++;
            dof++;
        }
        oldToNew[dof] = -1;
        dof++;
    }
    while (dof < GetNumRows())
    {
        oldToNew[dof] = dofCount;
        dofCount++;
        dof++;
    }
    
    // now, traverse all rows and renumber the entries
    int targetRow = 0;
    for(int sourceRow = 0; sourceRow < GetNumRows(); sourceRow++)
    {
        if (oldToNew[sourceRow] == -1)
        {
            //free(columnIndices[sourceRow]);
            //free(columnEntries[sourceRow]);
            continue;
        }
        
        int targetIndex = 0;
        for(int sourceIndex=0; sourceIndex<GetRowLength(sourceRow); sourceIndex++)
        {
            int oldIndex = columnIndices[sourceRow][sourceIndex];
            int newIndex = oldToNew[oldIndex];
            if (newIndex == -1)
                continue;
            columnIndices[sourceRow][targetIndex] = newIndex;
            columnEntries[sourceRow][targetIndex] = columnEntries[sourceRow][sourceIndex];
            ++targetIndex;
        }
        
        columnIndices[sourceRow].resize(targetIndex);
        columnEntries[sourceRow].resize(targetIndex);
        
        // replace target row with source row, fast`
        columnIndices[targetRow].swap(columnIndices[sourceRow]);
        columnEntries[targetRow].swap(columnEntries[sourceRow]);
        ++targetRow;
    }
    
    int newRowCount = targetRow;
    columnEntries.resize(newRowCount);
    columnIndices.resize(newRowCount);
}

void SparseMatrix::RemoveColumn(int index)
{
    // remove column 'index'
    for(int i=0; i<GetNumRows(); i++)
    {
        // traverse all rows
        int cachedRowLength = GetRowLength(i);
        for(int j=0; j<cachedRowLength; j++)
        {
            // seek for entry 'index'
            if (columnIndices[i][j] == index) // found
            {
                // shift all elements ahead one step back
                for(int k=j; k<cachedRowLength-1; k++)
                {
                    columnIndices[i][k] = columnIndices[i][k+1];
                    columnEntries[i][k] = columnEntries[i][k+1];
                } 
                --cachedRowLength;
            }
        }
        
        // decrease indices for DOFs above index
        for(int j=0; j<cachedRowLength; j++)
        {
            if(columnIndices[i][j] > index)     
            {
                // decrease index
                --columnIndices[i][j];
            }
        }   
    }
}

void SparseMatrix::RemoveColumns(int numRemovedColumns, int * removedColumns, int oneIndexed)
{
    // the removed dofs must be pre-sorted
    // build a map from old dofs to new ones
    int numColumns = GetNumColumns();
    
    // must increase numColumns to accommodate matrices with zero columns on the right
    for(int i=0; i<numRemovedColumns; i++)
    {
        int removedColumn0Indexed = removedColumns[i] - oneIndexed;
        int neededNumColumns = removedColumn0Indexed + 1;
        if (neededNumColumns > numColumns)
            numColumns = neededNumColumns;
    }
    
    vector<int> oldToNew(numColumns);
    int dof = 0;
    int dofCount = 0;
    for(int i=0; i<numRemovedColumns; i++)
    {
        while (dof < removedColumns[i] - oneIndexed)
        {
            oldToNew[dof] = dofCount;
            dofCount++;
            dof++;
        }
        oldToNew[dof] = -1;
        ++dof;
    }
    while (dof < numColumns)
    {
        oldToNew[dof] = dofCount;
        ++dofCount;
        ++dof;
    }
    
    // now, traverse all rows and renumber the entries
    for(int row = 0; row < GetNumRows(); row++)
    {
        int targetIndex = 0;
        int cachedRowLength = GetRowLength(row);
        for(int sourceIndex=0; sourceIndex<cachedRowLength; sourceIndex++)
        {
            int oldIndex = columnIndices[row][sourceIndex];
            int newIndex = oldToNew[oldIndex];
            if (newIndex == -1)
                continue;
            columnIndices[row][targetIndex] = newIndex;
            columnEntries[row][targetIndex] = columnEntries[row][sourceIndex];
            ++targetIndex;
        }
        
        size_t newLength = targetIndex;
        columnIndices[row].resize(newLength);
        columnEntries[row].resize(newLength);
    }
}

void SparseMatrix::RemoveColumnsSlow(int numColumns, int * columns, int oneIndexed)
{
    for(int i=0; i<numColumns; i++)
        RemoveColumn(columns[i]-i-oneIndexed);
}

void SparseMatrix::RemoveRow(int index)
{
    // remove row 'index'
    columnEntries.erase(columnEntries.begin()+index);
    columnIndices.erase(columnIndices.begin()+index);
}

void SparseMatrix::RemoveRowsSlow(int numRows, int * rows, int oneIndexed)
{
    for(int i=0; i<numRows; i++)
        RemoveRow(rows[i]-i-oneIndexed);
}

void SparseMatrix::RemoveRows(int numRemovedRows, int * removedRows, int oneIndexed)
{
    // the removed dofs must be pre-sorted
    // build a map from old dofs to new ones
    vector<int> oldToNew(GetNumRows());
    int dof = 0;
    int dofCount = 0;
    for(int i=0; i<numRemovedRows; i++)
    {
        while (dof < removedRows[i] - oneIndexed)
        {
            oldToNew[dof] = dofCount;
            dofCount++;
            dof++;
        }
        oldToNew[dof] = -1;
        ++dof;
    }
    while (dof < GetNumRows())
    {
        oldToNew[dof] = dofCount;
        ++dofCount;
        ++dof;
    }
    
    // now, traverse all rows and renumber the entries
    int targetRow = 0;
    for(int sourceRow = 0; sourceRow < GetNumRows(); sourceRow++)
    {
        if (oldToNew[sourceRow] == -1)
        {
            //free(columnIndices[sourceRow]);
            //free(columnEntries[sourceRow]);
            continue;
        }
        
        columnIndices[targetRow].swap(columnIndices[sourceRow]);
        columnEntries[targetRow].swap(columnEntries[sourceRow]);
        ++targetRow;
    }
    
    int newRowCount = targetRow;
    columnEntries.resize(newRowCount);
    columnIndices.resize(newRowCount);
}

double SparseMatrix::GetInfinityNorm() const
{
    double norm = 0.0;
    
    for(int i=0; i<GetNumRows(); i++)
    {
        double absRowSum = 0;
        
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            absRowSum += fabs(columnEntries[i][j]);
        }
        
        if (absRowSum > norm)
            norm = absRowSum;
    }
    
    return norm;
}

void SparseMatrix::DoOneGaussSeidelIteration(double * x, const double * b) const
{
    for(int i=0; i<GetNumRows(); i++)
    {
        double buffer = b[i];
        int diagIndex = -1;
        
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            int column = columnIndices[i][j];
            if (column != i)
                buffer -= columnEntries[i][j] * x[column];
            else
                diagIndex = j;
        }
        x[i] = buffer / columnEntries[i][diagIndex];
    }
}

void SparseMatrix::ComputeResidual(const double * x, const double * b, double * residual) const
{
    MultiplyVector(x,residual);
    for(int i=0; i<GetNumRows(); i++)
        residual[i] -= b[i];
}

double SparseMatrix::CheckLinearSystemSolution(const double * x, const double * b, int verbose, double * buffer) const
{
    double * bufferv = NULL;
    
    if (buffer == NULL)
    {
        bufferv = (double*) malloc (sizeof(double) * GetNumRows());
        buffer = bufferv;
    }
    
    MultiplyVector(x,buffer);
    
    double inftyNorm = 0;
    double inftyNorm_b = 0;
    for(int i=0; i<GetNumRows(); i++)
    {
        if (fabs(buffer[i] - b[i]) > inftyNorm)
            inftyNorm = fabs(buffer[i] - b[i]);
        
        if (fabs(b[i]) > inftyNorm_b)
            inftyNorm_b = fabs(b[i]);
    }
    
    if (verbose)
    {
        printf("Infinity residual norm ||Ax-b|| is %G. ||b|| is %G.\n", inftyNorm, inftyNorm_b);
        printf("Relative infinity residual norm ||Ax-b||/||b|| is %G.\n", inftyNorm / inftyNorm_b);
    }
    
    free(bufferv);
    
    return inftyNorm / inftyNorm_b;
}

void SparseMatrix::MakeDenseMatrix(double * denseMatrix) const
{
    memset(denseMatrix, 0, sizeof(double) * (GetNumRows() * GetNumColumns()));
    for(int i=0; i< GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
            denseMatrix[GetNumRows() * columnIndices[i][j] + i] = columnEntries[i][j];
    }
}

void SparseMatrix::MakeDenseMatrixTranspose(int numColumns, double * denseMatrix) const
{
    // note: we cannot use GetNumColumns() here because the rightmost columns of the sparse matrix can be zero and the GetNumColumns() will not be accurate
    memset(denseMatrix, 0, sizeof(double) * (GetNumRows() * numColumns));
    for(int i=0; i<GetNumRows(); i++)
    {
        int offset = i * numColumns;
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
            denseMatrix[offset + columnIndices[i][j]] = columnEntries[i][j];
    }
}

void SparseMatrix::MultiplyRow(int row, double scalar) // multiplies all elements in row 'row' with scalar 'scalar'
{
    int rowLength = GetRowLength(row);
    for(int j=0; j<rowLength; j++)
        columnEntries[row][j] *= scalar;
}

int SparseMatrix::GetNumColumns() const 
{
    int numColumns = -1;
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            numColumns = std::max(columnIndices[i][j], numColumns);
        }
    }
    return numColumns + 1;
}

void SparseMatrix::IncreaseNumRows(int numAddedRows)
{
    int newn = GetNumRows() + numAddedRows;
    
    columnIndices.resize(newn);
    columnEntries.resize(newn);
}

SparseMatrix SparseMatrix::ConjugateMatrix(SparseMatrix & U, int verbose)
{
    SparseMatrixOutline outline(U.GetNumColumns());
    
    for(int i=0; i<GetNumRows(); i++)
    {
        if (verbose)
        {
            if (i % 100 == 1)
                printf("Processing row %d / %d...\n", i, GetNumRows());
        }
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            int I = i;
            int J = columnIndices[i][j];
            double scalar = columnEntries[i][j];
            
            // compute tensor product of rows I and J of U
            int URowLengthI = U.GetRowLength(I);
            for(int k=0; k<URowLengthI; k++)
            {
                int URowLengthJ = U.GetRowLength(J);
                for(int l=0; l<URowLengthJ; l++)
                {
                    int K = U.columnIndices[I][k];
                    int L = U.columnIndices[J][l];
                    
                    // there is an entry at (I, K), and another entry at (J, L); compute their contribution to tensor product:
                    outline.AddEntry(K, L, scalar * U.columnEntries[I][k] * U.columnEntries[J][l]);
                }
            }
        }
    }
    
    if (verbose)
        printf("Creating sparse matrix from outline...\n");
    
    return SparseMatrix(&outline);
}

void SparseMatrix::BuildConjugationIndices(SparseMatrix & U, SparseMatrix & MTilde, precomputedIndicesType * precomputedIndices)
{
    typedef pair< pair<int,int>, pair<int, int> > fourTuple;
    typedef vector<fourTuple> listOfFourTuples;
    typedef map<int, listOfFourTuples> rowMap;
    vector<rowMap> rowMaps;
    for(int i=0; i<MTilde.GetNumRows(); i++)
        rowMaps.push_back(rowMap());
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
        {
            int I = i;
            int J = columnIndices[i][j];
            
            // compute tensor product of rows I and J of U
            int URowLengthI = U.GetRowLength(I);
            for(int k=0; k<URowLengthI; k++)
            {
                int URowLengthJ = U.GetRowLength(J);
                for(int l=0; l<URowLengthJ; l++)
                {
                    int K = U.columnIndices[I][k];
                    int L = U.columnIndices[J][l];
                    //outline.AddEntry(K, L, scalar * U.columnEntries[I][k] * U.columnEntries[J][l]);
                    
                    fourTuple tuple(make_pair(make_pair(i,j), make_pair(k,l)));
                    
                    rowMap::iterator iter = rowMaps[K].find(L);
                    if (iter == rowMaps[K].end())
                    {
                        listOfFourTuples singletonList;
                        singletonList.push_back(tuple);            
                        rowMaps[K].insert(make_pair(L, singletonList));
                    }
                    else
                    {
                        (iter->second).push_back(tuple);
                    }
                }
            }
        }
    }
    
    // copy map to precomputedIndices
    (*precomputedIndices) = (int***) malloc (sizeof(int**) * MTilde.GetNumRows());
    for(int i=0; i<MTilde.GetNumRows(); i++)
    {
        (*precomputedIndices)[i] = (int**) malloc (sizeof(int*) * rowMaps[i].size());
        int j = 0;
        for(rowMap::iterator iter = rowMaps[i].begin(); iter != rowMaps[i].end(); iter++)
        {
            (*precomputedIndices)[i][j] = (int*) malloc (sizeof(int) * (4 * ((iter->second).size()) + 1));
            ((*precomputedIndices)[i][j])[0] = (int)(iter->second).size();
            for(int k=0; k<((*precomputedIndices)[i][j])[0]; k++)
            {
                ((*precomputedIndices)[i][j])[1+4*k+0] = ((iter->second)[k]).first.first;
                ((*precomputedIndices)[i][j])[1+4*k+1] = ((iter->second)[k]).first.second;
                ((*precomputedIndices)[i][j])[1+4*k+2] = ((iter->second)[k]).second.first;
                ((*precomputedIndices)[i][j])[1+4*k+3] = ((iter->second)[k]).second.second;
            }
            j++;
        }
    }
}

void SparseMatrix::ConjugateMatrix(precomputedIndicesType precomputedIndices, SparseMatrix & U, SparseMatrix & MTilde)
{
    MTilde.ResetToZero();
    for(int row=0; row<MTilde.GetNumRows(); row++)
    {
        int ** rowIndices = precomputedIndices[row];
        int MTildeRowLength = MTilde.GetRowLength(row);
        for(int j=0; j<MTildeRowLength; j++)
        {
            int * entryIndices = rowIndices[j];
            int numSummationTerms = entryIndices[0];
            for(int k=0; k<numSummationTerms; k++)
            {
                int * entryIndex = &entryIndices[4*k+1];
                int rowOfM = entryIndex[0];
                int columnIndexOfM = entryIndex[1];
                int columnOfM = columnIndices[rowOfM][columnIndexOfM];
                int columnIndexofU_for_MTilde_row = entryIndex[2];
                int columnIndexofU_for_MTilde_column = entryIndex[3];
                (MTilde.columnEntries)[row][j] += columnEntries[rowOfM][columnIndexOfM] * U.columnEntries[row][columnIndexofU_for_MTilde_row] * U.columnEntries[columnOfM][columnIndexofU_for_MTilde_column];
            }
        }
    }
}

void SparseMatrix::ConjugateMatrix(double * U, int r, double * UTilde)
{
    double * MU = (double*) malloc (sizeof(double) * GetNumRows() * r);
    MultiplyMatrix(GetNumRows(), r, U, MU);
    
    // compute U^T * MU
    for(int i=0; i<r; i++)
        for(int j=0; j<r; j++)
        {
            double entry = 0.0;
            for(int k=0; k<GetNumRows(); k++)
                entry += U[i * GetNumRows() + k] * MU[j * GetNumRows() + k];
            UTilde[j * r + i] = entry;
        }
    
    free(MU);
}

SparseMatrix * SparseMatrix::Transpose(int numColumns)
{
    if (numColumns < 0)
        numColumns = GetNumColumns();
    
    SparseMatrixOutline outline(numColumns);
    
    for(int i=0; i<GetNumRows(); i++)
    {
        int rowLength = GetRowLength(i);
        for(int j=0; j<rowLength; j++)
            outline.AddEntry(columnIndices[i][j], i, columnEntries[i][j]);
    }
    
    return new SparseMatrix(&outline);;
}

void SparseMatrix::SetRows(SparseMatrix * source, int startRow, int startColumn) 
{
    for(int i=0; i<source->GetNumRows(); i++)
    {
        int row = startRow + i;
        if (row >= GetNumRows())
            return;
        
        int newRowLength = source->GetRowLength(i);
        columnIndices[row].resize(newRowLength);
        columnEntries[row].resize(newRowLength);
        for(int j=0; j<newRowLength; j++)
        {
            columnIndices[row][j] = startColumn + source->columnIndices[i][j];
            columnEntries[row][j] = source->columnEntries[i][j];
        }
    }
}

void SparseMatrix::AppendRowsColumns(SparseMatrix * source)
{
    vector<int> oldRowLength = GetRowLengths();
    vector<int> rowLength = GetRowLengths();
    
    int oldNumRows = GetNumRows();
    IncreaseNumRows(source->GetNumRows());
    SetRows(source, oldNumRows);
    
    // add transpose of rows:
    
    // first, establish new column lengths
    for(int row=0; row<source->GetNumRows(); row++)
    {
        for(int j=0; j<source->GetRowLength(row); j++)
        {
            int column = source->GetColumnIndex(row, j);
            ++rowLength[column];
        }
    }
    
    // extend size
    for(int row=0; row<oldNumRows; row++)
    {
        columnIndices[row].resize(rowLength[row]);
        columnEntries[row].resize(rowLength[row]);
    }
    
    // restore old row length
    for(int i=0; i<oldNumRows; i++)
        rowLength[i] = oldRowLength[i];
    
    // write entries into their place
    for(int row=0; row<source->GetNumRows(); row++)
    {
        for(int j=0; j<source->GetRowLength(row); j++)
        {
            int column = source->GetColumnIndex(row, j);
            columnIndices[column][rowLength[column]] = oldNumRows + row;
            columnEntries[column][rowLength[column]] = source->GetEntry(row, j);
            rowLength[column]++;
        }
    }
    
    // append zero diagonal in lower-right block (helps with some solvers)
    for(int row=0; row<source->GetNumRows(); row++)
    {
        rowLength[oldNumRows + row]++;
        columnIndices[oldNumRows + row].push_back(oldNumRows + row);
        columnEntries[oldNumRows + row].push_back(0.0);
    }
}

SparseMatrix * SparseMatrix::CreateIdentityMatrix(int numRows)
{
    SparseMatrixOutline * outline = new SparseMatrixOutline(numRows);
    for (int row=0; row<numRows; row++)
        outline->AddEntry(row, row, 1.0);
    SparseMatrix * mat = new SparseMatrix(outline);
    delete outline;
    return mat;
}

vector<int> SparseMatrix::GetRowLengths() const
{
    vector<int> result(GetNumRows());
    for (int i=0; i<GetNumRows(); i++)
    {
        result[i] = GetRowLength(i);
    }
    return result;
}

SparseMatrixOutline SparseMatrix::GenerateOutline() const
{
    SparseMatrixOutline outline(GetNumRows());
    for (int row=0; row<GetNumRows(); row++)
    {
        for (int entry=0; entry<GetRowLength(row); entry++)
        {
            outline.AddEntry(row, GetColumnIndex(row, entry));
        }
    }
    return outline;
}

