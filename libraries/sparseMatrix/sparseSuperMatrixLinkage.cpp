//
//  sparseSuperMatrixLinkage.cpp
//  VegaFEM
//
//  Created by Damian Stewart on 19/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include "sparseSuperMatrixLinkage.h"

SparseSuperMatrixLinkage::SparseSuperMatrixLinkage(shared_ptr<SparseMatrix> matrix_, const vector<int>& fixedRows, const vector<int>& fixedColumns, shared_ptr<SparseMatrix> superMatrix_, bool oneIndexed)
: matrix(matrix_), superMatrix(superMatrix_)/*, indexRemapper(superMatrix, matrix)*/
{
    BuildSuperMatrixIndices(fixedRows, fixedColumns, oneIndexed);
    
    //RemoveRowsColumnsFromIndexRemapper(fixedRows, fixedColumns);
}

static void BuildRenumberingVector(int nConstrained, int nSuper, const vector<int>& fixedDOFs, vector<int>& superDOFs, int oneIndexed)
{
    // superRows[i] is the row index in the super matrix corresponsing to row i of constrained matrix
    superDOFs.resize(nConstrained);
    int constrainedDOF = 0;
    int superDOF = 0;
    int numFixedDOFs = fixedDOFs.size();
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
            superDOFs[constrainedDOF] = superDOF;
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
        superDOFs[constrainedDOF] = superDOF;
        
        constrainedDOF++;
        superDOF++;
    }
}

void SparseSuperMatrixLinkage::BuildSuperMatrixIndices(const vector<int>& fixedRows, const vector<int>& fixedColumns, bool oneIndexed)
{
    int numFixedColumns = (int)fixedColumns.size();
    int numFixedRows = (int)fixedRows.size();
    int numSuperColumns = superMatrix->GetNumColumns();
    int numColumns = numSuperColumns - numFixedColumns;
    
    if ((matrix->GetNumRows() + numFixedRows != superMatrix->GetNumRows()) || (matrix->GetNumColumns() + numFixedColumns > numSuperColumns) )
    {
        printf("Error in BuildSuperMatrixIndices: number of constrained DOFs does not match the size of the two matrices.\n");
        printf("num rows: %d num fixed rows in super matrix: %d num rows in super matrix: %d\n", matrix->GetNumRows(), numFixedRows, superMatrix->GetNumRows());
        printf("num columns: %d num fixed columns in super matrix: %d num columns in super matrix: %d\n", numColumns, numFixedColumns, numSuperColumns);
        exit(1);
    }
    
    // build row renumbering function:
    BuildRenumberingVector(matrix->GetNumRows(), superMatrix->GetNumRows(), fixedRows, superRows, oneIndexed?1:0);
    // build column renumbering function:
    vector<int> superColumns_;
    BuildRenumberingVector(numColumns, numSuperColumns, fixedColumns, superColumns_, oneIndexed?1:0);
    
    // superRows[i] is the row index in the super matrix corresponsing to row i of constrained matrix
    // superColumns_[i] is the dense column index in the super matrix corresponsing to the dense column i of constrained matrix
    
    // build column indices
    superMatrixIndices.resize(matrix->GetNumRows());
    const auto& columnIndices = matrix->GetColumnIndices();
    for(int i=0; i < matrix->GetNumRows(); i++)
    {
        int rowLength = matrix->GetRowLength(i);
        superMatrixIndices[i].resize(rowLength);
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
    
    
    int numRows = superRows.size();
    printf("sparseSuperMatrixLinkage:");
    printf("subMatrix Row index -> superMatrix row index: superMatrixColumnIndices\n");
    for (int row=0; row<numRows; row++)
    {
        printf("%2i -> %2i:", row, superRows[row]);
        for (int j=0; j<superMatrixIndices[row].size(); j++)
        {
            printf("%2i ", superMatrixIndices[row][j]);
        }
        printf("\n");
    }
    
}

void SparseSuperMatrixLinkage::RemoveRowsColumnsFromIndexRemapper(SparseMatrixIndexRemapper& indexRemapper, const vector<int>& rowsToRemove, const vector<int>& columnsToRemove)
{
    auto rowsToRemoveSorted = rowsToRemove;
    std::sort(rowsToRemoveSorted.begin(), rowsToRemoveSorted.end());
    // iterate in reverse order to avoid having to compensate for already deleted earlier rows
    for (auto it = rowsToRemoveSorted.rbegin(); it != rowsToRemoveSorted.rend(); ++it)
    {
        indexRemapper.RemoveSuperRowFromSubMatrix(*it);
    }
    
    auto columnsToRemoveSorted = columnsToRemove;
    std::sort(columnsToRemoveSorted.begin(), columnsToRemoveSorted.end());
    for (auto it = columnsToRemoveSorted.rbegin(); it != columnsToRemoveSorted.rend(); ++it)
    {
        printf("removing column %i\n", *it);
        indexRemapper.RemoveSuperColumnFromSubMatrix(*it);
    }
}

void SparseSuperMatrixLinkage::AssignFromSuperMatrix()
{
    auto& columnEntries = matrix->GetDataHandle();
    const auto& superMatrixColumnEntries = superMatrix->GetDataHandle();
    for(int i=0; i<matrix->GetNumRows(); i++)
    {
        const vector<double>& row = superMatrixColumnEntries[superRows[i]];
        const vector<int>& indices = superMatrixIndices[i];
        int rowLength = matrix->GetRowLength(i);
        for(int j=0; j < rowLength; j++)
            columnEntries[i][j] = row[indices[j]];
    }
}

