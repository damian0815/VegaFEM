//
//  sparseSubMatrixLinkage.cpp
//  VegaFEM
//
//  Created by Damian Stewart on 19/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include "sparseSubMatrixLinkage.h"

SparseSubMatrixLinkage::SparseSubMatrixLinkage(shared_ptr<SparseMatrix> matrix, shared_ptr<SparseMatrix> subMatrix_)
: superMatrix(matrix), subMatrix(subMatrix_)
{
    BuildSubMatrixIndices();
}


void SparseSubMatrixLinkage::BuildSubMatrixIndices()
{
    subMatrixIndices.clear();
    subMatrixIndices.resize(superMatrix->GetNumRows());
    
    for(int i=0; i<superMatrix->GetNumRows(); i++)
    {
        int submatrixRowLength = subMatrix->GetRowLength(i);
        subMatrixIndices[i].resize(submatrixRowLength);
        const vector<int>& indices = subMatrix->GetColumnIndices()[i];
        
        for(int j=0; j < submatrixRowLength; j++)
        {
            // finds the position in row i of element with column index jDense
            // int GetInverseIndex(int i, int jDense);
            subMatrixIndices[i][j] = superMatrix->GetInverseIndex(i, indices[j]);
            if (subMatrixIndices[i][j] == -1)
            {
                printf("Error (BuildSubMatrixIndices): given matrix is not a submatrix of this matrix. The following index does not exist in this matrix: (%d,%d)\n", i, indices[j]);
                exit(1);
            }
        }
    }
}

void SparseSubMatrixLinkage::AddSubMatrixToSuperMatrix(double factor)
{
    for(int i=0; i<superMatrix->GetNumRows(); i++)
    {
        const auto& indices = subMatrixIndices[i];
        int subMatrixRowLength = subMatrix->GetRowLength(i);
        auto& superColumnEntries = superMatrix->GetDataHandle();
        const auto& subColumnEntries = subMatrix->GetDataHandle();
        for(int j=0; j < subMatrixRowLength; j++)
            superColumnEntries[i][indices[j]] += factor * subColumnEntries[i][j];
    }
}