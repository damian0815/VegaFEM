//
//  sparseSubMatrixLinkage.cpp
//  VegaFEM
//
//  Created by Damian Stewart on 19/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include "sparseSubMatrixLinkage.h"

SparseSubMatrixLinkage::SparseSubMatrixLinkage(shared_ptr<SparseMatrix> matrix, shared_ptr<SparseMatrix> subMatrix_, int denseRowColumnOffset)
: superMatrix(matrix), subMatrix(subMatrix_), indexRemapper(superMatrix, subMatrix, denseRowColumnOffset)
{
}


void SparseSubMatrixLinkage::AddSubMatrixToSuperMatrix(double factor)
{
    auto& superColumnEntries = superMatrix->GetDataHandle();
    const auto& subColumnEntries = subMatrix->GetDataHandle();
    
    for(int subRow=0; subRow<subMatrix->GetNumRows(); subRow++) {
        int superRow = indexRemapper.GetSuperMatrixRowForSubMatrixRow(subRow);
        int subMatrixRowLength = subMatrix->GetRowLength(subRow);
        for(int sparseSubJ=0; sparseSubJ < subMatrixRowLength; sparseSubJ++) {
            int sparseSuperJ = indexRemapper.GetSuperMatrixSparseColumnForSubMatrixSparseColumn_SubMatrixRow(subRow, sparseSubJ);
            superColumnEntries.at(superRow).at(sparseSuperJ) += factor * subColumnEntries.at(subRow).at(sparseSubJ);
        }
    }
}