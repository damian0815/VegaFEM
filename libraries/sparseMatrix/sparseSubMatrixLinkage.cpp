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
    //const auto& subToSuperIndicesAllRows = indexRemapper.GetsubMatrixSparseToSuperMatrixSparseColumnMaps();
    auto& superColumnEntries = superMatrix->GetDataHandle();
    const auto& subColumnEntries = subMatrix->GetDataHandle();
    
    for(int row=0; row<superMatrix->GetNumRows(); row++)
    {
        //const auto& subToSuperIndices = subToSuperIndicesAllRows[row];
        //int subMatrixRowLength = subMatrix->GetRowLength(row);
        int subMatrixRowLength = subMatrix->GetRowLength(row);
        for(int sparseSubJ=0; sparseSubJ < subMatrixRowLength; sparseSubJ++)
        {
            int sparseSuperJ = indexRemapper.GetSuperMatrixSparseColumnForSubMatrixSparseColumn(row, sparseSubJ);
            //int sparseSuperJ = subToSuperIndices[sparseSubJ];
            superColumnEntries.at(row).at(sparseSuperJ) += factor * subColumnEntries.at(row).at(sparseSubJ);
        }
    }
}