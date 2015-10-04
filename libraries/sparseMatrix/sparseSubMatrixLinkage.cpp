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


