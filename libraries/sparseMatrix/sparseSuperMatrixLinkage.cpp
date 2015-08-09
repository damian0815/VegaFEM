//
//  sparseSuperMatrixLinkage.cpp
//  VegaFEM
//
//  Created by Damian Stewart on 19/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include "sparseSuperMatrixLinkage.h"

SparseSuperMatrixLinkage::SparseSuperMatrixLinkage(shared_ptr<SparseMatrix> superMatrix, shared_ptr<SparseMatrix> subMatrix, int rowColumnOffset)
: indexRemapper(superMatrix, subMatrix, rowColumnOffset)
{
}


