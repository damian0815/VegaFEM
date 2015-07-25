//
//  sparseSuperMatrixLinkage.h
//  VegaFEM
//
//  Created by Damian Stewart on 19/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#ifndef __VegaFEM__sparseSuperMatrixLinkage__
#define __VegaFEM__sparseSuperMatrixLinkage__

#include <iostream>
#include "sparseMatrix.h"
#include "sparseMatrixIndexRemapper.h"

using std::shared_ptr;
using std::vector;

class SparseSuperMatrixLinkage
{
public:
    SparseSuperMatrixLinkage(shared_ptr<SparseMatrix> superMatrix, shared_ptr<SparseMatrix> subMatrix);
    
    shared_ptr<SparseMatrix> GetSuperMatrix() { return indexRemapper.GetSuperMatrix(); }
    shared_ptr<SparseMatrix> GetSubMatrix() { return indexRemapper.GetSubMatrix(); }
    
    void AssignSubMatrixFromSuperMatrix() { indexRemapper.AssignSubMatrixFromSuperMatrix(); }
    
    SparseMatrixIndexRemapper& GetIndexRemapper() { return indexRemapper; }
    
private:
    SparseMatrixIndexRemapper indexRemapper;
};

#endif /* defined(__VegaFEM__sparseSuperMatrixLinkage__) */
