//
//  sparseSubMatrixLinkage.h
//  VegaFEM
//
//  Created by Damian Stewart on 19/07/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#ifndef __VegaFEM__sparseSubMatrixLinkage__
#define __VegaFEM__sparseSubMatrixLinkage__

#include <iostream>
#include "sparseMatrix.h"
#include "sparseMatrixIndexRemapper.h"

using std::shared_ptr;

class SparseSubMatrixLinkage
{
public:
    SparseSubMatrixLinkage(shared_ptr<SparseMatrix> matrix, shared_ptr<SparseMatrix> subMatrix);
    
    shared_ptr<SparseMatrix> GetSuperMatrix() { return superMatrix; }
    shared_ptr<SparseMatrix> GetSubMatrix() { return subMatrix; }
    
    void AddSubMatrixToSuperMatrix(double factor);
    
private:
    void BuildSubMatrixIndices();
    
    
    shared_ptr<SparseMatrix> superMatrix;
    shared_ptr<SparseMatrix> subMatrix;
    
       /*
     length(subMatrixIndices) == number of rows
     */
    //vector<vector<int> > subMatrixIndices;
    SparseMatrixIndexRemapper indexRemapper;
    
};

#endif /* defined(__VegaFEM__sparseSubMatrixLinkage__) */
