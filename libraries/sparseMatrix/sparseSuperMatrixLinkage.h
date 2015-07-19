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

using std::shared_ptr;
using std::vector;

class SparseSuperMatrixLinkage
{
public:
    SparseSuperMatrixLinkage(shared_ptr<SparseMatrix> matrix, const vector<int>& fixedRows, const vector<int>& fixedColumns, shared_ptr<SparseMatrix> superMatrix, bool oneIndexed=false);
    
    shared_ptr<SparseMatrix> GetSuperMatrix() { return superMatrix; }
    shared_ptr<SparseMatrix> GetSubMatrix() { return matrix; }
    
    void AssignFromSuperMatrix();
    
private:
    void BuildSuperMatrixIndices(const vector<int>& fixedRows, const vector<int>& fixedColumns, bool oneIndexed);
    
    shared_ptr<SparseMatrix> matrix;
    shared_ptr<SparseMatrix> superMatrix;
    
       /*
     length(subMatrixIndices) == number of rows
     */
    vector<vector<int> > superMatrixIndices;
    vector<int> superRows;
};

#endif /* defined(__VegaFEM__sparseSuperMatrixLinkage__) */
