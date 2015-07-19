/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "integrator" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC     *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integratorBaseSparse.h"
#include "sparseMatrixOutline.h"

IntegratorBaseSparse::IntegratorBaseSparse(int r, double timestep, SparseMatrix * massMatrix_, ForceModel * forceModel_, int numConstrainedDOFs_, int * constrainedDOFs_, double dampingMassCoef, double dampingStiffnessCoef): IntegratorBase(r, timestep, dampingMassCoef, dampingStiffnessCoef), massMatrix(massMatrix_), forceModel(forceModel_), numConstrainedDOFs(numConstrainedDOFs_)
{
    systemSolveTime = 0.0;
    forceAssemblyTime = 0.0;
    
    constrainedDOFs = (int*) malloc (sizeof(int) * numConstrainedDOFs);
    memcpy(constrainedDOFs, constrainedDOFs_, sizeof(int) * numConstrainedDOFs);
    
    SparseMatrixOutline outline(r);
    auto outlinePtr = &outline;
    dampingMatrix = std::make_shared<SparseMatrix>(outlinePtr);
}

IntegratorBaseSparse::~IntegratorBaseSparse()
{
    free(constrainedDOFs);
}

void IntegratorBaseSparse::SetDampingMatrix(std::shared_ptr<SparseMatrix> dampingMatrix_)
{
    dampingMatrix = dampingMatrix_;
}

double IntegratorBaseSparse::GetKineticEnergy()
{
    return 0.5 * massMatrix->QuadraticForm(qvel);
}

double IntegratorBaseSparse::GetTotalMass()
{
    return massMatrix->SumEntries();
}

