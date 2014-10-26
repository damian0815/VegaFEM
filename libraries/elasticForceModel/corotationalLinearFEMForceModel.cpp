/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "forceModel" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC     *
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

#include "corotationalLinearFEMForceModel.h"

CorotationalLinearFEMForceModel::CorotationalLinearFEMForceModel(CorotationalLinearFEM * corotationalLinearFEM_, int warp_): corotationalLinearFEM(corotationalLinearFEM_), warp(warp_)
{
  r = 3 * corotationalLinearFEM->GetTetMesh()->getNumVertices();
}

CorotationalLinearFEMForceModel::~CorotationalLinearFEMForceModel() {}

void CorotationalLinearFEMForceModel::GetInternalForce(double * u, double * internalForces)
{
  corotationalLinearFEM->ComputeForceAndStiffnessMatrix(u, internalForces, NULL, warp);
}

void CorotationalLinearFEMForceModel::GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix)
{
  corotationalLinearFEM->GetStiffnessMatrixTopology(tangentStiffnessMatrix);
}

void CorotationalLinearFEMForceModel::GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix)
{
  corotationalLinearFEM->ComputeForceAndStiffnessMatrix(u, NULL, tangentStiffnessMatrix, warp);
} 

void CorotationalLinearFEMForceModel::GetForceAndMatrix(double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix)
{
  corotationalLinearFEM->ComputeForceAndStiffnessMatrix(u, internalForces, tangentStiffnessMatrix, warp);
}

