/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "massSpringSystem" library, Copyright (C) 2007 CMU, 2009 MIT,         *
 *                                           2014 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic, Daniel Schroeder                          *
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

#ifndef MASSSPRINGSYSTEMFROMTETMESHCONFIGFILE_H_
#define MASSSPRINGSYSTEMFROMTETMESHCONFIGFILE_H_

#include "massSpringSystem.h"

class MassSpringSystemTetMeshConfiguration
{
public:
  MassSpringSystemTetMeshConfiguration() : tetMeshFilename(NULL), surfaceMeshFilename(NULL) {}
  ~MassSpringSystemTetMeshConfiguration() { free(tetMeshFilename); free(surfaceMeshFilename); }

  char * tetMeshFilename;
  char * surfaceMeshFilename;
  double density;
  double tensileStiffness;
  double damping;
  int addGravity;
};

class MassSpringSystemFromTetMeshConfigFile
{
public:

  // generates a mass spring system from the given configuration file
  // if massSpringSystemConfiguration is not NULL, it also returns the mass spring parameters read from the file (massSpringSystemConfiguration is output parameter only)
  int GenerateMassSpringSystem(const char * configFilename, MassSpringSystem ** massSpringSystem, MassSpringSystemTetMeshConfiguration * massSpringSystemTetMeshConfiguration = NULL);

protected:
  char * DuplicateString(const char*);
};

#endif

