/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.0                               *
 *                                                                       *
 * "objMeshGPUDeformer" library , Copyright (C) 2007 CMU, 2009 MIT,      *
 *                                                        2013 USC       *
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

#ifndef _OBJMESHGPUDEFORMER_COARSETOFINE_FBO_H_
#define _OBJMESHGPUDEFORMER_COARSETOFINE_FBO_H_

#include "objMeshGPUDeformer_coarseToFine.h"

class ObjMeshGPUDeformer_coarseToFine_fbo : public ObjMeshGPUDeformer_coarseToFine
{
public:

  virtual ~ObjMeshGPUDeformer_coarseToFine_fbo();

protected:
  GLuint fbo;

  virtual void EnableRTT();
  virtual void DisableRTT();
  virtual void BindRT();
  virtual void UnbindRT();
  virtual int InitRTT();

  virtual void * GetDerivedData() { return &fbo; }
  virtual void SetDerivedData(void * data);

  void BindDynamicTextureToFBO();
  void CheckFramebufferStatus();
  int InitExtensions();     
};

#endif

