/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "sparseMatrix" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC   *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "sparseMatrix.h"
using namespace std;

SparseMatrixOutline::SparseMatrixOutline(int numRows_): numRows(numRows_)
{
    Allocate();
}

SparseMatrixOutline::~SparseMatrixOutline()
{
    // deallocate column entries
    for(int i=0; i<numRows; i++)
        columnEntries[i].clear();
    columnEntries.clear();
}

SparseMatrixOutline::SparseMatrixOutline(int numRows_, double diagonal): numRows(numRows_)
{
    Allocate();
    
    pair<int,double> entry;
    
    for(int i=0; i<numRows; i++)
    {
        entry.first = i;
        entry.second = diagonal;
        columnEntries[i].insert(entry);
    }
}

SparseMatrixOutline::SparseMatrixOutline(int numRows_, double * diagonal): numRows(numRows_)
{
    Allocate();
    pair<int,double> entry;
    for(int i=0; i<numRows; i++)
    {
        entry.first = i;
        entry.second = diagonal[i];
        columnEntries[i].insert(entry);
    }
}

SparseMatrixOutline::SparseMatrixOutline(const char * filename, int expand)
{
    if (expand <= 0)
    {
        printf("Error: invalid expand factor %d in SparseMatrixOutline constructor.\n", expand);
        throw 1;
    }
    
    FILE * inputMatrix = fopen(filename,"r");
    if (!inputMatrix)
    {
        printf("Error: couldn't open input sparse matrix file %s.\n", filename);
        throw 2;
    }
    
    // read input size
    int m1, n1;
    if (fscanf(inputMatrix, "%d\n%d\n", &m1, &n1) < 2)
    {
        printf("Error: could not read sparse matrix dimensions in file %s.\n", filename);
        throw 3;
    }
    
    numRows = expand * m1;
    
    printf("Loading matrix from %s... Size is %d x %d .\n", filename, numRows, expand * n1);fflush(NULL);
    
    Allocate();
    
    char s[4096];
    while (fgets(s,4096,inputMatrix) != NULL)
    {
        int i1,j1;
        double x;
        sscanf(s, "%d %d %lf\n", &i1, &j1, &x);
        for(int e=0; e<expand; e++)
            AddEntry(expand * i1 + e, expand * j1 + e, x);
    }
    
    fclose(inputMatrix);
}

void SparseMatrixOutline::Allocate()
{
    // allocate empty datastructure for row entries
    columnEntries.clear();
    map<int,double> emptyMap;
    for (int i=0; i<numRows; i++)
        columnEntries.push_back(emptyMap);
}

void SparseMatrixOutline::IncreaseNumRows(int numAddedRows)
{
    map<int,double> emptyMap;
    for(int i=0; i<numAddedRows; i++)
        columnEntries.push_back(emptyMap);
    
    numRows += numAddedRows;
}

void SparseMatrixOutline::AddEntry(int i, int j, double value)
{
    map<int,double>::iterator pos = columnEntries[i].find(j);
    if (pos != columnEntries[i].end())
        pos->second += value;
    else
    {
        pair<int,double> entry(j,value);
        columnEntries[i].insert(entry);
    }
}

// add a block matrix, starting at row i, and column j
void SparseMatrixOutline::AddBlockMatrix(int iStart, int jStart, const SparseMatrix * block, double scalarFactor)
{
    int nBlock = block->GetNumRows();
    for(int i=0; i<nBlock; i++)
    {
        int rowLength = block->GetRowLength(i);
        for(int j=0; j<rowLength; j++)
            AddEntry(iStart + i, jStart + block->GetColumnIndex(i,j), scalarFactor * block->GetEntry(i,j));
    }
}

void SparseMatrixOutline::MultiplyRow(int row, double scalar)
{
    for(map<int,double>::iterator iter = columnEntries[row].begin(); iter != columnEntries[row].end(); iter++)
        iter->second *= scalar;
}

void SparseMatrixOutline::AddBlock3x3Entry(int i, int j, double * matrix3x3)
{
    for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
            AddEntry(3*i+k,3*j+l,matrix3x3[3*k+l]);
}

double SparseMatrixOutline::GetEntry(int i, int j) const
{
    map<int,double>::const_iterator pos = columnEntries[i].find(j);
    if (pos != columnEntries[i].end())
        return (pos->second);
    else
        return 0;
}

int SparseMatrixOutline::GetNumColumns() const
{
    int numColumns = -1;
    for(int i=0; i<numRows; i++)
    {
        map<int,double>::const_iterator j1;
        // traverse all row entries
        for(j1=columnEntries[i].begin(); j1 != columnEntries[i].end(); j1++)
            if (j1->first > numColumns)
                numColumns = j1->first;
    }
    return numColumns + 1;
}

int SparseMatrixOutline::Save(const char * filename, int oneIndexed) const
{
    FILE * fout = fopen(filename, "w");
    if (!fout)
        return 1;
    
    fprintf(fout, "%d\n%d\n", numRows, GetNumColumns());
    for(int i=0; i<numRows; i++)
    {
        map<int,double>::const_iterator j1;
        // traverse all row entries
        for(j1=columnEntries[i].begin(); j1 != columnEntries[i].end(); ++j1)
            fprintf(fout,"%d %d %.15f\n",i,j1->first + oneIndexed, j1->second + oneIndexed);
        
    }
    fclose(fout);
    
    return 0;
}

void SparseMatrixOutline::Print() const
{
    for (int i=0; i<numRows; i++)
    {
        for (int j=0; j<numRows; j++)
            printf("%f ",GetEntry(i,j));
        printf("\n");
    }
}

int SparseMatrixOutline::GetNumEntries() const
{
    int num = 0;
    for(int i=0; i<numRows; i++)
        num += columnEntries[i].size();
    return num;
}
