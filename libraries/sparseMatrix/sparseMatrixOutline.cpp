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

SparseMatrixOutline::SparseMatrixOutline(int numRows_)
{
    Allocate(numRows_);
}

SparseMatrixOutline::~SparseMatrixOutline()
{
    columnEntries.clear();
}

SparseMatrixOutline::SparseMatrixOutline(int numRows_, double diagonal)
{
    Allocate(numRows_);
    
    pair<int,double> entry;
    
    for(int i=0; i<GetNumRows(); i++)
    {
        entry.first = i;
        entry.second = diagonal;
        columnEntries[i].insert(entry);
    }
}

SparseMatrixOutline::SparseMatrixOutline(int numRows_, double * diagonal)
{
    Allocate(numRows_);
    pair<int,double> entry;
    for(int i=0; i<GetNumRows(); i++)
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
    
    int numRows = expand * m1;
    
    printf("Loading matrix from %s... Size is %d x %d .\n", filename, numRows, expand * n1);fflush(NULL);
    
    Allocate(numRows);
    
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

void SparseMatrixOutline::Allocate(int numRows)
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
}

void SparseMatrixOutline::AddEntry(int i, int j, double value)
{
    map<int,double>::iterator pos = columnEntries.at(i).find(j);
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

void SparseMatrixOutline::AddBlock3x3Entry(int i, int j)
{
    for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
            AddEntry(3*i+k, 3*j+l);
}

void SparseMatrixOutline::AppendEntries(const SparseMatrixOutline& other)
{
    int offset = GetNumRows();
    IncreaseNumRows(other.GetNumRows());
    
    int otherNumColumns = other.GetNumColumns();
    for (int row=0; row<other.GetNumRows(); row++)
    {
        for (int column=0; column<otherNumColumns; column++)
        {
            if (other.HasEntry(row, column))
            {
                AddEntry(row+offset, column+offset);
            }
        }
    }
}


double SparseMatrixOutline::GetEntry(int i, int j) const
{
    map<int,double>::const_iterator pos = columnEntries[i].find(j);
    if (pos != columnEntries[i].end())
        return (pos->second);
    else
        return 0;
}

vector<pair<int, int>> SparseMatrixOutline::GetEntries() const // return all entries as pairs (i,j)
{
    vector<pair<int,int>> entries;
    for (int row=0; row<GetNumRows(); row++) {
        for (const auto& column: columnEntries[row]) {
            entries.push_back(make_pair(row, column.first));
        }
    }
    return entries;
}

bool SparseMatrixOutline::HasEntry(int i, int j) const
{
	auto pos = columnEntries[i].find(j);
	return (pos != columnEntries[i].end());
}

int SparseMatrixOutline::GetNumColumns() const
{
    int numColumns = -1;
    for(int i=0; i<GetNumRows(); i++)
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
    
    fprintf(fout, "%d\n%d\n", GetNumRows(), GetNumColumns());
    for(int i=0; i<GetNumRows(); i++)
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
    for (int i=0; i<GetNumRows(); i++)
    {
        for (int j=0; j<GetNumRows(); j++)
            printf("%f ",GetEntry(i,j));
        printf("\n");
    }
}

void SparseMatrixOutline::PrintSparse() const
{
    for (int i=0; i<GetNumRows(); i++)
    {
        const auto& row = columnEntries[i];
        if (row.size()>0)
        {
            printf("%3i: ", i);
            for (auto it: row)
            {
                printf("%3i,", it.first);
            }
            printf("\n");
        }
    }
}

static void PrintMultiple(const char* toPrint, int count)
{
	for (int i=0; i<count; i++)
	{
		printf(toPrint);
	}
}

void SparseMatrixOutline::PrintTopology(int clusterSize) const
{
	int numColumns = GetNumColumns();
	
	printf("+-");
	PrintMultiple("-", (numColumns*2)/clusterSize);
	printf("+");
	printf("\n");
	
	for (int i=0; i<GetNumRows(); i+=clusterSize)
	{
		printf("| ");
		for (int j=0; j<numColumns; j+=clusterSize)
		{
			if (HasEntry(i,j))
			{
				printf("* ");
			}
			else
			{
				printf("  ");
			}
		}
		printf("|");
		printf("\n");
	}
	
	printf("+-");
	PrintMultiple("-", (numColumns*2)/clusterSize);
	printf("+");
	
	printf("\n");
}

int SparseMatrixOutline::GetNumEntries() const
{
    int num = 0;
    for(int i=0; i<GetNumRows(); i++)
        num += columnEntries[i].size();
    return num;
}
