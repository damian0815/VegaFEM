//
//  main.cpp
//  testSparseMatrixIndexRemapper
//
//  Created by Damian Stewart on 19/08/15.
//  Copyright (c) 2015 Damian Stewart. All rights reserved.
//

#include <iostream>
#include <limits>

#include "sparseMatrix.h"
#include "sparseMatrixOutline.h"
#include "sparseMatrixIndexRemapper.h"
#include "sparseSubMatrixLinkage.h"


using namespace std;

static double getRandomDouble()
{
	return double(int(((double)arc4random()/numeric_limits<uint32_t>::max())*100));
}

static shared_ptr<SparseMatrix> buildRandomSparseMatrix(int minSize=3, int maxSize=25)
{
	int rowColumnCount = (arc4random_uniform(maxSize-minSize)+minSize)*3;
	SparseMatrixOutline outline(rowColumnCount);
	
	for (int i=0; i<arc4random_uniform(rowColumnCount/2)+3; i++)
	{
		int j = arc4random_uniform(rowColumnCount);
		int k = arc4random_uniform(rowColumnCount);
		if (!outline.HasEntry(j, k))
		{
			outline.AddEntry(j, k);
		}
	}
	
	auto matrix = make_shared<SparseMatrix>(outline);
	
	for (int i=0; i<rowColumnCount; i++)
	{
		for (int j=0; j<matrix->GetRowLength(i); j++)
		{
			matrix->SetEntry(i, j, getRandomDouble());
		}
	}
	
	return matrix;
}


static void insertRandomEntry(shared_ptr<SparseMatrix> matrix)
{
	auto range = matrix->GetNumRows();
	int j, k;
	while (true) {
		j = arc4random_uniform(range);
		if (matrix->GetRowLength(j) == 0) {
			continue;
		}
		k = arc4random_uniform(range);
		
		if (matrix->GetInverseIndex(j, k) == -1) {
			break;
		}
	}
	
	auto sparseK = matrix->InsertNewEntry(j, k);
	matrix->SetEntry(j, sparseK, -getRandomDouble());
}

struct SparseMatrixMapping
{
	SparseMatrixMapping(SparseMatrixIndexRemapper& remapper)
	{
		auto superMatrix = remapper.GetSuperMatrix();
		auto subMatrix = remapper.GetSubMatrix();
		mSuperDenseToSubDenseColumnMappings.resize(superMatrix->GetNumRows());
		
		for (int superRow=0; superRow < superMatrix->GetNumRows(); superRow++) {
			if (remapper.HasSubMatrixRowForSuperMatrixRow(superRow)) {
				int subRow = remapper.GetSubMatrixRowForSuperMatrixRow(superRow);
				mRowMapping[superRow] = subRow;
				
				for(int superSparseColumn=0; superSparseColumn < superMatrix->GetRowLength(superRow); superSparseColumn++) {
					if (remapper.HasSubMatrixSparseColumnForSuperMatrixSparseColumn(superRow, superSparseColumn)) {
						int subSparseColumn = remapper.GetSubMatrixSparseColumnForSuperMatrixSparseColumn(superRow, superSparseColumn);
						
						int superDenseColumn = superMatrix->GetColumnIndex(superRow, superSparseColumn);
						int subDenseColumn = subMatrix->GetColumnIndex(subRow, subSparseColumn);
						
						mSuperDenseToSubDenseColumnMappings[superRow][superDenseColumn] = subDenseColumn;
					}
				}
			}
		}
	}
	

	void Print() const
	{
		for (auto rowMap: mRowMapping) {
			int superRow = rowMap.first;
			int subRow = rowMap.second;
			if (mSuperDenseToSubDenseColumnMappings.at(superRow).size()) {
				printf("%3i->%3i: ", superRow, subRow);
				for (auto columnMap: mSuperDenseToSubDenseColumnMappings.at(superRow)) {
					printf("%3i:%3i ", columnMap.first, columnMap.second);
				}
				printf("\n");
			}
		}
	}
	
	map<int,int> mRowMapping;
	vector<map<int,int>> mSuperDenseToSubDenseColumnMappings; // one for each row in the super matrix
	
	bool operator==(const SparseMatrixMapping& other) const
	{
		return mRowMapping == other.mRowMapping && mSuperDenseToSubDenseColumnMappings == other.mSuperDenseToSubDenseColumnMappings;
	}
	
};


static void PrintLinkageInfo(shared_ptr<SparseMatrix> superMatrix, shared_ptr<SparseMatrix> subMatrix, const SparseMatrixMapping& mapping)
{
	superMatrix->Print(true);
	subMatrix->Print(true);
	mapping.Print();
}

static bool testBasicSubmatrix(bool noisy)
{
	auto matrix = buildRandomSparseMatrix();
	auto subMatrix = make_shared<SparseMatrix>(*matrix);
	
	auto linkage = matrix->AttachSubMatrix(subMatrix);
	
	SparseMatrixMapping preMapping(linkage->GetIndexRemapper());
	if (noisy) {
		PrintLinkageInfo(matrix, subMatrix, preMapping);
		cout << endl;
	}
	
	for (int i=0; i<arc4random_uniform(matrix->GetNumRows()/3); i++) {
		insertRandomEntry(matrix);
	}
	SparseMatrixMapping postMapping(linkage->GetIndexRemapper());
	
	if (noisy) {
		PrintLinkageInfo(matrix, subMatrix, postMapping);
	}
	
	return postMapping == preMapping;
}


static bool testOffsetSubmatrix(bool noisy)
{
	auto subMatrix = buildRandomSparseMatrix();
	auto superMatrix = make_shared<SparseMatrix>(subMatrix->GetNumRows()*2);
	
	auto offset = arc4random_uniform(subMatrix->GetNumRows());
	superMatrix->CreateEntriesIfNecessary(subMatrix->GetTopology(), offset);
	auto linkage = superMatrix->AttachSubMatrix(subMatrix, offset);
	
	SparseMatrixMapping preMapping(linkage->GetIndexRemapper());
	if (noisy) {
		PrintLinkageInfo(superMatrix, subMatrix, preMapping);
		cout << endl;
	}
	
	for (int i=0; i<arc4random_uniform(superMatrix->GetNumRows()/3); i++) {
		insertRandomEntry(superMatrix);
	}
	SparseMatrixMapping postMapping(linkage->GetIndexRemapper());
	
	if (noisy) {
		PrintLinkageInfo(superMatrix, subMatrix, postMapping);
		cout << endl;
		cout << endl;
	}
	
	return postMapping == preMapping;
}

static bool testAppendSubmatrix(bool noisy)
{
	auto superMatrix = make_shared<SparseMatrix>(0);
	
	vector<shared_ptr<SparseMatrix>> subMatrices;
	vector<shared_ptr<SparseSubMatrixLinkage>> linkages;
	vector<shared_ptr<SparseMatrixMapping>> preMappings;
	
	for (int i=0; i<arc4random_uniform(3)+3; i++) {
		auto subMatrix = buildRandomSparseMatrix(10,20);
		subMatrices.push_back(subMatrix);
		
		auto offset = superMatrix->GetNumRows();
		superMatrix->Append(subMatrix.get());
		
		auto linkage = superMatrix->AttachSubMatrix(subMatrix, offset);
		linkages.push_back(linkage);
	}
	
	for (auto linkage: linkages) {
		preMappings.push_back(make_shared<SparseMatrixMapping>(linkage->GetIndexRemapper()));
	}
	
	
	for (int i=0; i<arc4random_uniform(10); i++) {
		insertRandomEntry(superMatrix);
	}
	
	bool success = true;
	for (int i=0; i<linkages.size(); i++) {
		auto linkage = linkages[i];
		SparseMatrixMapping postMapping(linkage->GetIndexRemapper());
		
		if ( !(preMappings[i]->operator==(postMapping) ) )
		{
			cout << "mapping " << i << " test failed" << endl;
			success = false;
		}
	}
	
	return success;
}

int main(int argc, const char * argv[])
{
	cout << "************" << endl;
	cout << "basic" << endl;
	cout << "************" << endl;
	
	for (int i=0; i<10000; i++) {
		bool success = testBasicSubmatrix((i<5));
		if (!success) {
			cout << "testbasicSubmatrix failed!" << endl;
		}
	}
	
	cout << "************" << endl;
	cout << "offset" << endl;
	cout << "************" << endl;
	
	for (int i=0; i<10000; i++) {
		bool success = testOffsetSubmatrix((i<5));
		if (!success) {
			cout << "testBasicSubmatrix failed!" << endl;
		}
	}
	
	cout << "************" << endl;
	cout << "append" << endl;
	cout << "************" << endl;
	
	for (int i=0; i<10000; i++) {
		bool success = testAppendSubmatrix((i<5));
		if (!success) {
			cout << "testAppendSubmatrix failed!" << endl;
		}
	}
	
	
	return 0;
}
