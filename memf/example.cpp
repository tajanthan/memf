#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <climits>
#include <crtdbg.h>	// order important
#endif

#include <cassert>
#include <fstream>
#include <iostream>

#include "MEMF.h"
#include "util/util.h"

using namespace std;

typedef unsigned int nodeid;
typedef short labelid;
typedef int captype;

static nodeid width;
static nodeid height;
static labelid labels;
static nodeid numNodes;
static captype *unaryPotentialArray;
static captype *binaryWeightArray;
static captype *binaryFunctionArray;
static captype **newUnaryPotentialArray;
static nodeid changedUnaries = 0;

void readUnaryPotential(char* unaryF)
{
	ifstream fs(unaryF);
	assert(fs);

	char c;
	nodeid ii;
	labelid l;
	captype val;
	nodeid count = 0;
	while (fs >> ii >> c >> l >> c >> val) {	// read comma separated values		
		unaryPotentialArray[ii * labels + l] = val;
		++count;
	}
	assert(count == numNodes * labels);
	fs.close();
}

void readBinaryPotential(char* binaryWF, char* binaryF)
{
	ifstream fs(binaryWF);
	assert(fs);

	nodeid numEdges = width * (height - 1) + height * (width - 1);

	char c;
	nodeid i, j;
	captype val;
	nodeid count = 0;
	while (fs >> i >> c >> j >> c >> val) {	// read comma separated values
		if (j == i + 1) binaryWeightArray[i] = val;	// right neighbor
		if (j == i + width) binaryWeightArray[i + numNodes] = val;	// down neighbor
		++count;
	}
	assert(count == numEdges);
	fs.close();

	fs.open(binaryF);
	assert(fs);
	labelid li, lj;
	count = 0;
	while (fs >> li >> c >> lj >> c >> val) {	// read comma separated values
		binaryFunctionArray[li * labels + lj] = val;
		++count;
	}
	assert(count == labels * labels);
	fs.close();
}

captype unaryPotential(nodeid i, labelid li)
{
	return unaryPotentialArray[i * labels + li];
}

captype binaryWeights(nodeid i, nodeid j)
{
	if (i + 1 == j) return binaryWeightArray[i];	// right
	else if (i + width == j) return binaryWeightArray[i + numNodes];	// down
	else if (j + 1 == i) return binaryWeightArray[j];	// left
	else if (j + width == i) return binaryWeightArray[j + numNodes];	// up
	else assert(0);
	return 0;
}

captype binaryFunction(labelid li, labelid lj)
{
	return binaryFunctionArray[li * labels + lj];
}

void readNewUnaryPotential(char* newUnaryF)
{
	ifstream fs(newUnaryF);
	assert(fs);

	char c;
	nodeid ii;
	labelid l;
	captype val;
	nodeid count = 0;
	fs >> changedUnaries;
	newUnaryPotentialArray = new captype*[changedUnaries]();
	while (fs >> ii >> c >> l >> c >> val) {	// read comma separated values		
		newUnaryPotentialArray[count] = new captype[3]();
		newUnaryPotentialArray[count][0] = ii;
		newUnaryPotentialArray[count][1] = l;
		newUnaryPotentialArray[count][2] = val;
		++count;
	}
	assert(changedUnaries == count);
	fs.close();
}

void createLocalArrays(nodeid w, nodeid h, labelid l, char* unaryF, char* binaryWF, char* binaryF, char* newUnaryF)
{
	width = w;
	height = h;
	labels = l;	
	numNodes = width * height;
	unaryPotentialArray = new captype[numNodes*labels]();
	binaryWeightArray = new captype[numNodes * 2]();	// right, down neighbor
	binaryFunctionArray = new captype[labels*labels]();

	readUnaryPotential(unaryF);
	readBinaryPotential(binaryWF, binaryF);

	if (newUnaryF) readNewUnaryPotential(newUnaryF);
}

void deleteLocalArrays()
{
	delete[] unaryPotentialArray;
	delete[] binaryWeightArray;
	delete[] binaryFunctionArray;
	for (nodeid i = 0; i < changedUnaries; ++i) delete[] newUnaryPotentialArray[i];
	delete[] newUnaryPotentialArray;
}

void errorFunction(char* msg)
{
	cout << endl << "Error!\n " << msg << endl;
	exit(1);
}

captype newUnaryPotential(nodeid ii, labelid li)
{
	for (nodeid i = 0; i < changedUnaries; ++i) {
		if (ii == (nodeid)newUnaryPotentialArray[i][0] && li == (labelid)newUnaryPotentialArray[i][1]) return newUnaryPotentialArray[i][2];
	}
	return unaryPotentialArray[ii * labels + li];
}

captype minimumEnergy(MEMF<nodeid, labelid, captype> &memf, int riter = 0)
{
	captype energy = 0;
	for (nodeid ii = 0; ii < numNodes; ++ii) {
		if (!riter) energy += unaryPotential(ii, memf.getLabel(ii));
		else energy += newUnaryPotential(ii, memf.getLabel(ii));

		nodeid i = ii / width;
		nodeid j = ii % width;
		if (j < width - 1) {	// right
			nodeid jj = ii + 1;
			energy += binaryWeights(ii, jj) * binaryFunction(memf.getLabel(ii), memf.getLabel(jj));
		}
		if (i < height - 1) {	// down
			nodeid jj = ii + width;
			energy += binaryWeights(ii, jj) * binaryFunction(memf.getLabel(ii), memf.getLabel(jj));
		}
	}
	return energy;
}

captype minimumCut(MEMF<nodeid, labelid, captype> &memf, char* logf, int riter = 0)
{
	ofstream fout;
	fout.open(logf, ofstream::app);
	fout << "\n# Minimum cut\n";
	for (nodeid ii = 0; ii < numNodes; ++ii) {
		fout << ii << ": " << memf.getLabel(ii) << endl;
	}

	captype energy = minimumEnergy(memf, riter);
	cout << "\nMinimum energy: " << energy << endl << endl;
	fout << "\nMinimum energy: " << energy << endl << endl;
	return energy;
}

int main(int argc, char **argv)
{
#ifdef _WIN32
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);	// dump memory leaks
#endif
		
	char* usage = " [width] [height] [labels] [unary potential file (V*L)] [binary weights file (E)] [binary potential file (L*L)] [new-unary-file]";
	if (argc < 7) {
		cout << "Error!" << endl;
		cout << "Usage: " << argv[0] << usage << endl;
		return 1;
	}
	
	// for dynamic MRF - need more description!!
	char* newUnaryF = NULL;
	if (argc > 7) newUnaryF = argv[7];

	createLocalArrays(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), argv[4], argv[5], argv[6], newUnaryF);	

	captype flow = 0, energy = 0;
	char* logf = (char*)"memf.out";
	MEMF<nodeid, labelid, captype> memf(width, height, labels, logf, errorFunction);
	memf.setEnergy(unaryPotential, binaryWeights, binaryFunction);

	if (!newUnaryF) {
		flow = memf.optimize();
		energy = minimumCut(memf, logf);
	}
	else {
		flow = memf.optimize(true);
		energy = minimumCut(memf, logf);

		// iter = 2
		for (nodeid i = 0; i < changedUnaries; ++i) {
			nodeid ii = (nodeid)newUnaryPotentialArray[i][0];
			labelid li = (labelid)newUnaryPotentialArray[i][1];
			captype val = newUnaryPotentialArray[i][2];
			memf.updateUnary(ii, li, val - unaryPotential(ii, li));
		}
		flow = memf.optimize(true);
		energy = minimumCut(memf, logf, 1);
	}

	cout << "Max-flow = " << flow << endl;
	
	deleteLocalArrays();

	return 0;
}




