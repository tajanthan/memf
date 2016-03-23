
#include <sstream>

#include "MEMF.h"
#include "util/util.h"

template <typename nodeid, typename labelid, typename captype>
MEMF<nodeid, labelid, captype>::MEMF(nodeid width, nodeid height, labelid labels, char* logf, void(*err_function)(char *))
	: width(width), height(height), labels(labels), err_function(err_function), reuse(false), iter(0), optiter(0), titer(0), flow(0), grp(NULL)
{
	imagexy = width * height;
	totPathLength = 0;
	peakPathLength = 0;
	fout.open(logf);
	std::cout << "\n########## MEMF ##########\n";
	fout << "\n########## MEMF ##########\n";
	std::cout << width << " x " << height << " x " << labels << std::endl;
	fout << width << " x " << height << " x " << labels << std::endl;

	theta_i = new captype[imagexy * labels]();
	gamma_ij = new captype[2 * imagexy]();
	g = new captype[labels * labels]();
	gpp_ij = new captype[(labels + 1) * (labels + 1)]();
	gpp_ji = new captype[(labels + 1) * (labels + 1)]();

	M = new captype[2 * imagexy * 2 * labels]();	// initialize to zero

	_theta_ij = new captype[labels * labels]();
	_phi_ij = new captype[(labels + 1) * (labels + 1)]();
	_phi_ji = new captype[(labels + 1) * (labels + 1)]();
	_phi_i = new captype[labels]();
	_phi_j = new captype[labels]();
	_f_ji = new captype[labels - 1]();
	_f_ij = new captype[labels - 1]();

	labelling = new labelid[imagexy]();
	augPathFlag = new short[imagexy]();

	mem = sizeof(captype);
	size_t len = imagexy * labels + 4 * imagexy + labels * labels * 5 + 2 * imagexy * 2 * labels + 4 * labels;
	mem = mem * len;

	maxBlockCount = 0;
	maxArcCount = 0;
#ifdef MEMF_DEBUG
	augPathLenArray = new nodeid[2000]();
#endif
}

template <typename nodeid, typename labelid, typename captype>
MEMF<nodeid, labelid, captype>::~MEMF()
{
	fout.close();
	delete[] theta_i;
	delete[] gamma_ij;
	delete[] g;
	delete[] gpp_ij;
	delete[] gpp_ji;
	delete[] M;
	delete[] _theta_ij;
	delete[] _phi_ij;
	delete[] _phi_ji;
	delete[] _phi_i;
	delete[] _phi_j;
	delete[] _f_ji;
	delete[] _f_ij;
	delete[] labelling;
	delete[] augPathFlag;
	if (grp) delete grp;
#ifdef MEMF_DEBUG
	delete[] augPathLenArray;
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::setEnergy(captype(*unaryPotential)(nodeid, labelid),
	captype(*binaryWeights)(nodeid, nodeid), captype(*binaryFunction)(labelid, labelid))
{
	for (nodeid ii = 0; ii < imagexy; ++ii) {
		for (labelid xi = 0; xi < labels; ++xi) {
			theta_i[getUnaryId(ii, xi)] = unaryPotential(ii, xi);
		}
	}

	for (nodeid ii = 0; ii < imagexy; ++ii) {
		nodeid i = ii / width;
		nodeid j = ii % width;
		if (j < width - 1) {	// right
			nodeid jj = ii + 1;
			gamma_ij[getEdgeId(ii, jj)] = binaryWeights(ii, jj);
		}
		if (i < height - 1) {	// down
			nodeid jj = ii + width;
			gamma_ij[getEdgeId(ii, jj)] = binaryWeights(ii, jj);
		}
	}

	for (labelid xi = 0; xi < labels; ++xi) {
		for (labelid xj = 0; xj < labels; ++xj) {
			g[getLabelId(xi, xj)] = binaryFunction(xi, xj);
		}
	}

	for (labelid xi = 1; xi < labels; ++xi) {
		for (labelid xj = 1; xj < labels; ++xj) {
			captype e = binaryFunction(xi, xj - 1) + binaryFunction(xi - 1, xj) 
				- binaryFunction(xi - 1, xj - 1) - binaryFunction(xi, xj);
			if (xi < xj) {
				gpp_ij[getPhiId(xi, xj)] = e;
			}
			else if (xi == xj) {
				captype c = e / (captype)2;
				gpp_ij[getPhiId(xi, xj)] = c;
				gpp_ji[getPhiId(xj, xi)] = c;
			}
			else {
				gpp_ji[getPhiId(xj, xi)] = e;
			}
		}
	}
#if MEMF_DEBUG >= 5
	fout << "\ngpp_ij\n";
	printArray(gpp_ij, labels + 1, labels + 1);
	fout << "\ngpp_ji\n";
	printArray(gpp_ji, labels + 1, labels + 1);
	fout << "---------\n";
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::setEnergy(captype *unaryPotentialArray,
	captype *binaryWeightsArray, captype *binaryFunctionArray)
{
	for (nodeid ii = 0; ii < imagexy; ++ii) {
		for (labelid xi = 0; xi < labels; ++xi) {
			theta_i[getUnaryId(ii, xi)] = unaryPotentialArray[getUnaryId(ii, xi)];
		}
	}

	for (nodeid ii = 0; ii < imagexy; ++ii) {
		nodeid i = ii / width;
		nodeid j = ii % width;
		if (j < width - 1) {	// right
			nodeid jj = ii + 1;
			gamma_ij[getEdgeId(ii, jj)] = binaryWeightsArray[getEdgeId(ii, jj)];
		}
		if (i < height - 1) {	// down
			nodeid jj = ii + width;
			gamma_ij[getEdgeId(ii, jj)] = binaryWeightsArray[getEdgeId(ii, jj)];
		}
	}

	for (labelid xi = 0; xi < labels; ++xi) {
		for (labelid xj = 0; xj < labels; ++xj) {
			g[getLabelId(xi, xj)] = binaryFunctionArray[getLabelId(xi, xj)];
		}
	}

	for (labelid xi = 1; xi < labels; ++xi) {
		for (labelid xj = 1; xj < labels; ++xj) {
			captype e = binaryFunctionArray[getLabelId(xi, xj - 1)] + binaryFunctionArray[getLabelId(xi - 1, xj)]
				- binaryFunctionArray[getLabelId(xi - 1, xj - 1)] - binaryFunctionArray[getLabelId(xi, xj)];
			if (xi < xj) {
				gpp_ij[getPhiId(xi, xj)] = e;
			}
			else if (xi == xj) {
				captype c = e / (captype)2;
				gpp_ij[getPhiId(xi, xj)] = c;
				gpp_ji[getPhiId(xj, xi)] = c;
			}
			else {
				gpp_ji[getPhiId(xj, xi)] = e;
			}
		}
	}
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::makeEdgePositiveInit(nodeid ii, nodeid jj)
{	// initiall no need to call Ishikawa params and no swap required, only need to know a positive phi_ij or not!!!!
#ifdef MEMF_DEBUG
	assert(grp);
#endif
#if MEMF_DEBUG >= 5
	fout << std::endl << ii << ", " << jj << std::endl;
#endif

	//// clear arrays
	clearArray(_phi_ij, (labels + 1) * (labels + 1));
	clearArray(_phi_ji, (labels + 1) * (labels + 1));
	if (!grp->validBlockCount[ii]) clearArray(_phi_i, labels);
	if (!grp->validBlockCount[jj]) clearArray(_phi_j, labels);

	captype* phi_ij = _phi_ij;
	captype* phi_ji = _phi_ji;
	captype* phi_i = _phi_i;	
	captype* phi_j = _phi_j;	
	
	ishikawaParamsInit_ij(ii, jj, phi_ij, phi_ji);

#if MEMF_DEBUG >= 5
	fout << "\nphi_ij\n";
	printArray(phi_ij, labels + 1, labels + 1);
	fout << "\nphi_ji\n";
	printArray(phi_ji, labels + 1, labels + 1);
	fout << "---------\n";
#endif

#if MEMF_DEBUG >= 2
	assertUpperTriangular(phi_ij);
	assertUpperTriangular(phi_ji);
#endif

#if MEMF_DEBUG >= 2
	assertPositive(phi_ij);
	assertPositive(phi_ji);
#endif

	if (!grp->validBlockCount[ii]) addUnaryParams(ii, phi_i);
	if (!grp->validBlockCount[jj]) addUnaryParams(jj, phi_j);

#if MEMF_DEBUG >= 2
	if (grp->validBlockCount[ii] == 0) assertPositive(phi_i, labels);
	if (grp->validBlockCount[jj] == 0) assertPositive(phi_j, labels);
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::makeEdgePositiveDirect(nodeid ii, nodeid jj)
{
#ifdef MEMF_DEBUG
	assert(grp);
#endif
#if MEMF_DEBUG >= 5
	fout << std::endl << ii << ", " << jj << std::endl;
#endif

	captype* m_ji = M + getMessageId(ii, jj, 0);			// length: labels
	captype* m_ij = M + getMessageId(ii, jj, labels);
	if (allZero(m_ji, labels) && allZero(m_ij, labels)) {
		makeEdgePositiveInit(ii, jj);
		return;
	}

	//// clear arrays
	clearArray(_phi_ij, (labels + 1) * (labels + 1));
	clearArray(_phi_ji, (labels + 1) * (labels + 1));
	if (!grp->validBlockCount[ii]) clearArray(_phi_i, labels);
	if (!grp->validBlockCount[jj]) clearArray(_phi_j, labels);
	clearArray(_f_ji, labels - 1);
	clearArray(_f_ij, labels - 1);

	captype* phi_ij = _phi_ij;
	captype* phi_ji = _phi_ji;
	captype* phi_i = _phi_i;
	captype* phi_j = _phi_j;
	captype* f_ji = _f_ji;
	captype* f_ij = _f_ij;

	ishikawaParamsInit_ij(ii, jj, phi_ij, phi_ji);

	derivative(m_ji, f_ji);
	derivative(m_ij, f_ij);

#if MEMF_DEBUG >= 5
	fout << "\nStep-0\nphi_ij\n";
	printArray(phi_ij, labels + 1, labels + 1);
	fout << "\nphi_ji\n";
	printArray(phi_ji, labels + 1, labels + 1);
	fout << "\nm_ji\n";
	printArray(m_ji, labels);
	fout << "\nm_ij\n";
	printArray(m_ij, labels);
	fout << "\nf_ji\n";
	printArray(f_ji, labels - 1);
	fout << "\nf_ij\n";
	printArray(f_ij, labels - 1);
	fout << "---------\n";
#endif

	// right edges
	for (labelid xi = 1; xi < labels; ++xi) {
		labelid xi1 = xi - 1;
		if (f_ji[xi1] >= 0) continue;	// check for negative --> (right) outgoing flow!
		for (labelid xj = xi; xj < labels; ++xj) {
			labelid xj1 = xj - 1;
			nodeid lij = getPhiId(xi, xj);
			if (!phi_ij[lij]) continue;
			if (f_ij[xj1] <= 0) continue;	// check for positive --> incoming flow!
			captype f = std::min(-f_ji[xi1], f_ij[xj1]);
			f = std::min(f, phi_ij[lij]);
			
			phi_ij[lij] -= f;
			phi_ji[getPhiId(xj, xi)] += f;

			f_ji[xi1] += f;
			f_ij[xj1] -= f;
			if (!f_ji[xi1]) break;
		}
	}

#if MEMF_DEBUG >= 5
	fout << "\nStep-1\nphi_ij\n";
	printArray(phi_ij, labels + 1, labels + 1);
	fout << "\nphi_ji\n";
	printArray(phi_ji, labels + 1, labels + 1);
	fout << "\nf_ji\n";
	printArray(f_ji, labels - 1);
	fout << "\nf_ij\n";
	printArray(f_ij, labels - 1);
	fout << "---------\n";
	fout.flush();
#endif

#if MEMF_DEBUG >= 2
	assertPositive(f_ji, labels - 1);
#endif

	// left edges
	for (labelid xj = 1; xj < labels; ++xj) {
		labelid xj1 = xj - 1;
		if (f_ij[xj1] >= 0) continue;	// check for negative --> (left) outgoing flow!
		for (labelid xi = xj; xi < labels; ++xi) {
			labelid xi1 = xi - 1;
			nodeid lji = getPhiId(xj, xi);
			if (!phi_ji[lji]) continue;
			if (f_ji[xi1] <= 0) continue;	// check for positive --> incoming flow!
			captype f = std::min(-f_ij[xj1], f_ji[xi1]);
			f = std::min(f, phi_ji[lji]);

			phi_ji[lji] -= f;
			phi_ij[getPhiId(xi, xj)] += f;

			f_ij[xj1] += f;
			f_ji[xi1] -= f;
			if (!f_ij[xj1]) break;
		}
	}

	// TO DO
	/* Check for non-zero exit-flows f_ji, f_ij and push flow
		--> to make the reconstructed flow satisfy the exit-flows 
		*In the current implementation, this situation doesn't seem to occur --> but not sure why?
		--> updated exit-flows are all zero at this point */
	// checking code!!!
#if MEMF_DEBUG >= 2
	for (labelid i = 0; i < labels; ++i) {
		if(f_ji[i]) fout << "\nBUG:: f_ji[" << i << "] = " << f_ji[i] << " != 0 <== (iter: " << iter << ", ii: " << ii << ", jj: " << jj << ")\n";
		if(f_ij[i]) fout << "\nBUG:: f_ij[" << i << "] = " << f_ij[i] << " != 0 <== (iter: " << iter << ", ii: " << ii << ", jj: " << jj << ")\n";
	}
#endif
	// end


#if MEMF_DEBUG >= 5
	fout << "\nStep-2\nphi_ij\n";
	printArray(phi_ij, labels + 1, labels + 1);
	fout << "\nphi_ji\n";
	printArray(phi_ji, labels + 1, labels + 1);
	fout << "\nf_ji\n";
	printArray(f_ji, labels - 1);
	fout << "\nf_ij\n";
	printArray(f_ij, labels - 1);
	fout << "---------\n";
	fout.flush();
#endif

#if MEMF_DEBUG >= 2
	assertZero(f_ji, labels - 1);
	assertZero(f_ij, labels - 1);
#endif

#if MEMF_DEBUG >= 2
	assertPositive(phi_ij);
	assertPositive(phi_ji);
#endif

	if (!grp->validBlockCount[ii]) addUnaryParams(ii, phi_i);
	if (!grp->validBlockCount[jj]) addUnaryParams(jj, phi_j);

#if MEMF_DEBUG >= 2
	if (grp->validBlockCount[ii] == 0) assertPositive(phi_i, labels);
	if (grp->validBlockCount[jj] == 0) assertPositive(phi_j, labels);
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::reduceGraph(nodeid ii, nodeid jj)
{
	// _phi_ij, _phi_ji, _phi_i, _phi_j are updated in makeEdgePositive!!
	// edges cannot be repeated!!!

#ifdef MEMF_DEBUG
	assert(grp);	
#endif
	if (!grp->validBlockCount[ii]) {
#ifdef MEMF_DEBUG
		assert(*std::min_element(_phi_i, _phi_i + labels) == 0);
#endif
		addBlocks(ii, _phi_i);	// add blocks if not already added
	}

	if (!grp->validBlockCount[jj]) {
#ifdef MEMF_DEBUG
		assert(*std::min_element(_phi_j, _phi_j + labels) == 0);
#endif
		addBlocks(jj, _phi_j);
	}

	// not optimal to call twice?
	addBlockEdges(ii, jj, _phi_ij);
	addBlockEdges(jj, ii, _phi_ji);	
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::construct(bool changed)
{
	if (!changed) {
		for (nodeid ii = 0; ii < imagexy; ++ii) {
			nodeid i = ii / width;
			nodeid j = ii % width;
			if (j < width - 1) {	// right
				nodeid jj = ii + 1;

				makeEdgePositiveInit(ii, jj);

				reduceGraph(ii, jj);
			}
			if (i < height - 1) {	// down
				nodeid jj = ii + width;

				makeEdgePositiveInit(ii, jj);

				reduceGraph(ii, jj);
			}
		}
	}
	else {	// process changed list
#ifdef MEMF_DEBUG
		assert(grp);
		assert(!changedList.empty());
#endif
		for (typename std::vector<nodeid>::iterator it = changedList.begin(); it != changedList.end(); ++it) {
			nodeid ii = *it;
			if (augPathFlag[ii] != TO_PROCESS) continue;	// process only once
			augPathFlag[ii] = PROCESSED;

			// delete affected portion
			for (nodeid dir = 0; dir < 4; ++dir) {
				if (neighborExists(ii, dir)) {
					nodeid kk = getNeighborNode(ii, dir);
					if (!augPathFlag[kk]) {		// not part of the remaining augPath & not already processed
						grp->deleteEdges(kk, ii);	// delete edges kk --> ii
					}
				}
			}
			grp->deleteBlocks(ii);	// remove blocks in ii

			// reconstruct affected portion
			for (nodeid dir = 0; dir < 4; ++dir) {
				if (neighborExists(ii, dir)) {
					nodeid kk = getNeighborNode(ii, dir);
					if (augPathFlag[kk] != TO_PROCESS) {	// not part of the remaining augPath
						if (kk < ii) {
							makeEdgePositiveDirect(kk, ii);	// may not be needed to spearate?
							reduceGraph(kk, ii);
						}
						else {
							makeEdgePositiveDirect(ii, kk);
							reduceGraph(ii, kk);
						}
					}
				}
			}
		}
		grp->repairBFSTree();		
	}
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::construct(std::vector<RGblock*>& augPathList)
{
#ifdef MEMF_DEBUG
	assert(grp);
	assert(!augPathList.empty());
#endif

	for (typename std::vector<RGblock*>::iterator it = augPathList.begin(); it != augPathList.end(); ++it) {
		RGblock* iiB = *it;
		augPathFlag[iiB->node] = TO_PROCESS;
	}

	for (typename std::vector<RGblock*>::iterator it = augPathList.begin(); it != augPathList.end(); ++it) {
		nodeid ii = (*it)->node;
		if (augPathFlag[ii] != TO_PROCESS) continue;	// process only once
		augPathFlag[ii] = PROCESSED;

		// delete affected portion
#if MEMF_DEBUG >= 1
		stimer[2].resume();
#endif
		for (nodeid dir = 0; dir < 4; ++dir) {
			if (neighborExists(ii, dir)) {
				nodeid kk = getNeighborNode(ii, dir);
				if (!augPathFlag[kk]) {		// not part of the remaining augPath & not already processed
					grp->deleteEdges(kk, ii);	// delete edges kk --> ii
				}
			}
		}
		grp->deleteBlocks(ii);	// remove blocks in ii
#if MEMF_DEBUG >= 1
		stimer[2].pause();
#endif

		// reconstruct affected portion
		for (nodeid dir = 0; dir < 4; ++dir) {
			if (neighborExists(ii, dir)) {
				nodeid kk = getNeighborNode(ii, dir);
				if (augPathFlag[kk] != TO_PROCESS) {	// not part of the remaining augPath
					if (kk < ii) {
#if MEMF_DEBUG >= 1
						stimer[3].resume();
#endif
						makeEdgePositiveDirect(kk, ii);	// may not be needed to spearate?
#if MEMF_DEBUG >= 1
						stimer[3].pause();
						stimer[2].resume();
#endif
						reduceGraph(kk, ii);
#if MEMF_DEBUG >= 1
						stimer[2].pause();
#endif
					}
					else {
#if MEMF_DEBUG >= 1
						stimer[3].resume();
#endif
						makeEdgePositiveDirect(ii, kk);
#if MEMF_DEBUG >= 1
						stimer[3].pause();
						stimer[2].resume();
#endif
						reduceGraph(ii, kk);
#if MEMF_DEBUG >= 1
						stimer[2].pause();
#endif
					}
				}
			}
		}
	}
#if MEMF_DEBUG >= 1
	stimer[4].resume();
#endif
	grp->repairBFSTree();
#if MEMF_DEBUG >= 1
	stimer[4].pause();
#endif

	for (typename std::vector<RGblock*>::iterator it = augPathList.begin(); it != augPathList.end(); ++it) {
		nodeid ii = (*it)->node;
		augPathFlag[ii] = 0;
	}
#if MEMF_DEBUG >= 2
	for (nodeid i = 0; i < imagexy; ++i) assert(augPathFlag[i] == 0);
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::augment(std::vector<RGblock*>& augPathList)
{
	RGblock* iiBlock = augPathList.front();
#if MEMF_DEBUG >= 4
	fout << std::endl << "(" << iiBlock->node << ", " << iiBlock->lcheck << "," << iiBlock->lhat << ")";
#endif
	for (typename std::vector<RGblock*>::iterator it = std::next(augPathList.begin()); it != augPathList.end(); ++it) {
		RGblock* jjBlock = *it;

		if (iiBlock->node != jjBlock->node) {
			// forward reparametrization
			reparametrize(iiBlock->node, iiBlock->lcheck, jjBlock->node, jjBlock->lcheck);
#if MEMF_DEBUG >= 4
			fout << " --> (" << jjBlock->node << ", " << jjBlock->lcheck << "," << jjBlock->lhat << ")";
#endif
			//// backward reparametrization
			//backReparametrize(jjBlock->node, jjBlock->lcheck, iiBlock->node, iiBlock->lcheck);
		}

		iiBlock = jjBlock;
	}

	// subtract minimum from unary 
	captype* th_i = _phi_i;
	clearArray(th_i, labels);
	addUnaryParams(iiBlock->node, th_i);

	captype m = *std::min_element(th_i, th_i + labels);
	fout.flush();
#ifdef MEMF_DEBUG
	assert(m > 0);
#endif
	trivialAugmentation(iiBlock->node, m);

	// backward reparametrization to make the reduced graph sparse
	// this also ensures _phi_i & _phi_j are all zero after make_positive!!?
	nodeid lastNode = iiBlock->node;
#if MEMF_DEBUG >= 4
	fout << std::endl << "Backward: (" << iiBlock->node << ", " << iiBlock->lcheck << "," << iiBlock->lhat << ")";
#endif
	for (typename std::vector<RGblock*>::reverse_iterator rit = std::next(augPathList.rbegin()); rit != augPathList.rend(); ++rit) {
		RGblock* jjBlock = *rit;

		if (iiBlock->node != jjBlock->node) {
			// backward reparametrization
			backReparametrize(iiBlock->node, iiBlock->lcheck, jjBlock->node, jjBlock->lcheck);
#if MEMF_DEBUG >= 4
			fout << " --> (" << jjBlock->node << ", " << jjBlock->lcheck << "," << jjBlock->lhat << ")";
#endif
		}

		iiBlock = jjBlock;
	}

	clearArray(th_i, labels);
	addUnaryParams(lastNode, th_i);
	captype m1 = *std::min_element(th_i, th_i + labels);
	//assert(m1 == 0); //==> WHY?
	if (m1 > 0) trivialAugmentation(lastNode, m1);

#if MEMF_DEBUG >= 4
	fout << " ==> m1 = " << m1 << std::endl;
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::trimAugPath(std::vector<RGblock*>& augPathList, RGblock* startBlock)
{
#ifdef MEMF_DEBUG
	assert(grp);
#endif
	RGblock* iiBlock = startBlock;
#if MEMF_DEBUG >= 4
	fout << std::endl << "(" << iiBlock->node << ", " << iiBlock->lcheck << "," << iiBlock->lhat << ")";
#endif
	while (iiBlock) {		
		if (!iiBlock->dummy) {
			augPathList.push_back(iiBlock);
			// mark all above blocks dummy
			RGblock* column_i = grp->blockGrid + grp->getGridId(iiBlock->node, 0);
			for (labelid i = iiBlock->id + 1; i < grp->validBlockCount[iiBlock->node]; ++i) {
				column_i[i].dummy = 1;
			}
		}
		else {
			while (augPathList.back()->node != iiBlock->node) {
				RGblock* v = augPathList.back();
				RGblock* column_i = grp->blockGrid + grp->getGridId(v->node, 0);
				for (labelid i = v->id + 1; i < grp->validBlockCount[v->node]; ++i) {
					column_i[i].dummy = 0;
				}
				augPathList.pop_back();
			}
			RGblock* v = augPathList.back();
			for (nodeid aid = v->firstChild; aid != NONE; aid = grp->arcs[aid].nextChild) {
				RGarc* a = grp->arcs + aid;
				RGblock* w = grp->blockGrid + a->dest;
#ifdef MEMF_DEBUG
				assert(iiBlock->next != NONE);
#endif
				RGblock* jjB = grp->blockGrid + iiBlock->next;
#ifdef MEMF_DEBUG
				assert(jjB);
#endif
				if (jjB != w && jjB->node == w->node) {
					augPathList.push_back(w);	// hack
					break;
				}
			}
		}
		if (iiBlock->next == NONE) break;
		iiBlock = grp->blockGrid + iiBlock->next;
#if MEMF_DEBUG >= 4
		if (iiBlock) fout << " --> (" << iiBlock->node << ", " << iiBlock->lcheck << "," << iiBlock->lhat << ")";
#endif
		fout.flush();
	}
	resetFlag(startBlock);
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::resetFlag(RGblock* startBlock)
{
#ifdef MEMF_DEBUG
	assert(grp);
#endif
	RGblock* iiBlock = startBlock;
	while (iiBlock) {
		// reset dummy flag
		RGblock* column_i = grp->blockGrid + grp->getGridId(iiBlock->node, 0);
		for (labelid i = iiBlock->id + 1; i < grp->validBlockCount[iiBlock->node]; ++i) {
			column_i[i].dummy = 0;
		}
		if (iiBlock->next == NONE) break;
		iiBlock = grp->blockGrid + iiBlock->next;
	}
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::finalLabelling()
{
	// find labelling from reduced graph
	for (nodeid ii = 0; ii < imagexy; ++ii) {
		labelling[ii] = grp->getCutLabel(ii);
	}

#if MEMF_DEBUG >= 2
	for (nodeid ii = 0; ii < imagexy; ++ii) {
		assert(labelling[ii] < labels);
	}
#endif
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::updateUnary(nodeid ii, labelid li, captype val)
{
	if (!reuse) err_function((char*)"Function call is invalid, since reuse == false!");

	if (!augPathFlag[ii]) {
		augPathFlag[ii] = TO_PROCESS;
		changedList.push_back(ii);
	}

	theta_i[getUnaryId(ii, li)] += val;
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::init()
{
	if (grp) err_function((char*)"Invalid call to optimize, reuse must be true!");

	// subtract minimum from each node unary potential
	for (nodeid ii = 0; ii < imagexy; ++ii) {
		trivialAugmentation(ii);
	}
	{
		std::stringstream ss;
		ss << "optimize (flow = " << flow << ", iter = " << iter << ", titer = " << titer << ")";
		timer.printData(fout, ss.str());
	}

	grp = new ReducedGraph(imagexy, width, height, labels);
	construct();

#if MEMF_DEBUG >= 1
	std::cout << "\nEstimated memory: \nEMF: " << (float)mem / (1024.0 * 1024.0) << " MB, RG: "
		<< (float)grp->rmem / (1024.0 * 1024.0) << " MB, TOTAL: " << (float)(mem + grp->rmem) / (1024.0 * 1024.0) << " MB" << std::endl;
	fout << "\nEstimated memory: \nEMF: " << (float)mem / (1024.0 * 1024.0) << " MB, RG: "
		<< (float)grp->rmem / (1024.0 * 1024.0) << " MB, TOTAL: " << (float)(mem + grp->rmem) / (1024.0 * 1024.0) << " MB" << std::endl;
	fout << "captype: " << sizeof(captype) << std::endl;
	fout << "nodeid: " << sizeof(nodeid) << std::endl;
	fout << "labelid: " << sizeof(labelid) << std::endl;
	fout << "size_t: " << sizeof(size_t) << std::endl;
	fout << "unsigned short: " << sizeof(unsigned short) << std::endl;

	fout << "block: " << sizeof(RGblock) << std::endl;
	fout << "block*: " << sizeof(RGblock*) << std::endl;
	fout << "arc: " << sizeof(RGarc) << std::endl;
	fout << "arc*: " << sizeof(RGarc*) << std::endl;
#endif

	grp->initBFS();
}

template <typename nodeid, typename labelid, typename captype>
void MEMF<nodeid, labelid, captype>::initReuse()
{
#ifdef MEMF_DEBUG
	assert(reuse);
#endif
	iter = 0;
	//flow = 0;
	peakPathLength = 0;
	totPathLength = 0;
	
	if (changedList.empty()) return;	// nothing to do

	for (typename std::vector<nodeid>::iterator it = changedList.begin(); it != changedList.end(); ++it) {
		nodeid ii = *it;
		// subtract minimum from unary 
		captype* th_i = _phi_i;
		clearArray(th_i, labels);
		addUnaryParams(ii, th_i);

		captype m = *std::min_element(th_i, th_i + labels);
		trivialAugmentation(ii, m);
	}
	{
		std::stringstream ss;
		ss << "optimize (flow = " << flow << ", iter = " << iter << ")";
		timer.printData(fout, ss.str());
	}

	construct(true);

	for (typename std::vector<nodeid>::iterator it = changedList.begin(); it != changedList.end(); ++it) {
		augPathFlag[*it] = 0;
	}
	changedList.clear();
#if MEMF_DEBUG >= 2
	for (nodeid i = 0; i < imagexy; ++i) assert(augPathFlag[i] == 0);
#endif
}

template <typename nodeid, typename labelid, typename captype>
captype MEMF<nodeid, labelid, captype>::optimize(bool reuse)
{
	timer.startTimer();	

	this->reuse = reuse;
	if (!optiter) init();
	else if (reuse) initReuse();	
	++optiter;

#if MEMF_DEBUG >= 1
	for (int i = 0; i < 5; ++i) stimer[i].startTimer();
#endif

	std::vector<RGblock*> augPathList;
	while (true) {		
		++iter;

#if MEMF_DEBUG >= 1
		stimer[0].resume();
#endif
		RGblock* startBlock = grp->findAugPath();
#if MEMF_DEBUG >= 1
		stimer[0].pause();
#endif
		if (startBlock == NULL) break;

#if MEMF_DEBUG >= 1
		stimer[1].resume();
#endif
		trimAugPath(augPathList, startBlock);	// remove unusable blocks!
		augment(augPathList);	// do reparametrization
#if MEMF_DEBUG >= 1
		stimer[1].pause();
#endif

#ifdef MEMF_DEBUG
		nodeid len = (nodeid)augPathList.size();
		totPathLength += len;
		if (peakPathLength < len) {
			peakPathLength = len;
			augPathList.shrink_to_fit();
		}
		++augPathLenArray[len];
#if MEMF_DEBUG >= 3
		{
			std::stringstream ss;
			ss << "optimize (flow = " << flow << ", iter = " << iter << ")";
#if MEMF_DEBUG >= 1
			maxBlockCount = std::max(maxBlockCount, grp->blockCount);
			maxArcCount = std::max(maxArcCount, grp->arcCount);
			ss << "\nblockCount = " << grp->blockCount << ", arcCount = " << grp->arcCount;
#endif
			timer.printData(fout, ss.str());
		}
#else
		if(iter % 1000 == 0) {
			std::stringstream ss;
			ss << "optimize (flow = " << flow << ", iter = " << iter << ", avg-length = " << (totPathLength / iter) 
				<< ", peak-length = " << peakPathLength << ")";
#if MEMF_DEBUG >= 1
			maxBlockCount = std::max(maxBlockCount, grp->blockCount);
			maxArcCount = std::max(maxArcCount, grp->arcCount);
			ss << "\nblockCount = " << grp->blockCount << ", arcCount = " << grp->arcCount;
#endif
			timer.printData(fout, ss.str());
		}
#endif
#endif

		// delete aug-path portion and reconstruct
		construct(augPathList);
		augPathList.clear();
	}	
#if MEMF_DEBUG >= 2
	grp->testConsistency(2);
#endif

	{
		std::stringstream ss;
#ifdef MEMF_DEBUG
		ss << "optimize (flow = " << flow << ", iter = " << iter << ", avg-length = " << (totPathLength / iter) 
			<< ", peak-length = " << peakPathLength << ")";
#else
		ss <<  "optimize (flow = " << flow << ", iter = " << iter << ")";
#endif
#if MEMF_DEBUG >= 1
		maxBlockCount = std::max(maxBlockCount, grp->blockCount);
		maxArcCount = std::max(maxArcCount, grp->arcCount);
		ss << "\nmaxBlockCount = " << maxBlockCount << ", maxArcCount = " << maxArcCount;
		ss << "\nBlock % = " << (float)maxBlockCount / (imagexy * labels) * 100 << ", Arc % = " 
			<< (float)maxArcCount / (2 * imagexy * 2 * labels) * 100;
#endif
		timer.printData(fout, ss.str());
	}

#ifdef MEMF_DEBUG
	fout << "\n##AugPathLenArray\n";
	for (nodeid i = 1; i <= peakPathLength; ++i) fout << i << ": " << augPathLenArray[i] << std::endl;
#endif
	
	double totTime = timer.getTotalTime();
	timer.stopTimer();

	finalLabelling();

	std::cout << "\nEstimated memory: \nMEMF: " << (float)mem / (1024.0 * 1024.0) << " MB, RG: " 
		<< (float)grp->rmem / (1024.0 * 1024.0) << " MB, TOTAL: " << (float)(mem + grp->rmem) / (1024.0 * 1024.0) << " MB" << std::endl;
	fout << "\nEstimated memory: \nMEMF: " << (float)mem / (1024.0 * 1024.0) << " MB, RG: "
		<< (float)grp->rmem / (1024.0 * 1024.0) << " MB, TOTAL: " << (float)(mem + grp->rmem) / (1024.0 * 1024.0) << " MB" << std::endl;

	size_t memUsage = getMemoryUsage();
	float mb = (float)memUsage / (1024.0 * 1024.0);
	std::cout << "PeakPagefileUsage: " << mb << " MB\n";
	fout << "PeakPagefileUsage: " << mb << " MB\n";
	fout.flush();

#if MEMF_DEBUG >= 1
	fout << "\nSubroutine runtimes %:\n";
	fout << "augmenting_path: " << stimer[0].getTotalTime()*100/totTime << "\n";
	fout << "augment: " << stimer[1].getTotalTime()*100/totTime << "\n";
	fout << "simplify_graph: " << stimer[2].getTotalTime()*100/totTime << "\n";
	fout << "compute_edges: " << stimer[3].getTotalTime()*100/totTime << "\n";
	fout << "tree_recycling: " << stimer[4].getTotalTime()*100/totTime << "\n";
	fout.flush();
#endif
		
	return flow;
}

template <typename nodeid, typename labelid, typename captype>
captype MEMF<nodeid, labelid, captype>::optimize(int& iterations)
{
	captype f = optimize();
	iterations = iter;
	return f;
}

#include "memf_instances.inc"
