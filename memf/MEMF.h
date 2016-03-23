/*
	This code implements the MEMF algorithm
	described in the following paper 

	"Memory Efficient Max Flow for Multi-label Submodular MRFs", 
	Thalaiyasingam Ajanthan, Richard Hartley and Mathieu Salzmann,
	IEEE Conference on Computer Vision and Pattern Recognition,
	June 2016.

	*Code Assumptions:
		1. The MRF energy has the following form
			E(x) = \sum \theta_{i}(x_i) + \sum \theta_{ij} (x_i, x_j),
			where \theta_{ij} (x_i, x_j) = \gamma_{ij} g(x_i, x_j)
		2. The code currently supports MRF with 4-connected grid structure only
			nodes labelled from 0 --> width * height - 1, in a grid structure.
			E.g. width = 3, height = 2
			0 -- 1 -- 2
			|	 |	  |
			3 -- 4 -- 5

	If you use this code, please consider citing the aforementioned paper 
	in any resulting publication.

	This code is for research purposes only, if you want to use it for commercial purpose
	please contact us.
	
	*Contact: thalaiyasingam.ajanthan@nicta.com.au
*/

#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <new>
#include <queue>
#include <vector>

#include "util/util.h"
//#define MEMF_DEBUG 0

// Possible instances are in emf_instances.inc
template <typename nodeid, typename labelid, typename captype> class MEMF
{
public:

	/* Constructor
		width - image width, height - image height, labels - number of labels 
		err_function - called in case of an error */
	MEMF(nodeid width, nodeid height, labelid labels, char* logf, void(*err_function)(char *) = NULL);

	/* Destructor */
	~MEMF();

	/* Set the MRF energy
		unaryPotential(i, li)  - returns the unary potential \theta_{i}(li)
		binaryWeights(i, j)    - returns the constant binary weights \gamma_{ij}
		binaryFunction(li, lj) - returns the smoothness cost \theta(li, lj) */
	void setEnergy(captype(*unaryPotential)(nodeid, labelid), captype(*binaryWeights)(nodeid, nodeid),
		captype(*binaryFunction)(labelid, labelid));

	/* Set the MRF energy
	unaryPotentialArray    - unary potential \theta_{i}(li) --> index: getUnaryId(i, li)
	binaryWeightsArray     - binary weights \gamma_{ij} --> index: getEdgeId(i, j)
	binaryFunctionArray    - smoothness cost \theta(li, lj) --> index: getLabelId(li, lj) */
	void setEnergy(captype *unaryPotentialArray, captype *binaryWeightsArray,
		captype *binaryFunctionArray);

	/* Run algorithm 
		returns maximum flow */
	captype optimize(int& iter);

	/* Run algorithm
	returns maximum flow */
	captype optimize(bool reuse = false);

	/* Get the label
		returns the label of the given node */
	labelid getLabel(nodeid ii)
	{
		return labelling[ii];
	}

	/* Update unary potential, if reuse == true
		val	- new_value - old_value --> val = th'_{i}(li) - th_{i}(li) */
	void updateUnary(nodeid ii, labelid li, captype val);

private:

	static const unsigned short ZERO = 1;	// node 0
	static const unsigned short ONE = 2;	// node 1
	static const int NONE = -1;	// equivalent of NULL for index

	static const int TO_PROCESS = 1;	// flags used in augPathFlag[];
	static const int PROCESSED = 2;

	bool reuse;		// flag to indicate reuse-tree
	nodeid width, height;
	labelid labels;
	nodeid imagexy;		// image grid
	size_t iter, optiter, titer;
	std::ofstream fout;	// log file
	captype flow;		// total flow
	nodeid totPathLength;	// sum of aug-path lengths
	nodeid peakPathLength;	// peak aug-path length
	size_t mem;				// memory estimate

	captype* theta_i;	// \theta_{i}(x_i), length: imagexy*labels
	captype* gamma_ij;	// \gamma_{ij}, length: 2*images
	captype* g;			// g(x_i, x_j),	length: labels*labels
	captype* gpp_ij;	// initial ishikawa weights, length: (labels+1)*(labels+1)	// i --> j
	captype* gpp_ji;	// initial ishikawa weights, length: (labels+1)*(labels+1)	// i --> j

	captype* M;			// message vectors, length: 2*imagexy * 2*labels (2*labels per edge [m_ji; m_ij], 4-connected graph ==> atmost 2 edges (right, down) per node)
	
	captype* _theta_ij;	// multi-label graph parameters for the edge (i,j), length: labels*labels (clear before update)
	captype* _phi_ij;	// ishikawa parameters for the edge (i,j), length: (labels+1)*(labels+1) (clear before update)
	captype* _phi_ji;	// ishikawa parameters for the edge (j,i), length: (labels+1)*(labels+1) (clear before update)
	captype* _phi_i;	// ishikawa parameters for the node i, length: labels (clear before update)
	captype* _phi_j;	// ishikawa parameters for the node j, length: labels (clear before update)
	captype* _f_ji;		// flow vector for message m_ji, length: labels - 1 (clear before update)
	captype* _f_ij;		// flow vector for message m_ij, length: labels - 1 (clear before update)

	labelid* labelling;	// result labelling, length: imagexy
	short* augPathFlag;	// flag nodeids in aug-path, length: imagexy

	void(*err_function)(char *);

	nodeid maxBlockCount, maxArcCount;
	nodeid *augPathLenArray;	// length of aug-path --> frequency array, length: max-length (2000);

	Timer timer;
	Timer stimer[5];
	// inner class for the reduced graph
	////////////////////
	class ReducedGraph
	{
	public:

		struct arc;
		// reduced graph block
		struct block
		{	// order --> largest to smallest
			size_t t;				// time when marked discovered --> cycle avoidance!

			nodeid firstChild;		// first arc with this source - id
			nodeid firstParent;		// first arc with this dest - id
			nodeid next, prev;		// next and prev blocks in the path - ids

			nodeid node;			// nodeid in the image grid
			labelid id;				// block id
			labelid lcheck, lhat;	// low and high label

			unsigned short discovered : 1;	// currently in the Q or already visited
			unsigned short active : 1;		// currently in the queue --> (if active == true ==> discovered == true)
			unsigned short repair : 1;		// mark to be repaired
			unsigned short valid : 1;		// valid block
			unsigned short dummy : 1;		// mark dummy block --> lower node in path
			unsigned short terminal : 2;	// indicates source or sink			

			inline void init()
			{
				discovered = 0;
				active = 0;
				repair = 0;
				valid = 0;
				dummy = 0;
				terminal = 0;
				t = 0;
				node = 0;
				id = 0;
				lcheck = 0;
				lhat = 0;
				next = NONE;
				prev = NONE;
				firstChild = NONE;
				firstParent = NONE;
			}

			block() 
			{
				init();
			}

			inline void set(nodeid node, nodeid id, labelid lcheck, labelid lhat, unsigned short terminal)
			{
				init();
				this->valid = 1;
				this->node = node;
				this->id = id;
				this->lcheck = lcheck;
				this->lhat = lhat;
				this->terminal = terminal;
				//if (cap > 0) terminal = ZERO;		// connected to node 0
				//else if(cap < 0) terminal = ONE;	// connected to node 1
			}
		};

		struct arc
		{
			nodeid source;
			nodeid dest;
			nodeid nextChild;
			nodeid nextParent;

			arc() : source(NONE), dest(NONE), nextChild(NONE), nextParent(NONE)
			{};

			void set(block* sourceb, block* destb, nodeid source, nodeid dest, nodeid id)
			{
				this->source = source;
				this->dest = dest;
				this->nextChild = sourceb->firstChild;
				sourceb->firstChild = id;
				this->nextParent = destb->firstParent;
				destb->firstParent = id;
			}
		};

		block *blockGrid;				// array of block: length = imagexy * labels
		labelid *validBlockCount;		// stores the valid block count: length = imagexy
		arc *arcs;						// storce arc array: length = 2 * imagexy * labels * 2
		nodeid imagexy_;
		labelid labels_;
		nodeid width_;
		nodeid height_;
		size_t rmem;		// memory estimate
		nodeid blockCount, arcCount;

		size_t time;			// monotonically increasing counter
		
		std::queue<block*> Q;	// queue to store the frontier of BFS search tree --> may be changed later?
		std::vector<nodeid> augPathList;	// store the aug-path for repair
		std::list<block*> repairList;	// store the affected neighbor blocks to repair

		// constructor
		ReducedGraph(nodeid imagexy, nodeid width, nodeid height, labelid labels)
			: imagexy_(imagexy), width_(width), height_(height), labels_(labels), time(0)
		{
			blockGrid = new block[imagexy_ * labels_]();
			validBlockCount = new labelid[imagexy_]();
			arcs = new arc[2 * imagexy_ * 2 * labels_]();

			rmem = sizeof(block) * imagexy_ * labels_ +
				sizeof(labelid) * imagexy_ + sizeof(arc) * 2 * imagexy_ * 2 * labels_;
#if MEMF_DEBUG >= 1
			blockCount = 0;
			arcCount = 0;
#endif
		}

		// destructor
		~ReducedGraph()
		{
			delete[] blockGrid;
			delete[] validBlockCount;
			delete[] arcs;
		}

		// returns the unary id (ii,x_i)
		inline nodeid getGridId(nodeid ii, labelid xi)
		{
			return ii * labels_ + xi;
		}

		// returns the message id (ii, jj, xx) [mji, mij]
		inline nodeid getArcId(nodeid e, bool swap, labelid xx)
		{
			labelid ll = xx;
			if (swap) {
				ll = (xx >= labels_) ? xx - labels_ : xx + labels_;
			}
			return e * 2 * labels_ + ll;
		}

		// create and add block, given start block
		inline void addBlock(block* bp, nodeid node, labelid lcheck, labelid lhat, unsigned short terminal)
		{
#ifdef MEMF_DEBUG
			assert(lhat >= lcheck);
#endif
			nodeid id = validBlockCount[node];
			block* bi = bp + id;
			bi->set(node, id, lcheck, lhat, terminal);
			validBlockCount[node] = id + 1;
#if MEMF_DEBUG >= 1
			++blockCount;
#endif
		}

		// add an edge, given blocks and edge id
		inline void addEdge(block* bi, block* bj, nodeid bid, nodeid bjd, nodeid e, bool swap)
		{
#ifdef MEMF_DEBUG
			assert(bi->valid && bj->valid);
#endif

			nodeid id = getArcId(e, swap, bi->id);
			arc* a = arcs + id;
			a->set(bi, bj, bid, bjd, id);
#if MEMF_DEBUG >= 1
			++arcCount;
#endif
		}

		// find augpath, returns the starting block
		block* findAugPath()
		{
#if MEMF_DEBUG >= 2
			testConsistency(1);
#endif
			block* endBlock = BFS();
#if MEMF_DEBUG >= 2
			testConsistency(1);
#endif
			if (!endBlock) return NULL;	// no augmenting path possible

			block* startBlock = blockGrid + endBlock->prev;
			startBlock->next = (nodeid)(endBlock - blockGrid);
			while (startBlock->terminal != ZERO) {
				block* tmp = startBlock;
#ifdef MEMF_DEBUG
				assert(tmp->prev != NONE);
#endif
				startBlock = blockGrid + tmp->prev;
#ifdef MEMF_DEBUG
				assert(startBlock);
				assert(startBlock->discovered);
#endif
				startBlock->next = (nodeid)(tmp - blockGrid);
			}
#ifdef MEMF_DEBUG
			assert(startBlock->prev == NONE);
#endif
			return startBlock;
		}

		// initialize BFS
		inline void initBFS()
		{
			for (nodeid i = 0; i < imagexy_; ++i) {
				block* v = blockGrid + getGridId(i, validBlockCount[i] - 1);
				if (v->terminal == ZERO) {
					Q.push(v);
					v->discovered = 1;
					v->active = 1;
				}
			}
		}

		// do BFS
		block* BFS()
		{
			while (!Q.empty()) {
				block* v = Q.front();
				Q.pop();
				if (v->active) {	// process active blocks only
					
					if (v->terminal == ONE) return v;	// sink found

					v->active = false;	

					// grow tree
					for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
						arc* a = arcs + aid;
						block* w = blockGrid + a->dest;
#ifdef MEMF_DEBUG
						assert(w->valid);
#endif
						if (!w->discovered) {	// not discovered yet
							Q.push(w);
							w->discovered = 1;
							w->active = 1;
							w->prev = (nodeid)(v - blockGrid);
							w->t = ++time;
						}
					}
				}
			}
			return NULL;	// no augpath exists
		}
		
		// delete edges from kk-blocks to ii-blocks (kk --> ii)
		void deleteEdges(nodeid kk, nodeid ii)
		{
			block* column_k = blockGrid + getGridId(kk, 0);
			for (labelid lk = 0; lk < validBlockCount[kk]; ++lk) {
				block* v = column_k + lk;

#if MEMF_DEBUG >= 6
				std::cout << "\nchild (" << v->node << ", " << v->lcheck << "," << v->lhat << ")";
				for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
					arc* a = arcs + aid;
					block* w = blockGrid + a->dest;
					std::cout << " --> (" << w->node << ", " << w->lcheck << "," << w->lhat << ")";
				}
#endif
				arc* pa = NULL;	// prev-a
				for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
					arc* a = arcs + aid;
					block* w = blockGrid + a->dest;
					if (w->node == ii) {
						if (pa) pa->nextChild = a->nextChild;
						else v->firstChild = a->nextChild;
#if MEMF_DEBUG >= 1
						--arcCount;
#endif
						break;	// only one
					}
					pa = a;
				}
				
#if MEMF_DEBUG >= 6
				std::cout << "\nparent (" << v->node << ", " << v->lcheck << "," << v->lhat << ")";
				for (nodeid aid = v->firstParent; aid != NONE; aid = arcs[aid].nextParent) {
					arc* a = arcs + aid;
					block* w = a->source;
					std::cout << " --> (" << w->node << ", " << w->lcheck << "," << w->lhat << ")";
				}
#endif
				bool start = false;
				pa = NULL;	// prev-a
				nodeid ca = NONE;	// cur-a-id
				for (nodeid aid = v->firstParent; aid != NONE; aid = arcs[aid].nextParent) {
					arc* a = arcs + aid;
					block* w = blockGrid + a->source;
					if (w->node == ii) {
						if (!start) start = true;
					}
					else {
						if (start) {
							ca = aid;
							break;
						}
						pa = a;
					}
				}
				if (start) {
					if (pa) pa->nextParent = ca;
					else v->firstParent = ca;
				}
#if MEMF_DEBUG >= 2
#if MEMF_DEBUG >= 6
				std::cout << "\nchild (" << v->node << ", " << v->lcheck << "," << v->lhat << ")";
#endif
				for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
					arc* a = arcs + aid;
					block* w = blockGrid + a->dest;
#if MEMF_DEBUG >= 6
					std::cout << " --> (" << w->node << ", " << w->lcheck << "," << w->lhat << ")";
#endif
					assert(w->node != ii);
				}
#if MEMF_DEBUG >= 6
				std::cout << "\nparent (" << v->node << ", " << v->lcheck << "," << v->lhat << ")";
#endif
				for (nodeid aid = v->firstParent; aid != NONE; aid = arcs[aid].nextParent) {
					arc* a = arcs + aid;
					block* w = blockGrid + a->source;
#if MEMF_DEBUG >= 6
					std::cout << " --> (" << w->node << ", " << w->lcheck << "," << w->lhat << ")";
#endif
					assert(w->node != ii);
				}
#endif

				if (v->discovered) {
					// check for prev of v
					if (v->prev != NONE && blockGrid[v->prev].node == ii) {
						v->prev = NONE;
						if (!v->repair) {
							v->repair = 1;
							repairList.push_back(v);
						}
					}
				}
				if (v->discovered && !v->active && (v->prev != NONE || v->terminal == ZERO)) {
					if (!v->repair) {
						v->repair = 1;
						repairList.push_back(v);
					}
				}
			}
		}

		// delete blocks in ii
		void deleteBlocks(nodeid ii)
		{
			augPathList.push_back(ii);
			block* column_i = blockGrid + getGridId(ii, 0);
			for (labelid li = 0; li < validBlockCount[ii]; ++li) {
				block* v = column_i + li;
#if MEMF_DEBUG >= 1
				--blockCount;
				for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
					--arcCount;
				}
#endif
				v->init();
			}
			validBlockCount[ii] = 0;
		}

		// returns true if w can be a valid prev to v
		inline bool isValidPrev(block* v, block* w)
		{
#ifdef MEMF_DEBUG
			assert(v->valid && w->valid);
#endif
			bool res = false;
			if (w->terminal == ONE) res = false;
			else if (w->terminal == ZERO) res = true;
			else if (w->discovered/* && !w->active*/ && w->prev != NONE) {
				if (w->t <= v->t) res = true;
			}
			return res;
		}

		// find valid prev from parent list
		inline void findValidPrev(block* v)
		{
			for (nodeid aid = v->firstParent; aid != NONE; aid = arcs[aid].nextParent) {
				arc* a = arcs + aid;
				block* w = blockGrid + a->source;
#ifdef MEMF_DEBUG
				assert(w->valid);
#endif
				if (isValidPrev(v, w)) {
					v->prev = (nodeid)(w - blockGrid);
					v->repair = 0;
					return;
				}
			}
		}

		// repair neighbor, since v is no longer discovered, from children list
		inline void repairNeighbor(block* v)
		{
			for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
				arc* a = arcs + aid;
				block* w = blockGrid + a->dest;
#ifdef MEMF_DEBUG
				assert(w->valid);
#endif
				if (w->prev != NONE && (blockGrid + w->prev) == v) {
					w->prev = NONE;
					if (!w->repair) {
						w->repair = 1;
						repairList.push_back(w);	// repair w
					}
				}
			}
		}

		// repair affected neighbor blocks and augPath nodes
		void repairBFSTree()
		{
			for (typename std::vector<nodeid>::iterator it = augPathList.begin(); it != augPathList.end(); ++it) {
				nodeid ii = *it;
				block* v = blockGrid + getGridId(ii, validBlockCount[ii] - 1);
				if (v->terminal == ZERO) {
					v->discovered = 1;
					v->active = 1;
					Q.push(v);
				}
			}
			
			for (typename std::list<block*>::iterator it = repairList.begin(); it != repairList.end(); ++it) {
				block* v = *it;
#ifdef MEMF_DEBUG
				assert(v->repair);	
#endif
				if (v->discovered && !v->active && (v->prev != NONE || v->terminal == ZERO)) {// discovered but not active, with valid prev
#ifdef MEMF_DEBUG
					assert(v->discovered && !v->active);
#endif
					v->active = 1;
					v->repair = 0;
					Q.push(v);
				}	
				else {	// discovered, without valid prev
#ifdef MEMF_DEBUG
					assert(v->discovered && v->prev == NONE && v->terminal != ZERO);
#endif
					findValidPrev(v);
					if (!v->repair) {
						if (v->discovered && !v->active && v->terminal != ONE) {
							for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
								arc* a = arcs + aid;
								block* w = blockGrid + a->dest;
								if (!w->discovered) {
									v->active = 1;
									Q.push(v);
								}
							}
						}
					} else {	// not repaired!
						repairNeighbor(v);

						v->active = 0;		// need to remove from Q
						v->discovered = 0;
						v->repair = 0;

						for (nodeid aid = v->firstParent; aid != NONE; aid = arcs[aid].nextParent) {
							arc* a = arcs + aid;
							block* w = blockGrid + a->source;
							if (!w->repair && w->discovered && !w->active && w->terminal != ONE) {
								w->active = 1;
								Q.push(w);
							}
						}
					}
				}
			}

#if MEMF_DEBUG >= 2
			testConsistency(1);
#endif

			repairList.clear();
			augPathList.clear();
		}

#ifdef MEMF_DEBUG
		void testConsistency(int cond)
		{
			for (nodeid ii = 0; ii < imagexy_; ++ii) {
				for (labelid li = 0; li < labels_; ++li) {
					block* v = blockGrid + getGridId(ii, li);
					if (cond == 1) {
						if (!v->valid) {
							assert(v->firstChild == NONE && v->firstParent == NONE);
							break;
						}

						for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
							arc* a = arcs + aid;
							block* w = blockGrid + a->dest;
							bool found = false;
							for (nodeid bid = w->firstParent; bid != NONE; bid = arcs[bid].nextParent) {
								arc* b = arcs + bid;
								block* x = blockGrid + b->source;
								if (x == v) {
									found = true;
									break;
								}
							}
							assert(found);
						}

						if (v->terminal == ZERO) assert(v->discovered);
						if (v->discovered) {
							bool valid = (v->prev != NONE || v->terminal == ZERO);
							assert(valid);
						}
						if (!v->discovered) assert(v->prev == NONE);
						if (v->prev != NONE) {
							assert(blockGrid[v->prev].valid);
							assert(blockGrid[v->prev].discovered);
							block* x = blockGrid + v->prev;
							while (x->prev != NONE) {	// avoid cycle! 
								x = blockGrid + x->prev;
								assert(x != v);
							}
						}
						if (v->prev == NONE && v->terminal != ZERO) assert(!v->discovered);
						if (v->discovered && !v->active && v->terminal != ONE) {
							for (nodeid aid = v->firstChild; aid != NONE; aid = arcs[aid].nextChild) {
								arc* a = arcs + aid;	
								block* w = blockGrid + a->dest;
								assert(w->discovered);
							}
						}
					}
					else if (cond == 2) {
						if (!v->discovered) {
							for (nodeid aid = v->firstParent; aid != NONE; aid = arcs[aid].nextParent) {
								arc* a = arcs + aid;
								block* w = blockGrid + a->source;
								assert(!w->discovered);
							}
						}
					}
				}
			}
		}
#endif

		// return the cut point for the nodeid ii (change from non-visited to visited)
		labelid getCutLabel(nodeid ii) {
			block* column_i = blockGrid + getGridId(ii, 0);
			for (labelid li = 0; li < validBlockCount[ii]; ++li) {
				block* v = column_i + li;
				if (v->discovered)  return v->lcheck - 1;
			}
			return column_i[validBlockCount[ii] - 1].lhat;
		}
	};
	////////////////////

	typedef typename ReducedGraph::block RGblock;
	typedef typename ReducedGraph::arc RGarc;

	ReducedGraph* grp;
	std::vector<nodeid> changedList;

	// returns true if neighbor exists
	inline bool neighborExists(nodeid ii, nodeid dir)
	{
		nodeid i = ii / width;	// ver-id
		nodeid j = ii % width;	// hor-id
		switch (dir){
		case 0: return (j != width - 1);	// right
		case 1: return (i != height - 1);	// down
		case 2: return (j != 0);			// left
		case 3:	return (i != 0);			// up
		default: 
#ifdef MEMF_DEBUG
				assert(0);
#endif
				break;
		}
		return false;
	}

	// return neighbor id
	inline nodeid getNeighborNode(nodeid ii, nodeid dir)
	{
		// assumes i has the desired neighbor!
		switch (dir){
		case 0: return ii + 1;		// right
		case 1: return ii + width;	// down
		case 2: return ii - 1;		// left
		case 3: return ii - width;	// up
		default:  
#ifdef MEMF_DEBUG
				assert(0);
#endif
				break;
		}
		return 0;
	}

	// returns the unary id (ii,x_i)
	inline nodeid getUnaryId(nodeid ii, labelid xi)
	{
		return ii * labels + xi;
	}

	// returns the label id (xi, xj)
	inline nodeid getLabelId(labelid xi, labelid xj)
	{
		return xi * labels + xj;
	}

	inline nodeid getPhiId(labelid xi, labelid xj)
	{
		return xi * (labels + 1) + xj;
	}

	inline bool swapDirection(nodeid ii, nodeid jj)
	{
		return (jj < ii);
	}

	// returns the edge id (ii,jj)
	inline nodeid getEdgeId(nodeid ii, nodeid jj)
	{
		if (jj == ii + 1) return ii;
		else if (jj == ii + width) return ii + imagexy;
		else if (jj == ii - 1) return jj;
		else if (jj == ii - width) return jj + imagexy;
#ifdef MEMF_DEBUG
		assert(0);
#endif
		return 0;
	}

	// returns the message id (ii, jj, xx) [mji, mij]
	inline nodeid getMessageId(nodeid ii, nodeid jj, labelid xx)
	{
		nodeid e = getEdgeId(ii, jj);
		labelid ll = xx;
		if (swapDirection(ii, jj)) {
			ll = (xx >= labels) ? xx - labels : xx + labels;
		}
		return e * 2 * labels + ll;
	}

	// clear the given captype array
	inline void clearArray(captype* arr, nodeid len)
	{
		for (nodeid i = 0; i < len; ++i) {
			arr[i] = 0;
		}
	}

	// check if all are zero
	inline bool allZero(captype* arr, labelid len)
	{
		for (labelid i = 0; i < len; ++i) {
			if (arr[i]) return false;
		}
		return true;
	}

	// subtract minimum unary corresponding to ii, updates theta_i
	inline void trivialAugmentation(nodeid ii)
	{
		captype* th_i = theta_i + getUnaryId(ii, 0);
		captype m = *std::min_element(th_i, th_i + labels);
		if (m) {
			for (labelid li = 0; li < labels; ++li) {
				th_i[li] -= m;
			}
			flow += m;
			++titer;
		}
	}

	// return init ishikawa params
	inline void ishikawaParamsInit_ij(nodeid ii, nodeid jj, captype* phi_ij, captype* phi_ji)
	{
		bool swap = swapDirection(ii, jj);
		captype w = gamma_ij[getEdgeId(ii, jj)];

		for (labelid xi = 1; xi < labels; ++xi) {
			for (labelid xj = 1; xj < labels; ++xj) {
				nodeid ll_ij = getPhiId(xi, xj);
				nodeid ll_ji = getPhiId(xj, xi);
				if (xi < xj) {
					phi_ij[ll_ij] = swap ? w * gpp_ji[ll_ji] : w * gpp_ij[ll_ij];
				}
				else if (xi == xj) {
					phi_ij[ll_ij] = swap ? w * gpp_ji[ll_ji] : w * gpp_ij[ll_ij];
					phi_ji[ll_ji] = swap ? w * gpp_ij[ll_ij] : w * gpp_ji[ll_ji];
				}
				else {
					phi_ji[ll_ji] = swap ? w * gpp_ij[ll_ij] : w * gpp_ji[ll_ji];
				}
			}
		}
	}

	// calculate first order difference of m_pq
	inline void derivative(captype* m_pq, captype* f_pq)
	{
		for (labelid xq = 1; xq < labels; ++xq) {
			f_pq[xq - 1] = m_pq[xq] - m_pq[xq - 1];
		}
	}

	// add ishikawa unary params for pp
	inline void addUnaryParams(nodeid pp, captype* phi_p)
	{
		captype* th_p = theta_i + getUnaryId(pp, 0);
		captype *m_p[4];
		for (nodeid dir = 0; dir < 4; ++dir) {	// 4-connected
			if (neighborExists(pp, dir)) {
				nodeid qq = getNeighborNode(pp, dir);
				m_p[dir] = M + getMessageId(pp, qq, 0);	// m_qp
			}
			else {
				m_p[dir] = 0;
			}
		}

		for (labelid xp = 0; xp < labels; ++xp) {
			phi_p[xp] += th_p[xp];
			for (nodeid dir = 0; dir < 4; ++dir) {
				if (m_p[dir]) {
					phi_p[xp] += m_p[dir][xp];	// m_qp:xp
				}
			}
		}
	}

	// deduct m from all elements
	inline void trivialAugmentation(nodeid pp, captype m, captype* phi_p)
	{
		for (labelid lp = 0; lp < labels; ++lp) {
			phi_p[lp] -= m;
		}
		trivialAugmentation(pp, m);
	}

	// add blocks for the node pp
	inline void addBlocks(nodeid pp, captype* phi_p)
	{
#ifdef MEMF_DEBUG
		assert(grp);
#endif
		RGblock* bs = grp->blockGrid + grp->getGridId(pp, 0);
		labelid lcheck = 1;	// start with 1 not 0
		labelid lhat;
		labelid l1 = labels - 1;
		for (labelid lp = 1; lp < l1; ++lp) {	// 1 --> labels - 2
			if (phi_p[lp]) continue;
#ifdef MEMF_DEBUG
			assert(phi_p[lp] == 0);
#endif
			lhat = lp;
			unsigned short terminal = 0;
			if (lcheck == 1 && phi_p[0]) terminal = ONE;
			grp->addBlock(bs, pp, lcheck, lhat, terminal);
			lcheck = lhat + 1;
		}
		lhat = l1;
		unsigned short terminal = 0;
		if (lcheck == 1 && phi_p[0]) terminal = ONE;
		else if (phi_p[l1]) terminal = ZERO;
		grp->addBlock(bs, pp, lcheck, lhat, terminal);
	}

	// add edges between blocks for the edge pp --> qq
	inline void addBlockEdges(nodeid pp, nodeid qq, captype* phi_pq)
	{
		// no vertical edges in the reduced graph!!
		// not an efficient way? 
#ifdef MEMF_DEBUG
		assert(grp);
#endif
		nodeid e = getEdgeId(pp, qq);
		bool swap = swapDirection(pp, qq);
		labelid lastId = grp->validBlockCount[qq];
		for (labelid lg = grp->validBlockCount[pp] - 1; lg >= 0; --lg) {
			bool edgeAdded = false;			
			nodeid pid = grp->getGridId(pp, lg);
			RGblock* bp = grp->blockGrid + pid;
			for (labelid ld = 0; ld < grp->validBlockCount[qq]; ++ld) {
				nodeid qid = grp->getGridId(qq, ld);
				RGblock* bq = grp->blockGrid + qid;
				if (ld == lastId) {	// add an edge, possible to reach from block above (going up and across)
					grp->addEdge(bp, bq, pid, qid, e, swap);
					edgeAdded = true;
					break;
				}
				for (labelid lgg = bp->lcheck; lgg <= bp->lhat; ++lgg) {
					for (labelid ldd = bq->lcheck; ldd <= bq->lhat; ++ldd) {
						if (phi_pq[getPhiId(lgg, ldd)]) {
							grp->addEdge(bp, bq, pid, qid, e, swap);
							lastId = ld;
							edgeAdded = true;
							break;
						}
					}
					if (edgeAdded) break;
				}
				if (edgeAdded) break;
			}
		}
	}

	// subtract m from unary potentials corresponding to ii, updates theta_i
	inline void trivialAugmentation(nodeid ii, captype m)
	{
		captype* th_i = theta_i + getUnaryId(ii, 0);
		for (labelid li = 0; li < labels; ++li) {
			th_i[li] -= m;
		}
		flow += m;
	}

	// do forward reparametrization
	inline void reparametrize(nodeid ii, labelid li, nodeid jj, labelid lj)
	{
		nodeid e = getEdgeId(ii, jj);
		captype* m_ji = M + getMessageId(ii, jj, 0);			// length: labels
		captype* m_ij = M + getMessageId(ii, jj, labels);

		bool swap = swapDirection(ii, jj);
		captype w = gamma_ij[e];
		
		captype* th_i = theta_i + getUnaryId(ii, 0);
		captype *m_i[4];
		for (nodeid dir = 0; dir < 4; ++dir) {	// 4-connected
			if (neighborExists(ii, dir)) {
				nodeid kk = getNeighborNode(ii, dir);
				if (kk != jj) m_i[dir] = M + getMessageId(ii, kk, 0);	// m_ki
				else m_i[dir] = 0;
			}
			else {
				m_i[dir] = 0;
			}
		}
		for (labelid xi = li; xi < labels; ++xi) {
			m_ji[xi] = -th_i[xi];
			for (nodeid dir = 0; dir < 4; ++dir) {
				if (m_i[dir]) {
					m_ji[xi] -= m_i[dir][xi];	// m_ki:xi
				}
			}
		}

		// m'_ij:xj = \min_{xi} theta_ij:xixj - m'_ji:xi
		for (labelid xj = lj; xj < labels; ++xj) {
			m_ij[xj] = std::numeric_limits<captype>::max();
			for (labelid xi = 0; xi < labels; ++xi) {
				captype th_ij = swap ? w * g[getLabelId(xj, xi)] : w * g[getLabelId(xi, xj)];
				m_ij[xj] = std::min(m_ij[xj], (th_ij - m_ji[xi]));
			}
		}

//		// backward --> make-canonical
//		// m'_ji:xi = \min_{xj} theta_ij:xixj - m'_ij:xj
//		for (labelid xi = li; xi < labels; ++xi) {
//			m_ji[xi] = std::numeric_limits<captype>::max();
//			for (labelid xj = 0; xj < labels; ++xj) {
//				captype th_ij = swap ? w * g[getLabelId(xj, xi)] : w * g[getLabelId(xi, xj)];
//				m_ji[xi] = std::min(m_ji[xi], (th_ij - m_ij[xj]));
//			}
//		}
	}

	// do backward reparametrization
	inline void backReparametrize(nodeid ii, labelid li, nodeid jj, labelid lj)
	{
		nodeid e = getEdgeId(ii, jj);
		captype* m_ji = M + getMessageId(ii, jj, 0);			// length: labels
		captype* m_ij = M + getMessageId(ii, jj, labels);

		bool swap = swapDirection(ii, jj);
		captype w = gamma_ij[e];

		// m'_ij:xj = \min_{xi} theta_ij:xixj - m'_ji:xi
		for (labelid xj = lj; xj < labels; ++xj) {
			m_ij[xj] = std::numeric_limits<captype>::max();
			for (labelid xi = 0; xi < labels; ++xi) {
				captype th_ij = swap ? w * g[getLabelId(xj, xi)] : w * g[getLabelId(xi, xj)];
				m_ij[xj] = std::min(m_ij[xj], (th_ij - m_ji[xi]));
			}
		}
	}

#ifdef MEMF_DEBUG

	void assertUpperTriangular(captype* arr)
	{
		for (labelid i = 1; i <= labels; ++i) {
			for (labelid j = 0; j < i; ++j) {
				assert(arr[getPhiId(i, j)] == 0);
			}
		}
	}

	void assertPositive(captype* arr)
	{
		for (labelid i = 1; i <= labels; ++i) {
			for (labelid j = 0; j < labels; ++j) {
				assert(arr[getPhiId(i, j)] >= 0);
			}
		}
	}

	void assertPositive(captype* arr, labelid len)
	{
		for (labelid i = 0; i < len; ++i) {
			assert(arr[i] >= 0);
		}
	}	

	void assertZero(captype* arr, labelid len)
	{
		for (labelid i = 0; i < len; ++i) {
			assert(arr[i] == 0);
		}
	}

	void printArray(captype* arr, labelid len)
	{
		for (labelid i = 0; i < len; ++i) {
			fout << "\t" << arr[i];
		}
		fout << std::endl;
		fout.flush();
	}

	void printArray(captype* arr, labelid len1, labelid len2)
	{
		for (labelid i = 0; i < len1; ++i) {
			for (labelid j = 0; j < len2; ++j) {
				fout << "\t" << arr[i * len2 + j];
			}
			fout << std::endl;
		}
		fout << std::endl;
		fout.flush();
	}
#endif

	// construct the initial reduced graph
	void construct(bool changed = false);

	// construct the reduced graph given augpath
	void construct(std::vector<RGblock*>&);
		
	// make ishikawa edge params positive - init
	void makeEdgePositiveInit(nodeid, nodeid);

	// make ishikawa edge params positive - directly
	void makeEdgePositiveDirect(nodeid, nodeid);

	// construct reduced graph portion
	void reduceGraph(nodeid, nodeid);

	// do reparametrization starting from startBlock
	void augment(RGblock*);

	// do reparametrization on augPath
	void augment(std::vector<RGblock*>&);

	// remove unusable blocks from aug-path
	void trimAugPath(std::vector<RGblock*>&, RGblock*);

	// reset dummy flag
	void resetFlag(RGblock*);

	// find the final labelling
	void finalLabelling();

	// initialize optimization
	void init();

	// initialize reuse optimization (optiter > 0)
	void initReuse();
	
};

