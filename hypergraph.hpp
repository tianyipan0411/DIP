#ifndef _HYPERGRAPH_H_
#define _HYPERGRAPH_H_

#include "rwgraph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include <algorithm> 
#include "mappedHeap.hpp"
#include "HeapData.hpp"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

using namespace std;

/*
* building the hypergraph procedure which generates hyperedges following IC model
*/
void addHyperedge(vector<Graph> mtg, vector<HyperGraph> & hg, int t, long long num, bool lt, vector<double> speed_thres, 
					vector<double> cov_thres, vector<bool> speed_sample, vector<bool> cov_sample, double time_threshold, double dec)
{
	int numNodes = mtg[0].getSize();

	omp_set_num_threads(t);

	int base_sample = num - hg[0].getNumEdge_base();
	if(base_sample > 0)
	{
		int sample_per_thread = base_sample/t;
		vector<bool> graph_assignment (t,false);
		#pragma omp parallel
		{
			int graph_index = 0;
			double th_local = 0.0;
			double dec_local = 0.0;
			#pragma omp critical
			{
				for(unsigned int index =0; index < graph_assignment.size(); ++index){
					if(!graph_assignment[index])
					{
						graph_index = index;
						graph_assignment[index] = true;
						break;
					}
				}
				th_local = time_threshold;
				dec_local = dec;
			}

			Graph g = mtg[graph_index];
			vector<NodeTime> visit_mark(numNodes+1,{0,0.0});
			vector<bool> visit(numNodes+1,false);
			for (int i = 0; i < sample_per_thread; ++i){
				hg[graph_index].pollingIC_base(g,visit, visit_mark, th_local, dec_local);
			}

			//cout<<graph_index<<" finish!"<<endl;	
		}
		for(int i=1; i < t; ++i){
			hg[0].addBaseEdges(hg[i].getBaseEdges());
			hg[i].clearAllEdges();
		}
	}

	for(unsigned int sped_index = 0; sped_index < speed_thres.size(); ++sped_index)
	{
		//if sample of the speed index is needed
		if(speed_sample[sped_index])
		{
			long long speed_sample = num - hg[0].getNumEdge_sped(sped_index);
			vector<int> visit_mark(numNodes+1,0);
		    for (int i = 0; i < speed_sample; ++i){
				hg[0].pollingIC_speedup(visit_mark,speed_thres[sped_index],sped_index);
			}
		}
		hg[0].updateDeg_sped(sped_index);
	}

	for(unsigned int cov_index = 0; cov_index < cov_thres.size(); ++cov_index)
	{
		// if sample of the cov index is needed
		if(cov_sample[cov_index])
		{
			long long cov_sample = num - hg[0].getNumEdge_cov(cov_index);
			vector<int> visit_mark(numNodes+1,0);
			for (int i = 0; i < cov_sample; ++i){
				hg[0].pollingIC_coverage(visit_mark,cov_thres[cov_index],cov_index);
			}
		}
		hg[0].updateDeg_cov(cov_index);
	}
	//return hg.getNumEdge();
}

/*
* find seed nodes procedure using greedy algorithm
*/
void buildSeedSet(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, 
					vector<double> speed_ratio, vector<double> cov_ratio,  vector<int> &speed_deg, vector<int> &cov_deg)
{	
	unsigned int j,l,m,maxInd;
	vector<int> e, nList;

	vector<int> nodeDegree(n,0);
	// effective gain considered in the heap
	vector<vector<int>> nodeDegree_speed;
	vector<vector<int>> nodeDegree_cov;
	// number of edges not included in the seed set per node
	vector<vector<int>> inc_edge_speed;
	vector<vector<int>> inc_edge_cov;
	vector<int> indx(n,0);
	vector<vector<int>> speed_indx;
	vector<vector<int>> cov_indx;
	vector<int> speed_thres;
	vector<int> cov_thres;
	for(l=0;l<speed_ratio.size();++l){
		speed_thres.push_back(((double)hg.getNumEdge_sped(l))*speed_ratio[l]);
		nodeDegree_speed.push_back(vector<int>(n,0));
		inc_edge_speed.push_back(vector<int>(n,0));
		speed_indx.push_back(vector<int>(n,0));
	}
	
	for(l=0;l<cov_ratio.size();++l){
		cov_thres.push_back(((double)hg.getNumEdge_cov(l))*cov_ratio[l]);
		nodeDegree_cov.push_back(vector<int>(n,0));
		inc_edge_cov.push_back(vector<int>(n,0));
		cov_indx.push_back(vector<int>(n,0));
	}
	//cout<<"sped thre: "<<speed_thres[0]<<" cov thre: "<<cov_thres[0]<<endl;
	for (j = 0; j < n; ++j){
		indx[j] = j;
		for(l=0;l<speed_ratio.size();++l){
			speed_indx[l][j]=j;
			int size = hg.getNode_sped(j+1, l).size();
			int gain =  min(speed_thres[l], size);
			nodeDegree[j] += gain;
			nodeDegree_speed[l][j] = gain;
			inc_edge_speed[l][j] = size;
		}
		for(l=0;l<cov_ratio.size();++l){
			cov_indx[l][j]=j;
			int size = hg.getNode_cov(j+1,l).size();
			int gain = min(cov_thres[l],size);
			nodeDegree[j] += gain;
			nodeDegree_cov[l][j] = gain;
			inc_edge_cov[l][j] = size;
		}
	}
	InfCost<int> hd(&nodeDegree[0]);
	vector<InfCost<int>> hd_speed;
	vector<InfCost<int>> hd_cov;
	
	// two vector of heaps to calculate the max degree for each hyperedge set
	vector<MappedHeap<InfCost<int> >> heap_speed;
	vector<MappedHeap<InfCost<int> >> heap_cov;

	// check if an edge is removed
	vector<vector<bool>> edgeMark_speed;
	vector<vector<bool>> edgeMark_cov;	
	
	for(l=0;l<speed_ratio.size();++l){
		hd_speed.push_back(InfCost<int>(&nodeDegree_speed[l][0]));
		MappedHeap<InfCost<int> > tmp_heap(speed_indx[l],hd_speed[l]);
		heap_speed.push_back(tmp_heap);
		edgeMark_speed.push_back(vector<bool>(hg.getNumEdge_sped(l),false));
	}
	
	for(l=0;l<cov_ratio.size();++l){
		hd_cov.push_back(InfCost<int>(&nodeDegree_cov[l][0]));
		MappedHeap<InfCost<int> > tmp_heap(cov_indx[l],hd_cov[l]);
		heap_cov.push_back(tmp_heap);
		edgeMark_cov.push_back(vector<bool>(hg.getNumEdge_cov(l),false));
	}
	
	// node heap with adj. edges as gain. the main one, combine all edge sets and considering the ratios 
	// (bound the gain from each edge set)
	MappedHeap<InfCost<int> > heap(indx,hd);

/* 	int index1 = heap_speed[0].pop();
	int index2 = heap_speed[0].top();
	nodeDegree_speed[0][index2]=0;
	heap_speed[0].push(index1);
	heap_speed[0].heapify(index2);
	heap_speed[0].heapify(index1);
	cout<<"large"<<nodeDegree_speed[0][index1]<<"small"<<nodeDegree_speed[0][index2]<<"new "<<nodeDegree_speed[0][heap_speed[0].top()]<<endl; 
	//long long numEdge = hg.getNumEdge(); */


	//vector<bool> edgeMark(numEdge, false);
	vector<bool> nodeMark(n+1, true);
	
	double totalCost = 0;

	//i=1;

	// building each seed at a time
	while(totalCost < k && !heap.empty()){
		// index of heap one larger than nodeDegree
		maxInd = heap.pop()+1;
		nodeMark[maxInd] = false;

		totalCost++;
		//// degree is for storing greedy coverage after each iteration
		//degree[i] = degree[i-1]+nodeDegree[maxInd-1];
		//totalDegree = degree[i];
		seeds.push_back(maxInd);
		// index of hyperedges
		
		for(l=0;l<speed_ratio.size();++l){
			speed_deg[l] += nodeDegree_speed[l][maxInd-1];
			speed_thres[l] -= nodeDegree_speed[l][maxInd-1];
 /* 			if(speed_thres[0] < 0)
			{
				cout<<"node "<<maxInd<<" degree "<< nodeDegree_speed[l][maxInd-1]<<endl;
				cout<<"top node "<< heap_speed[l].top() <<" degree "<<nodeDegree_speed[l][heap_speed[l].top()]<<endl;
				cout <<"speed!!"<<speed_thres[0]<<endl;
			}  */
			e = hg.getNode_sped(maxInd,l);

			for (j = 0; j < e.size(); ++j){
				if (edgeMark_speed[l][e[j]]){
					continue;
				}
				// nList is a hyperedge
				nList = hg.getEdge_sped(e[j],l);
				for (m = 0; m < nList.size(); ++m){
					inc_edge_speed[l][nList[m]-1]--;
					if(inc_edge_speed[l][nList[m]-1]<nodeDegree_speed[l][nList[m]-1])
					{
						nodeDegree[nList[m]-1]--;
						nodeDegree_speed[l][nList[m]-1]--;
						heap_speed[l].heapify(nList[m]-1);
						if (nodeMark[nList[m]]){
							heap.heapify(nList[m]-1);
						}
					}

				}
				edgeMark_speed[l][e[j]] = true;
			}
		

			//alter gain for those whose inclusion will exceed the required threshold
			int ind = heap_speed[l].top();
			while(speed_thres[l] < nodeDegree_speed[l][ind])
			{
				nodeDegree[ind] -= nodeDegree_speed[l][ind] - speed_thres[l];
				nodeDegree_speed[l][ind] = speed_thres[l];
				if(nodeMark[ind+1]){
					heap_speed[l].heapify(ind);
					heap.heapify(ind);
				}
				else{
					heap_speed[l].pop();
					//cout<<"here?"<<endl;
				}
				//TODO: check why this is happening... perhaps index issues?
				if(nodeDegree_speed[l][heap.top()] > nodeDegree_speed[l][heap_speed[l].top()])
				{
					//cout<<"!!!!"<<endl;
					heap_speed[l].heapify(heap.top());
				}
				ind = heap_speed[l].top();
			}
/* 			 if(speed_thres[l] < nodeDegree_speed[l][heap.top()])
			{
				cout<<"cur thres: "<<speed_thres[l]<<"largest deg: "<<nodeDegree_speed[l][heap.top()]<<endl;
			}  */
		}
		
		for(l=0;l<cov_ratio.size();++l){
			cov_deg[l] += nodeDegree_cov[l][maxInd-1];
			cov_thres[l] -= nodeDegree_cov[l][maxInd-1];
/* 			if(cov_thres[0] < 0)
			{
				cout <<"cov!!"<<speed_thres[0]<<endl;
			} */
			e = hg.getNode_cov(maxInd,l);

			for (j = 0; j < e.size(); ++j){
				if (edgeMark_cov[l][e[j]]){
					continue;
				}
				// nList is a hyperedge
				nList = hg.getEdge_cov(e[j],l);
				for (m = 0; m < nList.size(); ++m){
					inc_edge_cov[l][nList[m]-1]--;
					// dec the heap related data here only if the number of new edges (inc_edge_cov) is less than what is considered
					if(inc_edge_cov[l][nList[m]-1]<nodeDegree_cov[l][nList[m]-1])
					{
						nodeDegree[nList[m]-1]--;
						nodeDegree_cov[l][nList[m]-1]--; 
						heap_cov[l].heapify(nList[m]-1);
						if (nodeMark[nList[m]]){
							heap.heapify(nList[m]-1);
						}						
					}
					
				}
				edgeMark_cov[l][e[j]] = true;
			}
			

  
			//alter gain for those whose inclusion will exceed the required threshold
			int ind = heap_cov[l].top();
			while(cov_thres[l] < nodeDegree_cov[l][ind])
			{
				nodeDegree[ind] -= nodeDegree_cov[l][ind] - cov_thres[l];
				nodeDegree_cov[l][ind] = cov_thres[l];
				if(nodeMark[ind+1]){
					heap_cov[l].heapify(ind);
					heap.heapify(ind);
				}
				else{
					heap_cov[l].pop();
				}
				ind = heap_cov[l].top();
			}

		}
		//i++;
	}
	for(l=0;l<speed_ratio.size();++l){
		vector<int>().swap(nodeDegree_speed[l]);
		vector<bool>().swap(edgeMark_speed[l]);
	}
	
	for(l=0;l<cov_ratio.size();++l){
		vector<int>().swap(nodeDegree_cov[l]);
		vector<bool>().swap(edgeMark_cov[l]);
	}
	vector<int>().swap(nodeDegree);
	vector<int>().swap(e);
	vector<int>().swap(nList);
	vector<int>().swap(indx);
	//vector<bool>().swap(edgeMark);
}

void CheckHeap()
{
	_CheckHeap heap;
	heap.Check();
}
#endif
