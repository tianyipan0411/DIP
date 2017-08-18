#include "rwgraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <algorithm>
#include "mappedHeap.hpp"
#include "HeapData.hpp"

using namespace std;


const vector<int> & Graph::operator [] (int u) const
{
	return adjList[u];
}


const vector<int> & Graph::operator [] (int u)
{
	return adjList[u];
}


const vector<UI> & Graph::getWeight (int u) const
{
        return weights[u];
}

const vector<UI> & Graph::getWeight (int u)
{
        return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const
{
	return adjList[u].size();
}

/*
* get the max degree
*/
int Graph::getMaxDegree() const
{
	return maxDegree;
}

/*
* get the number of nodes
*/
int Graph::getSize() const
{
	return numNodes;
}

/*
* get the number of edges
*/
int Graph::getEdge() const
{
	return numEdges;
}

/*
* read binary graph input for LT model
* difference between LT and IC is for LT we accumulate the weights for fast choosing a random node
*/
void Graph::readGraphLT(const char* filename)
{
   	FILE * pFile;
    	pFile = fopen(filename, "rb");
    	fread(&numNodes, sizeof(int), 1, pFile);
    	fread(&numEdges, sizeof(long long), 1, pFile);
    	node_deg=vector<int>(numNodes + 1);
    	fread(&node_deg[1], sizeof(int), numNodes, pFile);
        
	vector<int> a;
    	vector<UI> b;
    	adjList.push_back(a);
    	weights.push_back(b);
	
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);

                adjList.push_back(tmp);
        }

        for (unsigned int i = 1; i <= numNodes; ++i){
		vector<float> tmp(node_deg[i] + 1, 0);
                vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp[j] += tmp[j-1];
                        if (tmp[j] >= 1){
                                tmp1[j] = UI_MAX;
                        } else {
                                tmp1[j] = tmp[j]*UI_MAX;
                        }
                }

                weights.push_back(tmp1);
                node_deg[i]++;
        }
}

/*
* read input graph for IC model
*/
void Graph::readGraphIC(const char* filename)
{
    	FILE * pFile;
    	pFile = fopen(filename, "rb");
    	fread(&numNodes, sizeof(int), 1, pFile);
    	fread(&numEdges, sizeof(long long), 1, pFile);
    	node_deg=vector<int>(numNodes + 1);
    	fread(&node_deg[1], sizeof(int), numNodes, pFile);
	vector<int> a;
	vector<UI> b;
    	adjList.push_back(a);
    	weights.push_back(b);
	
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                adjList.push_back(tmp);
        }

        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<float> tmp(node_deg[i] + 1, 0);
		vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp1[j] = tmp[j]*UI_MAX;
                }

		if (tmp1[node_deg[i]] <= 0)
                        tmp1[node_deg[i]] = UI_MAX;
		
                weights.push_back(tmp1);
        }
		maxDegree = 0;
		for(unsigned i = 1; i <= numNodes; ++i)
		{
			if(maxDegree < getDegree(i)) maxDegree = getDegree(i);
		}
}

void Graph::writeToFile(const char * filename)
{/*
	ofstream output(filename);
	for (unsigned int i = 0; i < numNodes; ++i){
		for (unsigned int j = 0; j < adjList[i].size(); ++j){
			if (adjList[i][j] > i){
				output << adjList[i][j] << " " << i << " " << weights[i][j] << endl;
			}
		}
	}
	output.close();
*/	
}



bool HyperGraph::verifySeedSet(Graph &g, set<int> seedSet, double speed_rate, double speed_req, double cov_req, int &num_inf_nodes, double & sped_time)
{
	vector<double> inf_times(g.getSize(), 10);
	vector<int> indx(g.getSize(), 0);
	set<int> inf_nodes;
	for (int j = 0; j < g.getSize(); ++j){
		indx[j] = j;
		if(seedSet.count(j+1))
		{
			inf_times[j]=0;
			//cout<<"!!"<<endl;
			inf_nodes.insert(j+1);
		}
	}
	TimeInf<double> hd(&inf_times[0]);
	MappedHeap<TimeInf<double>> heap_inf_times(indx,hd);
	double delay_dec = 1;
	double max_time = 10;
	double cur_time = 0;
	int cur_inf_nodes = seedSet.size();
	bool speed_changed = false;
	set<int> poped_nodes;
	while(cur_time < max_time && !heap_inf_times.empty())
	{
		int cur_id = heap_inf_times.pop();
		poped_nodes.insert(cur_id);
		cur_time = inf_times[cur_id];
		if(!speed_changed && poped_nodes.size() >= speed_req)
		{
			delay_dec = 1 / speed_rate;
			speed_changed = true;
			sped_time = cur_time;
		}
		//cout<<"time:"<<cur_time<<" neighbors: "<<g[cur_id+1].size()<<"deg: "<<g.node_deg[cur_id+1]<<endl;
		// the edge weights
		const vector<UI> &w=g.getWeight(cur_id+1);
		// the neighbors
		const vector<int> &neigh = g[cur_id+1];
		for (int i = 0; i < g.node_deg[cur_id+1]; ++i){
			//cout<<(1.0/((double)g[neigh[i]].size()))<<endl;
			if ( sfmt_genrand_real1(&sfmtSeed) < (1.0/((double)g[neigh[i]].size()))){
				int neighbor_id = neigh[i];
				double delay = distribution(generator);
				double inf_time = delay * delay_dec + inf_times[cur_id];
				//cout<<"delay: "<<delay<<" inf_time: "<<inf_time<<endl;
				if (inf_time < max_time){
					if(!inf_nodes.count(neighbor_id))
					{
						cur_inf_nodes++;
						inf_nodes.insert(neighbor_id);
					}
					if(inf_times[neighbor_id-1]>inf_time)
					{
						inf_times[neighbor_id-1] = inf_time;
						//cout<<"inf_time: "<<inf_time<<endl;
						heap_inf_times.heapify(neighbor_id-1);
					}
				}
			}
		}
	}
	//cout<<"nodes: "<<cur_inf_nodes<<"total nodes: "<<g.getSize()<<endl;
	num_inf_nodes = cur_inf_nodes;
	return cur_inf_nodes  >= cov_req * (1-0.1);
}


// choose a random edge in LT model based on linear search
inline int HyperGraph::randIndex_lin(const vector<UI> &w, unsigned int si)
{
        UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ranNum > w[si - 1])
                return -1;

        for (unsigned int i = 1; i < si; ++i){
                if (ranNum <= w[i])
                        return i;
        }
        return -1;
}

// choose a random live edge in LT model based on binary search
inline int HyperGraph::randIndex_bin(const vector<UI> &w, unsigned int si)
{
	UI ran = sfmt_genrand_uint32(&sfmtSeed);
	if (si <= 1 || ran > w[si - 1])
                return -1;
        int left = 1;
        int right = si - 1;
        int prob;
        for (unsigned int i = 0; i < si; ++i){
                prob = (left + right)/2;
                if (w[prob - 1] > ran){
                        right = prob - 1;
                        continue;
                }
                if (w[prob] <= ran){
                        left = prob + 1;
                        continue;
                }
                break;
        }
        return prob;
}


HyperGraph::HyperGraph(unsigned int n, int speed_up_times, int cov_reqs)
{
	sfmt_init_gen_rand(&sfmtSeed, rand());
	//generator.seed (sfmt_genrand_uint32(&sfmtSeed));
	sped_node_edge = vector<vector<vector<int>>>(speed_up_times);
	for(unsigned int i = 0; i < sped_node_edge.size(); i++){
		sped_node_edge[i] = vector<vector<int>>(n+1);
	}
	sped_edge_node = vector<vector<vector<int>>>(speed_up_times);
	cov_node_edge = vector<vector<vector<int>>>(cov_reqs);
	for(unsigned int i = 0; i < cov_node_edge.size(); i++){
		cov_node_edge[i] = vector<vector<int>>(n+1);
	}
	cov_edge_node = vector<vector<vector<int>>>(cov_reqs);
	maxDegree = 0;
	numNodes = n;
	sped_curEdge=vector<int> (speed_up_times,0);
	cov_curEdge=vector<int> (cov_reqs,0);
	distribution = std::weibull_distribution<double>(1.0,4.0);
}

void HyperGraph::updateDeg_sped(int index){
	unsigned int num=sped_edge_node[index].size();
	for (unsigned int i = sped_curEdge[index]; i < num; ++i){
		unsigned int num2 = sped_edge_node[index][i].size();
		for (unsigned int j=0;j<num2;++j){
			sped_node_edge[index][sped_edge_node[index][i][j]].push_back(i);
		}
	}
	sped_curEdge[index] = sped_edge_node[index].size();
}

void HyperGraph::updateDeg_cov(int index){
	unsigned int num=cov_edge_node[index].size();
	for (unsigned int i = cov_curEdge[index]; i < num; ++i){
		unsigned int num2 = cov_edge_node[index][i].size();
		for (unsigned int j=0;j<num2;++j){
			cov_node_edge[index][cov_edge_node[index][i][j]].push_back(i);
		}
	}
	cov_curEdge[index] = cov_edge_node[index].size();
}
/* void HyperGraph::updateEdge(int index){
	sped_curEdge[index] = sped_edge_node[index].size();
} */

/* 
// Add a hyperedge into the hypergraph

void HyperGraph::addEdge(vector<int> & edge)
{
	edge_node.push_back(edge);
	unsigned int ind = edge_node.size() - 1;
	for (unsigned int i = 0; i < edge.size(); ++i)
		node_edge[edge[i]].push_back(ind);
}

// Add a hyperedge into the hypergraph while keeping track of the node with max degree

void HyperGraph::addEdgeD(vector<int> & edge)
{
        edge_node.push_back(edge);
        int ind = edge_node.size() - 1;
        for (unsigned int i = 0; i < edge.size(); ++i){
                node_edge[edge[i]].push_back(ind);
		if (node_edge[edge[i]].size() > maxDegree)
			maxDegree = node_edge[edge[i]].size();
	}
} */


const vector<NodeTime> & HyperGraph::getBaseEdge(int e){
	return base_edge_node[e];
}

/*
* get the number of hyperedges
*/
long long HyperGraph::getNumEdge_base() const
{
        return base_edge_node.size();
}

/*
* get an edge from the hypergraph
*/
const vector<int> & HyperGraph::getEdge_sped(int e, int index) const{
	return sped_edge_node[index][e];
}

const vector<int> & HyperGraph::getEdge_sped(int e, int index){
	return sped_edge_node[index][e];
}


/*
* get the list of hyperedges incident to node n
*/
const vector<int> & HyperGraph::getNode_sped(int n, int index) const{
	return sped_node_edge[index][n];
}

const vector<int> & HyperGraph::getNode_sped(int n, int index){
	return sped_node_edge[index][n];
}

/*
* get the number of hyperedges
*/
long long HyperGraph::getNumEdge_sped(int index) const
{
        return sped_edge_node[index].size();
}

/*
* get an edge from the hypergraph
*/
const vector<int> & HyperGraph::getEdge_cov(int e, int index) const{
	return cov_edge_node[index][e];
}

const vector<int> & HyperGraph::getEdge_cov(int e, int index){
	return cov_edge_node[index][e];
}


/*
* get the list of hyperedges incident to node n
*/
const vector<int> & HyperGraph::getNode_cov(int n, int index) const{
	return cov_node_edge[index][n];
}

const vector<int> & HyperGraph::getNode_cov(int n, int index){
	return cov_node_edge[index][n];
}

/*
* get the number of hyperedges
*/
long long HyperGraph::getNumEdge_cov(int index) const
{
        return cov_edge_node[index].size();
}


/*
* get the maximum degree
*/
int HyperGraph::getMaxDegree()
{
	return maxDegree;
}

/*
* remove sped/cov hyperedges
*/
void HyperGraph::clearEdges()
{
	for(unsigned int index = 0; index < sped_edge_node.size();index++){
		vector<vector<int>>().swap(sped_edge_node[index]);
		vector<vector<int>>().swap(sped_node_edge[index]);
		sped_node_edge[index] = vector<vector<int>>(numNodes+1);
		sped_curEdge[index] = 0;
	}
	for(unsigned int index = 0; index < cov_edge_node.size();index++){
		vector<vector<int>>().swap(cov_edge_node[index]);
		vector<vector<int>>().swap(cov_node_edge[index]);
		cov_node_edge[index] = vector<vector<int>>(numNodes+1);
		cov_curEdge[index] = 0;
	}
	//cout << "clear edges!" << endl;
	maxDegree = 0;
}


/*
* remove all the hyperedges
*/
void HyperGraph::clearAllEdges()
{
	vector<vector<NodeTime>>().swap(base_edge_node);
	for(unsigned int index = 0; index < sped_edge_node.size();index++){
		vector<vector<int>>().swap(sped_edge_node[index]);
		vector<vector<int>>().swap(sped_node_edge[index]);
		sped_node_edge[index] = vector<vector<int>>(numNodes+1);
		sped_curEdge[index] = 0;
	}
	for(unsigned int index = 0; index < cov_edge_node.size();index++){
		vector<vector<int>>().swap(cov_edge_node[index]);
		vector<vector<int>>().swap(cov_node_edge[index]);
		cov_node_edge[index] = vector<vector<int>>(numNodes+1);
		cov_curEdge[index] = 0;
	}
	//cout << "clear all edges!" << endl;
	maxDegree = 0;
}
// polling process under LT model

/* bool HyperGraph::pollingLT2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{	
	unsigned int i;
	bool t = false;
        unsigned int gSize = g.getSize();
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
        unsigned int num_marked = 0;
        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
		num_marked++;
		if (link[cur] < k)
			t=true;
                int ind;
	        if (g.weights[cur].size() >= 32)
        	        ind = randIndex_bin(g.weights[cur],g.node_deg[cur]);
                else
                	ind = randIndex_lin(g.weights[cur],g.node_deg[cur]);

                if (ind == -1)
                        break;
                cur = g.adjList[cur][ind - 1];
        }
	edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
	return t;
}

bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        unsigned int i;
        bool t = false;
        unsigned int gSize = g.getSize();
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
        unsigned int num_marked = 0;
        for (i = 0; i < gSize; ++i){
		if (link[cur] < k){
                        t=true;
			break;
                }
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
		num_marked++;
		int ind;
                if (g.weights[cur].size() >= 32)
                        ind = randIndex_bin(g.weights[cur],g.node_deg[cur]);
                else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur]);

                if (ind == -1)
                        break;
                cur = g.adjList[cur][ind - 1];
        }
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
        return t;
}


void HyperGraph::pollingLT1(Graph &g, vector<bool> &visit, vector<int> &visit_mark)
{
        unsigned int i;
        unsigned int gSize = g.getSize();
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize+1;
        unsigned int num_marked = 0;
        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
		num_marked++;
                const vector<int> &neigh = g[cur];
                int ind;
                if (g.weights[cur].size() >= 32)
                        ind = randIndex_bin(g.weights[cur],g.node_deg[cur]);
                else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur]);

                if (ind == -1)
                        break;
                cur = neigh[ind - 1];
        }
	edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
} */

// poll a sample from a base sample. The sample is used for calculating the seed set of ensuring the speed up with given index will happen (assuming all previous speed ups.). 

void HyperGraph::pollingIC_speedup(vector<int> &visit_mark, double modified_threshold, int target_speedup_index)
{
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(base_edge_node.size());
		vector<NodeTime> sample = base_edge_node[cur];

		int L = 0;
		int R = sample.size()-1;
		int curPos = (L+R)*0.5;
		
        while(L < R){
			NodeTime cur_node = sample[curPos];

			if(cur_node.influenceTime > modified_threshold)
				R = curPos - 1;
			else
				L = curPos + 1;
			
			curPos = (L+R)*0.5;
        }
		
		for(i =0; i <= curPos; ++i){
			visit_mark[i] = sample[i].nodeID;
		}
		sped_edge_node[target_speedup_index].push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+curPos+1));

}

// transform the threshold. Current one is assuming normal speed up. Change to the corresponding one with full speed up.
// speedup_time_ints is time interval, the delay_decs is the ratio of speed of certain time over max speed
double HyperGraph::calc_modified_threshold(double time_threshold, vector<double> speedup_time_ints, vector<double> delay_decs)
{
	double time_remain = time_threshold;
	int cur_index = 0;
	double new_threshold = 0;
	while(time_remain > speedup_time_ints[cur_index])
	{
		new_threshold += speedup_time_ints[cur_index] * delay_decs[cur_index];
		time_remain -= speedup_time_ints[cur_index];
		cur_index++;
	}
	new_threshold += time_remain * delay_decs[cur_index];
	return new_threshold;
}

// poll a sample from a base sample. The sample is used for calculating the seed set of ensuring the coverage when the speed up happened at time t. 
// Assuming single speed up for now. This method samples for first t time upder normal speed, remaining time under fast speed. 
void HyperGraph::pollingIC_coverage(vector<int> &visit_mark, double modified_threshold, int cov_index)
{
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(base_edge_node.size());
		vector<NodeTime> sample = base_edge_node[cur];
		
		int L = 0;
		int R = sample.size()-1;
		int curPos = (L+R)*0.5;
		
        // find the node that is just speedup_time away from the origin and include all nodes before it in the sample
        while(L < R){
			NodeTime cur_node = sample[curPos];
			if(cur_node.influenceTime > modified_threshold)
				R = curPos - 1;
			else
				L = curPos + 1;
			
			curPos = (L+R)*0.5;
        }
		
		for(i =0; i <= curPos; ++i){
			visit_mark[i] = sample[i].nodeID;
		}
		cov_edge_node[cov_index].push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+curPos+1));

}

// polling process under IC model: the base case which records time and is the only one that really traverse the graph.
// the delay_dec is applied at time 0. (< 1, which is the reduction of delay, reciprocal of speed up)
void HyperGraph::pollingIC_base(Graph &g, vector<bool> &visit, vector<NodeTime> &visit_mark, double time_threshold, double delay_dec)
{
	int i;
	unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize())+1;
	int curPos=0;
    int num_marked=1;
	visit[cur] = true;
	NodeTime init_node = {cur, 0.0};
	visit_mark[0] = init_node;
	while(curPos < num_marked){
		NodeTime cur_node = visit_mark[curPos];
		cur = cur_node.nodeID;
		curPos++;
		// the edge weights
		const vector<UI> &w=g.getWeight(cur);
		// the neighbors
		const vector<int> &neigh = g[cur];
		for (i = 0; i < g.node_deg[cur]; ++i){
			//cout<<(1.0/((double)neigh.size()))<<endl;
			if ( sfmt_genrand_real1(&sfmtSeed) <  (1.0/((double)neigh.size()))){
				double delay = distribution(generator);
				double inf_time = delay * delay_dec + cur_node.influenceTime;
				if (!visit[neigh[i]] && inf_time < time_threshold){
					visit[neigh[i]] = true;
					NodeTime neigh_node = {(unsigned int)neigh[i], inf_time};
					visit_mark[num_marked] = neigh_node;
					num_marked++;
				}
			}
		}
	}
	
	std::sort(visit_mark.begin(), visit_mark.begin()+num_marked, TimeComp);

	base_edge_node.push_back(vector<NodeTime>(visit_mark.begin(),visit_mark.begin()+num_marked));

	for(i = 0; i < num_marked;++i){
		visit[visit_mark[i].nodeID]=false;
	}
}
void HyperGraph::addBaseEdges(vector<vector<NodeTime> > samples)
{
	base_edge_node.insert(base_edge_node.end(), samples.begin(), samples.end());
}

vector<vector<NodeTime> > HyperGraph::getBaseEdges()
{
	return base_edge_node;
}
/* bool HyperGraph::pollingIC(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize())+1;
	int curPos=0;
        int num_marked=1;
        visit[cur] = true;
        visit_mark[0] = cur;
        bool t = false;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
		curPos++;
                if (link[cur] < k){
                        t=true;
			break;
		}
                const vector<UI> &w=g.getWeight(cur);
                const vector<int> &neigh = g[cur];
                for (i = 0; i < g.node_deg[cur]; ++i){
                        if (sfmt_genrand_uint32(&sfmtSeed) <  w[i+1]){
                        	if (!visit[neigh[i]]){
                                        visit[neigh[i]] = true;
                                        visit_mark[num_marked]=neigh[i];
					num_marked++;
                                }
                        }
                }
        }
        for(i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
        }
        return t;
}

void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark)
{
	int i;
	unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize())+1;
	int curPos=0;
	int num_marked=1;
	visit[cur] = true;
	visit_mark[0] = cur;
	while(curPos < num_marked){
			cur = visit_mark[curPos];
			curPos++;
			const vector<UI> &w=g.getWeight(cur);
			const vector<int> &neigh = g[cur];
			for (i = 0; i < g.node_deg[cur]; ++i){
					if (sfmt_genrand_uint32(&sfmtSeed) <  w[i+1]){
						if (!visit[neigh[i]]){
									visit[neigh[i]] = true;
									visit_mark[num_marked]=neigh[i];
				num_marked++;
							}
					}
			}
	}
	for(i = 0; i < num_marked;++i){
			visit[visit_mark[i]]=false;
	}
}
bool HyperGraph::pollingIC2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize())+1;
	int curPos=0;
        int num_marked=1;
        visit[cur] = true;
        visit_mark[0] = cur;
        bool t = false;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
		curPos++;
                if (link[cur] < k){
                        t=true;
			break;
		}
                const vector<UI> &w=g.getWeight(cur);
                const vector<int> &neigh = g[cur];
                for (i = 0; i < g.node_deg[cur]; ++i){
                        if (sfmt_genrand_uint32(&sfmtSeed) <  w[i+1]){
                        	if (!visit[neigh[i]]){
                                        visit[neigh[i]] = true;
                                        visit_mark[num_marked]=neigh[i];
					num_marked++;
                                }
                        }
                }
        }
        for(i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
        }
        return t;
} */
/*
* convert from an integer to a string
*/
string intToStr(int i) {
        stringstream ss;
        ss << i;
        return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s) {
        unsigned int i;
        istringstream myStream(s);

        if (myStream>>i) {
                return i;
        } else {
                cout << "String " << s << " is not a number." << endl;
                return atoi(s.c_str());
        }
        return i;
}

/*
* measure the consumed memory
*/
float getCurrentMemoryUsage() {

        string pid = intToStr(unsigned(getpid()));
        string outfile = "tmp_" + pid + ".txt";
        string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
        system(command.c_str());

        string mem_str;
        ifstream ifs(outfile.c_str());
        std::getline(ifs, mem_str);
        ifs.close();

        mem_str = mem_str.substr(0, mem_str.size()-1);
        float mem = (float)strToInt(mem_str);

	command = "rm " + outfile;
        system(command.c_str());

        return mem/1024;

        return 0;
}
