#ifndef _RWGRAPH_H
#define _RWGRAPH_H
#include <vector>
#include <random>
#include <set>
#include "sfmt/SFMT.h"

typedef uint32_t UI;
typedef uint64_t ULL;

// Used for storing time information in nodes of RR sets.
// Only useful for the "base collection" that has the fastest time and the most edges
struct NodeTime
{
	unsigned int nodeID;
	double influenceTime;
};

// Graph class defines fundamental operators on graph
class Graph
{
	friend class HyperGraph;
	private:
		UI UI_MAX = 4294967295U;
		ULL ULL_MAX = 18446744073709551615ULL;

		// number of nodes
		unsigned int numNodes;
		// number of edges
		unsigned int numEdges;
		int maxDegree;
		// adjacency list
		std::vector<std::vector<int> > adjList;
		std::vector<int> node_deg;
		std::vector<std::vector<UI> > weights;
		std::vector<unsigned int> ben_ind;
		std::vector<double> ben;
	
	public:
		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u) const;
		// return weights of neighbours of node u
		const std::vector<UI> & getWeight(int u) const;

		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u);
		// return weights of neighbours of node u
		const std::vector<UI> & getWeight(int u);

		// get degree of node u
		int getDegree(int u) const;
		int getMaxDegree() const;
		// get size of the graph
		int getSize() const;
		// get number of edges
		int getEdge() const;
		// return benenfit index

		// read graph from a file
		void readGraphLT(const char * filename);
		// read graph from a file
		void readGraphIC(const char * filename);
		// write the graph to file
		void writeToFile(const char * filename);
		
};

class HyperGraph
{
	private:
		// store base hyperedges: includes the nodes reachable with full speed
		std::vector<std::vector<NodeTime> > base_edge_node;
		
		// two main groups: sped is for those samples used for triggering speedup at given time
		// cov is for coverage requirements (possibly multiple groups. in the base case, only one)
		// store the edges that a node is incident to
		std::vector<std::vector<std::vector<int> > >  sped_node_edge;	
		std::vector<std::vector<std::vector<int> > >  cov_node_edge;
		
		// store hyperedges for speed up triggering
		std::vector<std::vector<std::vector<int> > > sped_edge_node;
		std::vector<std::vector<std::vector<int> > >  cov_edge_node;
		
		std::vector<int> sped_curEdge;
		std::vector<int>  cov_curEdge;
		
		unsigned int maxDegree;
		unsigned int numNodes;
		std::vector<UI> bDist;
		sfmt_t sfmtSeed;
		inline int randIndex_bin(const std::vector<UI> &w, unsigned int si);
		inline int randIndex_lin(const std::vector<UI> &w, unsigned int si);

	public:
		HyperGraph(unsigned int n, int speed_up_times, int cov_reqs);
		void updateDeg_sped(int index);
		void updateDeg_cov(int index);
		//void updateEdge();
		//void addEdge(std::vector<int> & edge);
		//void addEdgeD(std::vector<int> & edge);
		int getMaxDegree();
		const std::vector<NodeTime> & getBaseEdge(int e);
		
		const std::vector<int> & getEdge_sped(int e, int index) const;
		const std::vector<int> & getEdge_sped(int e, int index);
		const std::vector<int> & getNode_sped(int n, int index) const;
		const std::vector<int> & getNode_sped(int n, int index);
		
		const std::vector<int> & getEdge_cov(int e, int index) const;
		const std::vector<int> & getEdge_cov(int e, int index);
		const std::vector<int> & getNode_cov(int n, int index) const;
		const std::vector<int> & getNode_cov(int n, int index);
		
		// get number of edges
        long long getNumEdge_sped(int index) const;
		long long getNumEdge_cov(int index) const;
		long long getNumEdge_base() const;
		std::vector<std::vector<NodeTime> > getBaseEdges();
		void addBaseEdges(std::vector<std::vector<NodeTime> > samples);
		// only clean sped/cov edges (when the speed up time changes)
		void clearEdges();
		// clean all edges (for the aux hgs used in multithreading)
		void clearAllEdges();
/* 		void pollingLT1(Graph &g, std::vector<bool> &visit, std::vector<int> &mark_visit);
		bool pollingLT2(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit);
		bool pollingLT(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit);
		void pollingIC1(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark);
		bool pollingIC2(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark);
		bool pollingIC(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark); */
		void pollingIC_base(Graph &g, std::vector<bool> &visit, std::vector<NodeTime> &visit_mark, double time_threshold, double delay_dec);
		void pollingIC_speedup(std::vector<int> &visit_mark, double modified_threshold, int target_speedup_index);
		void pollingIC_coverage(std::vector<int> &visit_mark, double time_threshold, int cov_index);
		double calc_modified_threshold(double time_threshold, std::vector<double> speedup_time_ints, std::vector<double> delay_decs);
		// Generator of the weibull distribution
		std::default_random_engine generator;
		std::weibull_distribution<double> distribution;
		// utility function for sorting NodeTime by time
		static bool TimeComp (NodeTime nt1, NodeTime nt2) {return (nt1.influenceTime < nt2.influenceTime);}
		
		bool verifySeedSet(Graph &g, std::set<int> seedSet, double speed_rate, double speed_req, double cov_req, int &num_inf_nodes, double &sped_time);
};



float getCurrentMemoryUsage();

#endif
