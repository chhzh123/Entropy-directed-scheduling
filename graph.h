// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implement of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.
// Our work has been contributed to ICCD 2018.

// This head file contains the basic graph class.

#ifndef GRAPH_H
#define GRAPH_H

struct VNode
{
	int num;
	std::string name;
	std::string type;
	// the default delay of an operation is 1
	int delay = 1;
	// the number of incoming edges
	int incoming = 0;
	// temporary number of incoming edges (used for topo)
	int tempIncoming = 0;
	// for easily reading info of its predecessors and successors
	// the nodes are stored in a adjcent list
	std::vector<VNode*> pred;
	std::vector<VNode*> succ;
	VNode(int _num,std::string _name,std::string _type,int _delay = 1):
		num(_num),name(_name),type(_type),delay(_delay){};
};

class graph
{
public:
	graph() = default;
	// read from dot file
	graph(std::ifstream& infile);
	~graph();

	// output
	void printAdjlist();
	void mainScheduling();

	// EDS starting from the front
	void EDS();
	// EDS starting from the last
	void EDSrev();
	// EDS for resource-constrained scheduling problems
	void RCEDS();

	// List scheduling for resource-constrained problems
	void RCLS();
	// Force-directed scheduling for resource-constrained problems
	void FDS(){};

	// ILP formulation
	void generateRCSILP(std::ofstream& outfile);
	void generateTCSILP(std::ofstream& outfile);

	// set basic parameters
	inline void setLC(double _LC) { LC = _LC; };
	inline void setConstrainedLatency(int conlatency) { ConstrainedLatency = conlatency; };
	inline void setMODE(std::vector<int> _MODE) { MODE = _MODE; };
	inline void setMAXRESOURCE(int mul,int alu)
		{ MAXRESOURCE.first = mul; MAXRESOURCE.second = alu; };

private:
	// initialization
	void initialize();
	void clearMark();
	void setDegrees(); // in-degree or out-degree
	void addVertex(std::string name,std::string type);
	bool addEdge(std::string vFrom,std::string vTo);
	VNode* findVertex(std::string name);
	inline std::string mapResourceType(std::string type);

	// preparation
	void topologicalSortingDFS();
	void topologicalSortingKahn();
	void dfsASAP(VNode* node);
	void dfsALAP(VNode* node);

	// scheduling
	void placeCriticalPath();
	bool scheduleNodeStep(VNode* node,int step);
	bool scheduleNodeStepResource(VNode* node,int step);

	// output
	void countResource();
	
	int vertex = 0;
	int edge = 0;
	int typeNum = 0;
	int numScheduledOp = 0;

	// critical path delay
	int cdepth = 0;
	// maximum latency of scheduled operations
	int maxLatency = 0;

	// Use adjacent list to store the graph
	std::vector<VNode*> adjlist;
	// mark for DFS-based topological sorting and scheduled operations
	std::vector<int> mark;
	// topological ordering
	std::vector<VNode*> order;
	// ordering without operations on critical path
	std::vector<VNode*> edsOrder;

	// ASAP & ALAP value
	std::vector<std::pair<int,int>> axap;
	// N_r
	std::map<std::string,int> nr;
	// N_r(t)
	std::vector<std::map<std::string,int>> nrt;

	// final results
	std::map<int,int> schedule;

	// variables used for generating ILP
	std::vector<std::vector<int>> ilp;
	std::map<int,std::map<std::string,std::vector<int>>> rowResource; // step type ops

	// latency factor
	double LC = 1;
	int ConstrainedLatency;
	// for resource-constrained scheduling
	std::pair<int,int> MAXRESOURCE;

	// MODE[0]: 0 EDS            1 EDSrev        2 EDSResourceConstrained
	//			3 LS for RCS     4 ILP for TCS   5 ILP for RCS
	// MODE[1]: 0 DFS    1 Kahn
	std::vector<int> MODE;
};

// for std::string split
std::vector<std::string> split(const std::string& input, const std::string& regex) // split the string
{
	// passing -1 as the submatch index parameter performs splitting
	std::regex re(regex);
	std::sregex_token_iterator
		first{input.begin(), input.end(), re, -1},
		last;
	return {first, last};
}

#endif // GRAPH_H