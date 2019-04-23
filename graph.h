// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.
// Our work has been contributed to ICCD 2018 and TCAD.

// This head file contains the basic graph class.

#ifndef GRAPH_H
#define GRAPH_H

#include<vector>
#include<map>
#include<algorithm>

#define MAXINT_ 0x3f3f3f3f

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
	// ASAP & ALAP value
	int asap = 1;
	int alap = MAXINT_;
	// length of original time frame
	int length = 0;
	// final results
	int cstep = 0;
	bool criticalPath = false;
	VNode(int _num,std::string _name,std::string _type,int _delay = 1):
		num(_num),name(_name),type(_type),delay(_delay){};
	inline void setASAP(int _asap)
	{
		asap = std::max(asap,_asap);
	}
	inline void setALAP(int _alap)
	{
		alap = std::min(alap,_alap);
	}
	inline void extendALAP(int _alap)
	{
		alap = _alap;
	}
	inline void setLength() // this function should be called as soon as topo sorting has been implemented
	{
		length = alap - asap + 1;
	}
	inline int getLength()
	{
		return alap - asap + 1;
	}
	void schedule(int step)
	{
		cstep = step;
		for (auto pnode = succ.cbegin(); pnode != succ.cend(); ++pnode) // updown behavior
			(*pnode)->setASAP(step + delay);
	}
	void scheduleBackward(int step)
	{
		cstep = step;
		for (auto pnode = pred.cbegin(); pnode != pred.cend(); ++pnode)
			(*pnode)->setALAP(step - (*pnode)->delay);
	}
	void scheduleAll(int step) // mainly for FDS
	{
		cstep = step;
		iterativeSetASAP(step);
		iterativeSetALAP(step);
		setLength();
	}
	void iterativeSetASAP(int step)
	{
		for (auto pnode = succ.cbegin(); pnode != succ.cend(); ++pnode)
			(*pnode)->iterativeSetASAP(step + delay);
		setASAP(step);
	}
	void iterativeSetALAP(int step)
	{
		for (auto pnode = pred.cbegin(); pnode != pred.cend(); ++pnode)
			(*pnode)->iterativeSetALAP(step - (*pnode)->delay);
		setALAP(step);
	}
	bool testValid(int step)
	{
		for (auto pnode = succ.cbegin(); pnode != succ.cend(); ++pnode)
			if (step + delay >= (*pnode)->cstep)
				return false;
		for (auto pnode = pred.cbegin(); pnode != pred.cend(); ++pnode)
			if ((*pnode)->cstep + (*pnode)->delay >= step)
				return false;
		return true;
	}
};

class graph
{
public:
	graph() = default;
	~graph();

	// read from dot file
	void readFile(std::ifstream& infile);

	// output
	void printAdjlist() const;
	void printTimeFrame() const;
	void mainScheduling(int mode = 0);

	// EDS starting from the front
	void TC_EDS (int order_mode = 0);
	void TC_IEDS(int order_mode = 0);
	// EDS for resource-constrained scheduling problems
	void RC_EDS ();
	void RC_IEDS();

	// Force-directed scheduling for time-constrained problems
	void TC_FDS();
	// Force-directed scheduling for resource-constrained problems
	void RC_FDS();

	// List scheduling for time-constrained problems
	void TC_LS();
	void RC_LS();

	// test
	bool testFeasibleSchedule() const;

	// ILP formulation
	void generateRC_ILP(std::ofstream& outfile);
	void generateTC_ILP(std::ofstream& outfile);

	// set basic parameters
	inline void setLC(double _LC) { LC = _LC; };
	inline void setConstrainedLatency(int conlatency) { ConstrainedLatency = conlatency; };
	inline void setMODE(std::vector<int> _MODE) { MODE = _MODE; };
	inline void setMAXRESOURCE(const std::map<std::string,int> gr)
		{ MAXRESOURCE = gr; };
	inline void setPRINT(int mode) { if (mode == 0) PRINT = false; };
	inline double getLC() const {return LC;};
	inline int getMaxLatency() const {return maxLatency;};

private:
	// initialization
	void initialize();
	void clearMark();
	void setDegrees(); // in-degree or out-degree
	void addVertex(const std::string name,const std::string type);
	bool addEdge(const std::string vFrom,const std::string vTo);
	VNode* findVertex(const std::string name) const;
	inline std::string mapResourceType(const std::string type) const;

	// preparation
	void topologicalSortingDFS(bool aslap_order = false);
	void topologicalSortingKahn();
	void dfsASAP(VNode* const& node);
	void dfsALAP(VNode* const& node);

	// scheduling
	void placeCriticalPath();
	bool scheduleNodeStep(VNode* const& node,int step,int mode);
	bool newScheduleNodeStep(VNode* const& node,int step);
	bool scheduleNodeStepResource(VNode* const& node,int step,int mode);
	double calForce(int a,int b,int na,int nb,const std::vector<double>& DG,int delay) const;
	double calPredForce(VNode* const& v,int cstep,const std::map<std::string,std::vector<double>>& DG) const;
	double calSuccForce(VNode* const& v,int cstep,const std::map<std::string,std::vector<double>>& DG) const;

	// output
	void standardOutput() const;
	void simplifiedOutput() const;
	void printGanttGraph() const;
	void countResource() const;
	void countEachStepResource() const;
	void countTF();
	inline void print(const std::string str) const;
	
	int vertex = 0;
	int edge = 0;
	int typeNum = 0;
	int numScheduledOp = 0;

	// mul_delay
	int MUL_DELAY = 2;

	// critical path delay
	int cdepth = 0;
	// maximum latency of scheduled operations
	int maxLatency = 0;

	// Use adjacent list to store the graph
	std::vector<VNode*> adjlist;
	// mark for DFS-based topological sorting and list scheduling
	std::vector<int> mark;
	// topological ordering
	std::vector<VNode*> order;
	// ordering without operations on critical path
	std::vector<VNode*> edsOrder;

	// N_r
	std::map<std::string,int> nr;
	std::map<std::string,int> r_delay;
	std::map<std::string,std::vector<int>> TFcount;
	// N_r(t)
	std::vector<std::map<std::string,int>> nrt;
	// max N_r(t)
	std::map<std::string,int> maxNrt;

	// variables used for generating ILP
	std::vector<std::vector<int>> ilp;
	std::map<int,std::map<std::string,std::vector<int>>> rowResource; // step type ops

	// latency factor
	double LC = 1;
	int ConstrainedLatency;
	// for resource-constrained scheduling
	std::map<std::string,int> MAXRESOURCE;

	std::vector<int> MODE;
	bool PRINT = true;
};

#endif // GRAPH_H