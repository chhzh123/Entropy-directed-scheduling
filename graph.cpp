// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.
// Our work has been contributed to ICCD 2018 and TCAD.

// This head file contains the implementation of the basic graph operations and our EDS algorithm.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <utility> // pairs
#include <map>
#include <string>
#include <algorithm>
#include <regex> // regular expression for string split
#include <iterator>
#include <cmath>

#include "watch.h" // for high-accuracy time counting
#include "graph.h"
using namespace std;

#define MUL_DELAY 2 // default MUL delay

// for high-accuracy time counting
stop_watch watch;
// for string split
std::vector<std::string> split(const std::string& input, const std::string& regex);

graph::~graph()
{
	for (auto node : adjlist)
		delete node;
}

void graph::clearMark()
{
	mark.clear();
	for (int i = 0; i < vertex; ++i)
		mark.push_back(0);
}

void graph::initialize()
{
	print("Begin initializing...");
	clearMark();
	print("Initialized successfully!\n");
}

// read from dot file
void graph::readFile(ifstream& infile)
{
	string str;
	// The first two lines in dot file are useless info
	getline(infile,str);
	getline(infile,str);
	print("Begin parsing...");
	// operation nodes info
	while (getline(infile,str) && str.find("-") == std::string::npos)
	{
		vector<string> op = split(str," *\\[ *label *= *| *\\];| +"); // reg exp
		addVertex(op[1],op[2]); // op[0] = ""
	}
	// edges info
	do {
		vector<string> arc = split(str," *\\[ *name *= *| *\\];| *-> *| +"); // reg exp
		if (!addEdge(arc[1],arc[2]))
			cout << "Add edge wrong!" << endl;
	} while (getline(infile,str) && str.size() > 1);
	print("Parsed dot file successfully!\n");
	initialize();
}

void graph::addVertex(const string name,const string type)
{
	int delay = 1;
	// set MUL delay
	if (mapResourceType(type) == "MUL")
		delay = MUL_DELAY;
	// be careful of the numbers!!! start labeling from 0
	VNode* v = new VNode(vertex++,name,type,delay);
	adjlist.push_back(v);
	if (nr.find(mapResourceType(type)) != nr.end())
		nr[mapResourceType(type)]++;
	else
	{
		typeNum++;
		nr[mapResourceType(type)] = 1;
		maxNrt[mapResourceType(type)] = 0;
	}
}

bool graph::addEdge(const string vFrom,const string vTo)
{
	VNode* vf = findVertex(vFrom);
	VNode* vt = findVertex(vTo);
	if (vf == nullptr || vt == nullptr)
		return false;
	if ((MODE.size() == 2 && MODE[1] == 0) || (MODE.size() > 2 && MODE[2] == 1)) // top-down behavior
	{
		vf->succ.push_back(vt);
		vt->pred.push_back(vf);
	}
	else // bottom-up behavior
	{
		vt->succ.push_back(vf);
		vf->pred.push_back(vt);
	}
	edge++;
	vector<int> cons = {vf->num,vt->num,(-1)*vf->delay};
	ilp.push_back(cons);
	return true;
}

VNode* graph::findVertex(const string name) const
{
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		if ((*pnode)->name == name)
			return (*pnode);
	return nullptr;
}

void graph::printAdjlist() const
{
	cout << "Adjacent list:" << endl;
	cout << "[ Format: node num ( node name ) : successor num ( successor name ) ]" << endl;
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
	{
		cout << (*pnode)->num+1 << "( " << (*pnode)->name << " ): ";
		for (auto adjnode = (*pnode)->succ.cbegin(); adjnode != (*pnode)->succ.cend(); ++adjnode)
			cout << (*adjnode)->num+1 << "( " << (*adjnode)->name << " ) ";
		cout << endl;
	}
}

string graph::mapResourceType(const string type) const
{
	// return "ALL";
	if (type == "mul" || type == "MUL" || type == "div" || type == "DIV")
		return "MUL";
	if (type == "sub" || type == "add" || type == "SUB" || type == "ADD" ||
		type == "NEG" || type == "AND" || type == "les" || type == "LSR" || type == "ASR" ||
		type == "imp" || type == "exp" || type == "MemR" || type == "MemW" ||
		type == "STR" || type == "LOD" || type == "BNE" || type == "BGE" || type == "LSL")
		return "ALU";
	return type;
}

void graph::setDegrees()
{
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		(*pnode)->tempIncoming = (*pnode)->incoming = (*pnode)->pred.size();
}

void graph::dfsASAP(VNode* const& node)
{
	if (mark[node->num])
		return;
	if (!(node->pred.empty()))
		for (auto pprec = node->pred.cbegin(); pprec != node->pred.cend(); ++pprec)
		{
			dfsASAP(*pprec);
			node->setASAP((*pprec)->asap + (*pprec)->delay);
		}
	cdepth = max(node->asap + node->delay - 1,cdepth); // critical path delay
	if (MODE[0] < 10 || MODE[0] == 13)
		setConstrainedLatency(int(cdepth*LC));
	else
		setConstrainedLatency(MAXINT_);
	mark[node->num] = 1;
	order.push_back(node);
}

void graph::dfsALAP(VNode* const& node) // different from asap
{
	if (mark[node->num])
		return;
	if (node->succ.empty())
		node->setALAP(ConstrainedLatency - node->delay + 1); // ConstrainedLatency is used here, dfsasap must be done first
	else for (auto psucc = node->succ.cbegin(); psucc != node->succ.cend(); ++psucc)
	{
		dfsALAP(*psucc);
		node->setALAP((*psucc)->alap - node->delay);
	}
	node->setLength();
	mark[node->num] = 1;
}

void graph::topologicalSortingDFS()
{
	setDegrees();
	print("Begin topological sorting...");
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode) // asap
		if ((*pnode)->succ.empty() && !mark[(*pnode)->num]) // out-degree = 0
			dfsASAP(*pnode);
	clearMark();
	// countASAP();
	// countASAP_RCS();
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode) // alap
		if ((*pnode)->pred.empty() && !mark[(*pnode)->num]) // in-degree = 0
			dfsALAP(*pnode);
	// regenerate order
	// sort(order.begin(),order.end(),[](VNode* const& node1,VNode* const& node2)
	// 	{
	// 		if (node1->alap < node2->alap)
	// 			return true;
	// 		else if (node1->alap == node2->alap)
	// 				if (node1->asap < node2->asap)
	// 					return true;
	// 				else
	// 					return false;
	// 			else
	// 				return false;
	// 	});
	print("Topological sorting done!");
	// printTimeFrame();
}

void graph::topologicalSortingKahn()
{
	print("Begin topological sorting (Kahn)...");
	// -------- DFS part --------
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode) // asap
		if ((*pnode)->succ.empty() && !mark[(*pnode)->num]) // out-degree = 0
			dfsASAP(*pnode);
	clearMark();
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode) // alap
		if ((*pnode)->pred.empty() && !mark[(*pnode)->num]) // in-degree = 0
			dfsALAP(*pnode);
	// --------------------------
	order.clear();
	vector<VNode*> temp;
	setDegrees();
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		if ((*pnode)->pred.empty()) // in-degree = 0
			temp.push_back(*pnode);
	while (!temp.empty())
	{
		order.push_back(temp[0]);
		for (auto pnode = temp[0]->succ.cbegin(); pnode != temp[0]->succ.cend(); ++pnode)
		{
			(*pnode)->tempIncoming--;
			if ((*pnode)->tempIncoming == 0)
				temp.push_back(*pnode);
		}
		temp.erase(temp.begin());
	}
	print("Topological sorting (Kahn) done!");
	clearMark();
}

bool graph::newScheduleNodeStep(VNode* const& node,int step)
{
	cout << node->num << " " << step << endl;
	if (step + node->delay - 1 > ConstrainedLatency) // important to minus 1
	{
		cout << "Invalid schedule!" << endl;
		return false;
	}
	auto Rtype = mapResourceType(node->type);
	for (int i = node->cstep; i < node->cstep + node->delay; ++i)
		nrt[i][Rtype]--;
	for (int i = step; i < step + node->delay; ++i)
		nrt[i][Rtype]++;
	maxNrt[Rtype] = 0;
	for (auto x: nrt)
		maxNrt[Rtype] = max(maxNrt[Rtype],x[Rtype]);
	node->schedule(step);
	maxLatency = max(maxLatency,step + node->delay - 1);
	return true;
}

bool graph::scheduleNodeStep(VNode* const& node,int step,int mode = 0)
{
	if (MODE[0] == 3)
	{
		// cout << node->num+1 << " " << step << endl;
		for (auto pnode = node->succ.cbegin(); pnode != node->succ.cend(); ++pnode)
			if ((*pnode)->cstep == 0 && max((*pnode)->asap,step + node->delay) > (*pnode)->alap)
				return false;
		for (auto pnode = node->pred.cbegin(); pnode != node->pred.cend(); ++pnode)
			if ((*pnode)->cstep == 0 && min((*pnode)->alap,step - (*pnode)->delay) < (*pnode)->asap)
				return false;
	}
	if (step + node->delay - 1 > ConstrainedLatency) // important to minus 1
	{
		cout << "Invalid schedule!" << endl;
		return false;
	}
	for (int i = step; i < step + node->delay; ++i)
	{
		nrt[i][mapResourceType(node->type)]++;
		maxNrt[mapResourceType(node->type)] = max(maxNrt[mapResourceType(node->type)],nrt[i][mapResourceType(node->type)]);
	}
	switch (mode)
	{
		case 0: node->schedule(step);break;
		case 1: node->scheduleBackward(step);break;
		case 2: node->scheduleAll(step);break;
		default: cout << "Invaild schedule mode!" << endl;return false;
	}
	maxLatency = max(maxLatency,step + node->delay - 1);
	numScheduledOp++;
	return true;
}

bool graph::scheduleNodeStepResource(VNode* const& node,int step,int mode = 0)
{
	for (int i = step; i < step + node->delay; ++i)
	{
		nrt[i][mapResourceType(node->type)]++;
		maxNrt[mapResourceType(node->type)] = max(maxNrt[mapResourceType(node->type)],nrt[i][mapResourceType(node->type)]);
	}
	switch (mode)
	{
		case 0: node->schedule(step);break;
		case 1: node->scheduleBackward(step);break;
		case 2: node->scheduleAll(step);break;
		default: cout << "Invaild schedule mode!" << endl;return false;
	}
	maxLatency = max(maxLatency,step + node->delay - 1); // important to minus 1
	numScheduledOp++;
	return true;
}

void graph::placeCriticalPath()
{
	print("Begin placing critical path...");
	// int minL = MAXINT_;
	// for (auto node : order)
	// 	minL = min(minL,node->alap-node->asap);
	for (auto node : order)
	{
		if (node->asap == node->alap)
		{
			node->criticalPath = true;
			scheduleNodeStep(node,node->asap);
		}
		else
			edsOrder.push_back(node);
		// if (node->alap - node->asap == minL)
		// 	node->criticalPath = true;
		// if (node->asap == node->alap)
		// 	scheduleNodeStep(node,node->asap,2);
		// else
		// 	edsOrder.push_back(node);
	}
	print("Placing critical path done!");
}

void graph::TC_EDS(int sorting_mode)
{
	print("Begin EDS scheduling...\n");
	watch.restart();
	if (sorting_mode == 0)
		topologicalSortingDFS();
	else
		topologicalSortingKahn();
	// initialize N_r(t)
	map<string,int> temp;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		temp[pnr->first] = 0;
	for (int i = 0; i <= ConstrainedLatency; ++i) // number+1
		nrt.push_back(temp);

	// placing operations on critical path
	placeCriticalPath();
	print("Critical path time delay: " + to_string(cdepth));
	// cout << "Constrained latency: " << ConstrainedLatency << endl;

	// main part of scheduling
	print("Begin placing other nodes...");
	for (auto pnode = edsOrder.cbegin(); pnode != edsOrder.cend(); ++pnode)
	{
		int a = (*pnode)->asap, b = (*pnode)->alap;
		// because of topo order, it's pred must have been scheduled
		double minnrt = MAXINT_;
		int minstep = a, maxnrt = -MAXINT_, maxstep = a, flag = 0;
		for (int t = a; t <= b; ++t)
		{
			double sumNrt = 0;
			for (int d = 1; d <= (*pnode)->delay; ++d)
			{
				string tempType = mapResourceType((*pnode)->type);
				sumNrt += nrt[t+d-1][tempType];
			}
			if (sumNrt < minnrt) // leave freedom to remained ops
			{
				minnrt = sumNrt;
				minstep = t;
			}
		}
		// cout << (*pnode)->num+1 << " (" << (*pnode)->name << "): " << a << " " << b << " Step: " << minstep << endl;
		scheduleNodeStep(*pnode,minstep);
	}
	watch.stop();
	print("Placing other nodes done!\n");
	print("Test small steps.\n");
	int cnt = 0;
	while (cnt != vertex)
	{
		cnt = 0;
		for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); pnode++)
		{
			// cout << (*pnode)-> num << " " << (*pnode)->cstep << " " << cnt << endl;
			string tempType = mapResourceType((*pnode)->type);
			int t = (*pnode)->cstep, cntin = 0;
			for (int d = 1; d <= (*pnode)->delay; ++d)
				if (t+d-2>0 && nrt[t+d-1][tempType] > nrt[t+d-2][tempType] + 1)
					cntin++;
			if (cntin == (*pnode)->delay)
				if (t - 1 > 0 && (*pnode)->testValid(t-1)){
					newScheduleNodeStep(*pnode,t-1);
					continue;	
				}
			cntin = 0;
			for (int d = 1; d <= (*pnode)->delay; ++d)
				if (t+d <= ConstrainedLatency && nrt[t+d-1][tempType] > nrt[t+d][tempType] + 1)
					cntin++;
			if (cntin == (*pnode)->delay)
				if (t + (*pnode)->delay - 1 <= ConstrainedLatency && (*pnode)->testValid(t+1)){
					newScheduleNodeStep(*pnode,t+1);
					continue;
				}
			cnt++;
		}
	}
	print("Finish EDS scheduling!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::RC_EDS(int sorting_mode) // Resource-constrained EDS
{
	print("Begin resource-constrained entropy-directed scheduling (EDS)...\n");
	watch.restart();
	if (sorting_mode == 0)
		topologicalSortingDFS();
	else
		topologicalSortingKahn();
	// initialize N_r(t)
	map<string,int> temp,maxNr;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		if (pnr->first == "mul" || pnr->first == "MUL")
		{
			temp[pnr->first] = 0;
			maxNr[pnr->first] = MAXRESOURCE.first;
		}
		else
		{
			temp[pnr->first] = 0;
			maxNr[pnr->first] = MAXRESOURCE.second;
		}
	nrt.push_back(temp); // nrt[0]

	// NO nrt.push_back! NO placeCriticalPath!
	print("Begin placing operations...");
	for (auto pnode = order.cbegin(); pnode != order.cend(); ++pnode)
	{
		int a = (*pnode)->asap, b = (*pnode)->alap;
		// because of topo order, it's pred must have been scheduled
		int maxstep = a, maxnrt = -1;
		// cout << (*pnode)->name << " " << a << " " << b << endl;
		for (int t = a; t <= max(a,maxLatency) + MUL_DELAY; ++t)
		{
			int flag = 1, sumNrt = 0;
			for (int d = 1; d <= (*pnode)->delay; ++d)
			{
				string tempType = mapResourceType((*pnode)->type);
				if (t+d-1 >= nrt.size())
					nrt.push_back(temp); // important!
				if (nrt[t+d-1][tempType]+1 > maxNr[tempType])
					flag = 0;
				sumNrt += nrt[t+d-1][tempType];
			}
			// if (flag == 1)
			// {
			// 	maxstep = t;
			// 	break;
			// }
			if (flag == 1 && sumNrt > maxnrt) // leave freedom to remained ops
			{
				maxnrt = sumNrt;
				maxstep = t;
				break;
			}
		}
		scheduleNodeStepResource(*pnode,maxstep); // some differences
	}
	watch.stop();
	print("Placing operations done!\n");
	print("Finish entropy-directed scheduling!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::countResource() const
{
	ofstream out("./Resource_"+to_string(LC)+".out",ios::app);
	for (auto ptype = nr.crbegin(); ptype != nr.crend(); ++ptype)
	{
		cout << ptype->first << ": " << maxNrt.at(ptype->first) << endl;
		out << maxNrt.at(ptype->first) << " ";
		if (PRINT)
		// if (ptype->first == "MUL" || ptype->first == "mul")
			for (int i = 1; i <= maxLatency; ++i) // ConstrainedLatency
				cout << "Step " << i << ": " << nrt[i].at(mapResourceType(ptype->first)) << endl;
	}
	out << endl;
}

void graph::mainScheduling(int mode)
{
	switch (MODE[0])
	{
		case 0: TC_EDS(0);break;
		case 1: TC_EDS(1);break;
		case 3: TC_FDS();break;
		case 4: TC_LS();break;
		case 10: RC_EDS(0);break;
		case 11: RC_EDS(1);break;
		case 13: RC_FDS();break;
		case 14: RC_LS();break;
		default: cout << "Invaild mode!" << endl;return;
	}
	if (mode == 0)
		standardOutput();
	else
		simplifiedOutput();
}

void graph::standardOutput() const
{
	if (!testFeasibleSchedule())
	{
		cout << "\nInfeasible schedule!" << endl;
		// return;
	}
	else
		cout << "\nThe schedule is valid!" << endl;
	cout << "Output as follows:" << endl;
	printAdjlist();
	cout << "Topological order:" << endl;
	for (auto pnode = order.cbegin(); pnode != order.cend(); ++pnode)
		cout << (*pnode)->num+1 << ":" << (*pnode)->name << ((pnode-order.cbegin()+1)%5==0 ? "\n" : "   \t");
	cout << endl;
	if (MODE[0] < 10)
		printTimeFrame();
	cout << "Final schedule:" << endl;
	for (int i = 0; i < vertex; ++i)
		cout << i+1 << ": " << adjlist[i]->cstep << ((i+1)%5==0 ? "\n" : "\t");
	cout << endl;
	cout << "Gantt graph:" << endl;
	cout << "    ";
	for (int i = 1; i <= maxLatency; ++ i)
		cout << i % 10;
	cout << endl;
	for (int i = 0; i < vertex; ++i)
	{
		cout << setw(4) << std::left << i+1;
		for (int j = 1; j < adjlist[i]->cstep; ++j)
			cout << " ";
		for (int j = 1; j <= adjlist[i]->delay; ++j)
			cout << (adjlist[i]->delay > 1 ? "X" : "O");
		cout << endl;
	}
	cout << "Total latency: " << maxLatency << endl;
	if (MODE[0] >= 10)
		cout << "Constrained resource:\n"
				"MUL: " << MAXRESOURCE.first << endl <<
				"ALU: " << MAXRESOURCE.second << endl;
	cout << "Resource used:" << endl;
	countResource();
}

void graph::simplifiedOutput() const
{
	if (!testFeasibleSchedule())
	{
		cout << "\nInfeasible schedule!" << endl;
		return;
	}
	cout << "Total latency: " << maxLatency << endl;
	if (MODE[0] < 10)
	{
		cout << "Resource used:" << endl;
		countResource();
	}
	cout << endl;
}

bool graph::testFeasibleSchedule() const
{
	int flag = 0;
	for (int i = 0; i < vertex; ++i)
		for (auto pnode = adjlist[i]->succ.cbegin(); pnode != adjlist[i]->succ.cend(); ++pnode)
			if (adjlist[i]->cstep + adjlist[i]->delay - 1 >= (*pnode)->cstep)
			{
				flag = 1;
				cout << "Schedule conflicts with Node " << adjlist[i]->num+1 << " (" << adjlist[i]->name << ") "
					 << "and Node " << (*pnode)->num+1 << " (" << (*pnode)->name << ")." << endl;
			}
	if (flag == 1)
		return false;
	else
		return true;
}

// generated ILP in CPLEX form
void graph::generateTC_ILP(ofstream& outfile)
{
	topologicalSortingDFS();
	cout << "Time frame:" << endl;
	int cnt = 1;
	for (auto pnode : adjlist)
	 	cout << pnode->num+1 << ": [ " << pnode->asap << " , " << pnode->alap << " ]" << endl;
	cout << endl;
	cout << "Start generating ILP formulas for latency-constrained problems..." << endl;

	outfile << "Minimize" << endl;
	outfile << "M1 + M2" << endl;

	outfile << "Subject To" << endl;
	cnt = 0;
	// Time frame constraints
	for (auto pnode : adjlist)
	{
		for (int i = pnode->asap; i <= pnode->alap; ++i)
			outfile << "x" << cnt << "," << i << (i == pnode->alap ? " = 1\n" : " + ");
		cnt++;
	}
	cout << "Time frame constraints generated." << endl;

	// Resource constraints
	cnt = 0;
	for (cnt = 0; cnt < vertex; ++cnt)
		for (int i = adjlist[cnt]->asap; i <= adjlist[cnt]->alap + adjlist[cnt]->delay - 1; ++i)
			// cout << i << " " << adjlist[cnt]->type << endl;
			rowResource[i][mapResourceType(adjlist[cnt]->type)].push_back(cnt); // push delay
	cout << "Critical path delay: " << ConstrainedLatency << endl;
	for (int i = 1; i <= ConstrainedLatency; ++i)
		for (auto ptype = nr.cbegin(); ptype != nr.cend(); ++ptype)
		{
			if (rowResource[i][ptype->first].size() < 2)
				continue;
			for (int j = 0; j < rowResource[i][ptype->first].size(); ++j)
				for (int d = 0; d < adjlist[rowResource[i][ptype->first][j]]->delay; ++d)
					if (i-d >= 1)
						outfile << "x" << rowResource[i][ptype->first][j]
								<< "," << i-d << ((j == rowResource[i][ptype->first].size()-1 && (d == adjlist[rowResource[i][ptype->first][j]]->delay-1 || i-d == 1)) ? "" : " + ");
					else
						break;
			if (ptype->first == "MUL") // ptype->first == "mul" || 
				outfile << " - M1 <= 0" << endl;
			else
				outfile << " - M2 <= 0" << endl;
		}
	cout << "Resource constraints generated." << endl;

	// Precedence constraints
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		for (auto psucc = (*pnode)->succ.cbegin(); psucc != (*pnode)->succ.cend(); ++psucc)
		{
			for (int i = (*pnode)->asap; i <= (*pnode)->alap; ++i)
				outfile << i << " x" << (*pnode)->num << ","
						<< i << (i == (*pnode)->alap ? "" : " + ");
			outfile << " - ";
			for (int i = (*psucc)->asap; i <= (*psucc)->alap; ++i)
				outfile << i << " x" << (*psucc)->num << ","
						<< i << (i == (*psucc)->alap ? "" : " - ");
			outfile << " <= -" << (*pnode)->delay << endl;
		}
	cout << "Precedence constraints generated." << endl;

	// Bounds NO VARIABLES RHS!
	outfile << "Bounds" << endl;
	for (int i = 0; i < vertex; ++i)
		for (int j = adjlist[i]->asap; j <= adjlist[i]->alap; ++j)
			// outfile << "x" << i << "," << j << " >= 0" <<endl;
			outfile << "0 <= x" << i << "," << j << " <= 1" <<endl;
	outfile << "M1 >= 1" << endl;
	outfile << "M2 >= 1" << endl;
	cout << "Bounds generated." << endl;

	// Generals
	outfile << "Generals" << endl;
	for (int i = 0; i < vertex; ++i)
		for (int j = adjlist[i]->asap; j <= adjlist[i]->alap; ++j)
			outfile << "x" << i << "," << j << "\n";
	outfile << "M1\nM2" << endl;
	cout << "Generals generated." << endl;
	// *******SDC*******
	// cnt = 1;
	// for (auto cons : ilp)
	// 	outfile << "c" << cnt++ << ": x"
	// 			<< cons[0] << " - x"
	// 			<< cons[1] << " <= " << cons[2] << endl;
	// outfile << "Bounds" << endl;
	// for (int i = 0; i < vertex; ++i)
	// 	outfile << "x" << i << " >= 0" << endl;
	// outfile << "Generals" << endl;
	// for (int i = 0; i < vertex; ++i)
	// 	outfile << "x" << i << " ";
	// outfile << endl;
	outfile << "End" << endl;
	cout << "Finished ILP generation!" << endl;
}

// generated in CPLEX form
void graph::generateRC_ILP(ofstream& outfile)
{
	topologicalSortingDFS();
	cout << "Time frame:" << endl;
	int cnt = 1;
	for (auto pnode = adjlist.begin(); pnode != adjlist.end(); ++pnode)
	{
		(*pnode)->setALAP(vertex); // set upper bound
		cout << cnt++ << ": [ " << (*pnode)->asap << " , " << (*pnode)->alap << " ]" << endl;
	}
	cout << endl;
	cout << "Start generating ILP formulas for resource-constrained problems..." << endl;

	outfile << "Minimize" << endl;
	outfile << "L" << endl;

	outfile << "Subject To" << endl;
	cnt = 0;
	// Time frame constraints
	for (auto pnode : adjlist)
	{
		for (int i = pnode->asap; i <= pnode->alap; ++i)
			outfile << "x" << cnt << "," << i << (i == pnode->alap ? " = 1\n" : " + ");
		// (t+d-1) x
		for (int i = pnode->asap; i <= pnode->alap; ++i)
			outfile << (i + pnode->delay - 1) << " x" << cnt << "," << i << " - L <= 0" << endl;
		cnt++;
	}
	cout << "Time frame and upper latency constraints generated." << endl;

	// Resource constraints
	cnt = 0;
	for (cnt = 0; cnt < vertex; ++cnt)
		for (int i = adjlist[cnt]->asap; i <= adjlist[cnt]->alap + adjlist[cnt]->delay - 1; ++i)
			// cout << i << " " << adjlist[cnt]->type << endl;
			rowResource[i][mapResourceType(adjlist[cnt]->type)].push_back(cnt); // push delay
	// cout << "Critical path delay: " << ConstrainedLatency << endl;
	for (int i = 1; i <= vertex; ++i) // ConstrainedLatency
		for (auto ptype = nr.cbegin(); ptype != nr.cend(); ++ptype)
		{
			if (rowResource[i][ptype->first].size() < 2)
				continue;
			for (int j = 0; j < rowResource[i][ptype->first].size(); ++j)
				for (int d = 0; d < adjlist[rowResource[i][ptype->first][j]]->delay; ++d)
					if (i-d >= 1)
						outfile << "x" << rowResource[i][ptype->first][j]
								<< "," << i-d << ((j == rowResource[i][ptype->first].size()-1 && (d == adjlist[rowResource[i][ptype->first][j]]->delay-1 || i-d == 1)) ? "" : " + ");
					else
						break;
			if (ptype->first == "MUL") // ptype->first == "mul" || 
				outfile << " <= " << MAXRESOURCE.first << endl; // differenet
			else
				outfile << " <= " << MAXRESOURCE.second << endl;
		}
	cout << "Resource constraints generated." << endl;

	// Precedence constraints
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		for (auto psucc = (*pnode)->succ.cbegin(); psucc != (*pnode)->succ.cend(); ++psucc)
		{
			for (int i = (*pnode)->asap; i <= (*pnode)->alap; ++i)
				outfile << i << " x" << (*pnode)->num << ","
						<< i << (i == (*pnode)->alap ? "" : " + ");
			outfile << " - ";
			for (int i = (*psucc)->asap; i <= (*psucc)->alap; ++i)
				outfile << i << " x" << (*psucc)->num << ","
						<< i << (i == (*psucc)->alap ? "" : " - ");
			outfile << " <= -" << (*pnode)->delay << endl;
		}
	cout << "Precedence constraints generated." << endl;

	// Bounds NO VARIABLES RHS!
	outfile << "Bounds" << endl;
	for (int i = 0; i < vertex; ++i)
		for (int j = adjlist[i]->asap; j <= adjlist[i]->alap; ++j)
			// outfile << "x" << i << "," << j << " >= 0" <<endl;
			outfile << "0 <= x" << i << "," << j << " <= 1" <<endl;
	outfile << "L >= 1" << endl;
	cout << "Bounds generated." << endl;

	// Generals
	outfile << "Generals" << endl;
	for (int i = 0; i < vertex; ++i)
		for (int j = adjlist[i]->asap; j <= adjlist[i]->alap; ++j)
			outfile << "x" << i << "," << j << "\n";
	outfile << "L" << endl;
	cout << "Generals generated." << endl;
	outfile << "End" << endl;
	cout << "Finished ILP generation!" << endl;
}

void graph::print(const string str) const
{
	if (PRINT)
		cout << str << endl;
}

std::vector<std::string> split(const std::string& input, const std::string& regex) // split the string
{
	// passing -1 as the submatch index parameter performs splitting
	std::regex re(regex);
	std::sregex_token_iterator
		first{ input.begin(), input.end(), re, -1 },
		last;
	return { first, last };
}

void graph::countASAP() const
{
	vector<int> countasap(cdepth+1,0);
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode) // asap
		countasap[(*pnode)->asap] += 1;
	// cout << "(";
	// for (int i = 1; i < countasap.size(); ++i)
	// 	cout << countasap[i] << ", ";
	// cout << ")" << endl;
	vector<int> sumv(cdepth+1,0);
	for (int i = 1; i < countasap.size(); ++i)
		sumv[i] = sumv[i-1] + countasap[i];
	for (int i = 1; i < countasap.size(); ++i)
		if (sumv[i] > (float)(sumv[countasap.size()-1])/2.0)
		{
			if (i < (float)(countasap.size()-1)/2.0)
				cout << (float)i/(countasap.size()-1) << " bottom-up" << endl;
			else
				cout << (float)i/(countasap.size()-1) << " top-down" << endl;
			break;
		}
}

void graph::countASAP_RCS() const
{
	vector<int> countMULasap(cdepth+1,0);
	vector<int> countALUasap(cdepth+1,0);
	for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode) // asap
		if (mapResourceType((*pnode)->type) == "MUL")
			countMULasap[(*pnode)->asap] += 1;
		else
			countALUasap[(*pnode)->asap] += 1;
	cout << "MUL: (";
	for (int i = 1; i < countMULasap.size(); ++i)
		cout << countMULasap[i] << ", ";
	cout << ")" << endl;
	cout << "ALU: (";
	for (int i = 1; i < countALUasap.size(); ++i)
		cout << countALUasap[i] << ", ";
	cout << ")" << endl;
	vector<int> sumMULv(cdepth+1,0);
	vector<int> sumALUv(cdepth+1,0);
	for (int i = 1; i < countMULasap.size(); ++i)
	{
		sumMULv[i] = sumMULv[i-1] + countMULasap[i];
		sumALUv[i] = sumALUv[i-1] + countALUasap[i];
	}
	cout << "Total MUL:" << sumMULv[countMULasap.size()-1] << endl;
	cout << "Total ALU:" << sumALUv[countMULasap.size()-1] << endl;
	cout << "Total nodes: " << sumMULv[countMULasap.size()-1] + sumALUv[countMULasap.size()-1] << endl;
	for (int i = 1; i < countMULasap.size(); ++i)
		if (sumMULv[i] > (float)(sumMULv[countMULasap.size()-1])/2.0)
		{
			if (i >= (float)(countMULasap.size()-1)/2.0)
				cout << (float)i/(countMULasap.size()-1) << " MUL bottom-up" << endl;
			else
				cout << (float)i/(countMULasap.size()-1) << " MUL top-down" << endl;
			break;
		}
	for (int i = 1; i < countALUasap.size(); ++i)
		if (sumALUv[i] > (float)(sumALUv[countALUasap.size()-1])/2.0)
		{
			if (i >= (float)(countALUasap.size()-1)/2.0)
				cout << (float)i/(countALUasap.size()-1) << " ALU bottom-up" << endl;
			else
				cout << (float)i/(countALUasap.size()-1) << " ALU top-down" << endl;
			break;
		}
}

double graph::calForce(int a,int b,int na,int nb,const vector<double>& DG,int delay) const // [a,b]->[na,nb]
{
	if ((na > nb) || (a > b))
		return 0;
	double res = 0, sum = 0;
	for (int i = na; i <= nb+delay-1; ++i)
		sum += DG[i];
	res += sum/(double)(nb-na+1);
	res += (double)(nb-na)/(double)(3*(nb-na+1)); // look-ahead: temp_DG[i]=DG[i]+x(i)/3 => (h-1)^2/(3h^2)+\sum 1/(3h^2)=(h-1)/(3h)
	sum = 0;
	for (int i = a; i <= b+delay-1; ++i)
		sum += DG[i];
	res -= sum/(double)(b-a+1);
	res -= (double)(b-a)/(double)(3*(b-a+1)); // look-ahead
	return res;
}

double graph::calSuccForce(VNode* const& v,int cstep,const map<string,vector<double>>& DG) const
{
	double f = 0;
	for (auto pnode = v->succ.cbegin(); pnode != v->succ.cend(); ++pnode)
		if (mapResourceType((*pnode)->type) == mapResourceType(v->type) && (*pnode)->asap <= cstep && (*pnode)->alap >= cstep) // type should be same
		{
			f += calForce((*pnode)->asap,(*pnode)->alap,cstep+1,(*pnode)->alap,
					DG.at(mapResourceType((*pnode)->type)),(*pnode)->delay);
			if (cstep + 1 == (*pnode)->alap) // recursion
			{
				f += calForce((*pnode)->asap,(*pnode)->alap,cstep+1,cstep+1,DG.at(mapResourceType((*pnode)->type)),(*pnode)->delay);
				f += calSuccForce((*pnode),cstep+1,DG);
				f += calPredForce((*pnode),cstep+1,DG);
			}
		}
	return f;
}

double graph::calPredForce(VNode* const& v,int cstep,const map<string,vector<double>>& DG) const
{
	double f = 0;
	for (auto pnode = v->pred.cbegin(); pnode != v->pred.cend(); ++pnode)
		if (mapResourceType((*pnode)->type) == mapResourceType(v->type) && (*pnode)->asap <= cstep && (*pnode)->alap >= cstep)
		{
			f += calForce((*pnode)->asap,(*pnode)->alap,(*pnode)->asap,cstep-1,
					DG.at(mapResourceType((*pnode)->type)),(*pnode)->delay);
			if (cstep - 1 == (*pnode)->asap)
			{
				f += calForce((*pnode)->asap,(*pnode)->alap,cstep-1,cstep-1,DG.at(mapResourceType((*pnode)->type)),(*pnode)->delay);
				f += calSuccForce((*pnode),cstep-1,DG);
				f += calPredForce((*pnode),cstep-1,DG);
			}
		}
	return f;
}

void graph::TC_FDS() // Time-constrained Force-Directed Scheduling
{
	print("Begin time-constrained force-directed scheduling (FDS)...\n");
	watch.restart();
	topologicalSortingDFS();
	// initialize N_r(t)
	map<string,int> temp;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		temp[pnr->first] = 0;
	for (int i = 0; i <= ConstrainedLatency; ++i) // number+1
		nrt.push_back(temp);

	print("Begin placing operations...");
	clearMark();
	while (numScheduledOp < vertex)
	{
		double minF = MAXINT_;
		int bestop = 0, beststep = 1;

		// build distribution graph
		map<string,vector<double>> DG;// type step dg
		for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		{
			vector<double> temp(ConstrainedLatency + MUL_DELAY,0);
			DG[mapResourceType(pnr->first)] = temp;
		}
		for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
			for (int i = (*pnode)->asap; i <= (*pnode)->alap; ++i)
				for (int d = 0; d < (*pnode)->delay; ++d)
					DG[mapResourceType((*pnode)->type)][i + d] += 1.0 / (double)((*pnode)->getLength());

		// find the op and step with lowest force
		vector<pair<int,pair<int,int>>> fv; // force, op, step
		for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		{
			if ((*pnode)->cstep != 0)
				continue;
			for (int i = (*pnode)->asap; i <= (*pnode)->alap; ++i)
			{
				double f = calForce((*pnode)->asap, (*pnode)->alap, i, i, DG.at(mapResourceType((*pnode)->type)), (*pnode)->delay);
				f += calSuccForce(*pnode, i, DG);
				f += calPredForce(*pnode, i, DG);
				pair<int,int> opstep = {(*pnode)->num,i};
				pair<int,pair<int,int>> temp = {f,opstep};
				fv.push_back(temp);
				// if (f < minF)
				// {
				// 	minF = f;
				// 	bestop = (*pnode)->num;
				// 	beststep = i;
				// }
			}
		}
		std::sort(fv.begin(),fv.end(),
				[](const pair<int,pair<int,int>>& f1, const pair<int,pair<int,int>>& f2) // lambda
				{return f1.first < f2.first;});
		for (auto node : fv)
			if (scheduleNodeStep(adjlist[node.second.first],node.second.second,2))
				break;
		// scheduleNodeStep(adjlist[bestop],beststep,2);
	}
	watch.stop();
	print("Placing operations done!\n");
	print("Finish force-directed scheduling!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::RC_FDS() // Resource-constrained Force-Directed Scheduling
{
	print("Begin resource-constrained force-directed scheduling (FDS)...\n");
	watch.restart();
	topologicalSortingDFS();
	// initialize N_r(t)
	map<string,int> temp,maxNr;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		if (pnr->first == "mul" || pnr->first == "MUL")
		{
			temp[pnr->first] = 0;
			maxNr[pnr->first] = MAXRESOURCE.first;
		}
		else
		{
			temp[pnr->first] = 0;
			maxNr[pnr->first] = MAXRESOURCE.second;
		}
	nrt.push_back(temp); // nrt[0]

	print("Begin placing operations...");
	int cstep = 0;
	vector<VNode*> readyList;
	clearMark();
	while (numScheduledOp < vertex)
	{
		cstep++;
		// extend maximum latency
		if (cstep > ConstrainedLatency)
			for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
				if ((*pnode)->cstep == 0) // operations that have not been scheduled
					(*pnode)->setALAP((*pnode)->alap+1);

		// determine ready operations in cstep
		// (i.e. ops whose time frame intersects the current cstep)
		for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		{
			int flag = 1;
			for (auto ppred = (*pnode)->pred.cbegin(); ppred != (*pnode)->pred.cend(); ++ppred)
				if ((*ppred)->cstep == 0)
					flag = 0;
			if (mark[(*pnode)->num] == 0 && (*pnode)->asap <= cstep && flag == 1)
			{
				readyList.push_back(*pnode);
				mark[(*pnode)->num] = 1; // have been pushed into readyList
			}
		}

		// build distribution graph
		map<string,vector<double>> DG;// type step dg
		for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		{
			vector<double> temp(max(cstep,cdepth) + MUL_DELAY,0);
			DG[mapResourceType(pnr->first)] = temp;
		}
		for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
			for (int i = (*pnode)->asap; i <= (*pnode)->alap; ++i)
				for (int d = 0; d < (*pnode)->delay; ++d)
					DG[mapResourceType((*pnode)->type)][i+d] += 1/(double)((*pnode)->getLength()); // (*pnode)->length

		// sort the readyList by priority function (force) in decresing order
		std::sort(readyList.begin(),readyList.end(),
				[this,cstep,&DG](VNode* const& v1, VNode* const& v2) // lambda
				{
					// original force
					double f1 = calForce(v1->asap,v1->alap,cstep,cstep,DG.at(mapResourceType(v1->type)),v1->delay); // cannot use [], which is non-const
					// successor force
					f1 += calSuccForce(v1,cstep,DG);
					// predecessor force
					f1 += calPredForce(v1,cstep,DG);
					double f2 = calForce(v2->asap,v2->alap,cstep,cstep,DG.at(mapResourceType(v2->type)),v2->delay);
					f2 += calSuccForce(v2,cstep,DG);
					f2 += calPredForce(v2,cstep,DG);
					return (f1 > f2);
				});

		// test if the operations in readyList can be placed in this cstep
		for (int i = 0; i < readyList.size(); )
		{
			// cout << readyList[i]->num+1 << " " << readyList[i]->asap << " " << readyList[i]->alap << endl;
			int flag = 1;
			for (int d = 1; d <= readyList[i]->delay; ++d)
			{
				if (cstep+d-1 >= nrt.size())
					nrt.push_back(temp); // important!
				if (nrt[cstep+d-1][mapResourceType(readyList[i]->type)]+1 > maxNr[mapResourceType(readyList[i]->type)])
					flag = 0;
			}
			if (flag == 1)
			{
				scheduleNodeStepResource(readyList[i],cstep,2); // update time frame
				readyList.erase(readyList.begin()+i);
				i--;
			}
			i++;
		}
	}
	watch.stop();
	print("Placing operations done!\n");
	print("Finish force-directed scheduling!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::printTimeFrame() const // need to be printed before scheduling
{
	cout << "Time frame:" << endl;
	int cnt = 1;
	for (auto pnode : adjlist)
		cout << pnode->num+1 << ": [ " << pnode->asap << " , " << pnode->alap << " ]" << endl;
}

void graph::TC_LS() // Time-constrained List Scheduling
{
	print("Begin time-constrained list scheduling (LS)...\n");
	watch.restart();

	// obtain time frame (ASAP & ALAP)
	topologicalSortingDFS();

	// initialize N_r(t)
	map<string,int> temp;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		temp[pnr->first] = 0;
	for (int i = 0; i <= ConstrainedLatency; ++i) // number+1
		nrt.push_back(temp);

	print("Begin placing operations...");
	clearMark();
	setDegrees();

	std::sort(order.begin(),order.end(),
		[this](VNode* const& v1, VNode* const& v2) // lambda
		// { return (v1->alap - v1->asap < v2->alap - v2->asap); });
		{ return (v1->alap < v2->alap); });

	while (numScheduledOp < vertex)
	{
		vector<VNode*> readyList;
		for (auto pnode = order.cbegin(); pnode != order.cend(); ++pnode)
			if (mark[(*pnode)->num] == 0)
				readyList.push_back(*pnode);
		// printf("Scheduled %d ops. Len readyList: %d\n",numScheduledOp,readyList.size());
		map<string,int> maxNr;
		for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
			maxNr[mapResourceType(pnr->first)] = maxNrt[mapResourceType(pnr->first)] + 1;
		for (auto pnode = readyList.cbegin(); pnode != readyList.cend(); ++pnode)
			for (int step = (*pnode)->asap; step <= (*pnode)->alap; ++step)
			{
				bool flag_in = true;
				for (int delay = 0; delay < (*pnode)->delay; delay++)
					if (nrt[step+delay][mapResourceType((*pnode)->type)] + 1 > maxNr[mapResourceType((*pnode)->type)])
					{
						flag_in = false;
						break;
					}
				if (flag_in)
				{
					scheduleNodeStep(*pnode,step,2);
					// printf("Schedule node %d at %d\n",(*pnode)->num,step);
					mark[(*pnode)->num] = 1;
					break;
				}
			}
	}

	watch.stop();
	print("Placing operations done!\n");

	print("Finish list scheduling!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::RC_LS() // Resource-constrained List Scheduling
{
	print("Begin resource-constrained list scheduling (LS)...\n");
	watch.restart();

	// obtain time frame (ASAP & ALAP)
	topologicalSortingDFS();

	// initialize N_r(t)
	map<string,int> temp,maxNr;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
		if (pnr->first == "MUL")
		{
			temp[pnr->first] = 0;
			maxNr[pnr->first] = MAXRESOURCE.first;
		}
		else
		{
			temp[pnr->first] = 0;
			maxNr[pnr->first] = MAXRESOURCE.second;
		}
	for (int i = 0; i < vertex; ++i)
		nrt.push_back(temp);

	print("Begin placing operations...");
	clearMark();
	setDegrees();
	vector<VNode*> readyList;

	// sort by priority function
	std::sort(order.begin(),order.end(),
			[this](VNode* const& v1, VNode* const& v2) // lambda
			{ return ((v1->alap - v2->asap) < (v1->alap - v2->alap)); });

	// while there're unscheduled operations
	for (int cstep = 1; numScheduledOp < vertex; ++cstep)
	{
		// determine the ready operations
		for (auto pnode = order.cbegin(); pnode != order.cend(); ++pnode)
			if (mark[(*pnode)->num] == 0 && (*pnode)->tempIncoming == 0 && (*pnode)->asap <= cstep) // in-degree = 0
			{
				readyList.push_back(*pnode);
				mark[(*pnode)->num] = 1; // have been pushed into readyList
			}

		// determine eligible operations in c-step and schedule them
		for (int i = 0; i < readyList.size(); )
		{
			bool flag = true;
			for (int d = 1; d <= readyList[i]->delay; ++d)
				if (nrt[cstep+d-1][mapResourceType(readyList[i]->type)]+1 > maxNr.at(mapResourceType(readyList[i]->type)))
					flag = false;
			if (flag)
			{
				scheduleNodeStepResource(readyList[i],cstep,2);
				for (auto pnode = readyList[i]->succ.begin(); pnode != readyList[i]->succ.end(); ++pnode)
					(*pnode)->tempIncoming--;
				readyList.erase(readyList.begin()+i);
				i--;
			}
			i++;
		}
	}

	watch.stop();
	print("Placing operations done!\n");

	print("Finish list scheduling!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}