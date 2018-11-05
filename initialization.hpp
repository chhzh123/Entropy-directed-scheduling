// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file contains the implementation of the graph initialization part.

#include <regex> // regular expression for string split
#include <fstream>
using namespace std;

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
	if (nr.find(mapResourceType(type)) != nr.end()) // exist
		nr[mapResourceType(type)]++;
	else
	{
		typeNum++;
		nr[mapResourceType(type)] = 1;
		r_delay[mapResourceType(type)] = mapResourceType(type) == "MUL" ? MUL_DELAY : 1;
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

bool graph::newScheduleNodeStep(VNode* const& node,int step)
{
	// cout << node->num << " " << step << endl;
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

std::vector<std::string> split(const std::string& input, const std::string& regex) // split the string
{
	// passing -1 as the submatch index parameter performs splitting
	std::regex re(regex);
	std::sregex_token_iterator
		first{ input.begin(), input.end(), re, -1 },
		last;
	return { first, last };
}