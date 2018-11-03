// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file contains the implementation of the output part.

#include <fstream>
using namespace std;

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

void graph::countResource() const
{
	int sum_r = 0;
	ofstream out("./Resource_"+to_string(LC)+".out",ios::app);
	for (auto ptype = nr.crbegin(); ptype != nr.crend(); ++ptype)
	{
		cout << ptype->first << ": " << maxNrt.at(ptype->first) << endl;
		out << maxNrt.at(ptype->first) << " ";
		sum_r += maxNrt.at(ptype->first);
		if (PRINT)
		// if (ptype->first == "MUL" || ptype->first == "mul")
			for (int i = 1; i <= maxLatency; ++i) // ConstrainedLatency
				cout << "Step " << i << ": " << nrt[i].at(mapResourceType(ptype->first)) << endl;
	}
	out << "\t" << sum_r << endl;
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
	printGanttGraph();
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

void graph::printGanttGraph() const
{
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

void graph::print(const string str) const
{
	if (PRINT)
		cout << str << endl;
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

void graph::printTimeFrame() const // need to be printed before scheduling
{
	cout << "Time frame:" << endl;
	int cnt = 1;
	for (auto pnode : adjlist)
		cout << pnode->num+1 << ": [ " << pnode->asap << " , " << pnode->alap << " ]" << endl;
}