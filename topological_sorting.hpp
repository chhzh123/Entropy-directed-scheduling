// Copyright (c) 2018 Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file contains the implementation of the topological sorting.

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

void graph::topologicalSortingDFS(bool aslap_order)
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
	if (aslap_order)
		sort(order.begin(),order.end(),[](VNode* const& node1,VNode* const& node2)
			{
				if (node1->alap < node2->alap)
					return true;
				else if (node1->alap == node2->alap)
						if (node1->asap < node2->asap)
							return true;
						else
							return false;
					else
						return false;
			});
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