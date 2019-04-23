// Copyright (c) 2018 Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file contains the implementation of list scheduling (LS).

void graph::TC_LS() // Time-constrained List Scheduling
{
	print("Begin time-constrained list scheduling (LS)...\n");
	auto t1 = Clock::now();

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

	auto t2 = Clock::now();
	print("Placing operations done!\n");

	print("Finish list scheduling!\n");
	cout << "Total time used: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " ns" << endl;
}

void graph::RC_LS() // Resource-constrained List Scheduling
{
	print("Begin resource-constrained list scheduling (LS)...\n");
	auto t1 = Clock::now();

	// obtain time frame (ASAP & ALAP)
	topologicalSortingDFS();

	// initialize N_r(t)
	map<string,int> temp,maxNr;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
	{
		temp[pnr->first] = 0;
		maxNr[pnr->first] = MAXRESOURCE[pnr->first];
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

	auto t2 = Clock::now();
	print("Placing operations done!\n");

	print("Finish list scheduling!\n");
	cout << "Total time used: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " ns" << endl;
}