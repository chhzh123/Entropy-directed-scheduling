// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file contains the implementation of the main EDS algorithm.

void graph::TC_EDS(int order_mode)
{
	print("Begin EDS...\n");
	watch.restart();
	topologicalSortingDFS(order_mode);
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
	print("Finish EDS!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::TC_IEDS(int order_mode)
{
	print("Begin IEDS...\n");
	watch.restart();
	topologicalSortingDFS(order_mode);
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
	// first set virtual constrained resource numbers
	map<string,vector<int>> gr;
	for (auto pnr = nr.cbegin(); pnr != nr.cend(); ++pnr)
	{
		int totnr = pnr->second * r_delay[pnr->first];
		int num = (totnr % ConstrainedLatency == 0 ? totnr / ConstrainedLatency : totnr / ConstrainedLatency + 1); // ceiling function
		int q = totnr - ConstrainedLatency * (num - 1);
		vector<int> n(1,0);
		for (int i = 0; i < q; ++i)
			n.push_back(num);
		for (int i = 0; i < ConstrainedLatency - q; ++i)
			n.push_back(num-1);
		gr[pnr->first] = n;
	}
	print("Begin placing other nodes...");
	for (auto pnode = edsOrder.cbegin(); pnode != edsOrder.cend(); ++pnode)
	{
		int a = (*pnode)->asap, b = (*pnode)->alap;
		bool flag_out = false;
		// cout << (*pnode)->name << " " << a << " " << b << endl;
		for (int t = a; t <= b; ++t)
		{
			bool flag = true;
			for (int d = 1; d <= (*pnode)->delay; ++d)
			{
				string tempType = mapResourceType((*pnode)->type);
				if (nrt[t+d-1][tempType]+1 > gr[tempType][t+d-1])
					flag = 0;
			}
			if (flag)
			{
				scheduleNodeStep(*pnode,t,2);
				flag_out = true;
				break;
			}
		}
		if (!flag_out)
		{
			double minnrt = MAXINT_;
			int minstep = a;
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
			scheduleNodeStep(*pnode,minstep,2);
		}
	}
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
				if (t+d-2>0 && nrt[t+d-1][tempType] >= nrt[t+d-2][tempType] + 1)
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
	watch.stop();
	print("Finish IEDS!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
}

void graph::RC_EDS(int sorting_mode) // Resource-constrained EDS
{}

void graph::RC_IEDS(int order_mode) // Resource-constrained EDS
{
	print("Begin IEDS...\n");
	watch.restart();
	topologicalSortingDFS(order_mode);
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
	print("Finish IEDS!\n");
	cout << "Total time used: " << watch.elapsed() << " micro-seconds" << endl;
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