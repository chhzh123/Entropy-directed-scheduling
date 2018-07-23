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
					DG[mapResourceType((*pnode)->type)][i + d] += 1 / (double)((*pnode)->length);

		// find the op and step with lowest force
		for (auto pnode = adjlist.cbegin(); pnode != adjlist.cend(); ++pnode)
		{
			if ((*pnode)->cstep != 0)
				continue;
			for (int i = (*pnode)->asap; i <= (*pnode)->alap; ++i)
			{
				double f = calForce((*pnode)->asap, (*pnode)->alap, i, i, DG.at(mapResourceType((*pnode)->type)), (*pnode)->delay);
				f += calSuccForce(*pnode, i, DG);
				f += calPredForce(*pnode, i, DG);
				if (f < minF)
				{
					minF = f;
					bestop = (*pnode)->num;
					beststep = i;
				}
			}
		}
		scheduleNodeStep(adjlist[bestop],beststep,2);
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
					DG[mapResourceType((*pnode)->type)][i+d] += 1/(double)((*pnode)->length);

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