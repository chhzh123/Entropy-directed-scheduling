// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file contains the implementation of the integer linear programming (ILP).

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