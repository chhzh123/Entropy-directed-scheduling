// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This head file contains the implementation of the basic graph operations and our EDS algorithm.

#include <iostream>
#include <iomanip>
#include <vector>
#include <utility> // pairs
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "watch.h" // for high-accuracy time counting
// for high-accuracy time counting
stop_watch watch;

#include "graph.h"
#include "initialization.hpp"
#include "topological_sorting.hpp"
#include "output.hpp"
#include "ILP.hpp"
#include "FDS.hpp"
#include "LS.hpp"
#include "EDS.hpp"
using namespace std;

void graph::mainScheduling(int mode)
{
	switch (MODE[0])
	{
		case 0: TC_EDS(0);break;
		case 1: TC_IEDS(0);break;
		case 3: TC_FDS();break;
		case 4: TC_LS();break;
		case 10: RC_EDS(0);break;
		case 11: RC_IEDS(0);break;
		case 13: RC_FDS();break;
		case 14: RC_LS();break;
		default: cout << "Invaild mode!" << endl;return;
	}
	if (mode == 0)
		standardOutput();
	else
		simplifiedOutput();
}