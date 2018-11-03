// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file is the main function of EDS.

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "graph.h"
#include "graph.cpp"
using namespace std;

// Benchmarks for our experiments can be downloaded at
// https://www.ece.ucsb.edu/EXPRESS/benchmark/
const string dot_file[22] = {
	"",
	"hal",
	"horner_bezier_surf_dfg__12"/*c*/,
	"arf",
	"motion_vectors_dfg__7"/*c*/,
	"ewf",
	"fir2",
	"fir1",
	"h2v2_smooth_downsample_dfg__6",
	"feedback_points_dfg__7"/*c*/,
	"collapse_pyr_dfg__113"/*c*/,
	"cosine1",
	"cosine2",
	"write_bmp_header_dfg__7"/*c*/,
	"interpolate_aux_dfg__12"/*c*/,
	"matmul_dfg__3"/*c*/,
	"idctcol_dfg__3"/*c*/,
	"jpeg_idct_ifast_dfg__5"/*c*/,
	"jpeg_fdct_islow_dfg__6"/*c*/,
	"smooth_color_z_triangle_dfg__31"/*c*/,
	"invert_matrix_general_dfg__3"/*c*/,
	// "dag_500",
	// "dag_1000",
	"dag_1500"
};

// if you need to load from other path, please modify here
string path = "./Benchmarks/";

// default constraints for resource-constrained scheduling
// these fine-tuned constraints are first presented in the paper below
// -----------------------------------------------------------------
// @Article{
//  Author	=	{G. Wang and W. Gong and B. DeRenzi and R. Kastner},
//  title	=	{Ant Colony Optimizations for Resource- and Timing-Constrained Operation Scheduling},
//  journal	=	{IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD)},
//  year	=	{2007}
// }
// -----------------------------------------------------------------
const int RC[22][2] = {
	{0,0},/*null*/
	{2,1},/*1*/
	{2,1},/*2*/
	{3,1},/*3*/
	{3,4},/*4*/
	{1,2},/*5*/
	{2,3},/*6*/
	{2,3},/*7*/
	{1,3},/*8*/
	{3,3},/*9*/
	{3,5},/*10*/
	{4,5},/*11*/
	{5,8},/*12*/
	{1,9},/*13*/
	{9,8},/*14*/
	{9,8},/*15*/
	{5,6},/*16*/
	{10,9},/*17*/
	{5,7},/*18*/
	{8,9},/*19*/
	{15,11},/*20*/
	{25,20}/*21*/
};

void interactive()
{
	while (1)
	{
		cout << "\nPlease enter the file num." << endl;
		for (int i = 1; i < 21; ++i)
			cout << i << ": " << dot_file[i] << endl;
		int file_num;
		cin >> file_num;
		cout << "Begin reading dot file..." << endl;
		ifstream infile(path + dot_file[file_num] + ".dot");
		if (infile)
			cout << "Read in dot file successfully!\n" << endl;
		else
		{
			cout << "Error: No such files!" << endl;
			return;
		}
	
		graph gp;
		vector<int> MODE;
		cout << "\nPlease enter the scheduling mode:" << endl;
		cout << "Time-constrained(TC):\t0  EDS(DFS)\t1  EDS(Kahn)\t2  ILP\t3  FDS\t4  LS" << endl;
		cout << "Resource-constrained(RC):\t10 EDS(DFS)\t11 EDS(Kahn)\t12 ILP\t13 FDS\t 14 LS" << endl;
		int mode;
		cin >> mode;
		MODE.push_back(mode);
		if (MODE[0] < 10)
		{
			double lc;
			cout << "Please enter the latency factor." << endl;
			cin >> lc;
			gp.setLC(lc);
		}
		else
			gp.setMAXRESOURCE(RC[file_num][0],RC[file_num][1]);

		cout << "\nPlease enter the scheduling order:" << endl;
		cout << "0. Top-down\t 1.Bottom-up" << endl;
		cin >> mode;
		MODE.push_back(mode);
		gp.setMODE(MODE);
		gp.readFile(infile);

		if (MODE[0] == 2)
		{
			try{
				system("md TC_ILP");
			}
			catch(...){}
			char str[10];
			snprintf(str,sizeof(str),"%.1f",gp.getLC());
			ofstream outfile("./TC_ILP/"+dot_file[file_num]+"_"+string(str)+".lp");
			gp.generateTC_ILP(outfile);
			outfile.close();
		}
		else if (MODE[0] == 12)
		{
			try{
				system("md RC_ILP");
			}
			catch(...){}
			ofstream outfile("./RC_ILP/"+dot_file[file_num]+".lp");
			gp.generateRC_ILP(outfile);
			outfile.close();
		}
		else
			gp.mainScheduling();
	
		infile.close();
	}
}

// set these argv from cmd
// argv[0] default file path: needn't give
// argv[1] scheduling mode:
// 			time-constrained(TC):		0  EDS(DFS)  1  EDS(Kahn)    2  ILP    3  FDS   4  LS
//			resource-constrained(RC):	10 EDS(DFS)  11 EDS(Kahn)    12 ILP    13 FDS   14 LS
// ****** If the arguments below are not needed, you needn't type anything more. ******
// argv[2] latency factor (LC) or scheduling order
//                                0 top-down  1 bottom-up
void commandline(char *argv[])
{
	vector<int> MODE;
	MODE.push_back(stoi(string(argv[1]))); // scheduling mode
	switch (MODE[0])
	{
		case 0:
		case 1:
		case 3:
		case 4: MODE.push_back(stoi(string(argv[3])));break;
		case 10:
		case 11:
		case 13:
		case 14: MODE.push_back(stoi(string(argv[2])));break;
		case 2: MODE.push_back(stoi(string(argv[2])));break;
		case 12: MODE.push_back(stoi(string(argv[1])));break;
		default: cout << "Error: Mode wrong!" << endl;break;
	}

	for (int file_num = 1; file_num < 22; ++file_num)
	{
		ifstream infile(path + dot_file[file_num] + ".dot");
		if (!infile)
		{
			cout << "Error: No such files!" << endl;
			return;
		}

		graph gp;
		gp.setMODE(MODE);
		gp.setPRINT(0);
		gp.readFile(infile);
		if (MODE[0] >= 10)
			gp.setMAXRESOURCE(RC[file_num][0],RC[file_num][1]);
		else
			gp.setLC(stod(string(argv[2])));
		if (MODE[0] == 2)
		{
			try{
				system("md TC_ILP");
			}
			catch(...){}
			char str[10];
			snprintf(str,sizeof(str),"%.1f",gp.getLC());
			ofstream outfile("./TC_ILP/"+dot_file[file_num]+"_"+string(str)+".lp");
			gp.generateTC_ILP(outfile);
			outfile.close();
		}
		else if (MODE[0] == 12)
		{
			try{
				system("md RC_ILP");
			}
			catch(...){}
			ofstream outfile("./RC_ILP/"+dot_file[file_num]+".lp");
			gp.generateRC_ILP(outfile);
			outfile.close();
		}
		else
		{
			cout << "File # " << file_num << " (" << dot_file[file_num] << ") :" <<endl;
			gp.mainScheduling(1);
		}

		infile.close();
	}
}

int main(int argc,char *argv[])
{
	if (argc == 1) // interactive
		interactive();
	else // read from cmd
		commandline(argv);
	return 0;
}

// Some abbreviations appear in ExPRESS
// 1.MUL
// mul: multipy
// div: divide
// 2.ALU
// imp: import
// exp: export
// MemR: memory read
// MemW: memory write
// str: store register
// les: less
// lod: load data into memory?
// bge: ?
// neg: negetive
// add: add
// sub: subtract
// and: logical and
// lsr: logical shift right
// asr: arithmetic shift right
// bne: Branch if Not Equal