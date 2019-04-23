// Copyright (c) 2018 Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This file is the main function of EDS.

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <chrono> // timing

using Clock = std::chrono::high_resolution_clock;

#include "graph.h"
#include "graph.hpp"
using namespace std;

// Benchmarks for our experiments can be downloaded at
// https://www.ece.ucsb.edu/EXPRESS/benchmark/
const vector<string> dot_file = {
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
	"dag_500",
	"dag_1000",
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
// const vector<map<string,int>> RC = {
// 	{{"MUL",0 }, {"ALU",0 }},/*null*/
// 	{{"MUL",2 }, {"ALU",1 }},/*1*/
// 	{{"MUL",2 }, {"ALU",1 }},/*2*/
// 	{{"MUL",3 }, {"ALU",1 }},/*3*/
// 	{{"MUL",3 }, {"ALU",4 }},/*4*/
// 	{{"MUL",1 }, {"ALU",2 }},/*5*/
// 	{{"MUL",2 }, {"ALU",3 }},/*6*/
// 	{{"MUL",2 }, {"ALU",3 }},/*7*/
// 	{{"MUL",1 }, {"ALU",3 }},/*8*/
// 	{{"MUL",3 }, {"ALU",3 }},/*9*/
// 	{{"MUL",3 }, {"ALU",5 }},/*10*/
// 	{{"MUL",4 }, {"ALU",5 }},/*11*/
// 	{{"MUL",5 }, {"ALU",8 }},/*12*/
// 	{{"MUL",1 }, {"ALU",9 }},/*13*/
// 	{{"MUL",9 }, {"ALU",8 }},/*14*/
// 	{{"MUL",9 }, {"ALU",8 }},/*15*/
// 	{{"MUL",5 }, {"ALU",6 }},/*16*/
// 	{{"MUL",10}, {"ALU",9 }},/*17*/
// 	{{"MUL",5 }, {"ALU",7 }},/*18*/
// 	{{"MUL",8 }, {"ALU",9 }},/*19*/
// 	{{"MUL",15}, {"ALU",11}},/*20*/
// 	{{"MUL",5 }, {"ALU",9 }},/*21*/
// 	{{"MUL",6 }, {"ALU",12}},/*22*/
// 	{{"MUL",7 }, {"ALU",13}},/*23*/
// };

// Complex FU library
const vector<map<string,int>> RC = {
	{{"MUL",0 }, {"ALU",0 }},/*null*/
	{{"MUL",2 }, {"add",1 },{"sub",1},{"les",1}},/*1*/
	{{"MUL",1 }, {"ADD",1 },{"LOD",1},{"STR",1}},/*2*/
	{{"MUL",3 }, {"ADD",1 }},/*3*/
	{{"MUL",3 }, {"LOD",1 },{"ADD",2},{"STR",1}},/*4*/
	{{"MUL",1 }, {"ADD",2 }},/*5*/
	{{"MUL",2 }, {"add",1 },{"exp",1},{"imp",2}},/*6*/
	{{"MUL",2 }, {"ADD",2 },{"MemR",2},{"MemW",1}},/*7*/
	{{"MUL",1 }, {"ADD",2 },{"ASR",1},{"STR",1},{"LOD",1}},/*8*/
	{{"MUL",3 }, {"STR",2 },{"LOD",1},{"BGE",1},{"ADD",2}},/*9*/
	{{"MUL",3 }, {"ADD",3 },{"SUB",1},{"STR",3},{"LSL",1},{"LOD",3},{"ASR",1}},/*10*/
	{{"MUL",4 }, {"imp",6 },{"sub",1},{"exp",2},{"add",2}},/*11*/
	{{"MUL",4 }, {"add",1 },{"exp",2},{"imp",2},{"sub",2}},/*12*/
	{{"MUL",1 }, {"STR",3 },{"LSR",1},{"LOD",4},{"BNE",1},{"ASR",2},{"AND",2},{"ADD",4}},/*13*/
	{{"MUL",9 }, {"ADD",4 },{"SUB",2},{"STR",2},{"LOD",5}},/*14*/
	{{"MUL",8 }, {"STR",2 },{"LOD",3},{"ADD",3}},/*15*/
	{{"MUL",4 }, {"SUB",2 },{"STR",2},{"LSL",1},{"LOD",2},{"ASR",2},{"ADD",2}},/*16*/
	{{"MUL",4 }, {"SUB",1 },{"STR",2},{"LOD",4},{"ASR",1},{"ADD",4}},/*17*/
	{{"MUL",4 }, {"SUB",2 },{"STR",2},{"LOD",4},{"ASR",1},{"ADD",4}},/*18*/
	{{"MUL",8 }, {"SUB",3 },{"ADD",6},{"LOD",6}},/*19*/
	{{"MUL",14}, {"SUB",3 },{"STR",3},{"NEG",2},{"LOD",8},{"ADD",8}},/*20*/
	{{"MUL",5 }, {"add",9 }},/*21*/
	{{"MUL",6 }, {"add",12}},/*22*/
	{{"MUL",7 }, {"add",13}},/*23*/
};

void interactive()
{
	while (1)
	{
		cout << "\nPlease enter the file num." << endl;
		for (int file_num = 1; file_num < dot_file.size(); ++file_num)
			cout << file_num << ": " << dot_file[file_num] << endl;
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
		cout << "Time-constrained(TC):\t0  EDS\t1  IEDS\t2  ILP\t3  FDS\t4  LS" << endl;
		cout << "Resource-constrained(RC):\t10 EDS\t11 IEDS\t12 ILP\t13 FDS\t 14 LS" << endl;
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
			gp.setMAXRESOURCE(RC.at(file_num));

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
// 			time-constrained(TC):		0  EDS    1  IEDS    2  ILP    3  FDS   4  LS
//			resource-constrained(RC):	10 EDS    11 IEDS    12 ILP    13 FDS   14 LS
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

	for (int file_num = 1; file_num < dot_file.size(); ++file_num)
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
			gp.setMAXRESOURCE(RC.at(file_num));
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