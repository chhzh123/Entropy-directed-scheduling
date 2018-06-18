// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implement of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.
// Our work has been contributed to ICCD 2018.

// This file is the main function of our EDS.

#include <iostream>
#include <fstream>
#include <string>

#include "watch.h"
#include "graph.h"
#include "graph.cpp"
using namespace std;

// set these argv from cmd
// argv[0] default file path: needn't give
// argv[1] dot file num
// argv[2] scheduling mode:
// 			time-constrained(TC):		0  EDS		1 EDS(reverse)		2 LS 		3 FDS
//			resource-constrained(RC):	10 EDS 		11 LS 		12 FDS
//			ILP:						20 TC_ILP 	21 RC_ILP
// ****** If the arguments below are not needed, you needn't type anything more. ******
// argv[3] latency factor (LC)
// argv[4] topo mode: 0 DFS 1 Kahn

// Benchmarks for our experiments can be downloaded at
// https://www.ece.ucsb.edu/EXPRESS/benchmark/
const string dot_file[21] = {
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
	"invert_matrix_general_dfg__3"/*c*/
};

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
const int RC[21][2] = {
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
	{15,11}/*20*/
};

int main(int argc,char *argv[])
{

	cout << "Begin reading dot file..." << endl;
	// if you need to load from other path, please modify here
	string path = "./Benchmarks/";
	int file_num = stoi(string(argv[1]));
	ifstream infile(path + dot_file[file_num] + ".dot");
	if (infile)
		cout << "Read in dot file successfully!\n" << endl;
	else
	{
		cout << "Error: No such files!" << endl;
		return;
	}

	// initial the graph
	graph gp(infile);
	vector<int> MODE = {0,0};
	MODE[0] = stoi(string(argv[2])); // EDS mode
	if (MODE[0] == 0 || MODE[0] == 1)
		MODE[1] = stoi(string(argv[4])); // topo mode
	gp.setMODE(MODE);
	switch (MODE[0])
	{
		// time-constrained
		case 0: // EDS
		case 1: // EDS_rev
		case 2: // LS
		case 3: // FDS
				gp.setLC(stod(string(argv[3])));
				gp.mainScheduling();
				break;
		// resource-constrained
		case 10: // EDS
		case 11: // LS
		case 12: // FDS
				gp.setMAXRESOURCE(RC[file_num][0],RC[file_num][1]);
				gp.mainScheduling();
				break;
		case 20:{ // TC_ILP
					double lc = stod(string(argv[3]));
					gp.setLC(lc);
					ofstream outfile(dot_file[file_num]+"_"+to_string(lc)+".lp");
					gp.generateTC_ILP(outfile);
					outfile.close();
					break;
				}
		case 21:{ // RC_ILP
					ofstream outfile(dot_file[file_num]+".lp");
					gp.generateRC_ILP(outfile);
					outfile.close();
					break;
				}
		default: break;
	}

	infile.close();
	// getchar();
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