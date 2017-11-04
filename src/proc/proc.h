#ifndef PROC_H__
#define PROC_H__

#include "io.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>

using namespace std;

namespace g{
namespace proc{

//================== Sheng added for Transfer alignment between Renmin_style and Sheng_style ==================//

int Ali_To_Cor(vector <int> &ali2, vector <vector <int> > &AFP_Cor);
	
void Ali1_To_AliPair(int n1,int n2,vector <int> &ali1,
	vector<pair<int,int> > & alignment_out);

void Ali2_To_AliPair(int n1,int n2,vector <int> &ali2,
	vector<pair<int,int> > & alignment_out);

void AliPair_To_Ali1_Ali2(int n1,int n2,
	vector<pair<int,int> > & alignment_in,
	vector <int> &ali1,vector <int> &ali2);

void From_Align_Get_Bound(int moln1,int moln2,vector<pair<int,int> > &align,
	vector<pair<int,int> > &bound,int neib);

void Renmin_To_Sheng_align(int moln1,int moln2,
	vector<pair<int,int> > &align_in, vector<pair<int,int> > &align_out);

void Sheng_To_Renmin_align(int moln1,int moln2,
	vector<pair<int,int> > &align_in, vector<pair<int,int> > &align_out);

void Sheng_To_Renmin_bound(int moln1,int moln2,
	vector<pair<int,int> > &bound_in, vector<pair<int,int> > &bound_out);

//================================= sheng modify ======================//over

//-- Z normalization --//
void ZScoreNormalize(std::vector<double>& signals, double* avg = NULL, double* stdev = NULL);

//-- Peak pick --//
void LaplaceDiff(const std::vector<double>& raw, std::vector<double>& ldiff);
void Diff(const std::vector<double>& raw, std::vector<double>& diff);
void PeakPick(const std::vector<double>& raw, std::vector<std::pair<int, double> >& peaks);

//-- Dynamic Time Warping --//
double DynamicTimeWarping_global(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment);
double DynamicTimeWarping_local(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment);
double DynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment);

//-- constrained (bound) DTW --//
double BoundDynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment);
/** Left path is blocked*/
double BoundDynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment);

}

}


#endif

