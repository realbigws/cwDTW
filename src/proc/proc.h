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
//-> AfpCor relevant
long Ali1_To_AfpCor(long n1,long n2,vector <long> &ali1, 
	vector <vector <long> > &AFP_Cor, long thres=0);

long Ali2_To_AfpCor(long n1,long n2,vector <long> &ali2, 
	vector <vector <long> > &AFP_Cor, long thres=0);

void AfpCor_To_Ali1_Ali2(long n1,long n2,
	vector <vector <long> > &AFP_Cor,
	vector <long> &ali1,vector <long> &ali2);

//-> AliPair relevant
void Ali1_To_AliPair(long n1,long n2,vector <long> &ali1,
	vector<pair<long,long> > & alignment_out);

void Ali2_To_AliPair(long n1,long n2,vector <long> &ali2,
	vector<pair<long,long> > & alignment_out);

void AliPair_To_Ali1_Ali2(long n1,long n2,
	vector<pair<long,long> > & alignment_in,
	vector <long> &ali1,vector <long> &ali2);

//-> given alignment return boundary
void From_Align_Get_Bound(long moln1,long moln2,vector<pair<long,long> > &align,
	vector<pair<long,long> > &bound,long neib);

void Renmin_To_Sheng_align(long moln1,long moln2,
	vector<pair<long,long> > &align_in, vector<pair<long,long> > &align_out);

void Sheng_To_Renmin_align(long moln1,long moln2,
	vector<pair<long,long> > &align_in, vector<pair<long,long> > &align_out);

void Sheng_To_Renmin_bound(long moln1,long moln2,
	vector<pair<long,long> > &bound_in, vector<pair<long,long> > &bound_out);

//================================= sheng modify ======================//over

//-- Z normalization --//
void ZScoreNormalize(std::vector<double>& signals, double* avg = NULL, double* stdev = NULL);

//-- Peak pick --//
void LaplaceDiff(const std::vector<double>& raw, std::vector<double>& ldiff);
void Diff(const std::vector<double>& raw, std::vector<double>& diff);
void PeakPick(const std::vector<double>& raw, std::vector<std::pair<long, double> >& peaks);

//-- Dynamic Time Warping --//
double DynamicTimeWarping_global(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<long,long> >& alignment);
double DynamicTimeWarping_local(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<long,long> >& alignment);
double DynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<long,long> >& alignment);

//-- constrained (bound) DTW --//
double BoundDynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<long,long> >& bound, std::vector<std::pair<long,long> >& alignment);
/** Left path is blocked*/
double BoundDynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<long,long> >& bound, std::vector<std::pair<long,long> >& alignment);

}

}


#endif

