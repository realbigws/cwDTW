#include "proc.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <memory.h>
#include <climits>
#include <cfloat>
#define     EQN_EPS     1.e-9

using namespace std;

//================== Sheng added for Transfer alignment between Renmin_style and Sheng_style ==================//

//-------- Ali_To_Cor -------------//
long g::proc::Ali_To_Cor(vector <long> &ali2, vector <vector <long> > &AFP_Cor)
{
	long i,k;
	long num;
	long ii,jj;
	long count;
	long isFirst;
	long isLast;
	long type;
	long head1,head2;
	long index;

	//init
	count=-999999;
	num=0;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	type=0;
	ii=-1;
	jj=-1;
	long moln2=(long)ali2.size();
	long thres=0;
	AFP_Cor.clear();
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>=thres) 
				{
					vector <long> tmp_rec;
					tmp_rec.push_back(head1);
					tmp_rec.push_back(head2);
					tmp_rec.push_back(count);
					AFP_Cor.push_back(tmp_rec);
					num+=count;
				}
				count=0;
				isFirst=1;
			}
			continue;
		}

		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=ali2[i];
			type=1; // >0 mode
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
			continue;
		}
		if(i==jj+1&&ali2[i]==ii+1)
		{
			ii=ali2[i];
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>=thres) 
		{
			vector <long> tmp_rec;
			tmp_rec.push_back(head1);
			tmp_rec.push_back(head2);
			tmp_rec.push_back(count);
			AFP_Cor.push_back(tmp_rec);
			num+=count;
		}

		if(isLast==1)goto end;
		else goto ws_init;
	}

	if(count==999999)goto end;
	isLast=1;
	goto ws_end;
end:
	return num;
}

//------- given Ali1 return AliPair ------//
void g::proc::Ali1_To_AliPair(long n1,long n2,vector <long> &ali1,
	vector<pair<long,long> > & alignment_out)
{
	//init
	alignment_out.clear();
	//start
	long i,j;
	long ii,jj;
	long wlen;
	long pre_ii=0;
	long pre_jj=0;
	for(i=1;i<=n1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				alignment_out.push_back (pair<long,long>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment_out.push_back (pair<long,long>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment_out.push_back (pair<long,long>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=n1;i++)alignment_out.push_back (pair<long,long>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=n2;i++)alignment_out.push_back (pair<long,long>(-n1, i));  //Iy
}

//------- given Ali2 return AliPair ------//
void g::proc::Ali2_To_AliPair(long n1,long n2,vector <long> &ali2,
	vector<pair<long,long> > & alignment_out)
{
	//from ali1 to ali2
	vector <long> ali1(n1,-1);
	for(long i=0;i<n2;i++)
	{
		if(ali2[i]==-1)continue;
		long ii=ali2[i];
		ali1[ii]=i;
	}
	//proc
	Ali1_To_AliPair(n1,n2,ali1,alignment_out);
}

//------- given AliPair return Ali1 and Ali2 ------//
void g::proc::AliPair_To_Ali1_Ali2(long n1,long n2,
	vector<pair<long,long> > & alignment_in,
	vector <long> &ali1,vector <long> &ali2)
{
	//init
	ali1.resize(n1);
	ali2.resize(n2);
	for(long i=0;i<n1;i++)ali1[i]=-1;
	for(long i=0;i<n2;i++)ali2[i]=-1; 
	//proc
	for(long i=0;i<(long)alignment_in.size();i++)
	{
		long ii=alignment_in[i].first;
		long jj=alignment_in[i].second;
		if(ii>0 && jj>0)
		{
			ali1[ii-1]=jj-1;
			ali2[jj-1]=ii-1;
		}
	}
}


//---------- dynamic programming ----------//(bounded version)
// [motivation]: given an initial alignment during structural alignment,
// we could only search for a very narrow adjacent area around the alignment.
// That is to say, the searching space for DP could be restricted in O(n),
// as well as the searching time.
// [structure]: we use "bound" to record the additional data structure,
// length is n1, first is n1's correspondence, second is end position

//---- part 1. generate bound from a given alignment ------//
void g::proc::From_Align_Get_Bound(long moln1,long moln2,vector<pair<long,long> > &align,
	vector<pair<long,long> > &bound,long neib)
{
	long k;
	long ii,jj;
	long pre_ii,pre_jj;
	long size=(long)align.size();
	long first_=1;
	//[1]get real alignment
	bound.resize(moln1+1);
	bound[0].first=0;
	bound[0].second=0;
	pre_jj=0;
	for(k=0;k<size;k++)
	{
		ii=align[k].first;
		jj=align[k].second;
		if(ii>0&&jj>0)
		{
			bound[ii].first=jj;  //positive
			bound[ii].second=jj;
			pre_ii=ii+1;
			pre_jj=jj+1;
			//-> clear first
			first_=1;
		}
		else
		{
			//-> get next aligned position pre_jj_
			long pre_jj_=pre_jj;
			if(first_==1)
			{
				first_=0;
				long found=0;
				for(long k_=k+1;k_<size;k_++)
				{
					long ii_=align[k_].first;
					long jj_=align[k_].second;
					if(ii_>0&&jj_>0)
					{
						found=1;
						pre_jj_=jj_;
						break;
					}
				}
				if(found==0)
				{
					pre_jj_=moln2;
				}
			}
			//-> assing bound
			if(ii>0)
			{
				bound[ii].first=pre_jj;
				bound[ii].second=pre_jj_;
			}
			else
			{
				bound[abs(ii)].second++;
			}
		}
	}
	//[2]get horizontal expand
	for(k=0;k<=moln1;k++)
	{
		ii=bound[k].first;
		jj=bound[k].second;
		if(ii-neib<0)ii=0;
		else ii=ii-neib;
		if(jj+neib>moln2)jj=moln2;
		else jj=jj+neib;
		bound[k].first=ii;
		bound[k].second=jj;
	}
	//[3]get vertical expand
	long i;
	long first=1;
	long wscurr,wsprev;
	long wscurr_start=1;
	long wsprev_start=1;
	long ww1,ww2;
	wsprev=0;
	wscurr=0;
	for(k=0;k<size;k++)
	{
		ii=align[k].first;
		jj=align[k].second;
		if(ii>0&&jj>0)
		{
			if(first==1)
			{
				wscurr=ii;
				wscurr_start=jj;
				first=0;
				//record
				for(i=wsprev-neib;i<=wscurr+neib;i++)
				{
					if(i<1)continue;
					if(i>moln1)break;
					ww1=wsprev_start;
					if(ww1<1)ww1=1;
					ww2=wscurr_start;
					if(ww2>moln2)ww2=moln2;
					if(ww1<bound[i].first)bound[i].first=ww1;
					if(ww2>bound[i].second)bound[i].second=ww2;
				}
			}
		}
		else
		{
			if(first==0)
			{
				wsprev=abs(align[k-1].first)+1;
				wsprev_start=abs(align[k-1].second)+1;
				first=1;
			}
		}
	}
	//final record
	if(wsprev==moln1+1)
	{
		wscurr_start=moln2;
		for(i=wsprev-neib-1;i<=moln1;i++)
		{
			if(i<1)continue;
			if(i>moln1)break;
			ww1=wsprev_start;
			if(ww1<1)ww1=1;
			ww2=wscurr_start;
			if(ww2>moln2)ww2=moln2;
			if(ww1<bound[i].first)bound[i].first=ww1;
			if(ww2>bound[i].second)bound[i].second=ww2;
		}
	}
}

//-------------- Ren-min alignment to Sheng ---------------//
//[note1]: Ren-min's alignment starts from 0, and no negative value !!
//         we use 'first-come never-gone' strategy
//[note2]: moln1 MUST less than moln2
void g::proc::Renmin_To_Sheng_align(long moln1,long moln2,
	vector<pair<long,long> > &align_in, vector<pair<long,long> > &align_out)
{
	vector <long> ali1(moln1,-1);
	vector <long> ali2(moln2,-1);
	for(long i=0;i<(long)align_in.size();i++)
	{
		long pos1=align_in[i].first;
		long pos2=align_in[i].second;
		if(ali1[pos1]==-1 && ali2[pos2]==-1)
		{
			ali1[pos1]=pos2;
			ali2[pos2]=pos1;
		}
	}
	//-> ali1 to align
	Ali1_To_AliPair(moln1,moln2,ali1,align_out);
}

//-------------- Sheng alignment to Ren-min ---------------//
//[note1]: Sheng's alignment starts from 1, and have negative value !!
//         we use 'first-come never-gone' strategy
//[note2]: moln1 MUST less than moln2
void g::proc::Sheng_To_Renmin_align(long moln1,long moln2,
	vector<pair<long,long> > &align_in, vector<pair<long,long> > &align_out)
{
	vector <long> ali1(moln1,-1);
	vector <long> ali2(moln2,-1);
	AliPair_To_Ali1_Ali2(moln1,moln2,align_in,ali1,ali2);
	//-> ali1 to align
	vector<pair<long,long> > align_tmp;
	Ali1_To_AliPair(moln1,moln2,ali1,align_tmp);
	//-> get start and end
	long start=-1;
	for(long i=0;i<(long)align_tmp.size();i++)
	{
		long pos1=align_tmp[i].first;
		long pos2=align_tmp[i].second;
		if(pos1>0 && pos2>0)
		{
			if(pos1==1 || pos2==1)continue;
			start=i;
			break;
		}	
	}
	long end=-1;
	for(long i=(long)align_tmp.size()-1;i>=0;i--)
	{
		long pos1=align_tmp[i].first;
		long pos2=align_tmp[i].second;
		if(pos1>0 && pos2>0)
		{
			if(pos1==moln1-1 || pos2==moln2-1)continue;
			end=i;
			break;
		}	
	}
	if(start==-1 || end==-1)
	{
		align_out.clear();
		return;
	}
	//assign
	align_out.clear();
	align_out.push_back(pair<long,long>(0,0));
	for(long i=start;i<=end;i++)
	{
		long pos1=align_tmp[i].first-1;
		long pos2=align_tmp[i].second-1;
		align_out.push_back(pair<long,long>(pos1,pos2));
	}
	align_out.push_back(pair<long,long>(moln1-1,moln2-1));
}

//-------------- Sheng bound to Ren-min ---------------//
void g::proc::Sheng_To_Renmin_bound(long moln1,long moln2,
	vector<pair<long,long> > &bound_in, vector<pair<long,long> > &bound_out)
{
	bound_out.resize(moln1);
	for(long i=1;i<=moln1;i++)
	{
		long pos1=bound_in[i].first-1;
		long pos2=bound_in[i].second-1;
		if(pos1<0)pos1=0;
		if(pos2<0)pos2=0;
		bound_out[i-1].first=pos1;
		bound_out[i-1].second=pos2;
	}
}

//================================= sheng modify ======================//over

//--------- Z normalization ------------//
void g::proc::ZScoreNormalize(std::vector< double >& signals, double* avg, double* stdev)
{
	double sum = std::accumulate(signals.begin(), signals.end(), 0.0);
	double mean =  sum / signals.size();

	double acc = 0.0;
	for(long i = signals.size(); i--;){
		signals[i] = signals[i]-mean;
		acc += signals[i]*signals[i];
	}

	double deviation = std::sqrt(acc/signals.size());
	
	for(long i = signals.size(); i--;){
		signals[i] /= deviation;
	}
	
	if(avg){*avg = mean;}
	if(stdev){*stdev = deviation;}
}

//--------- Peak picking -------------//
void g::proc::Diff(const std::vector<double>& raw, std::vector<double>& diff)
{
	diff.resize(raw.size());
	diff[0] = 0;
	for(long i = 1; i < raw.size(); i++){
		diff[i] = raw[i]-raw[i-1];
	}
}

void g::proc::LaplaceDiff(const std::vector<double>& raw, std::vector<double>& ldiff)
{
	ldiff.resize(raw.size());
	double tmplt[] = {-1,2,-1};
	long tsize = 3;
	long h_tsize = tsize/2;
	long tbegin = h_tsize;
	long tend = raw.size()-h_tsize;
	
	memset(&ldiff[0], 0, sizeof(double)*ldiff.size());
	
	for(long i = tbegin; i < tend; i++){
		for(long j = 0; j < tsize; j++){
			ldiff[i] += tmplt[j]*raw[j+i-h_tsize];
		}
	}
}

void g::proc::PeakPick(const std::vector<double>& raw, std::vector<std::pair<long, double> >& peaks)
{
	peaks.clear();
	std::vector<double> diff;
	Diff(raw, diff);

	peaks.push_back(std::make_pair(0, raw[0]));
	
	for(long i = 0; i < diff.size()-1; i++){
		if((diff[i] > EQN_EPS && diff[i+1] < -1*EQN_EPS) || (diff[i] < -1*EQN_EPS && diff[i+1] > EQN_EPS)){		//local max min
			peaks.push_back(std::make_pair(i, raw[i]));
		}
	}
	
	peaks.push_back(std::make_pair(raw.size()-1, raw[raw.size()-1]));
}


//--------------- Dynamic Time Warping --------------//(global)
double g::proc::DynamicTimeWarping_global(
	const std::vector<double>& seq1, const std::vector<double>& seq2, 
	std::vector<std::pair<long,long> >& alignment)
{
	//-- create score matrix --//
	double* score[seq1.size()];
	for(long i = 0; i < seq1.size(); i++){
		score[i] = new double[seq2.size()];
	}

	//-- initialize score matrix --//
	for(long i = 0; i < seq1.size(); i++){
		for(long j = 0; j < seq2.size(); j++){
			score[i][j] = fabs(seq1[i]-seq2[j]);
		}
	}

	//-- initialize X-axix --//
	for(long i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}

	//-- initialize Y-axix --//
	for(long j = 1; j < seq2.size(); j++){
		score[0][j] += score[0][j-1];
	}

	//-- fill-up score matrix --//
	for(long i = 1; i < seq1.size(); i++){
		for(long j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(std::min(score[i-1][j], score[i][j-1]), score[i-1][j-1]);
		}
	}

	//-- obtain maximal score --//
	double diff = score[seq1.size()-1][seq2.size()-1];

	//-- generate alignment --//
	alignment.clear();
	long i = seq1.size()-1, j = seq2.size()-1;
	while(true){
		alignment.push_back(std::make_pair(i,j));
		long ipre = i-1 < 0 ? 0 : i-1;
		long jpre = j-1 < 0 ? 0 : j-1;
		
		double premin = std::min(std::min(score[ipre][j], score[i][jpre]), score[ipre][jpre]);
		
		if(premin == score[ipre][jpre]){
			i = ipre; j = jpre;
		}
		if(premin == score[ipre][j]){
			i = ipre;
		}
		if(premin == score[i][jpre]){
			j = jpre;
		}
		if(i == 0 && j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	std::reverse(alignment.begin(), alignment.end());

	//-- delete created matrix --//
	for(long i = 0; i < seq1.size(); i++){
		delete [] score[i];
	}

	//-- return distance --//
	return diff;
}

//--------------- Dynamic Time Warping --------------//(local)
//-> we let X (i.e., seq1) be the smaller input
double g::proc::DynamicTimeWarping_local(
	const std::vector<double>& seq1, const std::vector<double>& seq2, 
	std::vector<std::pair<long,long> >& alignment)
{
	//-- create score matrix --//
	double* score[seq1.size()];
	for(long i = 0; i < seq1.size(); i++){
		score[i] = new double[seq2.size()];
	}

	//-- initialize score matrix --//
	for(long i = 0; i < seq1.size(); i++){
		for(long j = 0; j < seq2.size(); j++){
			score[i][j] = fabs(seq1[i]-seq2[j]);
		}
	}

	//-- initialize X-axix --//
	for(long i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}

	//-- initialize Y-axix --//
//	for(long j = 1; j < seq2.size(); j++){
//		score[0][j] += score[0][j-1];
//	}

	//-- fill-up score matrix --//
	for(long i = 1; i < seq1.size(); i++){
		for(long j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(std::min(score[i-1][j], score[i][j-1]), score[i-1][j-1]);
		}
	}

	//-- obtain maximal score --//
	long minval=DBL_MAX;
	long min_j=0;
	for(long j = 0; j < seq2.size(); j++){
		if(score[seq1.size()-1][j]<minval){
			minval=score[seq1.size()-1][j];
			min_j=j;
		}
	}	
	double diff = score[seq1.size()-1][min_j];

	//-- generate alignment --//
	alignment.clear();
	long i = seq1.size()-1, j = min_j;
	while(true){
		alignment.push_back(std::make_pair(i,j));
		long ipre = i-1 < 0 ? 0 : i-1;
		long jpre = j-1 < 0 ? 0 : j-1;
		
		double premin = std::min(std::min(score[ipre][j], score[i][jpre]), score[ipre][jpre]);
		
		if(premin == score[ipre][jpre]){
			i = ipre; j = jpre;
		}
		if(premin == score[ipre][j]){
			i = ipre;
		}
		if(premin == score[i][jpre]){
			j = jpre;
		}
//		if(i == 0 /&& j == 0){
		if( i == 0 ){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	std::reverse(alignment.begin(), alignment.end());

	//-- delete created matrix --//
	for(long i = 0; i < seq1.size(); i++){
		delete [] score[i];
	}

	//-- return distance --//
	return diff;
}


//------- if not specified, then DTW indicates global DTW ---------//
double g::proc::DynamicTimeWarping(
	const std::vector<double>& seq1, const std::vector<double>& seq2, 
	std::vector<std::pair<long,long> >& alignment)
{	
	return DynamicTimeWarping_global(seq1,seq2,alignment);
}


//---------------------- constrained (bound) Dynamnic Time Warping ------------------------//
static inline double& SCORE(long i, long j, 
	std::vector<std::vector<double> > &score, const std::vector<std::pair<long,long> >& bound)
{
	static double invalid = DBL_MAX;
	if(bound[i].first <= j && j <= bound[i].second){	
		return score[i][j-bound[i].first];
	}
	else{
		return invalid;
	}
}

double g::proc::BoundDynamicTimeWarping(
	const std::vector<double>& sequence1, const std::vector<double>& sequence2, 
	const std::vector<std::pair<long,long> >& bound, std::vector<std::pair<long,long> >& alignment)
{
	/*check order*/
	const std::vector<double>* seq1, * seq2;
	seq1 = &sequence1;
	seq2 = &sequence2;

	//--- init score ---//	
	std::vector<std::vector<double> > score;
	score.resize(seq1->size());
	for(long i = 0; i < seq1->size(); i++){
		for(long j = bound[i].first; j <= bound[i].second; j++){
			score[i].push_back(std::fabs((*seq1)[i]-(*seq2)[j]));
		}
	}

	for(long i = 1; i < seq1->size(); i++){
		if(SCORE(i, 0, score, bound) == DBL_MAX){
			break;
		}
		
		SCORE(i, 0, score, bound) += SCORE(i-1, 0, score, bound);
	}
	
	for(long j = 1; j <= bound[0].second; j++){
		score[0][j] += score[0][j-1];
	}
	
	for(long i = 1; i < seq1->size(); i++){

printf("cur=%d\r",i);

		for(long j = bound[i].first; j <= bound[i].second; j++){
			double acc = std::min(std::min(SCORE(i-1, j, score, bound), SCORE(i, j-1, score, bound)), SCORE(i-1, j-1, score, bound));
			if(acc == DBL_MAX){
				score[i][j-bound[i].first] = acc;
			}
			else{
				score[i][j-bound[i].first] += acc;
			}
		}
	}
	
	double diff = SCORE(seq1->size()-1, seq2->size()-1, score, bound);
	
	
	long i = seq1->size()-1, j = seq2->size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		long ipre = i-1 < 0 ? 0 : i-1;
		long jpre = j-1 < 0 ? 0 : j-1;
		
		double scipjp = SCORE(ipre, jpre, score, bound);
		double scipj = SCORE(ipre, j, score, bound);
		double scijp = SCORE(i, jpre, score, bound);
		
		double premin = std::min(std::min(scipj, scijp), scipjp);
		
		if(premin == scipjp){
			i = ipre; j = jpre;
		}
		if(premin == scipj){
			i = ipre;
		}
		if(premin == scijp){
			j = jpre;
		}
		if(i == 0 && j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	
	std::reverse(alignment.begin(), alignment.end());

	return diff;
}

//--------- restricted cDTW (for nanopore analysis ONLY) ----------------//
double g::proc::BoundDynamicTimeWarpingR(
	const std::vector<double>& seq1, const std::vector<double>& seq2, 
	const std::vector<std::pair<long,long> >& bound, std::vector<std::pair<long,long> >& alignment)
{
	//--- init score ---//
	std::vector<std::vector<double> > score;
	score.resize(seq1.size());
	for(long i = 0; i < seq1.size(); i++){
		for(long j = bound[i].first; j <= bound[i].second; j++){
			score[i].push_back(std::fabs(seq1[i]-seq2[j]));
		}
	}
	
	for(long i = 1; i < seq1.size(); i++){
		if(SCORE(i, 0, score, bound) == DBL_MAX){
			break;
		}
		
		SCORE(i, 0, score, bound) += SCORE(i-1, 0, score, bound);
	}
	
	for(long j = 1; j <= bound[0].second; j++){
		score[0][j] += score[0][j-1];
	}
	
	for(long i = 1; i < seq1.size(); i++){
		for(long j = bound[i].first; j <= bound[i].second; j++){
			double acc = std::min(SCORE(i, j-1, score, bound), SCORE(i-1, j-1, score, bound));//std::min(std::min(SCORE(i-1, j, score, bound), SCORE(i, j-1, score, bound)), SCORE(i-1, j-1, score, bound));
			if(acc == DBL_MAX){
				score[i][j-bound[i].first] = acc;
			}
			else{
				score[i][j-bound[i].first] += acc;
			}
		}
	}
	
	double diff = SCORE(seq1.size()-1, seq2.size()-1, score, bound);
	
	
	long i = seq1.size()-1, j = seq2.size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		long ipre = i-1 < 0 ? 0 : i-1;
		long jpre = j-1 < 0 ? 0 : j-1;
		
		double scipjp = SCORE(ipre, jpre, score, bound);
// 		double scipj = SCORE(ipre, j, score, bound);
		double scijp = SCORE(i, jpre, score, bound);
		
		double premin = std::min(scijp, scipjp);//std::min(std::min(scipj, scijp), scipjp);
		
		if(premin == scipjp){
			i = ipre; j = jpre;
		}
// 		if(premin == scipj){
// 			i = ipre;
// 		}
		if(premin == scijp){
			j = jpre;
		}
		if(i == 0 || j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	
	std::reverse(alignment.begin(), alignment.end());
	
	return diff;
}
