#include "proc.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <memory.h>
#include <climits>
#include <cfloat>

using namespace std;

//================== Sheng added for Transfer alignment between Renmin_style and Sheng_style ==================//
//------- given Ali1 return AliPair ------//
void g::proc::Ali1_To_AliPair(int n1,int n2,vector <int> &ali1,
	vector<pair<int,int> > & alignment_out)
{
	//init
	alignment_out.clear();
	//start
	int i,j;
	int ii,jj;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
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
				alignment_out.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment_out.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment_out.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=n1;i++)alignment_out.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=n2;i++)alignment_out.push_back (pair<int,int>(-n1, i));  //Iy
}

//------- given Ali2 return AliPair ------//
void g::proc::Ali2_To_AliPair(int n1,int n2,vector <int> &ali2,
	vector<pair<int,int> > & alignment_out)
{
	//from ali1 to ali2
	vector <int> ali1(n1,-1);
	for(int i=0;i<n2;i++)
	{
		if(ali2[i]==-1)continue;
		int ii=ali2[i];
		ali1[ii]=i;
	}
	//proc
	Ali1_To_AliPair(n1,n2,ali1,alignment_out);
}

//------- given AliPair return Ali1 and Ali2 ------//
void g::proc::AliPair_To_Ali1_Ali2(int n1,int n2,
	vector<pair<int,int> > & alignment_in,
	vector <int> &ali1,vector <int> &ali2)
{
	//init
	ali1.resize(n1);
	ali2.resize(n2);
	for(int i=0;i<n1;i++)ali1[i]=-1;
	for(int i=0;i<n2;i++)ali2[i]=-1; 
	//proc
	for(int i=0;i<(int)alignment_in.size();i++)
	{
		int ii=alignment_in[i].first;
		int jj=alignment_in[i].second;
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
void g::proc::From_Align_Get_Bound(int moln1,int moln2,vector<pair<int,int> > &align,
	vector<pair<int,int> > &bound,int neib)
{
	int k;
	int ii,jj;
	int pre_ii,pre_jj;
	int size=(int)align.size();
	int first_=1;
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
			int pre_jj_=pre_jj;
			if(first_==1)
			{
				first_=0;
				int found=0;
				for(int k_=k+1;k_<size;k_++)
				{
					int ii_=align[k_].first;
					int jj_=align[k_].second;
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
	int i;
	int first=1;
	int wscurr,wsprev;
	int wscurr_start=1;
	int wsprev_start=1;
	int ww1,ww2;
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
void g::proc::Renmin_To_Sheng_align(int moln1,int moln2,
	vector<pair<int,int> > &align_in, vector<pair<int,int> > &align_out)
{
	vector <int> ali1(moln1,-1);
	vector <int> ali2(moln2,-1);
	for(int i=0;i<(int)align_in.size();i++)
	{
		int pos1=align_in[i].first;
		int pos2=align_in[i].second;
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
void g::proc::Sheng_To_Renmin_align(int moln1,int moln2,
	vector<pair<int,int> > &align_in, vector<pair<int,int> > &align_out)
{
	vector <int> ali1(moln1,-1);
	vector <int> ali2(moln2,-1);
	AliPair_To_Ali1_Ali2(moln1,moln2,align_in,ali1,ali2);
	//-> ali1 to align
	vector<pair<int,int> > align_tmp;
	Ali1_To_AliPair(moln1,moln2,ali1,align_tmp);
	//-> get start and end
	int start=-1;
	for(int i=0;i<(int)align_tmp.size();i++)
	{
		int pos1=align_tmp[i].first;
		int pos2=align_tmp[i].second;
		if(pos1>0 && pos2>0)
		{
			if(pos1==1 || pos2==1)continue;
			start=i;
			break;
		}	
	}
	int end=-1;
	for(int i=(int)align_tmp.size()-1;i>=0;i--)
	{
		int pos1=align_tmp[i].first;
		int pos2=align_tmp[i].second;
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
	align_out.push_back(pair<int,int>(0,0));
	for(int i=start;i<=end;i++)
	{
		int pos1=align_tmp[i].first-1;
		int pos2=align_tmp[i].second-1;
		align_out.push_back(pair<int,int>(pos1,pos2));
	}
	align_out.push_back(pair<int,int>(moln1-1,moln2-1));
}

//-------------- Sheng bound to Ren-min ---------------//
void g::proc::Sheng_To_Renmin_bound(int moln1,int moln2,
	vector<pair<int,int> > &bound_in, vector<pair<int,int> > &bound_out)
{
	bound_out.resize(moln1);
	for(int i=1;i<=moln1;i++)
	{
		int pos1=bound_in[i].first-1;
		int pos2=bound_in[i].second-1;
		if(pos1<0)pos1=0;
		if(pos2<0)pos2=0;
		bound_out[i-1].first=pos1;
		bound_out[i-1].second=pos2;
	}
}

//================================= sheng modify ======================//over


void g::proc::ZScoreNormalize(std::vector< double >& signals, double* avg, double* stdev)
{
	double sum = std::accumulate(signals.begin(), signals.end(), 0.0);
	double mean =  sum / signals.size();

	double acc = 0.0;
	for(size_t i = signals.size(); i--;){
		signals[i] = signals[i]-mean;
		acc += signals[i]*signals[i];
	}

	double deviation = std::sqrt(acc/signals.size());
	
	for(size_t i = signals.size(); i--;){
		signals[i] /= deviation;
	}
	
	if(avg){*avg = mean;}
	if(stdev){*stdev = deviation;}
}

void g::proc::RemoveHotSpot(std::vector< double >& signals, int thre)
{
	double mean =  std::accumulate(signals.begin(), signals.end(), 0.0) / signals.size();

	double acc = 0.0;
	for(size_t i = 0; i < signals.size(); i++){
		acc += signals[i]*signals[i];
	}

	double deviation = std::sqrt(acc/signals.size());

	double threshold = mean+thre*deviation;
	
	for(size_t i = 0; i < signals.size(); i++){
		if(signals[i] > threshold){
			signals[i] = threshold;
		}
		else if(signals[i] < -threshold){
			signals[i] = -threshold;
		}
	}
}

//used in WaveDenoise
static double GetThre(double* coeff, int size)  
{  
    double thr = 0.0;  
    double sigma = 0.0;  
	double copy[size];
	
    for(int i = 0; i < size; i++){  
        copy[i] = std::fabs(coeff[i]);
	}
  
    std::sort(copy, copy+size);  
  
    if (size % 2 == 0 && size >= 2){  
        sigma = (copy[size/2-1]+copy[size/2])/2/0.6745; 
	}
    else{
        sigma = copy[size/2]/0.6745;  
	}
  
    double N = size;  
    thr = sigma *sqrt(2.0*log(N));    

    return thr;  
}

// cutoff
static void WThresh(double* coeff, double thre, int size, bool soft)  
{   
	
    if (!soft){  //hard threshold
        for(int i = 0; i < size; i++){  
            if(std::fabs(coeff[i]) < thre){
                coeff[i] = 0.0;  
			}
        }  
    }  
    else{   //soft threshold
        for(int i = 0; i < size; i++){  
            if(std::fabs(coeff[i]) < thre){  
                coeff[i] = 0.0;  
            }  
            else{  
                if(coeff[i] < 0.0){  
                    coeff[i] = thre - std::fabs(coeff[i]);  
				}
                else{  
                    coeff[i] = std::fabs(coeff[i]) - thre;    
				}
            }  
        }  
	}  
} 

static double absmax(double *array, int N){
	double max;
	int i;

	max = 0.0;
	for (i = 0; i < N; ++i) {
		if (std::fabs(array[i]) >= max) {
			max = std::fabs(array[i]);
		}
	}

	return max;
}

/*
static void PyrWThre(wt_object wt, int lvl_thre){
	int J = wt->J;
	int t = wt->length[0];
	
	for(int i = 0; i < J; ++i){
		if(i+1 >= lvl_thre){
			int idx = t;
			int lvlen = wt->length[i+1];
			double thre = GetThre(wt->output+idx, lvlen);
			WThresh(wt->output+idx, thre, lvlen, true);  
		}
// 		printf("Level %d Access : output[%d] Length : %d \n", i + 1,t,wt->length[i+1]);
		t += wt->length[i+1];
	}
}
*/

/*
void g::proc::WaveDenoise(std::vector< double >& signals, bool soft)
{
	double* sigs = &signals[0];		//sst_nino3.dat

	size_t N = signals.size();
	int npyr = 10; 			// Total Number of scales

	wave_object obj = wave_init("db4");
	
	wt_object wt = wt_init(obj, "dwt", N, npyr);// Initialize the wavelet transform object
	setDWTExtension(wt, "sym");// Options are "per" and "sym". Symmetric is the default option
	setWTConv(wt, "direct");
	
	dwt(wt, sigs);// Perform DWT

	PyrWThre(wt, npyr-1);
	
	std::vector<double> output(wt->outlength);
	double* outs = &output[0];
	
	idwt(wt, outs); // Inverse Discrete Wavelet Packet Transform

// 	std::vector<double> diff;
// 	for(int i = 0; i < N; ++i){
// 		diff[i] = (sigs[i] - outs[i]);// /sigs[i];
// 	}
// 	wt_summary(wt); // Tree Summary
// 	printf("\n MAX %g \n", absmax(&diff[0], wt->siglength)); // If Reconstruction succeeded then the output should be a small value.

	wave_free(obj);
	wt_free(wt);
	
	signals = output;
}
*/

void g::proc::MedianFilter(std::vector<double>& signals, int width)
{
	double tmp[width];
	int width_2 = width*.5;
	std::vector<double> nsignl = signals;
	
	for(size_t i = width_2; i < signals.size()-width_2-1; i++){
		memcpy(tmp, &(signals[i-width_2]), sizeof(double)*width);
		std::sort(tmp, tmp+width);
		nsignl[i] = tmp[width_2];
	}
	
	signals = nsignl;
}

void g::proc::Diff(const std::vector<double>& raw, std::vector<double>& diff)
{
	diff.resize(raw.size());
	
	diff[0] = 0;
	
	for(int i = 1; i < raw.size(); i++){
		diff[i] = raw[i]-raw[i-1];
	}
}

void g::proc::LaplaceDiff(const std::vector<double>& raw, std::vector<double>& ldiff)
{
	ldiff.resize(raw.size());
	double tmplt[] = {-1,2,-1};
	int tsize = 3;
	int h_tsize = tsize/2;
	int tbegin = h_tsize;
	int tend = raw.size()-h_tsize;
	
	memset(&ldiff[0], 0, sizeof(double)*ldiff.size());
	
	for(int i = tbegin; i < tend; i++){
		for(int j = 0; j < tsize; j++){
			ldiff[i] += tmplt[j]*raw[j+i-h_tsize];
		}
	}
}

void g::proc::PeakPick(const std::vector<double>& raw, std::vector<std::pair<int, double> >& peaks)
{
	peaks.clear();
	std::vector<double> diff;
	Diff(raw, diff);

	peaks.push_back(std::make_pair(0, raw[0]));
	
	for(int i = 0; i < diff.size()-1; i++){
		if((diff[i] > 0 && diff[i+1] < 0) || (diff[i] < 0 && diff[i+1] > 0)){		//local max min
			peaks.push_back(std::make_pair(i, raw[i]));
		}
	}
	
	peaks.push_back(std::make_pair(raw.size()-1, raw[raw.size()-1]));
}

double g::proc::DynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment)
{	
	double* score[seq1.size()];
	double diff;
	
	for(int i = 0; i < seq1.size(); i++){
		score[i] = new double[seq2.size()];
	}
	
	for(int i = 0; i < seq1.size(); i++){
		for(int j = 0; j < seq2.size(); j++){
			score[i][j] = abs(seq1[i]-seq2[j]);
		}
	}
	
	for(int i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}
	
	for(int j = 1; j < seq2.size(); j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1.size(); i++){
		for(int j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(std::min(score[i-1][j], score[i][j-1]), score[i-1][j-1]);
		}
	}
	
	diff = score[seq1.size()-1][seq2.size()-1];
	
	int i = seq1.size()-1, j = seq2.size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		int ipre = i-1 < 0 ? 0 : i-1;
		int jpre = j-1 < 0 ? 0 : j-1;
		
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
	
	for(int i = 0; i < seq1.size(); i++){
		delete [] score[i];
	}
	
	return diff;
}

static inline double& SCORE(int i, int j, std::vector<double>* score, const std::vector<std::pair<int,int> >& bound){
	static double invalid = DBL_MAX;
	if(bound[i].first <= j && j <= bound[i].second){	
		return score[i][j-bound[i].first];
	}
	else{
		return invalid;
	}
}

double g::proc::BoundDynamicTimeWarping(const std::vector<double>& sequence1, const std::vector<double>& sequence2, 
	const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment)
{
	/*check order*/
	bool firstorder = true;
	const std::vector<double>* seq1, * seq2;
	
//	if(sequence1.size() > sequence2.size()){
//		firstorder = false;
//	}
	
	if(!firstorder){		//genome first
		seq1 = &sequence2;
		seq2 = &sequence1;
	}
	else{
		seq1 = &sequence1;
		seq2 = &sequence2;
	}
	
	std::vector<double> score[seq1->size()];
	double diff;
	
	for(int i = 0; i < seq1->size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			score[i].push_back(std::fabs((*seq1)[i]-(*seq2)[j]));
		}
	}
	
	for(int i = 1; i < seq1->size(); i++){
		if(SCORE(i, 0, score, bound) == DBL_MAX){
			break;
		}
		
		SCORE(i, 0, score, bound) += SCORE(i-1, 0, score, bound);
	}
	
	for(int j = 1; j <= bound[0].second; j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1->size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			double acc = std::min(std::min(SCORE(i-1, j, score, bound), SCORE(i, j-1, score, bound)), SCORE(i-1, j-1, score, bound));
			if(acc == DBL_MAX){
				score[i][j-bound[i].first] = acc;
			}
			else{
				score[i][j-bound[i].first] += acc;
			}
		}
	}
	
	diff = SCORE(seq1->size()-1, seq2->size()-1, score, bound);
	
	
	int i = seq1->size()-1, j = seq2->size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		int ipre = i-1 < 0 ? 0 : i-1;
		int jpre = j-1 < 0 ? 0 : j-1;
		
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
	
	if(!firstorder){
		for(int i = alignment.size(); i--;){
			std::swap(alignment[i].first, alignment[i].second);
		}
	}
	
	return diff;
}

double g::proc::BoundDynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment)
{	
	std::vector<double> score[seq1.size()];
	double diff;
	
	for(int i = 0; i < seq1.size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			score[i].push_back(std::fabs(seq1[i]-seq2[j]));
		}
	}
	
	for(int i = 1; i < seq1.size(); i++){
		if(SCORE(i, 0, score, bound) == DBL_MAX){
			break;
		}
		
		SCORE(i, 0, score, bound) += SCORE(i-1, 0, score, bound);
	}
	
	for(int j = 1; j <= bound[0].second; j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1.size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			double acc = std::min(SCORE(i, j-1, score, bound), SCORE(i-1, j-1, score, bound));//std::min(std::min(SCORE(i-1, j, score, bound), SCORE(i, j-1, score, bound)), SCORE(i-1, j-1, score, bound));
			if(acc == DBL_MAX){
				score[i][j-bound[i].first] = acc;
			}
			else{
				score[i][j-bound[i].first] += acc;
			}
		}
	}
	
	diff = SCORE(seq1.size()-1, seq2.size()-1, score, bound);
	
	
	int i = seq1.size()-1, j = seq2.size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		int ipre = i-1 < 0 ? 0 : i-1;
		int jpre = j-1 < 0 ? 0 : j-1;
		
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
