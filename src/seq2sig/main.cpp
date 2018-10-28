#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <proc/io.h>
#include "util/exception.h"
#include "kmer/kmer_index.h"

//---------- genome to signal (DNA) -------------//
bool Genomes2SignalSequence(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals,
	std::vector<std::string>& kmer_rec, int scale, 
	int FIVE_or_SIX, int ZSCO_or_NOT, int WANT_TAIL)
{
/*
	//------- simple version -----------//
	{
		long bound=genomes.size();
		index.assign(bound,0);
		signals.assign(bound,0);
		for(long i=0;i<bound;i++)
		{
			char c=genomes[i];
			double val=0;
			if(c=='A')val=1;
			if(c=='T')val=2;
			if(c=='C')val=-1;
			if(c=='G')val=-2;
			signals[i]=val;
		}
		return true;
	}
*/

	long bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		signals.assign(bound*scale,0);
		kmer_rec.resize(bound*scale);
		long cur=0;
		for(long i = 0; i < bound; i++){
			//-> get kmer
			std::string kmer(genomes.begin()+i,genomes.begin()+i+5);
			//-> get signal
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::AvgSignalAt_5mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208351)/12.832660;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			//-> assign
			for(int c = scale; c--;){
				signals[cur]=sigval;
				kmer_rec[cur]=kmer;
				cur++;
			}
		}
	}
	else               //-> 6mer model
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		signals.assign(bound*scale,0);
		kmer_rec.resize(bound*scale);
		long cur=0;
		for(long i = 0; i < bound; i++){
			//-> get kmer
			std::string kmer(genomes.begin()+i,genomes.begin()+i+6);
			//-> get signal
			double sigval;
			if(index[i]<0)sigval=0;
			else{
			 	sigval = g::Mer2Signal::AvgSignalAt_6mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208199)/12.868652;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			//-> assign
			for(int c = scale; c--;){
				signals[cur]=sigval;
				kmer_rec[cur]=kmer;
				cur++;
			}
		}
	}


	//---- tail k_mer ------//
	if(WANT_TAIL==1)
	{
		for(long i = bound; i < genomes.size(); i++){
			for(int c = scale; c--;){
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					signals.push_back(0);
				}
				else
				{
					signals.push_back(5.7*100+14);
				}
			}
		}
	}

}

//---------- genome to signal (RNA) -------------//
bool Genomes2SignalSequence_RNA(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals, 
	std::vector<std::string>& kmer_rec, int scale, 
	int mv200_or_mv180, int ZSCO_or_NOT, int WANT_TAIL)
{
/*
	//------- simple version -----------//
	{
		long bound=genomes.size();
		index.assign(bound,0);
		signals.assign(bound,0);
		for(long i=0;i<bound;i++)
		{
			char c=genomes[i];
			double val=0;
			if(c=='A')val=1;
			if(c=='T')val=2;
			if(c=='U')val=2;  //-> for RNA only
			if(c=='C')val=-1;
			if(c=='G')val=-2;
			signals[i]=val;
		}
		return true;
	}
*/

	long bound;
	if(mv200_or_mv180==1) //-> 200mv model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		signals.assign(bound*scale,0);
		kmer_rec.resize(bound*scale);
		long cur=0;
		for(long i = 0; i < bound; i++){
			//-> get kmer
			std::string kmer(genomes.begin()+i,genomes.begin()+i+5);
			//-> get signal
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::RnaSignalAt_5mer_200mv(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-87.240242)/14.171575;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			//-> assign
			for(int c = scale; c--;){
				signals[cur]=sigval;
				kmer_rec[cur]=kmer;
				cur++;
			}
		}
	}
	else                  //-> 180mv model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		signals.assign(bound*scale,0);
		long cur=0;
		for(long i = 0; i < bound; i++){
			//-> get kmer
			std::string kmer(genomes.begin()+i,genomes.begin()+i+5);
			//-> get signal
			double sigval;
			if(index[i]<0)sigval=0;
			else{
			 	sigval = g::Mer2Signal::RnaSignalAt_5mer_180mv(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-94.524117)/15.939267;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			//-> assign
			for(int c = scale; c--;){
				signals[cur]=sigval;
				kmer_rec[cur]=kmer;
				cur++;
			}
		}
	}


	//---- tail k_mer ------//
	if(WANT_TAIL==1)
	{
		for(long i = bound; i < genomes.size(); i++){
			for(int c = scale; c--;){
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					signals.push_back(0);
				}
				else
				{
					signals.push_back(5.7*100+14);
				}
			}
		}
	}

}


//---------------------- main ---------------------//
int main(int argc, char **argv)
{
	std::string input="";
	std::string output="";

	struct options opts;
	opts.scale=1;
	opts.kmer=1;
	opts.zsco=0;
	opts.rna=0;
	opts.name=0;
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	input=opts.input;
	output=opts.output;
	if(input=="" || output=="")
	{
		fprintf(stderr,"input or output is NULL \n");
		return -1;
	}

	std::vector<char> genomes;
	std::vector<int> index;
	std::vector<double> signals;
	g::io::ReadATCG(opts.input, genomes);

//printf("read done \n");

	std::vector<std::string> kmer_rec;
	if(opts.rna==0)   //-> DNA pore model
	{
		Genomes2SignalSequence(genomes, index, signals, kmer_rec,
			opts.scale, opts.kmer,opts.zsco,0);
	}
	else              //-> RNA pore model
	{
		Genomes2SignalSequence_RNA(genomes, index, signals, kmer_rec,
			opts.scale, opts.rna,opts.zsco,0);
	}

//printf("proc done \n");

	if(opts.zsco==1)
	{
		if(opts.name==1)
		{
			g::io::WriteSignalSequence_withName(opts.output, signals, kmer_rec);
		}
		else
		{
			g::io::WriteSignalSequence(opts.output, signals);
		}
	}
	else
	{
		std::vector<int> signals_ (signals.begin(), signals.end());
		if(opts.name==1)
		{
			g::io::WriteSignalSequence_int_withName(opts.output, signals_, kmer_rec);
		}
		else
		{
			g::io::WriteSignalSequence_int(opts.output, signals_);
		}
	}

//printf("write done \n");

	return 0;
}

