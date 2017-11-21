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
#include "5mer/5mer_index.h"

bool Genomes2SignalSequence(const std::vector<char>& genomes, 
	std::vector<int>& index, std::vector<double>& signals, 
	int scale, int FIVE_or_SIX, int ZSCO_or_NOT)
{
	size_t bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		signals.assign(bound*scale,0);
		int cur=0;
		for(size_t i = 0; i < bound; i++){
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
			for(int c = scale; c--;){
				signals[cur]=sigval;
				cur++;
			}
		}
	}
	else
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-6;//genomes.size()%5;
		signals.assign(bound*scale,0);
		int cur=0;
		for(size_t i = 0; i < bound; i++){
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
			for(int c = scale; c--;){
				signals[cur]=sigval;
				cur++;
			}
		}
	}

	//---- tail five_mer ------//
//	for(size_t i = bound; i < genomes.size(); i++){
//		for(int c = scale; c--;){
//			signals.push_back(100);
//		}
//	}
}

//---------------------- main ---------------------//
int main(int argc, char **argv)
{
	std::string input="";
	std::string output="";

	struct options opts;
	opts.scale=1;
	opts.kmer=0;
	opts.zsco=0;
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

	Genomes2SignalSequence(genomes, index, signals, 
		opts.scale, opts.kmer,opts.zsco);

//printf("proc done \n");

	if(opts.zsco==1){
		g::io::WriteSignalSequence(opts.output, signals);
	}
	else{
		std::vector<int> signals_ (signals.begin(), signals.end());
		g::io::WriteSignalSequence_int(opts.output, signals_);
	}

//printf("write done \n");

	return 0;
}

