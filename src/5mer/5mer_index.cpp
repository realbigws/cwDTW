#include <iostream>
#include "5mer_index.h"
#include "5mer_index_table.rc"
#include "5mer_index_official.rc"
#include "6mer_index_official.rc"


//---------- Nucleotide_to_Int ----------//
/* A = 0; C = 1; G = 2; T = 3 */
int g::Mer2Signal::Nucleotide_to_Int(char c,int &mark)
{
	int tag;
	mark=0;
	switch(c)
	{
		case 'A':
			tag = 0;
			break;
		case 'C':
			tag = 1;
			break;
		case 'G':
			tag = 2;
			break;
		case 'T':
			tag = 3;
			break;
		default:
			tag = 0;  //-> set to A
			mark = 1; //-> mark this position
//		default:
//			fprintf(stderr,"BAD CODE HERE !! %c \n",c);
//			exit(-1);
	}
	return tag;
}


//------------ given genome, transfer to 5-mer index ---------------//
void g::Mer2Signal::Genome2Index_5mer(const std::vector <char> &input, std::vector <int> &index)
{
	int div=256;
	long size=input.size();
	int mark=0;
	index.assign(size,-1);

	short idx = 0;
	for(int i = 0; i < 5; i++)
	{
		short tag=Nucleotide_to_Int(input[i],mark);
		tag <<= (4-i)*2;
		idx |= tag;
	}
	//-> assign mark
	if(mark){
		index[0]=-1;
	}
	else {
		index[0]=idx;
	}


	for(long i=5;i<size;i++)
	{
		idx=(idx%div)*4;
		short tag=Nucleotide_to_Int(input[i],mark);
		idx += tag;
		//-> assign mark
		if(mark){
			index[i-4]=-1;
		}
		else {
			index[i-4]=idx;
		}
	}
}

//------------ given genome, transfer to 6-mer index ---------------//
void g::Mer2Signal::Genome2Index_6mer(const std::vector <char> &input, std::vector <int> &index)
{
	int div=1024;
	long size=input.size();
	int mark=0;
	index.assign(size,-1);
	
	short idx = 0;
	for(int i = 0; i < 6; i++)
	{
		short tag=Nucleotide_to_Int(input[i],mark);
		tag <<= (5-i)*2;
		idx |= tag;
	}
	//-> assign mark
	if(mark){
		index[0]=-1;
	}
	else {
		index[0]=idx;
	}
	
	
	for(long i=6;i<size;i++)
	{
		idx=(idx%div)*4;
		short tag=Nucleotide_to_Int(input[i],mark);
		idx += tag;
		//-> assign mark
		if(mark){
			index[i-4]=-1;
		}
		else {
			index[i-5]=idx;
		}
	}
}

//----------------- FiveMer2Index ---------------------//
//-> self_model is based on 5mer pore model
double g::Mer2Signal::AvgSignalAt_SelfModel(int index)
{
	if(index<0)return 0;
	return index_table_self[index][0];
}

//-> 5mer case
double g::Mer2Signal::AvgSignalAt_5mer(int index)
{
	if(index<0)return 0;
	return index_table_5mer[index];
}

//-> 6mer case
double g::Mer2Signal::AvgSignalAt_6mer(int index)
{
	if(index<0)return 0;
	return index_table_6mer[index];
}

