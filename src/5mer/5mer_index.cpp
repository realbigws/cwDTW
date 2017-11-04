#include "5mer_index.h"
#include <iostream>

#include "5mer_index_table.rc"

/* A = 0; C = 1; G = 2; T = 3 */
int g::Mer2Signal::FiveMer2Index(char const* fivemer)
{
	const short atag = 0, ctag = 1, gtag = 2, ttag = 3;
	short idx = 0;
	
	for(int i = 0; i < 5; i++){
		short tag;
		switch(fivemer[i]){
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
			fprintf(stderr,"BAD CODE HERE !! %c \n",fivemer[i]);
			exit(-1);
		}
		tag <<= (4-i)*2;
		idx |= tag;
	}
	
	return idx;
}

int g::Mer2Signal::FiveMer2Index(char g0, char g1, char g2, char g3, char g4)
{
	char tmp[5] = {g0, g1, g2, g3, g4};
	return FiveMer2Index(tmp);
}

int g::Mer2Signal::FiveMer2Index(const std::vector<char>& fivemer)
{
	char tmp[5] = {fivemer[0], fivemer[1], fivemer[2], fivemer[3],  fivemer[4]};
	return FiveMer2Index(tmp);
}

int g::Mer2Signal::FiveMer2Index(const std::string& fivemer)
{
	return FiveMer2Index(fivemer.c_str());
}

double g::Mer2Signal::AvgSignalAt(int index)
{
	return index_table[index][0];
}
