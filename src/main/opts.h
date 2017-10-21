#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C" {
#include <getopt.h>
}
#include "util/exception.h"

struct options {
	char input[65532];
	char peer[65532];
	char output[65532];
	int radius;
	int level;
	float scale0;
	int mode;
	int dp_mode;
	int verbose;
	int test;
};

inline int GetOpts(int argc, char **argv, options* opts_) {

    static struct option longopts[] = {
	{ "help",            no_argument,            NULL,              'h' },
	{ "input",     required_argument,      NULL,    			  'i' },
	{ "peer",    	     required_argument,      NULL,      	'p' },
	{ "output",      	required_argument,      NULL,              'o' },
	{ "radius",      required_argument,      NULL,              'r' },
	{ "level",      required_argument,      NULL,              'l' },
	{ "scale",      required_argument,      NULL,              's' },
	{ "mode",       required_argument,      NULL,              'm' },
	{ "dp_mode",    required_argument,      NULL,              'M' },
	{ "verbose",         no_argument,            NULL,              'v' },
	{ "test",            no_argument,            NULL,              't' },
	{ NULL,              0,                      NULL,               0 }
    };

    if((argc !=5 && argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15 && argc != 17 && argc != 19 && argc!= 21) && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1) {
	EX_TRACE("WaveletDTW \n");
	EX_TRACE("required:\n");
        EX_TRACE("[-i SHORT SIGNAL][-p LONG SIGNAL][-o OUTPUT] \n");
	EX_TRACE("optional:\n");
	EX_TRACE("([-r NEIGHBOUR RADIUS])([-l LEVEL])([-s SCALE])([-m mode])([-M dp_mode])([-v verbose])([-t test])\n");
	EX_TRACE("[note]: by default, r=50, l=2, s=sqrt(2), m=0, M=0, v=0, t=0 \n");
        return -1;
    }

    int ch;
    while((ch = getopt_long(argc, argv, "hi:p:o:r:l:s:m:M:v:t:", longopts, NULL))!= -1) {
        switch(ch) {

        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;

        case 'h':
            EX_TRACE("WaveletDTW \n"
                     "required:\n"
                     "[-i SHORT SIGNAL][-p LONG SIGNAL][-o OUTPUT]\n"
                     "optional:\n"
                     "([-r NEIGHBOUR RADIUS])([-l LEVEL])([-s SCALE])([-m mode])([-M dp_mode])\n"
                     "([-v verbose])([-t test])\n"
                     "-------------------------------------------------------------\n"
                     "SHORT SIGNALE: (reference) sequence signal i.e. ATCG...;\n"
                     "LONG SIGNAL: (nanopore) signal that will be transformed to align with reference;\n"
                     "OUTPUT: signal alignment result;\n"
                     "NEIGHBOUR RADIUS: warp search radius (default 50);\n"
                     "LEVEL: sampling level in continous wavelet (default 2);\n"
                     "SCALE: base scale in continous wavelet (default sqrt(2));\n"
                     "mode: cDTW radius mode. 0 for 'set' and 1 for 'adapt' (default 0)\n"
                     "dp_mode: final cDTW mode. 0 for 'normal' and 1 for 'restrict' (default 0) \n"
                     "verbose: 0 for NO screenout message, and 1 for screenout (default 0);\n"
                     "test: test mode. 0 for not_use; 1 for equal_ave; 2 for peak_ave (default 0) \n");
            return -1;

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'p':
        {
            std::istringstream iss(optarg);
            iss >> opts_->peer;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
	case 'r':
        {
            std::istringstream iss(optarg);
            iss >> opts_->radius;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
	case 'l':
        {
            std::istringstream iss(optarg);
            iss >> opts_->level;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
	case 's':
        {
            std::istringstream iss(optarg);
            iss >> opts_->scale0;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

	case 'm':
	{
            std::istringstream iss(optarg);
            iss >> opts_->mode;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
	}
	break;

	case 'M':
	{
            std::istringstream iss(optarg);
            iss >> opts_->dp_mode;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
	}
	break;

	case 'v':
	{
             std::istringstream iss(optarg);
             iss >> opts_->verbose;
             if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
	}
	break;

	case 't':
	{
             std::istringstream iss(optarg);
             iss >> opts_->test;
             if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
	}
	break;

        case 0:
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif

