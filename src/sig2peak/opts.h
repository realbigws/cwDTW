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
    char output[65532];
    double scale0;
    int ZorNOT;
};

inline int GetOpts(int argc, char **argv, options* opts_){

    static struct option longopts[] = {
        { "help",            no_argument,            NULL,              'h' },
        { "input",           required_argument,      NULL,              'i' },
        { "output",          required_argument,      NULL,              'o' },
        { "scale",           required_argument,      NULL,              's' },
        { "ZorNOT",          no_argument,            NULL,              'z' },
        { NULL,              0,                      NULL,               0  }
    };
	
    if((argc != 7 && argc != 9 ) && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1){
		EX_TRACE("[-i SIGNAL INPUT][-o CWT OUTPUT][-s SCALE]([-z ZNORM_or_NOT=0])\n");
		EX_TRACE("[note]: if -z is set to 0, then NO Znormalize will be performd \n");
		return -1;
    }
    
    int ch;
    while((ch = getopt_long(argc, argv, "hi:o:s:z:", longopts, NULL))!= -1){
        switch (ch) {

        case '?':
        {
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;
        }

        case ':':
        {
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;
        }

        case 'h':
        {
            EX_TRACE("[-i GENOME INPUT][-o SIGNAL OUTPUT]([-z ZNORM_or_NOT=0])\n");
            EX_TRACE("[note]: if -z is set to 0, then NO Znormalize will be performd \n");
            return -1;
        }

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
        
        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 's':
        {
            std::istringstream iss(optarg);
            iss >> opts_->scale0;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'z':
        {
            std::istringstream iss(optarg);
            iss >> opts_->ZorNOT;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
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
