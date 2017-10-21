#ifndef EXCEPTION_H__
#define EXCEPTION_H__

#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <stdarg.h>
#include <cstring>
#include <cstdio>
#include <ctime>

#define DO_TRACE_
#define DO_NOTE_
#define DO_CLOCK_

#define EX_PRINT( MSG,args...)  printf(MSG, ##args);  fflush(stdout);
#define EX_ERROR( MSG,args...)  printf(MSG, ##args);  fflush(stderr);

#ifdef DO_TRACE_
#define EX_TRACE( MSG,args...)  printf(MSG, ##args);  fflush(stdout);
#define EX_CUR_WORKSPACE()	printf("%s", get_current_dir_name()); fflush(stdout);
#else 
#define EX_TRACE( MSG,args...)
#define EX_CUR_WORKSPACE()
#endif

#ifdef DO_NOTE_
#define EX_NOTE(X)			X
#else 
#define EX_NOTE(X)
#endif

#ifdef DO_CLOCK_
#define _DASH	"------------------------------------------------------------------------------\n"
#define _WAVE	"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
static const char* wday12310[] = {"Sun","Mon","Tue","Wed","Thu"," Fri","Sat"};
#define EX_TIME_BEGIN( MSG,args...) {time_t time_begin0112;clock_t clock_begin2290;{printf(MSG, ##args); \
	    time(&time_begin0112);clock_begin2290 = clock(); struct tm *p; p=localtime(&time_begin0112);\
	    printf("[Beginning] %d/%d/%d ", (1900+p->tm_year),(1+p->tm_mon), p->tm_mday); \
	    printf("%s %d:%d:%d\n", wday12310[p->tm_wday],p->tm_hour, p->tm_min, p->tm_sec);}  fflush(stdout);	    
#define EX_TIME_END( MSG,args...) time_t time_end8013;clock_t clock_end2520;{printf(MSG, ##args); \
	    time(&time_end8013);clock_end2520 = clock(); struct tm *p;p=localtime(&time_end8013);	\
	    printf("[Finished] %d/%d/%d ", (1900+p->tm_year),(1+p->tm_mon), p->tm_mday); \
	    printf("%s %d:%d:%d (time elapse: %ldms)\n", wday12310[p->tm_wday],p->tm_hour, p->tm_min, \
	    p->tm_sec, (clock_end2520-clock_begin2290)*1000/CLOCKS_PER_SEC);}}    fflush(stdout);
#define EX_BEGIN_CLOCK()	clock_t clock_begin2290;{clock_begin2290 = clock();}  
#define EX_END_CLOCK()		clock_t clock_end2520;{clock_end2520 = clock();}
#define EX_ELAPSE()		(clock_end2520-clock_begin2290)*1000/CLOCKS_PER_SEC  
#else
#define EX_TIME_BEGIN( MSG,args...)
#define EX_TIME_END( MSG,args...)
#define EX_BEGIN_CLOCK()
#define EX_END_CLOCK()
#define EX_ELAPSE()			0
#endif

namespace ex{
    
class Exception
{
private:
    std::string msg;
public:
    Exception(std::string _msg):msg(_msg){}
    const char* Msg(){ return msg.c_str();}
};

void EX_THROW(const char* x);

//ex::Puts("", MSG, ##args, "\0");
// void Puts(const char* msg, ...)
// {
//     va_list args;
//     char* str;
//     va_start(args, msg);
//     
//     while(1){
// 	str = va_arg(args, char*);
// 	if(strcmp(str, "\0")==0){
// 	    break;
// 	}
// 	std::cout<<str;
//     }
//     std::cout<<std::endl;
//     std::cout.flush();
//     va_end(args);
//     return;
// }

}

#endif