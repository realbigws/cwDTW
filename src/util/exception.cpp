#include "exception.h"

void ex::EX_THROW(const char* x)
{                                             
    std::ostringstream oss;                  
    oss << x;                                
    throw Exception(oss.str());              
}