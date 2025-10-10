#ifndef __utilities_h__
#define __utilities_h__

#include <algorithm>
#include <bitset>
#include <cblas.h>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <lapack.h>
#include <limits.h>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <tuple>
#include <type_traits>
#include <unistd.h>
#include <vector>

///Types of the states;
//For 4 sites or less
//typedef uint8_t sType;

//For 8 sites or less
//typedef uint16_t sType;

//For 16 sites or less
//typedef uint32_t sType;

//For more than 16 sites
typedef uint64_t sType;		

//Type of the vector 
typedef double vType;



typedef uint8_t uChar;
typedef uint16_t uShort;
typedef uint32_t uInt;
typedef uint64_t uLong;

//Declaration of external variables
extern int verbose;

#endif
