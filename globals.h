#pragma once

#include<string>
#include<vector>
#include <cstdint>
#include <immintrin.h>  //Header file for SSE instruction set.

//default value
extern int DBA;
extern int DBA_2;
extern int step;
extern int lengthOfAlignmentUnit;
extern int window;
extern double error_rate;

extern long long length_sum;
extern bool cigar;