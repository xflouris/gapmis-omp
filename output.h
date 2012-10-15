/**
    GapMis: a tool for pairwise sequence alignment with a single gap.
    Copyright (C) 2011 Solon P. Pissis, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifndef OUTPUT_H
#define OUTPUT_H

#include "types.h"

#define LINE_LNG 50

void usage ( void );

unsigned int results_many_to_many_scr ( const char * filename, struct TSeq * p, unsigned int cur_alloc_query, struct TSeq * t, unsigned int cur_alloc_target, struct Tout* out );

unsigned int results_many_to_many ( const char * filename, struct TSeq * p, unsigned int cur_alloc_query, struct TSeq * t, unsigned int cur_alloc_target, const struct Tin* in, struct Tout* out );

unsigned int print_alignment ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, char* seq_gap, char* mark_mis, const struct Tin* in, struct Tout* out );

void print_header ( FILE * out, const char * filename, const struct Tin* in );

void wrap ( const char * s1, const char * s1_header, const char * s2, const char * s2_header, const char * diff, int len, FILE * output );

void print_line ( const char * s, int start, int stop, int * nr_gaps, int end, int diff, FILE * output, const char * header );

#endif
