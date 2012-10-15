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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "types.h"
#include "output.h"

/* Usage of the tool */
void usage ( void )
 {
   fprintf ( stdout, "  ____             __  __ _\n" );
   fprintf ( stdout, " / ___| __ _ _ __ |  \\/  (_)___\n" );
   fprintf ( stdout, "| |  _ / _` | '_ \\| |\\/| | / __|\n" );
   fprintf ( stdout, "| |_| | (_| | |_) | |  | | \\__ \\\n" );
   fprintf ( stdout, " \\____|\\__,_| .__/|_|  |_|_|___/\n" );
   fprintf ( stdout, "            |_|\n\n" );

   fprintf ( stdout, "Usage: gapmis <options>\n" );
   fprintf ( stdout, "Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --sequences-a          <str>    Query sequences of length m filename.\n" );
   fprintf ( stdout, "  -b, --sequences-b          <str>    Target sequences of length n>=m filename.\n" );
   fprintf ( stdout, "Optional:\n" );
   fprintf ( stdout, "  -g, --gap-open-penalty    <float>   The gap open penalty is the score taken\n"
                     "                                      away when a gap is created. The best\n"
                     "                                      value depends on the choice of comparison\n"
                     "                                      matrix.   The   default   value   assumes\n"
                     "                                      you  are  using  the  EBLOSUM62  matrix\n"
                     "                                      protein  sequences,  and  the  EDNAFULL\n"
                     "                                      matrix for nucleotide sequences. Floating\n"
                     "                                      point number from 1.0 to 100.0. (default:\n"
                     "                                      10.0)\n" );


   fprintf ( stdout, "  -e, --gap-extend-penalty  <float>   The gap extension penalty is added to\n"
                     "                                      the standard gap penalty for each base or\n"
                     "                                      residue in the gap. This is how long gaps\n"
                     "                                      are penalized. Floating point number from\n"
                     "                                      0.0  to  10.0.  (default:  0.5)\n" );

   fprintf ( stdout, "  -o, --output-file         <str>     Output   alignment   filename   (default:\n"
                     "                                      gapmis.out)\n" );
   fprintf ( stdout, "  -f, --output-format       <int>     Output alignment format. 0 for outputting\n"
                     "                                      only  the  scores  and  1  for EMBOSS-like\n"
                     "                                      output. (default:  0)\n" );

   fprintf ( stdout, "  -d, --data-file           <str>     This  is  the  scoring  matrix  used  when\n"
                     "                                      comparing  sequences.  It  can  be  either\n"
                     "                                      `EBLOSUM62'  (for  protein  sequences)  or\n" 
                     "                                      `EDNAFULL'   (for  nucleotide  sequences).\n"  
                     "                                      (default: EDNAFULL)\n" );
   fprintf ( stdout, "  -t, --threads             <int>     Number of threads to be used  (default: 1)\n");
   fprintf ( stdout, "  -m, --max-gap             <int>     Limit the maximum gap size to this value\n"
                     "                                      (default: length of the longest sequence\n"
                     "                                      minus  1)\n\n" );

 }

/* Prints the header in the output file */
void print_header ( FILE * out, const char * filename, const struct Tin* in )
 {
   time_t               t;
   time ( &t );

   fprintf ( out, "####################################\n" );
   fprintf ( out, "# Program: GapMis\n" );
   fprintf ( out, "# Rundate: %s", ctime ( &t ) );
   fprintf ( out, "# Report file: %s\n", filename );
   fprintf ( out, "# Matrix: %s\n", ( in -> scoring_matrix ? "BLOSUM62" : "EDNAFULL" ) );
   fprintf ( out, "# Gap penalty: %.3f\n", in -> gap_open_pen );
   fprintf ( out, "# Extend penalty: %.3f\n", in -> gap_extend_pen );
   fprintf ( out, "####################################\n\n" );
 }

/* Creates the output file with the many to many scores */
unsigned int results_many_to_many_scr ( const char * filename, struct TSeq * p, unsigned int cur_alloc_query, struct TSeq * t, unsigned int cur_alloc_target, struct Tout* out )
 {

   FILE          * output;
   unsigned int t_i;
   unsigned int p_i;
   unsigned int o_i;
 
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      return ( 0 );
    }

   for ( t_i = 0, o_i = 0; t_i < cur_alloc_target; ++ t_i )
     {
       for ( p_i = 0; p_i < cur_alloc_query; ++ p_i, ++ o_i )
         {
           fprintf ( output, "%lf\n", out[o_i] . max_score );
         }
     }

  if ( fclose ( output ) ) 
   {
     return ( 0 );
   }
     
   return ( 1 );	
 }

/* Creates the output file with the many to many alignments (verbose) */
unsigned int results_many_to_many ( const char * filename, struct TSeq * p, unsigned int cur_alloc_query, struct TSeq * t, unsigned int cur_alloc_target, const struct Tin* in, struct Tout* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
   unsigned int t_i;
   unsigned int p_i;
   unsigned int o_i;
 
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      return ( 0 );
    }

   n = strlen ( t[0] . data );
   m = strlen ( p[0] . data );

   /* Dynamic memory allocation for seq_gap */
   if ( ! ( seq_gap = ( char * ) calloc ( n + 1 + ( in -> max_gap ) + 1, sizeof ( char ) ) ) )
     {
       return ( 0 );
     } 
   
   if ( ! ( mark_mis = ( char* ) calloc ( n + ( in -> max_gap ) + 1, sizeof( char ) ) ) )
     {
       return ( 0 );
     } 

   for ( t_i = 0, o_i = 0; t_i < cur_alloc_target; ++ t_i )
     {
       for ( p_i = 0; p_i < cur_alloc_query; ++ p_i, ++ o_i )
         {
           seq_gap[0]='\0';
           mark_mis[0]='\0';
           print_header ( output, filename, in );
   
           if ( out[o_i] . where == 1 ) //gap is in the text
            {
              print_alignment ( t[t_i] . data, n, p[p_i] . data, m, seq_gap, mark_mis, in, &out[o_i] );
              wrap ( p[p_i] . data, p[p_i] . header, seq_gap, t[t_i] . header, mark_mis, LINE_LNG, output ); 
            }
           else                     //gap is in the pattern
            {
              print_alignment ( p[p_i] . data, m, t[t_i] . data, n, seq_gap, mark_mis, in, &out[o_i] );
              wrap ( seq_gap, p[p_i] . header, t[t_i] . data, t[t_i] . header, mark_mis, LINE_LNG, output ); 
            }
   
           fprintf ( output, "\n" );
           fprintf ( output, "Alignment score: %lf\n", out[o_i] . max_score );
           fprintf ( output, "Number of mismatches: %d\n", out[o_i] . num_mis );
           fprintf ( output, "Length of gap: %d\n", out[o_i] . min_gap );
   
           if( out[o_i] . min_gap > 0 )
            {
              fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out[o_i] . where == 1 ) ? t[t_i] . header : p[p_i] . header, out[o_i] . gap_pos );
            } 
   
           fprintf ( output, "\n\n" );
         }
     }
  free ( mark_mis );
  free ( seq_gap ); 

  if ( fclose ( output ) ) 
   {
     return ( 0 );
   }
     
   return ( 1 );	
 }

/*
Creates seq_gap and mark_mis, and computes min_mis
*/
unsigned int print_alignment ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, char* seq_gap, char* mark_mis, const struct Tin* in, struct Tout* out )
 {
   unsigned int i, j;

   if ( out -> min_gap > 0 )
    {

      for ( i = 0; i < out -> gap_pos; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )	
          {
            mark_mis[i] = '.';
            out -> num_mis = out -> num_mis + 1;
          }
         else				
           mark_mis[i] = '|';
       }

      for ( j = 0; j < out -> min_gap; ++ j )
       {
         seq_gap[ j + i ] = '-'; 
         mark_mis[ j + i ] = ' ';
       }

      for ( ; i < seqb_len - out -> min_gap && i < seqa_len ; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         if ( seqa[i] != seqb[i + out -> min_gap] )	
          {
            mark_mis[ j + i ] = '.';
            out -> num_mis = out -> num_mis + 1;
          }
         else
           mark_mis[j + i] = '|';
       }
      
      for ( ; i < seqa_len; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         mark_mis[ j + i ] = '|';
       }
    }
   else
    {
      for ( i = 0; i < seqa_len; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )
          {
            mark_mis[i] = '.';
            out -> num_mis = out -> num_mis + 1;
          }
         else			
           mark_mis[i] = '|';
       }
    }	
   
   return ( 1 );
 }

void print_line ( const char * s, int start, int stop, int * nr_gaps, int end, int diff, FILE * output, const char * header )
 {
   int                  k;

   if ( start == stop ) return;

   if ( diff )
    {
      fprintf ( output, "%25s", "" );
    }
   else
    {
     if ( header )
       fprintf ( output, "%-13.13s %10d ", header, start + 1 - *nr_gaps );
     else
       fprintf ( output, "%-13.13s %10d ", "", start + 1 - *nr_gaps );
    }

   for ( ; start < stop; ++ start )
    {
      fputc ( s[start], output );
      if ( s[start] == '-' && ! diff ) ++ ( *nr_gaps );
    }

   if ( stop != end )
    {
      for ( k = stop; k < end; ++ k )
       {
         fputc ( ' ', output );
       }
    }
   if ( ! diff )  fprintf ( output, " %-10d", start - *nr_gaps );
   fprintf ( output, "\n" );
 }

/* Wrap two sequences s1 and s2 including the differences (diff) so that the line width is at most len */
void wrap ( const char * s1, const char * s1_header, const char * s2, const char * s2_header, const char * diff, int len, FILE * output )
 {
   int                  m, n, i, j;
   int                  nr_gaps_a;
   int                  nr_gaps_b;
   int                  nr_lines;

   if ( ! len ) 
     return;

   m = strlen ( s1 );
   n = strlen ( s2 );

   if ( ! n && ! m ) 
     return;

   i         = 0;
   j         = 0;
   nr_gaps_a = 0;
   nr_gaps_b = 0;

   //nr_lines = m / len;
   nr_lines = ( n > m ? m : n ) / len;
   for ( i = 0; i < nr_lines; ++ i )
    {
      /* Sequence s1 */
      print_line ( s1 , i * len, ( i + 1 ) * len, &nr_gaps_a, ( i + 1 ) * len, 0, output, s1_header );

      /* Difference */
      print_line ( diff, i * len, ( i + 1 ) * len, NULL, ( i + 1 ) * len, 1, output, NULL );

      /* Sequence s2 */
      print_line ( s2 , i * len, ( i + 1 ) * len, &nr_gaps_b, ( i + 1 ) * len, 0, output, s2_header );
      fprintf ( output, "\n" );
    }

   /* Last line of first sequence and difference */
   j = i * len;
   if ( j < m || j < n ) 
    {
      print_line ( s1, i * len, min ( m, n ), &nr_gaps_a, ( i + 1 ) * len, 0, output, s1_header );
      print_line ( diff, i * len, ( m < n ) ? m : n, NULL, ( i + 1 ) * len, 1, output, NULL );
      print_line ( s2, i * len, min ( n, m), &nr_gaps_b, ( i + 1 ) * len, 0, output, s2_header );
    }

 }
