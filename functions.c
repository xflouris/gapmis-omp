/**
    GapMis: a tool for pairwise sequence alignment with a single gap.
    Copyright (C) 2011 Solon P. Pissis, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <getopt.h>
#include "functions.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

/* long options for command-line switches */
static struct option long_options[] =
 {
   { "sequences-a",         required_argument, NULL, 'a' },
   { "sequences-b",         required_argument, NULL, 'b' },
   { "gap-open-penalty",    required_argument, NULL, 'g' },
   { "gap-extend-penalty",  required_argument, NULL, 'e' },
   { "output-file",         required_argument, NULL, 'o' },
   { "output-format",       required_argument, NULL, 'f' },
   { "data-file",           required_argument, NULL, 'd' },
   { "help",                no_argument,       NULL, 'h' },
   { "max-gap",             required_argument, NULL, 'm' },
   { "threads",             required_argument, NULL, 't' },
   { NULL,                  0,                 NULL, 0   }
 };

/* Decode the input switches */
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;

   /* initialisation */
   sw -> seq_a          = NULL;
   sw -> seq_b          = NULL;
   sw -> gap_open_pen   = 10.0;
   sw -> gap_extend_pen = 0.5;
   sw -> max_gap        = -1;
   sw -> out_file       = ( char * ) malloc ( 15 * sizeof ( char ) );
   sw -> matrix         = ( char * ) malloc ( 15 * sizeof ( char ) );
   sw -> threads        = 1;
   sw -> format         = 0;
   strcpy ( sw -> out_file, "gapmis.out" );
   strcpy ( sw -> matrix, "EDNAFULL" );

   while ( ( opt = getopt_long ( argc, argv, "a:b:g:e:o:f:d:m:t:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> seq_a = optarg;
           break;
         
         case 'b':
           sw -> seq_b = optarg;
           break;

         case 'f':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> format = val;
           break;
         
         case 't':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> threads = val;
           break;

         case 'o':
           free ( sw -> out_file );
           sw -> out_file = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> out_file, optarg );
           break;

         case 'd':
           free ( sw -> matrix );
           sw -> matrix = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> matrix, optarg );
           break;

         case 'g':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> gap_open_pen = val;
           break;
         
         case 'h':
           return ( 0 );

         case 'e':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> gap_extend_pen = val;
           break;

         case 'm':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> max_gap = val;
           break;

       }
    }

   return ( optind );
 }

/* Read fasta files with the sequences */
unsigned int read_fasta_files ( const char* query_fa,  struct TSeq ** p, unsigned int* cur_alloc_query, const char* target_fa, struct TSeq ** t, unsigned int* cur_alloc_target )
 {
   FILE*        fd_query;
   FILE*        fd_target;
   unsigned int max_alloc_target = 0;
   unsigned int max_alloc_query = 0;
   char         read[BUFFER_SIZE];     // Buffer for storing the read 
   char         readId[BUFFER_SIZE];   // Buffer for storing the Id of the read 

   *cur_alloc_target = *cur_alloc_query = 0;
   
   /* open file 1 (patterns) */
   if ( ! ( fd_query = fopen ( query_fa, "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file %s\n", query_fa );
      return ( 0 ); 
    }

   /* read file 1 (patterns) */
   while ( ! feof ( fd_query ) )
    {
      if ( fgetc ( fd_query ) && fgets ( readId, BUFFER_SIZE, fd_query ) && fgets ( read, BUFFER_SIZE, fd_query ) )
       {
          if ( *cur_alloc_query >= max_alloc_query )
           {
             *p = ( struct TSeq * ) realloc ( *p,   ( max_alloc_query + ALLOC_SIZE ) * sizeof ( struct TSeq ) ); 
             max_alloc_query += ALLOC_SIZE;
           }
          read[ strlen ( read ) - 1] = 0;
          readId[ strlen ( readId ) - 1] = 0;
          ( *p ) [ *cur_alloc_query ] . data   = strdup ( read );
          ( *p ) [ *cur_alloc_query ] . header = strdup ( readId );

          ++ ( *cur_alloc_query );
       }
    }
   fclose ( fd_query );

   /* open file 2 (texts) */
   if ( ! ( fd_target = fopen ( target_fa, "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file %s\n", target_fa );
      return ( 0 );
    }

   /* read file 2 (texts) */
   while ( ! feof ( fd_target ) )
    {
      if ( fgetc ( fd_target ) && fgets ( readId, BUFFER_SIZE, fd_target ) && fgets ( read, BUFFER_SIZE, fd_target ) )
       {
         if ( *cur_alloc_target >= max_alloc_target )
          {
            *t = ( struct TSeq * ) realloc ( *t, ( max_alloc_target + ALLOC_SIZE ) * sizeof ( struct TSeq ) );
            max_alloc_target += ALLOC_SIZE;
          }
         
         read[ strlen ( read ) - 1] = 0;
         readId[ strlen ( readId ) - 1] = 0;
         ( *t ) [ *cur_alloc_target ] . data   = strdup ( read );
         ( *t ) [ *cur_alloc_target ] . header = strdup ( readId );

         ++ ( *cur_alloc_target );
       }
    }
   fclose ( fd_target );

   /* adjust query and target sizes */
   *t = ( struct TSeq * ) realloc ( *t , ( *cur_alloc_target + 1 ) * sizeof ( struct TSeq ) );
   ( *t ) [ *cur_alloc_target ] . data   = NULL;
   ( *t ) [ *cur_alloc_target ] . header = NULL;

   *p = ( struct TSeq * ) realloc ( *p , ( *cur_alloc_query + 1 ) * sizeof ( struct TSeq ) );
   ( *p ) [ *cur_alloc_query ] . data   = NULL;
   ( *p ) [ *cur_alloc_query ] . header = NULL;

   return ( 1 );
 }


/*
The dynamic programming algorithm for calculating matrices G and H
*/
unsigned int dp_algorithm ( double ** G, unsigned int ** H, char * t, unsigned int n, char * p, unsigned int m, unsigned int matrix, unsigned int MAXgap )
{
	double gap;
	double mis;
	unsigned int i;
	unsigned int j;
	double matching_score;
	unsigned int j_min;
	unsigned int j_max;
	unsigned int valM;
   	unsigned int 	i_max;

   	i_max = min ( n, m + MAXgap );

   	for( i = 0; i < n + 1 ; i++ )	H[i][0] = i;
	for( j = 0; j < m + 1 ; j++ )	H[0][j] = j;

	for( i = 1; i < i_max + 1; i++ )		
        {
                j_min = max ( 1, (int) ( i - MAXgap ));
                j_max = min ( m, (int) ( i + MAXgap ));
		for( j = j_min; j <= j_max; j++ )
		{

			matching_score = ( matrix ? (double) pro_delta( t[i - 1], p[j - 1] ) : (double) nuc_delta( t[i - 1], p[j - 1] ) ) ;
			if ( matching_score == ERR )
				return 0;

		        mis = G[i - 1][j - 1] + matching_score;
	   		gap = G[j][j];
	   		valM = i - j;

	   		if( j > i )	
	     		{
	       			gap = G[i][i];
               			valM = j - i;
	     		}

           		if( gap > mis )		H[i][j] = valM;
	   		if( i == j )		gap = mis - 1;

           		G[i][j] = max ( mis, gap );
		}
        }
	return 1;
}

/* Returns the score for matching character a and b based on EDNAFULL matrix */
int nuc_delta ( char a, char b )
 {
   unsigned int index_a = nuc_char_to_index ( a );
   unsigned int index_b = nuc_char_to_index ( b );

   if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) )
     return ( EDNAFULL_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
int pro_delta ( char a, char b )
 {
   unsigned int index_a = pro_char_to_index( a );
   unsigned int index_b = pro_char_to_index( b );

   if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
     return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the index of char a in EDNAFULL matrix */
unsigned int nuc_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'T':
        index = 1; break;

      case 'G':
        index = 2; break;

      case 'C':
        index = 3; break;

      case 'S':
        index = 4; break;

      case 'W':
        index = 5; break;

      case 'R':
        index = 6; break;

      case 'Y':
        index = 7; break;

      case 'K':
        index = 8; break;

      case 'M':
        index = 9; break;

      case 'B':
        index = 10; break;

      case 'V':
        index = 11; break;

      case 'H':
        index = 12; break;

      case 'D':
        index = 13; break;

      case 'N':
        index = 14; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the nucleotide sequences!!!\n" );
        index = ERR; break;
    }
   
   return ( index );
 }

/* Returns the index of char a in EBLOSUM62 matrix */
unsigned int pro_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'R':
        index = 1; break;

      case 'N':
        index = 2; break;

      case 'D':
        index = 3; break;

      case 'C':
        index = 4; break;

      case 'Q':
        index = 5; break;

      case 'E':
        index = 6; break;

      case 'G':
        index = 7; break;

      case 'H':
        index = 8; break;

      case 'I':
        index = 9; break;

      case 'L':
        index = 10; break;

      case 'K':
        index = 11; break;

      case 'M':
        index = 12; break;

      case 'F':
        index = 13; break;

      case 'P':
        index = 14; break;

      case 'S':
        index = 15; break;

      case 'T':
        index = 16; break;

      case 'W':
        index = 17; break;

      case 'Y':
        index = 18; break;

      case 'V':
        index = 19; break;

      case 'B':
        index = 20; break;

      case 'Z':
        index = 21; break;

      case 'X':
        index = 22; break;

      case '*':
        index = 23; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the protein sequences!!!\n" );
        index = ERR; break;
    }
   return ( index );
 }

/* Computes the limits of the i-th coordinate for the matrix G in constant time */
unsigned int i_limits( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int MAXgap )
 {
   if ( (int) m - (int) MAXgap < 0 )    (* up )    = 0;
   else                                 (* up )    = m - MAXgap;
   if ( m + MAXgap > n )                (* down )  = n;
   else                                 (* down )  = m + MAXgap;
   return ( 0 );
 }

/* Computes the limits of the j-th coordinate for matrix G and H in constant time */
unsigned int j_limits ( unsigned int i, unsigned int m, unsigned int * left, unsigned int * right, unsigned int MAXgap )
 {
   if ( (int) i - (int) MAXgap > 0 )    (* left )   = i - MAXgap;
   else                                 (* left )   = 1;
   if ( i + MAXgap > m )                (* right )  = m;
   else                                 (* right )  = i + MAXgap;
   return ( 0 );
 }

/*
Computes the optimal alignment using matrix G in O(2*MAXgap+1) time
Note:	double gap_open_penalty, double gap_extend_penalty, double gap_open_offset_penalty are arguments given by the user to represent the gap penalty.
*/
unsigned int opt_solution ( 	double** G, 
				unsigned int n, 
				unsigned int m, 
				unsigned int MAXgap, 
				double gap_open_penalty, 
				double gap_extend_penalty, 
				double* MAXscore, 
				unsigned int* MINgap, 
				unsigned int* where, 
				unsigned int* start 
)
{
	double score = -DBL_MAX;
	unsigned int i, j;
			
	unsigned int up = 0;
	unsigned int down = 0;
	i_limits( n, m, &up, &down, MAXgap );			// computes the i coordinates for matrix G for the last column

	for ( i = up ; i <= down ; i++ )
	{
		double temp_score = 0.0;
		if ( i < m )
		{
			if ( m - i <= MAXgap )
			{
				temp_score = total_scoring ( m - i, G[i][m], gap_open_penalty, gap_extend_penalty );
				if ( temp_score > score )
				{
					score = temp_score;
					( *MAXscore ) = score; 
					( *MINgap ) = m - i;
					( *where ) = 1;		//where: gap is in the text and start backtracing from the last column
					( *start ) = i;		//backtrace from cell G[start,m]
				}
			}
		}
		else if ( i > m )
		{
			if ( i - m <= MAXgap )
			{
				temp_score = total_scoring( i - m, G[i][m], gap_open_penalty, gap_extend_penalty );
				if (  temp_score > score )
				{
					score = temp_score;
					( *MAXscore ) = score; 
					( *MINgap ) = i - m;
					( *where ) = 2;		//where: gap is in the pattern and start backtracing from last column
					( *start ) = i;		//backtrace from cell G[start,m]
				}
			}
		}
		else if ( i == m )
		{
			temp_score = total_scoring( 0, G[i][m], gap_open_penalty, gap_extend_penalty );
			if (  temp_score > score ) // mgap = 0 
			{
				score = temp_score;
				( *MAXscore ) = score; 
				( *MINgap ) = 0;
				( *where ) = 0;		//there is no gap
				( *start ) = m;		//no need to backtrace
			}
		}
	}

	unsigned int left = 0;
	unsigned int right = 0;
	j_limits ( n, m, &left, &right, MAXgap );	// computes the j coordinates for matrix G for the last row

	for ( j = left ; j < right ; j++ )
	{
		double temp_score = 0;
		if ( n - j <= MAXgap )
		{
			temp_score = total_scoring( n - j, G[n][j], gap_open_penalty, gap_extend_penalty );
			if (  temp_score > score )
			{
				score = temp_score;
				( *MAXscore ) = score; 
				( *MINgap ) = n - j;
				( *where ) = 3;		//where: gap is in the pattern and start backtracing from last row
				( *start ) = j;		//backtrace from cell G[n,start]
			}
		}
	}

	return 1;
}

/* Gives the position of the gap in O(m) time */
unsigned int backtracing ( unsigned int** H, unsigned int m, unsigned int n, unsigned int start, unsigned int where, unsigned int* gap_pos )
{
	unsigned int i, j;
	( *gap_pos ) = 0;

        if ( where == 1 || where == 2 )
	{
		i = start; j = m; 	//we start backtracing from the last column
        }
        else
	{
		i = n; j = start;	//we start backtracing from the last row
	}
	while ( i > 0 && j > 0)
	{
		if ( H[i][j] == 0 )
		{
			--i; --j;
		}
		else				
		{
			if ( i > j )	
				( *gap_pos ) = j;	
			else		
				( *gap_pos ) = i;
			break;	
		}
	}
	return 1;
}

/*
Gives the total score of an alignment in constant time
Note: double matrix_score is the value of G[i][m], i.e. the score of an alignment WITHOUT the gap penalties
*/
double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty )
 {
   return ( matrix_score + ( ( gap > 0 ) ? ( gap - 1 ) * gap_extend_penalty + gap_open_penalty : 0 ) );
 }

/* Swaps the text and the pattern in case m > n */
unsigned int swap_txt_pat ( struct TSeq ** seqa, unsigned int * n, struct TSeq ** seqb, unsigned int * m )
 {
   struct TSeq * tmp;

   tmp   = *seqa;
   *seqa = *seqb;
   *seqb = tmp;
   
   SWAP ( *n, *m );
   
   return ( 1 );
 }

