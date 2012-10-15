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
#include <omp.h>
#include "functions.h"
#include "output.h"

int main ( int argc, char ** argv)
 {
   struct       TSwitch  sw;
   
   unsigned int MAXgap;			//input arguments		
   double       gap_open_pen;
   double       gap_extend_pen;
   unsigned int scoring_matrix;
   char       * out_file;
   unsigned int threads;
   unsigned int format;
   struct       Tin    in;
     
   struct       TSeq * t = NULL;	//texts t
   unsigned int n; 			//length of text
   unsigned int t_i;
   struct       TSeq * p = NULL;	//patterns p
   unsigned int m; 			//length of pattern
   unsigned int p_i;
   
   struct       Tout * out = NULL;
   
   unsigned int i;

   unsigned int cur_alloc_target;
   unsigned int cur_alloc_query;
   
   /* checks the arguments */
   i = decode_switches ( argc, argv, &sw );

   if ( i < 5 || ! sw . seq_a || ! sw . seq_b ) 
    {
      usage ();
      return ( 1 );
    }
   else 
    {
      gap_open_pen   = - sw . gap_open_pen;	//the penalties should have a negative value
      gap_extend_pen = - sw . gap_extend_pen;
      out_file       =   sw . out_file;

      if ( ! strcmp ( "EDNAFULL", sw . matrix ) )       scoring_matrix = 0;
      else if ( ! strcmp ( "EBLOSUM62", sw . matrix ) ) scoring_matrix = 1;
      else
       {
         fprintf ( stderr, "Error: scoring matrix argument should be `EDNAFULL' for nucleotide sequences or `EBLOSUM62' for protein sequences!!!\n" );
         return ( 1 );
       }

      threads       =   sw . threads;
      format        =   sw . format;
    }

   /* reads the input data */
   read_fasta_files ( sw . seq_a,  &p, &cur_alloc_query, sw . seq_b, &t, &cur_alloc_target );
   
   /* allocate the space for the output : | querys | x | targets | */
   out = ( struct Tout * ) calloc ( cur_alloc_target * cur_alloc_query, sizeof ( struct Tout ) );
   if( out == NULL ) 
   {
      fprintf ( stderr, "Error: space for output could not be allocated!!!\n" );
      return ( 1 );
   }
   
   /* calculate text's and pattern's length */
   n = strlen ( t[0] . data );
   m = strlen ( p[0] . data );
   if( m > n )
    {
      fprintf ( stderr, "Error: the length of the patterns should be less than the length of the texts!!!\n" );
      return ( 1 );
    }

   /* checks the max gap length allowed: MAXgap < n */
   MAXgap =  ( sw . max_gap <= -1 ) ?  n - 1 : sw . max_gap;   
   if( MAXgap >= n )
    {
      fprintf ( stderr, "Error: the max gap length should be less than the length of the text!!!\n" );
      return ( 1 );
    }

   /* set the num of threads to be used */
   omp_set_num_threads( threads );
      
   /* many-to-many pattern matching */
   for ( t_i = 0; t_i < cur_alloc_target; ++ t_i )
     {
       #pragma omp parallel for private ( p_i )  
       for ( p_i = 0; p_i < cur_alloc_query; ++ p_i )
         {
           /* Allocate DP matrices G and H by the threads */
           double       ** G;
           unsigned int ** H;
           unsigned int    j;
           G = ( double ** ) malloc ( ( n + 1 ) * sizeof ( double * ) );
           G[0] = ( double * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof ( double ) );
           for ( j = 1; j < n + 1; ++ j )
             G[j] = ( void * ) G[0] + j * ( m + 1 ) * sizeof ( double );
           H = ( unsigned int ** ) malloc ( ( n + 1 ) * sizeof ( unsigned int * ) );
           H[0] = ( unsigned int * ) calloc ( ( n + 1 ) * ( m + 1 ) , sizeof ( unsigned int ) );
           for ( j = 1 ; j < n + 1 ; ++ j )
             H[j] = ( void * ) H[0] + j * ( m + 1 ) * sizeof ( unsigned int );

   	   unsigned int start    = 0;		//where to start backtracing
           unsigned int gap_pos  = 0;		//position of the gap
           unsigned int where    = 0;		//where is the gap: text or pattern
           unsigned int MINgap   = 0;		//to be computed
           double       MAXscore = 0;           //the score
           unsigned int o_i      = t_i * cur_alloc_query + p_i;

           /* dynamic programming algorithm */
           dp_algorithm( G, H, t[t_i] . data, n, p[p_i] . data, m, scoring_matrix, MAXgap );
              
           /* computes the optimal alignment based on the matrix score and the gap function */
           opt_solution ( G, n, m, MAXgap, gap_open_pen, gap_extend_pen, &MAXscore, &MINgap, &where, &start );
     
           /* computes the position of the gap */
           if ( MINgap > 0 ) backtracing ( H, m, n, start, where, &gap_pos );
           else gap_pos = 0;

           /* store the output */
	   out[o_i] . max_score = MAXscore;
   	   out[o_i] . min_gap   = MINgap;
           out[o_i] . where     = where;
           out[o_i] . gap_pos   = gap_pos;

           /* Deallocation by the threads*/
           free ( G[0] );
           free ( H[0] );
           free ( G );
           free ( H );
         }
     }

   /* output the results sequentially */
   in . scoring_matrix = scoring_matrix;
   in . gap_open_pen   = gap_open_pen;
   in . gap_extend_pen = gap_extend_pen;
   in . max_gap        = MAXgap;
   if ( ! format )
    {
      if ( ! ( results_many_to_many_scr ( out_file, p, cur_alloc_query, t, cur_alloc_target, out ) ) ) 
       {
         fprintf(stderr, "Error: results_many_to_many_scr() failed!!!\n");
         return ( 1 );	
       }
    } 
   else
    {
      if ( ! ( results_many_to_many ( out_file, p, cur_alloc_query, t, cur_alloc_target, &in, out ) ) ) 
       {
         fprintf(stderr, "Error: results_many_to_many() failed!!!\n");
         return ( 1 );	
       }
    }

   /* Dealloaction */
   for ( p_i = 0; p_i < cur_alloc_query; ++ p_i )
    {
      free ( ( void * ) p[p_i] . data );
      free ( ( void * ) p[p_i] . header );
    }
   free ( p );

   for ( t_i = 0; t_i < cur_alloc_target; ++ t_i )
    {
      free ( ( void * ) t[t_i] . data );
      free ( ( void * ) t[t_i] . header );
    }
   free ( t );
   
   free ( out );
   free ( sw . out_file );
   free ( sw . matrix );

   return ( 0 );
 }



