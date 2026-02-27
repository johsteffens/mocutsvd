/// Author and Copyright 2026 Johannes Bernhard Steffens

 /** This simple example program demonstrates the usage of mocutsvd.
  *  You may freely use code snippets from this file for your own integration of mocutsvd.
  *
  *  This example code is not designed for very large matrices.
  *  Use test.c for a comprehensive evaluation of mocutsvd on large matrices.
  *
  *  Compilation:
  *     gcc -o mocutsvd_example example.c mocutsvd.c -fopenmp -march=native -O3 -lm
  *     (-fopenmp is optional)
  *
  *  Usage:
  *     mocutsvd_example <rows> <cols>
  *
  *  What it does:
  *    - Generates a (<rows> x <cols>) matrix M of random values.
         <rows>, <cols> are taken from the command line.
  *    - Displays M
  *    - Performs the SVD on that matrix
  *    - Displays the decomposition U*, Σ, V*
  *    - Reconstructs M = U·Σ·V*
  *    - Displays the reconstructed M
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <memory.h>

#include "mocutsvd.h"

/**********************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------------

/** Full-cycle linear congruential random generator random generator with good lattice structure.
 *  (Not suitable for cryptographic purposes)
 */
uint64_t lcg_random( uint64_t v ) { return v * ( 505810533149048947ull * 4 + 1 ) + 1; }

// ---------------------------------------------------------------------------------------------------------------------

void mat_randomize( mocut_mat_s* m, uint64_t seed )
{
    uint64_t rval = seed;
    for( size_t i = 0; i < m->rows; i++ )
        for( size_t j = 0; j < m->cols; j++ )
            mocut_mat_s_set( m, i, j, ( ( double )( rval = lcg_random( rval ) ) * ( 1.0 / ( 1ull << 63 ) ) ) - 1.0 );
}

// ---------------------------------------------------------------------------------------------------------------------

void mat_print( const mocut_mat_s* m )
{
    if( m->rows > 1 ) printf( "(%zu x %zu)\n", m->rows, m->cols );
    for( size_t i = 0; i < m->rows; i++ )
    {
        for( size_t j = 0; j < m->cols; j++ ) printf( "% 6.3f ", mocut_mat_s_get( m, i, j ) );
        printf( "\n" );
    }
}

// ---------------------------------------------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    if( argc < 3 ) { fprintf( stderr, "Usage: example <rows> <cols>\n" ); exit( 1 ); }
    size_t rows = atoi( argv[ 1 ] );
    size_t cols = atoi( argv[ 2 ] );

    size_t m = rows;
    size_t n = cols;
    size_t k = m < n ? m : n;

    mocut_mat_s* mat_a = mocut_mat_s_create_alloc( m, n );
    mocut_mat_s* mat_u = mocut_mat_s_create_alloc( k, m );
    mocut_mat_s* mat_v = mocut_mat_s_create_alloc( k, n );

    // We fill matrix a with random values within range [-1.0, +1.0]
    mat_randomize( mat_a, 1234 );

    printf( "Original Matrix\n" );
    mat_print( mat_a );

    // SVD
    int err = 0;
    if( ( err = mocut_svd( mat_a, mat_u, mat_v ) ) )
    {
        if( err == MOCUT_WRN_CONVERGENCE )
        {
            printf( "SVD did not converge\n" );
        }
        else
        {
            fprintf( stderr, "%s\n", mocut_err_text( err ) );
            exit( 1 );
        }
    }

    printf( "\nMatrix Σ with diagonal singular values\n" );
    mat_print( mat_a );

    printf( "\nU*: Left singular vectors are row-vectors of U*\n" );
    mat_print( mat_u );

    printf( "\nV*: Right singular vectors are row vectors of V*\n" );
    mat_print( mat_v );

    // We now reconstruct the original matrix from singular values and vectors
    mocut_mat_s* mat_m = mocut_mat_s_create_alloc( m, n ); // this function also zeroes all elements

    for( size_t l = 0; l < k; l++ )
        for( size_t i = 0; i < m; i++ )
            for( size_t j = 0; j < n; j++ )
                *mocut_mat_s_ptr( mat_m, i, j ) +=
                    *mocut_mat_s_ptr( mat_u, l, i ) *
                    *mocut_mat_s_ptr( mat_a, l, l ) *
                    *mocut_mat_s_ptr( mat_v, l, j );

    printf( "\nReconstructed Matrix. This should look like the original.\n" );
    mat_print( mat_m );

    mocut_mat_s_discard( mat_a );
    mocut_mat_s_discard( mat_u );
    mocut_mat_s_discard( mat_v );
    mocut_mat_s_discard( mat_m );

    return 0;
}

// ---------------------------------------------------------------------------------------------------------------------

/**********************************************************************************************************************/

