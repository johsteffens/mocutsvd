/// Author and Copyright 2026 Johannes Bernhard Steffens

 /** This simple example program demonstrates the usage of mocutsvd.
  *  You may feely use code snippets from this file for your own integration of mocutsvd.
  *
  *  This example code is not designed for very large matrices.
  *  Use test.c for a comprehensive evaluation of mocutsvd on large matrices.
  *
  *  Compilation:
  *     gcc -o example example.c soputsvd.c -fopenmp -march=native -O3 -lm
  *     (-fopenmp is optional)
  *
  *  Usage:
  *     example <rows> <cols>
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

#include "soputsvd.h"

/**********************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------------

/** Full-cycle linear congruential random generator random generator with good lattice structure.
 *  (Not suitable for cryptographic purposes)
 */
uint64_t lcg_random( uint64_t v ) { return v * ( 505810533149048947ull * 4 + 1 ) + 1; }

// ---------------------------------------------------------------------------------------------------------------------

double* allocate_matrix( size_t rows, size_t cols )
{
    double* mat = NULL;
    size_t bytes = rows * cols * sizeof( double );
    if( bytes == 0 ) return NULL;

    // on large matrices try aligned_alloc( 0x40, bytes ) for a small performance boost on some platforms
    if( !( mat = malloc( bytes ) ) )
    {
        fprintf( stderr, "Failed to allocate %zu bytes", bytes );
        exit( 1 );
    }
    return mat;
}

// ---------------------------------------------------------------------------------------------------------------------

void free_matrix( double* m )
{
    if( !m ) return;
    free( m );
}

// ---------------------------------------------------------------------------------------------------------------------

void print_matrix( const double* mat, size_t rows, size_t cols )
{
    if( rows > 1 ) printf( "(%zu x %zu)\n", rows, cols );
    for( size_t i = 0; i < rows; i++ )
    {
        for( size_t j = 0; j < cols; j++ ) printf( "% 6.3f ", mat[ i * cols + j ] );
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

    double* buf_a = allocate_matrix( m, n ); // space for matrix -> singular values
    double* mat_u = allocate_matrix( k, m ); // space for left singular vectors
    double* mat_v = allocate_matrix( k, n ); // space for right singular vectors

    // We fill matrix a with random values within range [-1.0, +1.0]
    uint64_t rval = 1234;
    for( size_t i = 0; i < n * m; i++ )
    {
        buf_a[ i ] = ( ( double )( rval = lcg_random( rval ) ) * ( 1.0 / ( 1ull << 63 ) ) ) - 1.0;
    }

    printf( "Original Matrix\n" );
    print_matrix( buf_a, rows, cols );

    // SVD
    if( soput_svd( rows, cols, buf_a, mat_u, mat_v ) != 0 ) printf( "SVD did not converge\n" );

    printf( "\nSingular Values (Diagonal values of Σ)\n" );
    print_matrix( buf_a, 1, k );

    printf( "\nU*: Left singular vectors are row-vectors of U*\n" );
    print_matrix( mat_u, k, m );

    printf( "\nV*: Right singular vectors are row vectors of V*\n" );
    print_matrix( mat_v, k, n );

    // We now reconstruct the original matrix from singular values and vectors
    double* mat_m = allocate_matrix( rows, cols );

    // Reconstruction: M = U·Σ·V*
    memset( mat_m, 0, rows * cols * sizeof( double ) ); // init matrix with 0
    for( size_t l = 0; l < k; l++ )
    {
        for( size_t i = 0; i < m; i++ )
        {
            for( size_t j = 0; j < n; j++ )
            {
                mat_m[ i * cols + j ] += buf_a[ l ] * mat_u[ l * m + i ] * mat_v[ l * n + j ];
            }
        }
    }

    printf( "\nReconstructed Matrix. This should look like the original.\n" );
    print_matrix( mat_m, rows, cols );

    free_matrix( buf_a );
    free_matrix( mat_u );
    free_matrix( mat_v );
    free_matrix( mat_m );

    return 0;
}

// ---------------------------------------------------------------------------------------------------------------------

/**********************************************************************************************************************/

