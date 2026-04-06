/** Author and Copyright 2026 Johannes Bernhard Steffens
 *  https://github.com/johsteffens/mocutsvd
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

 /** This test program evaluates the performance of mocutsvd on arbitrary matrix sizes.
  *
  *  Compilation (compiler options are suggestions for best results):
  *     gcc -o mocutsvd_test test.c mocutsvd.c -fopenmp -march=native -O3 -lm
  *
  *  Usage:
  *     mocutsvd_test [UuVvOoRr] <rows> <cols>
  *     Run without arguments to get detailed help.
  *
  *  What it does:
  *     - Generates a (<rows> x <cols>) matrix M of random values.
  *     - Performs the SVD on M and measures the computation time.
  *     - Reconstructs R = U·Σ·V*.
  *     - Computes the RMS error of the reconstruction R from the original: RMS(R-M)
  *     - Tests U* and V* for orthonormality.
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "mocutsvd.h"

/**********************************************************************************************************************/

// ---------------------------------------------------------------------------------------------------------------------

// Measures the absolute computation time of 'expression' and stores it in time_var
#define ABS_TIME_OF( expression, time_var ) \
{ \
    struct timespec t0, t1; \
    clock_gettime( CLOCK_MONOTONIC, &t0 ); \
    expression; \
    clock_gettime( CLOCK_MONOTONIC, &t1 ); \
    time_var = t1.tv_sec - t0.tv_sec; \
    time_var += ( t1.tv_nsec - t0.tv_nsec ) * 1E-9; \
} \

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

/// RMS error between two arrays
double square( double v ) { return v * v; }
size_t size_min( size_t v1, size_t v2 ) { return v1 < v2 ? v1 : v2; }

double mat_rms_dev( const mocut_mat_s* m1, const mocut_mat_s* m2 )
{
    assert( m1->rows == m2->rows && m1->cols == m2->cols );
    double sum = 0;
    for( size_t i = 0; i < m1->rows; i++ )
    {
        double sum1 = 0;
        for( size_t j = 0; j < m1->cols; j++ )
        {
            sum1 += square( mocut_mat_s_get( m1, i, j ) - mocut_mat_s_get( m2, i, j ) );
        }
        sum += sum1;
    }
    return ( ( m1->rows * m1->cols ) > 0 ) ? sqrt( sum / ( m1->rows * m1->cols ) ) : 0;
}

// ---------------------------------------------------------------------------------------------------------------------

// reconstructs matrix M from U*, Σ, V*
void reconstruct( size_t rows, size_t cols, const double* d, const double* u, size_t u_stride, const double* v, size_t v_stride, double* r, size_t r_stride )
{
    size_t m = rows;
    size_t n = cols;
    size_t k = m < n ? m : n;
    const size_t nb = 32; // Data-local partitioning (L1 cache friendly)

    memset( r, 0, rows * r_stride * sizeof( double ) );

    // the outer loop runs in parallel
    #pragma omp parallel for
    for( size_t ib = 0; ib < m; ib += nb )
        for( size_t jb = 0; jb < n; jb += nb )
            for( size_t lb = 0; lb < n; lb += nb )
                for( size_t i = ib; i < size_min( ib + nb, m ); i++ )
                    for( size_t l = lb; l < size_min( lb + nb, k ); l++ )
                        for( size_t j = jb; j < size_min( jb + nb, n ); j++ )
                            r[ i * r_stride + j ] += u[ l * u_stride + i ] * d[ l ] * v[ l * v_stride + j ];

}

// ---------------------------------------------------------------------------------------------------------------------

/** Tests if all rows of M are orthonormal by returning RMS( E - M·M* )
 *  The return value should be close to zero on a successful test.
 */
double rms_dev_orthonormal( size_t rows, size_t cols, const double* m, size_t m_stride )
{
    mocut_mat_s* mat_r = mocut_mat_s_create_alloc( rows, rows );
    assert( mat_r );

    double* r = mat_r->data;
    size_t r_stride = mat_r->stride;

    const size_t nb = 32; // Data-local partitioning (L1 cache friendly)

    // the outer loop runs in parallel
    #pragma omp parallel for
    for( size_t ib = 0; ib < rows; ib += nb )
        for( size_t jb = 0; jb < rows; jb += nb )
            for( size_t lb = 0; lb < cols; lb += nb )
                for( size_t i = ib; i < size_min( ib + nb, rows ); i++ )
                    for( size_t l = lb; l < size_min( lb + nb, cols ); l++ )
                        for( size_t j = jb; j < size_min( jb + nb, rows ); j++ )
                            r[ i * r_stride + j ] += m[ i * m_stride + l ] * m[ j * m_stride + l ];

    double sum = 0;
    for( size_t i = 0; i < rows; i++ )
    {
        double sum_i = 0;
        for( size_t j = 0; j < rows; j++ ) sum_i += square( r[ i * r_stride + j ] - ( ( i == j ) ? 1 : 0 ) );
        sum += sum_i;
    }

    mocut_mat_s_discard( mat_r );

    return ( rows > 0 ) ? sqrt( sum / ( rows * rows ) ) : 0;
}

// ---------------------------------------------------------------------------------------------------------------------

void printhelp( void )
{
    printf
    (
        "This program evaluates the performance of mocutsvd as follows:\n"
        "- Generates a (<rows> x <cols>) matrix M of random values."
        "- Performs the SVD on M and measures the computation time."
        "- Reconstructs R = U·Σ·V*."
        "- Computes the RMS error of the reconstruction R from the original: RMS(R-M)"
        "- Tests U*, V* for orthonormality."
        "\n"
        "Usage: test <rows> <cols> [OPTION]\n"
        "\n"
        "Arguments:\n"
        "  rows: number of rows of the matrix\n"
        "  cols: number of columns of the matrix\n"
        "\n"
        "  OPTION:\n"
        "    U    (default) Computes singular vectors U*.\n"
        "    u    U* is not computed. No reconstruction.\n"
        "    V    (default) Computes singular vectors V*.\n"
        "    v    V* is not computed. No reconstruction.\n"
        "    O    (default) Tests U*, V* for orthonormality.\n"
        "    o    No orthonormality test.\n"
        "    R    (default) Given U,V: Reconstructs U·Σ·V* and computes error from the original.\n"
        "    r    No reconstruction.\n"
        "\n"
        "Examples:\n"
        "    test 4000 3000 \n"
        "    test 4000 3000 or\n"
        "    test 4000 3000 ou\n"
        "    test 4000 3000 ouv\n"
        "\n"
    );
}

// ---------------------------------------------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    int create_u = 1;
    int create_v = 1;
    int run_reconstruction = 1;
    int test_orthonormality = 1;
    size_t rows = 0;
    size_t cols = 0;

    if( argc < 3 || argc > 4 ) { printhelp();  exit( 0 ); }
    rows = atoi( argv[ 1 ] );
    cols = atoi( argv[ 2 ] );
    if( argc == 4 )
    {
        for( const char* opt = argv[ 3 ]; *opt != 0; opt++ )
        {
            switch( *opt )
            {
                case 'U': create_u = 1; break;
                case 'u': create_u = 0; break;
                case 'V': create_v = 1; break;
                case 'v': create_v = 0; break;
                case 'O': test_orthonormality = 1; break;
                case 'o': test_orthonormality = 0; break;
                case 'R': run_reconstruction = 1; break;
                case 'r': run_reconstruction = 0; break;
                default: { printhelp();  exit( 0 ); }
            }
        }
    }

    // maximum allowed error for a correct result
    double dev_limit = 1E-6;

    printf( "M = Randomized matrix of size (%zu x %zu).\n", rows, cols );

    size_t m = rows;
    size_t n = cols;
    size_t k = m < n ? m : n;

    mocut_mat_s* mat_m = mocut_mat_s_create_alloc( m, n );
    mocut_mat_s* mat_a = mocut_mat_s_create_alloc( m, n );
    mocut_mat_s* mat_u = create_u ? mocut_mat_s_create_alloc( k, m ) : NULL;
    mocut_mat_s* mat_v = create_v ? mocut_mat_s_create_alloc( k, n ) : NULL;

    mat_randomize( mat_m, 1234 );

    /// mat_a <- mat_m
    mocut_mat_s_copy( mat_a, mat_m );

    // SVD

    printf( "Running SVD M -> " );
    if( mat_u ) printf( "U*, " );
    printf( "Σ" );
    if( mat_v ) printf( ", V*" );
    printf( " : " );
    fflush( stdout );

    double time = 0;
    int svd_return = 0;
    ABS_TIME_OF( svd_return = mocut_svd( mat_a, mat_u, mat_v ), time )
    if( svd_return )
    {
        if( svd_return == MOCUT_WRN_CONVERGENCE )
        {
            printf( "SVD did not converge\n" );
        }
        else
        {
            fprintf( stderr, "%s\n", mocut_err_text( svd_return ) );
            exit( 1 );
        }
    }

    printf( "(%.2f sec)\n", time );
    if( mocut_mat_s_contains_nan( mat_a ) ) { fprintf( stderr, "mat_a: NAN detected." ); exit( 1 ); }
    if( mat_u && mocut_mat_s_contains_nan( mat_u ) ) { fprintf( stderr, "mat_u: NAN detected." ); exit( 1 ); }
    if( mat_v && mocut_mat_s_contains_nan( mat_v ) ) { fprintf( stderr, "mat_v: NAN detected." ); exit( 1 ); }

    int exit_status = 0;

    if( run_reconstruction && mat_u && mat_v )
    {
        printf( "Testing Reconstruction: R = U·Σ·V* : " );
        fflush( stdout );

        // we save singular values in vec_s and then set all values in mat_a zero to use it for reconstruction
        mocut_mat_s* vec_s = mocut_mat_s_create_alloc( 1, k );
        assert( vec_s != NULL );
        for( size_t i = 0; i < k; i++ ) mocut_mat_s_set( vec_s, 0, i, mocut_mat_s_get( mat_a, i, i ) );
        mocut_mat_s_set_zero( mat_a );

        ABS_TIME_OF( reconstruct( rows, cols, vec_s->data, mat_u->data, mat_u->stride, mat_v->data, mat_v->stride, mat_a->data, mat_a->stride ), time )
        printf( "(%.2f sec) ", time );

        double dev = mat_rms_dev( mat_a, mat_m );
        printf( "RMS(R-M)   : %.3g ", dev );
        if( dev <= dev_limit )
        {
            printf( "-> Successful\n" );
        }
        else
        {
            printf( "-> Failed\n" );
            exit_status = 1;
        }

        mocut_mat_s_discard( vec_s );
    }

    mocut_mat_s_discard( mat_a );
    mocut_mat_s_discard( mat_m );

    if( test_orthonormality && mat_u )
    {
        printf( "Testing Orthonormality of U*       : " );
        fflush( stdout );
        double dev = 0;

        ABS_TIME_OF( dev = rms_dev_orthonormal( mat_u->rows, mat_u->cols, mat_u->data, mat_u->stride ), time )
        printf( "(%.2f sec) ", time );
        printf( "RMS(U*·U-I): %.3g ", dev );
        if( dev <= dev_limit )
        {
            printf( "-> Successful\n" );
        }
        else
        {
            printf( "-> Failed\n" );
            exit_status = 1;
        }
    }

    if( test_orthonormality && mat_v )
    {
        printf( "Testing Orthonormality of V*       : " );
        fflush( stdout );
        double dev = 0;

        ABS_TIME_OF( dev = rms_dev_orthonormal( mat_v->rows, mat_v->cols, mat_v->data, mat_v->stride ), time )
        printf( "(%.2f sec) ", time );
        printf( "RMS(V*·V-I): %.3g ", dev );
        if( dev <= dev_limit )
        {
            printf( "-> Successful\n" );
        }
        else
        {
            printf( "-> Failed\n" );
            exit_status = 1;
        }
    }

    mocut_mat_s_discard( mat_u );
    mocut_mat_s_discard( mat_v );

    return exit_status;
}

// ---------------------------------------------------------------------------------------------------------------------

/**********************************************************************************************************************/

