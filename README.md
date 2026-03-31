# MocUT SVD - Singular Value Decomposition

## What it is

MocUT SVD is a stable, fast and easy-to-use implementation of the Singular Value Decomposition for C or C++ programs. It was completely redesigned from scratch and contains several improvements compared to traditional implementations. It is specially optimized for modern CPU architectures.

### Benefits

* **Stable**: Same or better stability as the widely used Golub-Kahan-Reinsch SVD approach.
* **Fast**: Efficient [true-scalable](doc/true_scalability.md) computation for small up to very large matrices.
* **Parallel**: By default all available CPU-cores are used.
* **In-Place**: Operates within the matrix-provided memory space and allocates no extra heap memory.
* **Stand-Alone**: Only depends on the C standard library. No third-party library required.
* **Portable**: Platform-agnostic design.

### Quick Start
* [Evaluation](#quick-test)
* [Example](#quick-example)
* [Performance Charts](#performance-charts)

## In a Nutshell

You only need the two files `mocutsvd.h` and `mocutsvd.c`. 

Include `mocutsvd.h` in your program (C or C++ code) and compile `mocutsvd.c` among your sources.

#### Example Code

``` C
#include "mocutsvd.h"

....

// Creates (20 x 10)-matrix A
mocut_mat_s* a = mocut_mat_s_create_alloc( 20, 10 );

// Fills matrix A with data ...
for( size_t i = 0; i < a->rows; i++ )
	for( size_t j = 0; j < a->cols; j++ )
		mocut_mat_s_set( a, i, j, /* your value at position [i][j] */  );

mocut_mat_s* u = mocut_mat_s_create(); // Matrix U*
mocut_mat_s* v = mocut_mat_s_create(); // Matrix V*

// SVD: A -> U*, Σ, V*; computing singular vectors in U* or V* is optional; pass NULL where not needed.
mocut_svd( a, u, v ); 

// At this point: 
//  - 'a' represents Σ. The diagonal elements of 'a' contain singular values in descending order.
//  - 'u', 'v', represent the matrices U*, V*:
//       - U*: (10 x 20); its row vectors represent the left singular vectors
//       - V*: (10 x 10); its row vectors represent the right singular vectors

// Access any matrix element via function mocut_mat_s_get or mocut_mat_s_ptr

....
    
// memory cleanup
mocut_mat_s_discard( a );
mocut_mat_s_discard( u );
mocut_mat_s_discard( v );

// For more details: See example.c
// For a comprehensive performance test: See test.c
```

#### Build

* Compile with gcc and flags below or use compatible flags of your preferred compiler-tool-chain.
* Compiler/Linker Flags (gcc):
  * Required: `-lm`. 
  * Recommended for general optimization: `-O3`
  * Recommended for advanced CPU specific optimization: `-march=native`. 
  * Recommended to unlock outer parallelity: `-fopenmp` (both: compiler and linker need this flag)

* **Example**: ```bash $ gcc -o mocutsvd_example example.c mocutsvd.c -fopenmp -march=native -O3 -lm```

* There is also a simple [makefile](makefile).

### Quick Test

```bash
git clone https://github.com/johsteffens/mocutsvd
cd mocutsvd
make mocutsvd_test

# Runs SVD on a random matrix of specified size; tests accuracy and measures comutation time.
./mocutsvd_test 5000 5000 
```

Possible output: (Test on a workstation with 16 HT cores)

```
M = Randomized matrix of size (5000 x 5000).
Running SVD M -> U*, Σ, V* : (7.37 sec)
Testing Reconstruction: R = U·Σ·V* : (1.00 sec) RMS(R-M)   : 9.49e-15 -> Successful
Testing Orthonormality of U*       : (1.41 sec) RMS(U*·U-I): 1.67e-16 -> Successful
Testing Orthonormality of V*       : (1.42 sec) RMS(V*·V-I): 1.73e-16 -> Successful
```

[`test.c`](test.c) is intended for performance testing. 
It creates a random matrix, decomposes it with a time measurement, reconstructs the original form the decomposed factors and computes the RMS-error of the reconstruction.

Usage: `./mocutsvd_test <rows> <columns>`
For a list of all options, run `./mocutsvd_test` without arguments.

For a quick example-based introduction, use [mocutsvd_example](#quick-example).

### Quick Example

```bash
git clone https://github.com/johsteffens/mocutsvd
cd mocutsvd
make mocutsvd_example

# Runs MocUT SVD on a (5x8) random matrix. Outputs all Matrix values to stdout.
./mocutsvd_example 5 8 
```
[`example.c`](example.c) contains a short example application intended to quickly learn the mocut-matrix-format and mocut-svd usage.
It outputs all matrix values to stdout. Use this program for smaller matrices.
To test the svd-performance on large matrices, use [mocutsvd_test](#quick-test).

___________________________________



## Detailed Description

Singular Value Decomposition is the method of finding the three components $U, \Sigma, V$ to an arbitrary (m x n) matrix $M$, which satisfy the equation $M = U \cdot \Sigma \cdot V^\ast$ . Where $\Sigma$ is a diagonal matrix containing singular values and $U$, $V$ are unitary matrices containing left and right singular vectors. The SVD algorithm is of critical importance in science and engineering.

MocUT SVD is an algorithm for singular value decomposition.  Given a matrix $M$, it calculates the matrices $U^\ast, \Sigma, V^\ast$ . It is derived from the Golub-Kahan-Reinsch approach and inherits its proven stability. It was developed from ground up with many performance-critical improvements.

#### About the Name

**MocUT** is a shortcut for Monoclinic Unitary Transformation, which is a special kind of recurring transformation pattern I designed for this SVD solution. More details can be found in this whitepaper: [MocUT SVD: Singular Value Decomposition via Monoclinic Unitary Transformations](doc/mocutsvd.md).

#### Platform Support

MocUT SVD is designed to utilize modern CPU architectures. Specifically: Multiple-cores, multi-layered caching, vectorization and hyper-threading.

At the same time the code maintains high portability: It only requires compliance to the C11 (or later) standard. Hence, the code is compilable by todays compiler-tool-chains for most operating systems and platforms.

Portability is achieved by utilizing [generic coding paradigms](doc/true_scalability.md) that help the compiler applying platform specific optimizations. Outer parallelity is achieved via [Open MP](https://en.wikipedia.org/wiki/OpenMP).

#### Transposed $U$, $V$

MocUT SVD generates the matrices of singular vectors in their transposed form: $U^\ast$, $V^\ast$, where singular vectors are row-vectors. If desired, you can convert the matrix back to the traditional form via function ```mocut_mat_s_copy_transposed```.

### Error Handling

An error condition is normally communicated by an integer error code as return value of a function. Returning `0` means `No Error`. Other conditions could be
* Incorrect usage of parameters. 
* Not enough memory to allocate a matrix.
* (Warning) SVD function did not converge.

The following function converts an error value to a descriptive text.

*  **```const char* mocut_err_text( int err_code )```**
  * Converts an error code to a descriptive 0-terminated text string.

### Matrix-Interface

The matrix is represented by structure ```mocut_mat_s```, which contains the following user-readable elements:

* ```size_t rows``` : Number of rows of the matrix
* ```size_t cols``` : Number of columns of the matrix
* ```size_t stride``` : Number of elements between the start of two successive rows. It is ```stride >= cols```.
* ```double* data```: Location of matrix data in memory.

The matrix uses a 'strided row-major' data-layout. This means that the element ```[i][j]``` is accessed as ```data[ i * stride + j ]```;  $i \in \{ 0, ..., \text{rows}-1 \},  j \in \{ 0, ..., \text{cols}-1 \}$.

The stride value is set automatically for optimal alignment during matrix-allocation. You may also setup a matrix manually, referencing external data via `matrix-setup` function.

#### Matrix-Functions

*  **```mocut_mat_s* mocut_mat_s_create( void );```**
   * Constructs an empty instance. For destruction, see ```mocut_mat_s_discard```.
   
   * **Return:** 
      * Pointer to constructed matrix.
      * NULL in case matrix could not be constructed.
   
* **```void mocut_mat_s_discard( mocut_mat_s* );```**
   * Destructs a matrix. Frees all memory this matrix owns. For construction, see ```mocut_mat_s_create```.
  
*  **```void mocut_mat_s_init( mocut_mat_s*  );```**
   * Initializes an empty matrix instance that was created on the stack.
   * To tear it down, call `mocut_mat_s_down`.
   * Do not use this function when you create a matrix via `mocut_mat_s_create`.
  
* **```void mocut_mat_s_down( mocut_mat_s* );```**
   * Tears down a matrix that was initialized by `mocut_mat_s_init`. 
   * Do not use this function when the matrix was created via `mocut_mat_s_create`.

* **```int mocut_mat_s_alloc( mocut_mat_s* o, size_t rows, size_t cols );```**
   * Allocates a (rows x cols )-matrix and initializes all values to zero.
   * This function takes care of proper data-alignment for optimal SVD-performance.
   * **Return:** 
      * `0`: Success    
      * `>0`: [Error Code](#error-handling) (Allocation failed).
   
* **```int mocut_mat_s_set( mocut_mat_s* o, size_t row, size_t col, double value );```**
   * Sets element value at position `[row][col]`.
   * **Return:**       
      * `0`: Success
      * `>0`: [Error Code](#error-handling) (Index out of rage).

* **```double mocut_mat_s_get( const mocut_mat_s* o, size_t row, size_t col );```**
   * Returns element value at position `[row][col]`. If `row` or `col` is out of range, `0` is returned.
   
* **```double* mocut_mat_s_ptr( const mocut_mat_s* o, size_t row, size_t col );```**
   * Returns a pointer to element value at position `[row][col]` for fast, compiler-optimizable read or write access. 
   * Same as: ```o->data[ i * o->stride + j ]```
   * **Note:** This function has no boundary check for `row`, `col`!

* **```void mocut_mat_s_clear( mocut_mat_s* o );```**
   * Clears memory from the latest allocation without destroying the matrix.
   
* **```int mocut_mat_s_setup( mocut_mat_s* o, size_t rows, size_t cols, size_t stride, double* data )```**
   * Makes the matrix use external data. The ABI of the external matrix must be compatible (strided row-major layout). The external matrix data is neither copied nor owned by `mocut_mat_s`. It will not be freed by `mocut_mat_s_discard`. The external data must stay alive and valid during the lifetime of the `mocut_mat_s` instance.
   * The purpose of this function is to simplify the integration of `mocutsvd` into a codebase that uses its own matrix representation.
   * **Note:** Using `mocut_mat_s_setup` with unaligned external data (`stride=cols`) is allowed. Unaligned rows can reduce SVD efficiency slightly.
   * **Return:** 
      * `0`: Success
      * `>0`: [Error Code](#error-handling)

* **```mocut_mat_s* mocut_mat_s_create_alloc( size_t rows, size_t cols )```**  
   * Convenience function: Combination of `mocut_mat_s_create` and `mocut_mat_s_alloc`. Returns NULL in case of error.
  
* **```mocut_mat_s* mocut_mat_s_create_setup( size_t rows, size_t cols, size_t stride, double* data )```**
   * Convenience function: Combination of `mocut_mat_s_create` and `mocut_mat_s_setup`. Returns NULL in case of error.
  
* **```int mocut_mat_s_copy( mocut_mat_s* o, const mocut_mat_s* m )```**  
   * Copies the matrix data from `m` to `o`. 
   * Both matrices must be allocated to the same size: `(o.rows==m.rows) && (o.cols==m.cols)`
   * **Return:** 
      * `0`: Success
      * `>0`: [Error Code](#error-handling)
  
* **```int mocut_mat_s_copy_transposed( mocut_mat_s* o, const mocut_mat_s* m )```**
   * Copies the transposed matrix data from `m` to `o`. 
   * Both matrices must be allocated to the respective transposed size: `(o.rows==m.cols) && (o.cols==m.rows)`
   * In case of a square matrix `(rows == cols)`, `o` and `m` may reference the same matrix: In-place transposition.
   * **Return:** 
      * `0`: Success
      * `>0`: [Error Code](#error-handling)

### SVD-Function

* **```int mocut_svd( mocut_mat_s* a, mocut_mat_s* u, mocut_mat_s* v );```**

```mocut_svd``` performs the thin SVD on a (m x n)-Matrix: $M \rightarrow U^\ast, \Sigma, V^\ast$ .

The matrices are being modified during execution. Matrix `a` must be initialized as $M$ before execution. After execution it represents $\Sigma$ : The diagonal elements represents the singular values, all other elements are set to zero.

Arguments `u`, `v` are optional. They represent $U^\ast$ and $V^\ast$ respectively containing the singular vectors as row vectors. Pass `NULL` as argument when not needed.

If $U^\ast$ or $V^\ast$ is needed, you can pass either an empty instance, or you can pass a pre-sized matrix (via `matrix-alloc` or `matrix-setup`). 

If it is empty, `mocut_svd` will allocate it to the correct size. If it is non-empty, the size must be `(k x m)` for `u` and `(k x n)` for `v` with `k = min(m,n)`.

#### Convergence

The final phase of the SVD is an iterative process that converges to the correct result. In the rare case of non-convrgence, SVD returns the warning `MOCUT_WRN_CONVERGENCE`. It indicates that the reconstruction  $U \cdot \Sigma \cdot V^\ast$ might deviate numerically from $M$.

#### Return Value

* `0`: Success
* MOCUT_WRN_CONVERGENCE: Insufficient convergence
* `Other >0`: [Error Code](#error-handling) 

## Thread Safety

MocUT-functions are thread-safe, provided the data passed as argument is not shared with other threads while the function is running.

Although `mocut_svd` may spawn its own threads, those are completely shielded from any multi-threaded environment you might be using in your application.


## Memory and Data Alignment
MocUT SVD offers some customization around memory handling.

Function `mocut_svd` uses only the memory space provided by the matrices. Per default the matrix interface uses `stdlib` functions `aligned_alloc` and `free` for allocation and destruction. Alternatively, you can declare your custom memory functions or supply external memory to a matrix via function `mocut_mat_s_setup`. You can also control memory alignment.

This opens possibilities for environments with limited or non-standard memory management:

### Custom Memory Management

If you wish to use a specific memory manager, define macros `MOCUT_MEM_ALLOC( size )` and `MOCUT_MEM_FREE( data )` before including `mocutsvd.h`. 

MOCUT_MEM_ALLOC( size ) should represent a dynamic allocation  function similar to stdlib's `malloc`.

MOCUT_MEM_FREE( data ) should represent a memory release function similar to stdlib's `free`.

In the example below, the memory manager [TBMAN](https://github.com/johsteffens/tbman) would be used:

``` C
#define MOCUT_MEM_ALLOC( size ) tbman_malloc( size )
#define MOCUT_MEM_FREE( data ) tbman_free( data )
#include "mocutsvd.h"
```

### No Memory Management

If you wish to prevent MocUT SVD from using any memory management at all, define `MOCUT_NO_MEM_ALLOC` before including `mocutsvd.h`.

``` C
#define MOCUT_NO_MEM_ALLOC
#include "mocutsvd.h"
```

In this case you cannot use functions `mocut_mat_s_create`, `mocut_mat_s_discard`, `mocut_mat_s_alloc`.

Instead, place the matrix instance on the stack and use functions `mocut_mat_s_init`, `mocut_mat_down`, `mocut_mat_setup`.

``` C
mocut_mat_s a; // instance 'a' placed on the stack
mocut_mat_s_init( &a ); // instance 'a' is initialized
mocut_mat_s_setup( &a, rows, cols, stride, data ); // assigning an external matrix data area to 'a'

... // working with 'a'

mocut_mat_s_down( &a ); // instance 'a' is cleaned up
```

### Custom Data Alignment

MocUT SVD [aligns](doc/true_scalability.md#data-alignment) matrix data to improve on inner parallelity and cache usage for most processors. If you wish to experiment with your own custom alignment, define MOCUT_VAL_ALIGN with a constant indicating the alignment in bytes before including `mucutsvd.h`:

``` C
#define MOCUT_VAL_ALIGN 32 // aligns matrix rows to multiple of 32 double values (32*8 bytes)
#include "mocutsvd.h"
```

## Side Effects and Remedies

#### Threads, Open MP
By default MocUT SVD spawns multiple threads according to the number of logical CPUs on the platform. 

* If you wish to run it only in the caller thread, avoid using the flag `-fopenmp` with compiler and linker. 

* If you want to prevent MocUT SVD using OpenMP at all, comment out the compiler directives​ `#pragma omp parallel for` in `mocutsvd.h` and `mocutsvd.c`

* If you wish more specific thread controls, look up the Open MP documentation: [https://www.openmp.org](https://www.openmp.org)

#### High CPU Load
On very large matrices, the function `mocut_svd` will put prolonged load on the CPU, running it near its rated power limits. Correctly configured machines can handle this type of load. Some platforms permit overriding safety limits at the users own risk (e.g. via overclocking, overvolting). Be aware that inadequately configured safety measures might cause malfunction or damage under high CPU load.

#### Valgrind
Advanced debugging tools like `valgrind` analyze the instructions and memory usage of a program:

* Valgrind might not recognize newer native CPU instructions, which occur with certain native compiler optimizations (such as `-march=native` ).

* Valgrind might not be compatible to a compilers integration of OpenMP  (`-fopenmp`).

To analyze code with `valgrind`, do not use these compiler flags.

## Performance Charts

![](doc/image/duration_tests.png)

Absolute processing time in seconds for a full decomposition of a ( $n$ x $n$ ) matrix: $M \rightarrow U^\ast, \Sigma, V^\ast$.

The following CPUs were used:
* `TR 7960`: AMD Ryzen™ Threadripper™ 7960X, containing 24 cores
* `RZ 7950`: AMD Ryzen™ 9 7950X, containing 16 cores
* `RZ 5800`: AMD Ryzen™ 7 5800X, containing 8 cores
* `i7 6700`: Intel® Core™ i7-6700, containing 4 cores

