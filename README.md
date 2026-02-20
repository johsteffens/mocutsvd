# MOCUT SVD - Singular Value Decomposition
# What it is

MOCUT SVD is an efficient and easy-to-use stand-alone-implementation of the thin Singular Value Decomposition for C or C++ programs.

It implements a specific SVD solution designed to run efficiently and mostly in parallel on modern CPUs. It makes minimal assumptions about the underlying hardware and runtime environment, thus providing a highly portable solution.

# How to use it

* Include the header `mocutsvd.h` in your program (C or C++ code) and compile `mocutsvd.c` among your sources. 

* Compile with following compiler and/or linker flags (or use compatible flags of your favorite Toolchain).  
  * Required: `-lm`. Optional for maximum performance: `-O3`, `-fopenmp` and `-march=native`. 

**Example:** ` gcc -o example example.c mocutsvd.c -fopenmp -march=native -O3 -lm`

* **Matrix ABI:**

A matrix is represented by a simple array of `double` precision values. The data arrangement is **row-major**. All rows of a matrix are placed subsequently in contiguous memory without gaps. The memory of all desired matrices must be properly allocated. 

* **SVD Function:**


```C
int mocut_svd( size_t m /*rows*/, size_t n /*cols*/, double* a, double* u, double* v );
```

It computes from the (m x n) matrix $M$ its singular values and optionally matrices $U^\ast$ and $V^\ast$ of singular row-vectors. The decomposition satisfies the equation $M = U \cdot \Sigma \cdot V^\ast$ .

Function `mocut_svd` operates within the memory space given by `a`, `u` and `v`. It does not change any allocation - only the data. Before execution, variable `a` points to the input matrix $M$; after execution it represents a vector of length `k` of singular values. If so desired, variables `u` and `v` represent matrices $U^\ast$ and $V^\ast$. (Note that both are thin and transposed.) Use `NULL` as argument for any unwanted matrix.

* **Thin SVD:**

MOCUT SVD computes the `k = min(m,n)` singular values and associated singular vectors according to the 🔗 [Thin SVD](https://en.wikipedia.org/wiki/Singular_value_decomposition#Reduced_SVDs). With $M$ being a `(m x n)` matrix, $U^\ast$ must be pre-allocated to a `(k x m)` matrix and $V^\ast$ must be pre-allocated to a `(k x n)` matrix.

# Benefits

* **Stable**: MOCUT SVD has the convergence-stability of Golub-Kahan-Reinsch SVD or better.
* **Fast and True-Scalable**: Cache-efficient parallel computation for small up to very large matrices.
* **Memory Efficient**: Operates within the user-provided memory space and allocates no extra heap memory.
* **Parallel**: Time-critical computation is (optionally) parallelized. By default all available CPU-cores and CPU-threads are used.

