# True Scalability: Paradigms for efficient generic algorithms in numerical analysis.

**Johannes B. Steffens**

#### Technical Whitepaper, March 2026



## Introduction

Numerical analysis is a fundamental discipline in sciences and engineering. The technical requirements for an efficient algorithm have shifted over the past decades. In early days the effort to perform elementary floating point operations by far dominated most other aspects of computation. It inspired the [big O notation](https://en.wikipedia.org/wiki/Big_O_notation) , which represents the asymptotic numerical complexity for one or more parameters. This notation is still used as indicator for the scalability of an algorithmic solution to a numeric problem.

Over time, the situation has shifted: Workstations, servers and even consumer- or embedded computing devices have one or more CPUs, each with with multiple independent processing cores. Each core typically offers a dedicated high-speed vectorized FPU, hyper-threading and multi-level caching. Processing speed today significantly depends on how well the algorithm is adapted to those architectural features. This specifically applies to [data-locality](https://en.wikipedia.org/wiki/Locality_of_reference) and [parallelity](https://en.wikipedia.org/wiki/Parallel_algorithm). While the traditional way of assessing numerical efficiency remains important, it may rank after above paradigms when it comes to assessing the expected computation time.

Hardware-specific optimization or adaptation optimal efficiency is less popular because it thwarts a developer's aim toward general purpose and platform-agnostic code. Luckily, the ability of compilers to generate optimized machine code has improved significantly as well. With [Open MP](https://www.openmp.org/), a standard was formed for parallelization, which has been adopted by many compiler tool-chains. 

Platform agnosticism and hardware specific efficiency need no longer be a chasm of a contradiction. It is bridged by a set of coding paradigms with which efficient, general purpose and future-proof code can be created.

Intuitively, we relate the execution-time of an algorithm solving a numeric problem to $\text{ numerical complexity } \over \text{ computational power }$. "*Computational power*" would be a function of : Number of CPU (-cores), CPU clock frequency, cache size & speed, RAM speed, etc. 

Therefore, in this document, I'd like to coin the term ***True-Scalable***, by describing five most important generic coding paradigms it involves. These are:

* **Data-Layout**
* **Inner Parallelity**
* **Outer Parallelity**
* **Data-Locality** 
* **Data-Alignment**

## Baseline

An algorithm shall be **true-scalable** when it is hardware-agnostic and its computation time is a linear function of  $\text{ numerical complexity } \over \text{ computational power }$. 

I'll describe a set of important generic paradigms to achieve true scalability. I use the simple multiplication of two n x n matrices as baseline example to illustrate how each paradigm is employed. Note that this baseline is not to be understood as generic recipe for improving an arbitrary algorithm in general. It is intended to demonstrate what each paradigm actually aims at and how the algorithm need to be restructured in order to approach the goal.

Let $a, b, c$ be (n x n) matrices. $c$ shall be initialized with zeros. We compute the matrix-matrix product $c = ab$.

<a id="example_1"></a>
**Example 1:** Baseline algorithm (Standard Formula)

``` C
for( int i = 0; i < n; i++ )
	for( int j = 0; j < n; j++ )
		for( int k = 0; k < n; k++ )
            c[i][j] += a[i][k] * b[k][j];
```

Although matrix-matrix multiplication has a numeric complexity of $O(n^3)$,  The time scalability of this textbook implementation is $O(n^x)$ with $x$ often significantly larger than $3$. The main reason is that the innermost loops access large sections and distant elements of memory. With small $n$, the entire matrix fits in low level cache and the algorithm runs fast, as $n$ rises, cache misses increase throughout all cache levels.

## Data Layout

We assume that the matrix is located in contiguous memory.  For $c_{ij}$ we used the generic notation`c[i][j]`, which in a strict sense is only correct C-code when the matrix has a fixed size and is located on the stack. The dynamic solution places the matrix data on the heap as contiguous one dimensional array. 

The [row-major](https://en.wikipedia.org/wiki/Row-_and_column-major_order) ordering, elements in a row have the same order in memory. All rows are stored as [strided array](https://en.wikipedia.org/wiki/Stride_of_an_array) with a fixed distance, called `stride`. It is  `stride >= columns`. In this ordering is $c_{ij}$ = `c[ i * stride + j ]`. For tightly packed matrices, `stride` would be equal to `columns`. We will show in section [Data Alignment](#data-alignment) that `stride` should be a multiple of a fixed block size.

The [column-major](https://en.wikipedia.org/wiki/Row-_and_column-major_order) ordering is the transposed counterpart with $c_{ij}$ = `c[ j * stride + i ]`.

There are other ordering methods but these two are most commonly used.

Modern [DRAM](https://en.wikipedia.org/wiki/Dynamic_random-access_memory) is accessed in [bursts](https://en.wikipedia.org/wiki/Burst_mode_(computing)) where multiple adjacent values are loaded at once into the cache. This happens even if a program needs only one or a few of them. It is therefore a good strategy to design the code such that the the data within a burst is mostly used such that the number of individual bursts is kept minimal. [Example 1](#example_1) should run faster with a row-major-order, because for two out of three matrices, rows are accesses sequentially in inner loops.

Going forward, we will assume that all matrices have row-major order.

## Inner Parallelity

The inner parallelity of a program addresses the synchronous parallel components of a processor. These are typically SIMD units. They use the same (or related) numerical operations but on different data. Key aspects is that each operation is fast (often just one CPU clock cycle). All parallel operations for a given instruction finish at the same time and all involved data is updated at the same time. These is no synchronization needed. 

Code that achieves this form of parallelity optimizes its innermost code-sections (typically innermost loops) for inner parallelity.

Modern compilers can detect code with inner parallelity and use the appropriate SIMD instructions automatically. This means that the developer can write readable and portable code and yet utilize the hardware advantages of the CPU when the code was compiled.

Code with inner parallelity is designed such that the most time consuming sections contain elementary operations which order of execution is permutable without affecting the result.

Example 1 has no inner parallelity because the innermost loop is not permutable.

<a id="example_2"></a>
**Example 2: With inner parallelity**

``` C
for( int i = 0; i < n; i++ )
	for( int k = 0; k < n; k++ )  // k-loop and j-loop are swapped
        for( int j = 0; j < n; j++ )  
            c[i][j] += a[i][k] * b[k][j]; 
			// The innermost loop is permutable: 
			// The compiler can optimize it (vectorization)
```
In [example 2](#example_2) the two inner loops are swapped, which makes the inner loop permutable. Many compilers can detect the inner parallelity and apply fast SIMD multiply-accumulate operations. Example 2 should therefor execute faster than example 1 when compiled with optimizations enabled (gcc: ```-O3 -march=native```).

The inner loop in example 2 has also better [data-locality](#data-locality) than example 1. It should therefore be faster even without any vectorization.

## Outer Parallelity

Parallelity can also be achieved by dividing a program into sections that can be executed in dedicated (independent) threads simultaneously. These threads typically run on different CPU cores, thus achieving a speed advantage over executing the program sections sequentially. If there are mutually independent code-sections, each can run in a dedicated thread. In this respect threads provide much flexibility. The drawback is that threads run asynchronously: Processing of thread-results can only proceed when all these threads have finished, hence threads need monitoring. This requires extra overhead and slows down overall processing. Multi-threaded programs run efficiently when the number of synchronization cycles is kept at a minimum. 

Therefore threaded sections are preferably long-running. The subdivision of a program into threads happens ideally in outer layers of the program. In case of nested loops, this would be the outer loop(s). Ensuring independence of code sections is in the providence of the developer. Compilers cannot do it reliably. However, once such sections are declared independent, a modern compiler can take care of parallelizing them. 

A platform agnostic solution is the [Open MP](https://www.openmp.org/) standard. It tells the compiler via [pragma-directive](https://en.wikipedia.org/wiki/Directive_(programming)) that all cycles of the outermost loop are independent and therefore can run in any order. A compiler supporting the Open MP standard will then attempt to distribute all cycles across multiple threads.

<a id="example_3"></a>
**Example 3: With inner and outer parallelity**

``` C
// outermost loop-cycles are independent and therefore can be parallelized
#pragma omp parallel for // this parallelizes the outermost loop
for( int i = 0; i < n; i++ )
	for( int k = 0; k < n; k++ )
        for( int j = 0; j < n; j++ )  
            c[i][j] += a[i][k] * b[k][j];
```

In this case the outermost loop in [example 3](#example_3) needs no adaptation. It is already independent because each cycle modifies a different block of data without overlaps. So, we can add a directive to the compiler to spread these cycles across multiple threads.

The developer must nevertheless carefully design and check his/her code for this property. The compiler will not be able to do it automaticallly, neither will the compiler be able to reliably detect harmful interdependencies.

## Data Locality

Modern CPUs can operate many thousand times faster than their early predecessors. However, memory latency has improved at a slower rate and there are physical limitations with respect to future developments. Information can travel at most at light speed but mass storage occupies space and needs to be accessible by all cores. It therefore cannot be placed in the immediate vicinity of a CPU core. This limitation will stay with us even when considering the ongoing progress in hardware packaging and miniaturization.

To mitigate the latency-problem, processors employ a layered memory structure: Small but fast memory, called cache, near the processing core. More dense but slower memory further away in memory banks as regular RAM. The cache mirrors a small portion of RAM and acts as fast temporary storage. Advanced architectures uses a hierarchy of multiple cache-layers.

Re-using cache and minimizing cache-RAM synchronizations is an important property of modern programs. The underlying design-criterion is called **[Data-Locality](https://en.wikipedia.org/wiki/Locality_of_reference)**. Its importance exceeds that of mere numerical efficiency. Since routines observing data locality can run hundreds of times faster than their inefficient counterparts, using more operations can be acceptable when the change yields a better data-locality.

Example 3 is not data-local because for each of the $n$ outer cycles, the entire matrix c is accessed and the accessed elements in one of the other matrices is spread out in memory  (in case of row-major layout, this would be $a$). As $n$ gets larger, the cache-efficiency is diminishing. Ultimately, the program will run extremely slow despite parallelism.

The principal goal for achieving data-locality is to limit the computational effort to a small data-area as far as possible. In case of matrix-multiplication this is typically done via [block-partitioning](https://en.wikipedia.org/wiki/Block_matrix) and [divide and conquer](https://en.wikipedia.org/wiki/Matrix_multiplication_algorithm#Divide-and-conquer_algorithm). 

<a id="example_4"></a>
**Example 4: With inner and outer parallelity and data-locality**

``` C
int min( int x, int y ) { return x < y ? x : y; }

// nb = block size; choose nb such that 3*nb*nb values fit well into the L1 cache
const int nb = 32; // (nb should be a power of 2) nb = 32 works well for most architectures 

#pragma omp parallel for
for( int ib = 0; ib < n; ib += nb )
    for( int kb = 0; kb < n; kb += nb )
        for( int jb = 0; jb < n; jb += nb )
            for( int i = ib; i < min( n, ib + nb ); i++ )
                for( int k = kb; k < min( n, kb + nb ); k++ )
                    for( int j = jb; j < min( n, jb + nb ); j++ )
                        c[i][j] += a[i][k] * b[k][j];
```

[Example 4](#example_4) shows block-partitioning. In this approach we have concentrated the main workload well inside the inner loops.

## Data Alignment
Hardware design favors address-alignment between cache-lines, vector-registers and DRAM. Therefore reading bulk-data from an aligned DRAM address is typically more efficient. For that reason it is beneficial to combine block-processing with data-alignment.

In case of a matrix representation, a useful strategy is aligning each matrix row accordingly in memory. This can be achieved by using function [`aligned_alloc`](https://en.cppreference.com/w/c/memory/aligned_alloc) for memory allocation and specifying the [`stride`](#data-layout) value as a multiple of the block-size used in data-blocking.

<a id="example_5"></a>
**Example 5: Aligning all matrices**

``` C
// nb = block size; choose nb such that 3*nb*nb values fit well into the L1 cache
const int nb = 32; // (nb should be a power of 2) nb = 32 works well for most architectures 

// stride is set to the smallest multiple of nb larger or equal n
int stride = n + ( ( n % nb ) > 0 ) ? ( nb - ( n % nb ) ) : 0;
double* a = aligned_alloc( nb * sizeof( double ), n * stride );
double* b = aligned_alloc( nb * sizeof( double ), n * stride );
double* c = aligned_alloc( nb * sizeof( double ), n * stride );

.... // fill a,b with data; set c to zero

int min( int x, int y ) { return x < y ? x : y; }

// time critical part ...
#pragma omp parallel for
for( int ib = 0; ib < n; ib += nb )
    for( int kb = 0; kb < n; kb += nb )
        for( int jb = 0; jb < n; jb += nb )
            for( int i = ib; i < min( n, ib + nb ); i++ )
                for( int k = kb; k < min( n, kb + nb ); k++ )
                    for( int j = jb; j < min( n, jb + nb ); j++ )
                        c[ i * stride + j ] += a[ i * stride + k ] * b[ k * stride + j ];
```

In [example 5](#example_5) we spell out the strided access for dynamically allocated matrices. The allocation address and specific `stride` value ensure proper alignment of all rows. With `n` <= `stride` < `n + nb` , only an insignificant amount of memory is wasted.

## Performance

The table below demonstrates the baseline performance of all discussed paradigms.

**Table:** Test on a homogeneous, hyper-threaded 24-core processor. The code was compiled with gcc-options `-O3`, `-march-native` and `-fopenmp`.

| Example                 | time/sec (n=3333) | time/sec (n=7173) |
| ----------------------- | ----------------- | ----------------- |
| 1 (Standard)            | 86.5              | 1133              |
| 2 (+ Inner Parallelity) | 13.2              | 144               |
| 3 (+ Outer Parallelity) | 0.564             | 5.72              |
| 4 (+ Data-Locality)     | 0.277             | 2.55              |
| 5 (+ Alignment)         | 0.184             | 1.65              |

The two values of n are chosen to be coprime to $n_b$, to show the effect of alignment. Their ratio is chosen such that the numerical complexity is around 10: $(7173/3333)^3 \approx 10$, to make the timing values easier comparable.

We can observe that the two data-local solutions show the best time-scaling behavior. The last two appear even slightly better than expected. The reason is that outer parallelity works better on larger matrices.

We can clearly observe the profound incremental effect on performance each paradigm provides.

## Scope and Outlook

#### Other Strategies

* The blocking-size can be further tweaked for better overall performance. In some situation it makes sense to consider using a variable size, dependent on other parameters of the algorithm.

* Hierarchical models of data-locality [1] can yield additional gain. We did not dive into this topic, because of potential trade-offs  with other paradigms.
* This document is intended as summary of important platform-agnostic paradigms. I therefore did not consider strategies, which are tied to more specialized hardware, such as the GPU. 

#### Strassen Algorithm

For sake of completeness, I point out that the [Strassen Algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm) has better numerical complexity of $O( n^{log_2 7})$ than the standard matrix-matrix multiplication. However, it requires more memory and is more difficult to parallelize. Therefore, it does not lend itself as well to a true-scalable optimization as the standard algorithm.

#### Other Computational Problems

The purpose of using the matrix-multiplication example is to illustrate the paradigm and help the reader to generalize across other computational problems. It is not intended to offer the best matrix-multiplication algorithm in circulation.

Finding a true-scalable algorithm can be difficult. Deciding whether such an algorithm exists is an np-hard problem. However, in many specific real-world scenarios these paradigms can help improving the overall computational efficiency.

## Conclusion

I described five of the most important platform-agnostic programing paradigms to achieve true scalability.

I have demonstrated that each has a large effect on computational performance. The differences in computational speed span three orders of magnitude. 

I conclude that mere numerical complexity is no sufficient assessment of expected computation time for an algorithm. 

The baseline algorithm, used here, is rather simple. More complex numerical tasks used in linear algebra, such as matrix decomposition, require significantly more elaborate algorithmic solutions to achieve true scalability. 

## Web References
The following references were used as inline links thoughout the document:

https://en.wikipedia.org/wiki/Big_O_notation

https://en.wikipedia.org/wiki/Locality_of_reference

https://en.wikipedia.org/wiki/Parallel_algorithm

https://www.openmp.org/

https://en.wikipedia.org/wiki/Row-_and_column-major_order

https://en.wikipedia.org/wiki/Stride_of_an_array

https://en.wikipedia.org/wiki/Dynamic_random-access_memory

https://en.wikipedia.org/wiki/Burst_mode_(computing)

https://en.wikipedia.org/wiki/Directive_(programming)

https://en.wikipedia.org/wiki/Block_matrix

https://en.wikipedia.org/wiki/Matrix_multiplication_algorithm#Divide-and-conquer_algorithm

https://en.cppreference.com/w/c/memory/aligned_alloc

https://en.wikipedia.org/wiki/Strassen_algorithm

## Literature

[1] Chengliang Zhang et al, A Hierarchical Model of Data Locality, POPL '06 Proceedings, Pages 16 - 29, https://dl.acm.org/doi/proceedings/10.1145/1111037


____

&copy; 2026 Johannes B. Steffens
