
# MOCUT SVD: Singular Value Decomposition via Monoclinic Unitary Transformations

#### Technical Whitepaper by Johannes B. Steffens

#### February 2026

## 1 Introduction

[Singular Value Decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) (SVD) is a fundamentally important discipline in linear algebra. Fast algorithms that perform the SVD reliably on any matrix are most sought after in sciences and engineering.

The SVD decomposes a *m x n* Matrix *M* into Matrices $ U $, $\Sigma$ , $ V $ such that $\Sigma$ is [diagonal](https://en.wikipedia.org/wiki/Diagonal_matrix) and  $ U $, $ V $ are both [unitary](https://en.wikipedia.org/wiki/Unitary_matrix)  and $M = U \cdot \Sigma \cdot V^\ast$. 

While the SVD is theoretically always possible with any matrix, developing an algorithm that accomplishes this task fast and numerically stable on any matrix is technically challenging.

The technical requirements for an efficient algorithm have shifted over the past decades: In early years, elementary numerical operations on floating point numbers were the most time consuming part. The time needed for fetching data from or storing to memory was of of secondary importance. Also workstations had a single processor (-core), hence multi threading was at that time not seen as a tool for speed improvement but rather a method of convenience for multitasking.

Today, the situation has significantly changed: Workstations, servers and even consumer- or embedded computing devices have one or more CPUs, each with with multiple independent processing cores. Each core typically offers a dedicated high-speed vectorized FPU, hyper-threading and multi-level caching. Processing speed depends on how well the algorithm is adapted to that architecture. Specifically, how [data-local](https://en.wikipedia.org/wiki/Locality_of_reference) and how [parallel](https://en.wikipedia.org/wiki/Parallel_algorithm) it can operate. Numerical efficiency remains important, but under certain circumstances ranks after above two paradigms.

This paper describes a general purpose and highly portable SVD algorithm, which I redesigned with all of above performance in appropriate relationship for modern general purpose CPUs as well as their evolution in the foreseeable future.

The name MOCUT is an acronym for the specific pattern, ordering and partitioning of unitary transformations described in this paper. It gives a performance edge over earlier SVD algorithms.

## 2 True Scalability

Numerical complexity is the amount of elementary numerical operations needed. These are operations like multiplications and additions.

In the early days of computing, the hardware was simple enough to permit assessment of computation time from the numerical complexity of an algorithm. Therefore numerical complexity became the most relevant descriptive property of an algorithm. It inspired the [big O notation](https://en.wikipedia.org/wiki/Big_O_notation) , which represents the asymptotic numerical complexity for one or more parameters. This notation is still used as indicator for the scalability of an algorithmic solution to a numeric problem.

Computing hardware keeps improving over time. Although all programs likely run faster on newer hardware, many no longer scale according to that simple metric given by numerical complexity. Newer hardware specific paradigms need proper consideration as well. On the other hand, developers strive to write platform-independent code to maximize portability. 

Intuitively, we relate the execution-time of an algorithm solving a numeric problem to $ \text{ numerical complexity } \over \text{ computational power } $. "*Computational power*" would be a function of : Number of CPU (-cores), CPU clock frequency, cache size & speed, RAM speed, etc. 

Therefore, I'd like to coin the term ***True-Scalable***:  An algorithm shall be **true-scalable** when it is hardware-agnostic and its computation time is a linear function of  $ \text{ numerical complexity } \over \text{ computational power } $.

This section discusses the most important efficiency paradigms to achieve true scalability. These are: **Data-Layout**, **Parallelity**, **Data-Locality** and **Data-Alignment**.

To illustrate this, we use the multiplication of two n x n matrices as example. Let $ a, b, c $ be (n x n) matrices. $ c $ shall be initialized with zeros. We compute the matrix-matrix product $ c = ab $.

**Example 1**: Standard formula

``` C
for( int i = 0; i < n; i++ )
	for( int j = 0; j < n; j++ )
		for( int k = 0; k < n; k++ )
            c[i][j] += a[i][k] * b[k][j];
```

#### 2.0 Data Layout



#### 2.1 Inner Parallelity

The inner parallelity of a program addresses the synchronous parallel components of a processor. These are typically SIMD units. They use the same (or related) numerical operations but on different data. Key aspects is that each operation is fast (often just one CPU clock cycle). All parallel operations for a given instruction finish at the same time and all involved data is updated at the same time. These is no synchronization needed. 

Code that achieves this form of parallelity optimizes its innermost code-sections (typically innermost loops) for inner parallelity.

Modern compilers can detect code with inner parallelity and use the appropriate SIMD instructions automatically. This means that the developer can write readable and portable code and yet utilize the hardware advantages of the CPU when the code was compiled.

Code with inner parallelity is designed such that the most time consuming sections contain elementary operations which order of execution is permutable without affecting the result.

Example 1 has no inner parallelity because the innermost loop is not permutable.

**Example 2: With inner parallelity**

``` C
for( int i = 0; i < n; i++ )
	for( int k = 0; k < n; k++ )  // k-loop and j-loop are swapped
        for( int j = 0; j < n; j++ )  
            c[i][j] += a[i][k] * b[k][j]; 
			// The innermost loop is permutable: 
			// The compiler can optimize it (vectorization)
```
In example 2 the two inner loops are swapped, which makes the inner loop permutable. Many compilers can detect the inner parallelity and apply fast SIMD multiply-accumulate operations. Example 2 should therefor execute faster than example 1 when compiled with optimizations enabled (gcc: ```-O3 -march=native```).

#### 2.2 Outer Parallelity

Parallelity can also be achieved by dividing a program into sections that can be executed in dedicated (independent) threads simultaneously. These threads typically run on different CPU cores, thus achieving a speed advantage over executing the program sections sequentially. If there are mutually independent code-sections, each can run in a dedicated thread. In this respect threads provide much flexibility. The drawback is that threads run asynchronously: Processing of thread-results can only proceed when all these threads have finished, hence threads need monitoring. This requires extra overhead and slows down overall processing. Multi-threaded programs run efficiently when the number of synchronization cycles is kept at a minimum. 

Therefore threaded sections are preferably long-running. The subdivision of a program into threads happens ideally in outer layers of the program. In case of nested loops, this would be the outer loop(s). Ensuring independence of code sections is in the providence of the developer. Compilers cannot do it reliably. However, once such sections are declared independent, a modern compiler can take care of parallelizing them. 

A platform agnostic solution is the [Open MP](https://www.openmp.org/) standard. It tells the compiler via [pragma-directive](https://en.wikipedia.org/wiki/Directive_(programming)) that all cycles of the outermost loop are independent and therefore can run in any order. A compiler supporting the Open MP standard will then attempt to distribute all cycles across multiple threads.

**Example 3: With inner and outer parallelity**

``` C
// outermost loop-cycles are independent and therefore can be parallelized
#pragma omp parallel for // this parallelizes the outermost loop
for( int i = 0; i < n; i++ )
	for( int k = 0; k < n; k++ )
        for( int j = 0; j < n; j++ )  
            c[i][j] += a[i][k] * b[k][j];
```

The outermost loop from example 2 is already independent because the data being modified in each cycle is not shared across cycles. Read-only data can be shared across cycles. Nevertheless the developer must carefully check this before allowing outer parallelity.

### 2.3 Data Locality

Modern CPUs can operate many thousand times faster than their early predecessors. However, memory latency has improved at a slower rate and there are physical limitations with respect to future developments. Information can travel at most at light speed but mass storage occupies space and needs to be accessible by all cores. It therefore cannot be placed in the immediate vicinity of a CPU core. This limitation will stay with us even when considering the ongoing progress in hardware packaging and miniaturization.

To mitigate the latency-problem, processors employ a layered memory structure: Small but fast memory, called cache, near the processing core. More dense but slower memory further away in memory banks as regular RAM. The cache mirrors a small portion of RAM and acts as fast temporary storage. Advanced architectures uses a hierarchy of multiple cache-layers.

Re-using cache and minimizing cache-RAM synchronizations is an important property of modern programs. The underlying design-criterion is called **[Data-Locality](https://en.wikipedia.org/wiki/Locality_of_reference)**. Its importance exceeds that of mere numerical efficiency. Since routines observing data locality can run hundreds of times faster than their inefficient counterparts, using more operations can be acceptable when the change yields a better data-locality.

Example 3 is not data-local because for each of the $ n $ outer cycles, the entire matrix c is accessed and the accessed elements in one of the other matrices is spread out in memory  (in case of row-major layout, this would be $ a $). As $ n $ gets larger, the cache-efficiency is diminishing. Ultimately, the program will run extremely slow despite parallelism.

The principal goal for achieving data-locality is to limit the computational effort to a small data-area as far as possible. In case of matrix-multiplication this is typically done via [block-partitioning](https://en.wikipedia.org/wiki/Block_matrix) and [divide and conquer](https://en.wikipedia.org/wiki/Matrix_multiplication_algorithm#Divide-and-conquer_algorithm). 

**Example 4: With inner and outer parallelity and data-locality**

``` C
int min( int x, int y ) { return x < y ? x : y; }

// nb = block size; choose nb such that 3*nb*nb values fit well into the L1 cache
const int nb = 32; // nb = 32 works well for most architectures

#pragma omp parallel for
for( int ib = 0; ib < n; ib += nb )
    for( int kb = 0; kb < n; kb += nb )
        for( int jb = 0; jb < n; jb += nb )
            for( int i = ib; i < min( n, ib + nb ); i++ )
                for( int k = kb; k < min( n, kb + nb ); k++ )
                    for( int j = jb; j < min( n, jb + nb ); j++ )
                        c[i][j] += a[i][k] * b[k][j];
```

Example 4 shows block-partitioning. In this approach we have concentrated the main workload well inside the inner loops.

### 2.4 Data Alignment

### 2.4 Conclusion

Although Numerical complexity is the most obvious paradigm for computational efficiency. It is sometimes necessary to divert from minimizing numerical complexity in favor of the other paradigms.

The standard matrix multiplication has not the best numerical complexity. The [Strassen Algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm) has a better numerical complexity $ O( n^{log_2 7}) $, but it requires more memory and is more difficult to parallelize. It might therefore have even worse time-complexity for practical range of $ n $.





**Table:** Test on a 16-core $\text{Ryzen}^{ \text{TM}} $ 9 7950x platform for examples 1, ..., 4; $ n $ = 1k and $ n $ = 10k.

| Example                 | time (n=1000) | time (n=10000) | time (n=9999) align=1 | time (n=9999) align=16 |
| ----------------------- | ------------- | -------------- | --------------------- | ---------------------- |
| 1 (Standard)            | 85 ms         | 805 s          |                       |                        |
| 2 (+ Inner Parallelity) | 68 ms         | 233 s          |                       |                        |
| 3 (+ Outer Parallelity) | 9.4 ms        | 22 s           |                       |                        |
| 4 (+ Data-Locality)     | 7.8 ms        | 10 s           |                       |                        |


By that definition, **example 4** would be close enough to be called true-scalable while **examples 1, 2 and 3** are not. 


## 3 Singular Value Decomposition

Given a mathematical Problem: Whether a true-scalable algorithm exists is not always decidable. Finding a true-scalable algorithm can be significantly more difficult than finding just any algorithmic solution.

In the following, I'd like to present a true-scalable algorithm for the Singular Value Decomposition.

### 3.1 Incremental Unitary Transformations

Lets express matrix $ M $ in by a generic decomposition $ U $, $ A $, $ V $ of which $ U $, $ V $ are unitary. $ A $ can be any matrix. 

(3.1) 	$M = UAV^*$

A trivial decomposition would be: $ U = V = \underline{1} $ and $ A = M $.

The product of a unitary matrix with its adjunct $(P^\ast P)$ is the unity and the product of two unitary matrices is unitary. This is used to convert one decomposition into another:

(3.2) 	$ UAV^\ast = U(P^\ast P)A(Q^\ast Q) V^\ast = (UP^\ast)(PAQ^\ast)(QV^\ast) = U_{new} A_{new} V_{new}^\ast$

Equation (3.2) represents an incremental step.

Most SVD algorithms diagonalize A incrementally. First by starting with the trivial decomposition and then fining a suitable sequence $ P_i $, $ Q_j $. 

There are many possible different sequences accomplishing this goal. An algorithm could be sufficiently described by the specific sequence it uses.

There are two classes of incremental unitary transformations, which are typically used: The Householder Reflection and the Givens Rotation.
#### 3.1.1 Householder Reflection (HR)

The [Householder Reflection](https://en.wikipedia.org/wiki/Householder_transformation) is determined by a normalized vector $ w $ and expressed in Matrix form as follows:

$ H_w = \underline{1} - 2 ww^* $; $ w^*w = 1 $

It can be easily verified that $ H_w $ is self-adjoint. It is also unitary because

$ H_w^*H_w = (\underline{1} - 2 ww^*)(\underline{1} - 2 ww^*) = \underline{1} - 4ww^* + 4ww^*ww^* = \underline{1} $

The HR is numerically efficient because applied to a vector $ v $ it can be implemented as $ H_w(v) = \underline{1}-2w(w^*v) $, which has a complexity of $ O(n) $ ($ n = dim( v ) $). 

The HR can be configured to set $ n-1 $ values to zero of a specified vector. This is typically a column or row of a matrix or a smaller section thereof.

#### 3.1.2 Givens Rotation (GR)

The [Givens Rotation](https://en.wikipedia.org/wiki/Givens_rotation) was the method of choice before the householder reflection became more popular. It is a 2D operation ($ n=2 $) and is typically used to set 1 value in a 2-vector to zero. To zero $ n-1 $ elements in an n-vector, one chains $ n-1 $ rotations together. The order of complexity $ O(n) $ is the same as with  the HR. But the actual number of operations needed is by a factor 1.5 higher for $ n \ggg 2 $. For $ n = 2 $ both GR and HR require the same computational effort.

### 3.2 The Golub-Reinsch Algorithm

The Golub-Reinsch Algorithm is a classic and popular SVD algorithm. It consists of two phases:

1. **Bi-Diagonalizing $ A $ via alternating left and right HR.**
2. **Diagonalizing $ A $ via alternating left and right GR.**

Bi-Diagonalizing means: Zeroing all elements in $ A $ except the main diagonal and one immediate sub-diagonal of $ A $.  This is often the upper sub-diagonal. 

Phase 1 zeros alternatingly the leftmost non-zero column under the main-diagonal via left UT $ P_i $ , then the uppermost non-zero row right from the sub-diagonal via $ Q_i $. This is done via householder transformation (s. Figure). 

$ P_i $ and $ Q_i $ are tightly coupled: To determine $ P_i $, the $ i $-th column of $ (A_{i-1}Q^*_{i-1}) $ must be known. To determine $ Q_i $, the $ i $-th row of  $ (P_i A_{i-1}) $ must be known. Consequently all of the residual not yet bi-diagonalized portion of A must be accessed before the next UT can be computed. This thwarts data-locality. The outermost loop is fairly long and not independent. This severely limits the effectiveness of outer parallelity because threads need frequent synchronizations.

Hence, phase 1 in the Golub-Reinsch Algorithm is not true-scalable.

Phase 2 uses left an right Givens Rotations to gradually eliminate sub-diagonal elements. I will show further down that phase 2 can be made true-scalable by using the MOCUT Approach.

### 3.3 The Band-Diagonal Approach

The band-diagonal approach is a method to overcome the limitations of the bi-diagonalization phase by splitting it into two stages**:**

1. **Band-Diagonalizing A**.
2. **Bi-Diagonalizing a band-diagonal A**.

#### 3.3.1 Band-Diagonalizing

Band-Diagonalizing means: Zeroing all elements in $ A $ except the main diagonal and a band of $ n_b $ immediate sub-diagonals of $ A $. This is done by alternating zeroing a block of $ n_b $ left columns and $ n_b $ upper rows. (s. Figure)

Within a single block $ P_i $ and $ Q_i $ are decoupled: To compute $ P_i $ only those rows of $ A $ need be known upfront, which are to be zeroed. In a transposed manner the same applies to $ Q_i $. 

Computing $ P_i $, $ Q_i $ has better data-locality. An accumulated bundle of one-sided householder reflections can be converted into a data-local matrix-matrix multiplication with good outer parallelity. The method is known as the WY-representation of accumulated Householder Reflections [2], [3]. 

I will show further down that the MOCUT Approach achieves true-scalability without requiring the WY-representation.

#### 3.3.2 Bi-Diagonalizing

Bi-diagonalizing a band-diagonal $ A $ is done by 









Hence, at any given state of the decomposition process, we can describe the matrices involved by using the product of unitary matrices $ P_i, Q_j $:

(1.2)	$ U_k = \prod_{i=1}^k P_i $

(1.3)	$ V_l = \prod_{j=1}^l Q_j $

(1.4)	$ A_{kl} = U_k M P_l^* $

We see the process as an iterative operation in which the physical representation of matrices $ A $, $ U $, $ V $ gradually change. $ A $ starts as $ M $ and ends as $\Sigma$. U, V start as unity ($E$) and end as matrix of singular vectors.

There are many different possible sequences of $ P_i $, $ Q_j $ achieving the final decomposition. Normally, they are chosen to set one or more elements in A to zero. For most elements, this can be done in a closed form (requiring a predictably finite number of operations).  The remainder is zeroed in an iterative fashion that converges to the desired result. 





## 2.2 Numerically efficient UT

There are two classes of UT, which are typically used. 

#### 2.2.1 Householder Reflection (HR)

The [Householder Reflection](https://en.wikipedia.org/wiki/Householder_transformation) is determined by a normalized vector $ w $ and expressed in Matrix form as follows:

$ H_w = \underline{1} - 2 ww^* $; $ w^*w = 1 $

It can be easily verified that $ H_w $ is self-adjoint. It is also unitary because

$ H_w^*H_w = (\underline{1} - 2 ww^*)(\underline{1} - 2 ww^*) = \underline{1} - 4ww^* + 4ww^*ww^* = \underline{1} $

The HR is numerically efficient because applied to a vector $ v $ it can be implemented as $ H_w(v) = \underline{1}-2w(w^*v) $, which has a complexity of $ O(n) $ ($ n = dim( v ) $). 

The HR can be configured to set $ n-1 $ values to zero of a specified vector. This is typically a column or row of a matrix or a smaller section thereof.

#### 2.2.2 Givens Rotation (GR)

The [Givens Rotation](https://en.wikipedia.org/wiki/Givens_rotation) was the method of choice before the householder reflection became more popular. It is a 2D operation ($ n=2 $) and is typically used to set 1 value in a 2-vector to zero. To zero $ n-1 $ elements in an n-vector, one chains $ n-1 $ rotations together. The order of complexity $ O(n) $ is the same as with  the HR. But the actual number of operations needed is by a factor 1.5 higher for $ n \ggg 2 $. For $ n = 2 $ both GR and HR require the same computational effort.

### 3 The Golub Reinsch SVD (GR_SVD) 

The [Golub-Reinsch-SVD](https://people.inf.ethz.ch/gander/talks/Vortrag2022.pdf) is today still the most widely used SVD algorithm. It consists of two phases:

1. **Bi-Diagonalizing $ A $ via alternating left and right HR.**
2. **Diagonalizing $ A $ via alternating left and right GR.**

Bi-Diagonalizing means: Zeroing all elements in $ A $ except the main diagonal and one immediate sub-diagonal of $ A $ (often the upper sub-diagonal). 





## 5 Locality of Matrix Operations

We have seen in Eq. (1) ... (4) that SVD is principally an algorithm that iteratively updates a matrix representation. The locality premise compels us to conceive a sequence of unitary operations that preserve locality in the operands.

To illustrate the problem, we fist describe phase 1 of the still widely use GR-SVD, which does not preserve locality. Without limiting generality, we consider a square ($ n $ x $ n $) matrix $ M $. Operations $ P_i $, $ Q_j $ are applied alternatingly ($P_1, Q_1, P_2, Q_2, ... $ ) . $ P_i $ introduces zeros in row $ i $: Elements $ \{(i+1), ..., n\} $. $ Q_j $ introduces zeros in column $ j $: Elements $ \{j+2, ..., n\} $. To determine either unitary transformation, the corresponding row or column of A must be known. That, however, requires the previous unitary transformation to have been executed. (See Figure).

Consequently before $ P_i $ or $ Q_j $ can be determined, the previous operation must be completed on $ A $. Since each changes a large portion of A, each operation likely thrashes the fastest cache layer(s). This operation is not data local and therefore executes slowly on large $ A $. 

A helpful intermediate stage is band-diagonalization: Here a band on $ n_b $ sub-diagonals are left non-zero.

So, effectively we get a 3-Phase SVD:

1. Band-diagonalizing A with a specified off-diagonal band-width $ n_b $.
2. Band- to Bi-Diagonalizing A
3. Diagonalizing A

This approach helps allows bundling left and right transformations in phase 1 (e.g. for $ n_b = 4 $:  $P_1, P_2, P_3, P_4, Q_1, Q_2, Q_3, Q_4, P_5, P_6, ... $   ). 

Let's call the bundles $ P_{b} $ and $ Q_{b} $ :

(4.1) 	$ P_{b_a} = \prod_{l = a}^{a + n_b} P_l $, $ a \in \{ 1, n_b, 2n_b, ... \}$

(4.2) 	$ Q_{b_a} = \prod_{l = a}^{a + n_b} Q_l $, $ a \in \{ 1, n_b, 2n_b, ... \}$

An accumulated bundle of one-sided householder reflections can be converted into a data-local matrix-matrix multiplication. The method is known as the WY representation of accumulated Householder Reflections [2], [3]. 

While I use a 3-phase SVD, I do not apply the WY representation. Instead, I conceived a different approach that I consider more flexible an can be applied efficiently to all 3 phases.

## 6 The MOCUT Approach

Lets assume we have a left sided transformation (-bundle) as described in Eq. (4.1).

We replace $ P_i $ by a different sequence $ \prod_j {\tilde{P}_{ij}} $ achieving the same objective, such as for example zeroing a set of columns below the diagonal. Note, we do not require strict equality:  $ P_i \neq  \prod_j {\tilde{P}_{ij}} $ . Only the objective of zeroing certain values in A shall be achieved. $ \tilde{P_{ij} } $ are chosen such that $ \tilde{P_{ij}}A $ modifies only a few consecutive rows of $ A $ beginning from the last rows up to the row at the diagonal. If $ \tilde{P_{ij} } $ modifies k rows in can zero $ k-1 $ consecutive values in column $ i $ of A. This also means that the affected rows have to overlap by exactly one row to zero out all values below the diagonal in column $ i $ of A. $ \tilde{P_{ij} } $ need not be a Householder Reflection, it could just as well be realized by multiple consecutive Givens Rotations at a slightly larger computational expense. We can now reformulate equation (4.1) for the entire bundle:

(5.1)	$ \tilde{P}_{b_a} = \prod_{l = a}^{a + n_b} \prod_{j} \tilde{ P }_{lj}$

Next, we require that $ \tilde{P}_{ik} $ and $ \tilde{P}_{lj} $ are commutative for $ i < l \and k > j $ . This can be achieved by organizing them in an oblique pattern: 

If $ \tilde{P}_{i,k} $ affects rows $ \alpha, ..., \alpha + \delta $ then $ \tilde{P}_{(i+1),k} $ affects rows $ \alpha + 1, ..., min( m, \alpha + \delta + 1) $ .

This gives us freedom to permute the factors in eq. (5.1):

(5.2)	$ \tilde{P}_{b_a} = \prod_{j} \prod_{l = a}^{a + n_b} \tilde{ P }_{lj}$

If we depict the inner product of (5.2) in the block of matrix $ A $ that is to be zeroed, we find that it has the shape of an oblique parallelogram. 

(5.3)	$ P_{t_{{a,j}}}  = \prod_{l = a}^{a + n_b} \tilde{ P }_{lj} $	(MOCUT-Tile)

If we cover A with all required transformations for band-diagonalization, we get a form of tiling. Each tile $ P_{t_{{a,j}}} $ overlaps with the one above or below it by a thin stripe of elements. 

A tile represents a transformation that can be easily implemented in a data-local fashion, because it affects only a sub-block of A of consecutive rows. Since $ P_{t_{{a,j}}} $ is applied to independent column-vectors, we can further partition the block into multiple sub-blocks of adjacent rows. For all values $ l $ in a specified tile, each $ \tilde{ P }_{lj} $ affects overlapping rectangular areas in that sub-block.

From these observations, I formulated the specific name of the approach:

MOCUT = ***Mo**no**c**linic **U**nitary **T**ransformation*

The name is inspired from crystallography, where the term *monoclinic crystal system* describes a special form of oblique crystal cell tiling. 
## Web References

TODO

## Literature

TODO



&copy; Johannes B. Steffens
